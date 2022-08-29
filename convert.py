"""
Convert input data from various sources to samples format.
"""
import argparse
import subprocess
import os
import sys
import math

import allel
import numpy as np
import tsinfer
import attr
import cyvcf2
import pysam
import tqdm
import pandas as pd
import multiprocessing

import tskit

GENERATION_TIME = 25


@attr.s()
class Site(object):
    position = attr.ib(None)
    alleles = attr.ib(None)
    genotypes = attr.ib(None)
    metadata = attr.ib({})
    inference = attr.ib(None)


def run_multiprocessing(args, function):
    """
    Run multiprocessing of sampledata files.
    We use multiple threads by splitting the VCF file into chunks and using the
    vcf_subset function of cyvcf2.
    """
    vcf_fn = args.data_file
    num_processes = args.num_threads
    if num_processes > 1:
        # Split the VCF into chunks
        callset = allel.read_vcf(vcf_fn, fields=["variants/CHROM", "variants/POS"])
        pos_list = callset["variants/POS"]
        chroms = callset["variants/CHROM"]
        assert np.all(chroms == chroms[0])
        chrom = str(chroms[0])

        def get_chromosome_chunks(lst, num_processes):
            length = len(lst)
            n = math.ceil(length / num_processes)
            chunks = list()
            for index, i in enumerate(range(0, length, n)):
                if index != num_processes - 1:
                    chunks.append(
                        (
                            args,
                            args.output_file + str(index),
                            (chrom + ":" + str(lst[i]) + "-" + str(lst[i + n])),
                        )
                    )
                else:
                    chunks.append(
                        (
                            args,
                            args.output_file + str(index),
                            (chrom + ":" + str(lst[i]) + "-" + str(lst[-1])),
                        )
                    )
            return chunks

        chunks = get_chromosome_chunks(pos_list, num_processes)
        chunks_iter = iter(chunks)
        reports = list()
        completed_files = list()
        with multiprocessing.Pool(processes=num_processes, maxtasksperchild=10) as pool:
            for index, row in enumerate(pool.map(function, chunks_iter)):
                reports.append(row)
                print(
                    "Processed Chunk {}: {} with {} sites added.".format(
                        index, chunks[index][2], row["num_sites"]
                    )
                )
                if row["num_sites"] > 0:
                    completed_files.append(index)
                else:
                    os.remove(args.output_file + str(index) + "-lock")

        # Combine reports and print
        master_report = reports[0]
        for report in reports[1:]:
            for var_type, val in report.items():
                master_report[var_type] += val
        print(master_report)

        # Combine sampledata files
        filenames = completed_files
        all_samples = []
        for name in filenames:
            all_samples.append(tsinfer.load(args.output_file + str(name)))
            os.remove(args.output_file + str(name))

        samples = all_samples[0].copy(args.output_file)
        samples.append_sites(*all_samples[1:])
        samples.finalise()
        assert np.all(np.diff(samples.sites_position[:]) > 0)

    else:
        raise ValueError


def make_sampledata(args):
    if isinstance(args, tuple):
        vcf_subset = args[2]
        args[0].output_file = str(args[1])
        args = args[0]
    else:
        vcf_subset = None
    data_provenance = {
        "ancestral_states_url": args.ancestral_states_url,
        "reference_name": args.reference_name,
    }

    # Get the ancestral states.
    fasta = pysam.FastaFile(args.ancestral_states_file)
    # NB! We put in an extra character at the start to convert to 1 based coords.
    ancestral_states = "X" + fasta.fetch(reference=fasta.references[0])
    # The largest possible site position is len(ancestral_states). Positions must
    # be strictly less than sequence_length, so we add 1.
    sequence_length = len(ancestral_states) + 1

    converter_class = {
        "1kg": ThousandGenomesConverter,
        "1kg_chrY": ThousandGenomesChrYConverter, 
    }
    try:
        with tsinfer.SampleData(
            path=args.output_file, num_flush_threads=1, sequence_length=sequence_length
        ) as samples:
            converter = converter_class[args.source](
                args.data_file, ancestral_states, samples, args.target_samples
            )
            if args.metadata_file:
                converter.process_metadata(args.metadata_file, args.progress)
            else:
                converter.process_metadata(args.progress)
            if vcf_subset is not None:
                report = converter.process_sites(
                    vcf_subset=vcf_subset,
                    show_progress=args.progress,
                    max_sites=args.max_variants,
                )
            else:
                report = converter.process_sites(
                    show_progress=args.progress, max_sites=args.max_variants
                )
            samples.record_provenance(
                command=sys.argv[0],
                args=sys.argv[1:],
                data=data_provenance,
            )
            assert np.all(np.diff(samples.sites_position[:]) > 0)
    except Exception as e:
        os.unlink(args.output_file)
        if report["num_sites"] == 0:
            return report
        raise e
    if report["num_sites"] == 0:
        os.unlink(args.output_file)
    return report


def filter_duplicates_target(vcf, target_sites_pos=None):
    """
    Returns the variants from this VCF with duplicate sites filtered
    out. If any site position appears more than once, throw all variants away.
    If target_sites_pos is not None, only returns variants from this VCF which
    are present in the target sampledata file.
    """
    if target_sites_pos is not None:

        def site_in_target(site):
            return site in target_sites_pos

    else:

        def site_in_target(site):
            return True

    row = next(vcf, None)
    bad_pos = -1
    for next_row in vcf:
        if bad_pos == -1 and next_row.POS != row.POS:
            if site_in_target(row.POS):
                yield row
        else:
            if bad_pos == -1:
                bad_pos = row.POS
            elif bad_pos != next_row.POS:
                bad_pos = -1
        row = next_row
    if row is not None and bad_pos != -1 and site_in_target(row.POS):
        yield row


class Converter(object):
    """
    Superclass of converters.
    """

    def __init__(self, data_file, ancestral_states, samples, target_samples=None):
        self.data_file = data_file
        self.ancestral_states = ancestral_states
        self.samples = samples
        if target_samples is not None:
            self.target_sites_pos = set(tsinfer.load(target_samples).sites_position[:])
        else:
            self.target_sites_pos = None
        self.num_samples = -1
        self.num_sites = 0
        # ancestral states counters.
        self.num_no_ancestral_state = 0
        self.num_low_confidence_ancestral_state = 0
        # Counters for genotypes and sites.
        self.num_unphased = 0
        self.num_missing_data = 0
        self.num_invariant = 0
        self.num_indels = 0
        self.num_non_biallelic = 0
        self.num_singletons = 0
        # (n - 1)-tons
        self.num_nmo_tons = 0

    def report(self):
        report_dict = {}
        report_dict["num_sites"] = self.num_sites
        report_dict["unphased"] = self.num_unphased
        report_dict["missing_data"] = self.num_missing_data
        report_dict["invariant"] = self.num_invariant
        report_dict["num_indels"] = self.num_indels
        report_dict["non_biallelic"] = self.num_non_biallelic
        report_dict["no_ancestral_state"] = self.num_no_ancestral_state
        report_dict[
            "low_confidence_ancestral_state"
        ] = self.num_low_confidence_ancestral_state
        report_dict["num_singletons"] = self.num_singletons
        report_dict["num_(n - 1)_tons"] = self.num_nmo_tons
        return report_dict

    def process_metadata(self, metadata_file):
        pass

    def get_ancestral_state(self, position):
        # From the ancestral states README:
        # The convention for the sequence is:
        #    ACTG : high-confidence call, ancestral state supported by other 2 sequences
        #    actg : low-confindence call, ancestral state supported by one sequence only
        #    N    : failure, the ancestral state is not supported by any other sequence
        #    -    : the extant species contains an insertion at this postion
        #    .    : no coverage in the alignment

        ret = None
        # NB: we assume that this array is modified so that the 1-indexed coordinates
        # work correctly!
        ancestral_state = self.ancestral_states[position]
        if ancestral_state in [".", "N", "-"]:
            self.num_no_ancestral_state += 1
            ret = None
            inference = False
        elif ancestral_state in ["a", "c", "t", "g"]:
            self.num_low_confidence_ancestral_state += 1
            ret = ancestral_state.upper()
            inference = True
        else:
            assert ancestral_state in ["A", "C", "T", "G"]
            ret = ancestral_state
            inference = True
        return ret, inference


class VcfConverter(Converter):
    def convert_genotypes(self, row, ancestral_state, inference):
        def return_genotype(allele, ancestral_state):
            if allele == ".":
                return tskit.MISSING_DATA
            else:
                return allele != ancestral_state

        ret = None
        num_diploids = self.num_samples // 2
        a = np.zeros(self.num_samples, dtype=np.int8)
        # Use upper case version of ancestral state (keep original
        # for checking low-confidence ancestral state)
        if ancestral_state is not None:
            all_alleles = set([ancestral_state])
        else:
            all_alleles = set([])
        # Fill in a with genotypes.
        bases = np.array(row.gt_bases)
        for j in range(num_diploids):
            missing = False
            if "|" in bases[j]:
                alleles = bases[j].split("|")
            else: 
                self.num_unphased += 1
                alleles = bases[j]
                #print(row, bases[j])
                #alleles = bases[j].split("/")
            if len(alleles) != 2:
                break
            for allele in alleles:
                if allele == ".":
                    self.num_missing_data += 1
                    missing = True
                else:
                    all_alleles.add(allele)
            if ancestral_state is None:
                ancestral_state = alleles[0]
            a[2 * j] = alleles[0] != ancestral_state
            a[2 * j + 1] = alleles[1] != ancestral_state

            if missing:
                if alleles[0] == ".":
                    a[2 * j] = tskit.MISSING_DATA
                if alleles[1] == ".":
                    a[2 * j + 1] = tskit.MISSING_DATA
        else:
            freq = np.sum(a == 1)
            
            if freq == self.num_samples or freq == 0:
                self.num_invariant += 1
            elif freq == self.num_samples - 1:
                self.num_nmo_tons += 1
            else:
                if len(all_alleles) > 2:
                    self.num_non_biallelic += 1
                    inference = False
                if any(len(allele) != 1 for allele in all_alleles):
                    self.num_indels += 1
                    inference = False
                metadata = {"ID": row.ID, "REF": row.REF}
                if freq == 1:
                    self.num_singletons += 1
                all_alleles.remove(ancestral_state)
                alleles = [ancestral_state, all_alleles.pop()]
                ret = Site(
                    position=row.POS, alleles=alleles, genotypes=a, metadata=metadata, inference=inference
                )
        return ret

    def process_sites(self, vcf_subset=None, show_progress=False, max_sites=None):
        num_data_sites = int(
            subprocess.check_output(
                ["bcftools", "index", "--nrecords", self.data_file]
            )
        )

        progress = tqdm.tqdm(total=num_data_sites, disable=not show_progress)
        self.num_sites = 0
        if vcf_subset is None:
            vcf = cyvcf2.VCF(self.data_file)
        else:
            vcf = cyvcf2.VCF(self.data_file)(vcf_subset)
        exclude_sites = []
        for row in filter_duplicates_target(vcf, self.target_sites_pos):
            ancestral_state, inference = self.get_ancestral_state(row.POS)
            #if ancestral_state is not None:
            site = self.convert_genotypes(row, ancestral_state, inference)
            if site is not None:
                if site.inference is False:
                    exclude_sites.append(site.position)
                self.samples.add_site(
                    position=site.position,
                    genotypes=site.genotypes,
                    alleles=site.alleles,
                    metadata=site.metadata,
                )
                #else:
                #    self.samples.add_site(
                #        position=site.position,
                #        genotypes=site.genotypes,
                #        alleles=site.alleles,
                #        metadata=site.metadata,
                #    )

                progress.set_postfix(used=str(self.num_sites))
                self.num_sites += 1
                if self.num_sites == max_sites:
                    break
            progress.update()
        progress.close()
        np.savetxt(self.data_file.split(".")[0] + ".excluded_sites.csv", np.array(exclude_sites).astype(int), fmt='%i')
        report_dict = self.report()
        return report_dict


class ThousandGenomesConverter(VcfConverter):
    """
    Converts data for the 1000 Genomes.
    """

    def process_metadata(self, metadata_file, show_progress=False):
        """
        Adds the 1000 genomes populations metadata.
        """
        # Based on
        # http://www.internationalgenome.org/faq/which-populations-are-part-your-study/
        populations = [
            ["CHB", "Han Chinese in Beijing, China", "EAS"],
            ["JPT", "Japanese in Tokyo, Japan", "EAS"],
            ["CHS", "Southern Han Chinese", "EAS"],
            ["CDX", "Chinese Dai in Xishuangbanna, China", "EAS"],
            ["KHV", "Kinh in Ho Chi Minh City, Vietnam", "EAS"],
            [
                "CEU",
                "Utah Residents (CEPH) with Northern and Western European Ancestry",
                "EUR",
            ],
            ["TSI", "Toscani in Italia", "EUR"],
            ["FIN", "Finnish in Finland", "EUR"],
            ["GBR", "British in England and Scotland", "EUR"],
            ["IBS", "Iberian Population in Spain", "EUR"],
            ["YRI", "Yoruba in Ibadan, Nigeria", "AFR"],
            ["LWK", "Luhya in Webuye, Kenya", "AFR"],
            ["GWD", "Gambian in Western Divisions in the Gambia", "AFR"],
            ["MSL", "Mende in Sierra Leone", "AFR"],
            ["ESN", "Esan in Nigeria", "AFR"],
            ["ASW", "Americans of African Ancestry in SW USA", "AFR"],
            ["ACB", "African Caribbeans in Barbados", "AFR"],
            ["MXL", "Mexican Ancestry from Los Angeles USA", "AMR"],
            ["PUR", "Puerto Ricans from Puerto Rico", "AMR"],
            ["CLM", "Colombians from Medellin, Colombia", "AMR"],
            ["PEL", "Peruvians from Lima, Peru", "AMR"],
            ["GIH", "Gujarati Indian from Houston, Texas", "SAS"],
            ["PJL", "Punjabi from Lahore, Pakistan", "SAS"],
            ["BEB", "Bengali from Bangladesh", "SAS"],
            ["STU", "Sri Lankan Tamil from the UK", "SAS"],
            ["ITU", "Indian Telugu from the UK", "SAS"],
        ]

        population_id_map = {}
        for pop in populations:
            pop_id = self.samples.add_population(
                dict(zip(["name", "description", "super_population"], pop))
            )
            population_id_map[pop[0]] = pop_id

        with open(metadata_file, "r") as ped_file:
            # Parse the individual metadata out of the ped file.
            columns = next(ped_file).split("\t")
            sane_names = [col.replace(" ", "_").lower().strip() for col in columns]
            metadata = {}
            populations = {}
            for line in ped_file:
                row = dict(zip(sane_names, line.strip().split("\t")))
                name = row["individual_id"]
                population_name = row.pop("population")
                populations[name] = population_id_map[population_name]
                # The value '0' seems to be used to encode missing, so insert None
                # instead to be more useful.
                nulled = {}
                for key, value in row.items():
                    if value == "0":
                        value = None
                    nulled[key] = value
                metadata[name] = nulled

        vcf = cyvcf2.VCF(self.data_file)
        individual_names = list(vcf.samples)
        vcf.close()
        self.num_samples = len(individual_names) * 2
        # Add in the metadata rows in the order of the VCF.
        for index, name in enumerate(individual_names):
            self.samples.add_individual(
                metadata=metadata[name], population=populations[name], ploidy=2
            )


def main():
    parser = argparse.ArgumentParser(
        description="Script to convert VCF files into tsinfer input."
    )
    parser.add_argument(
        "source",
        choices=["1kg"],
        help="The source of the input data.",
    )
    parser.add_argument("data_file", help="The input data file pattern.")
    parser.add_argument(
        "ancestral_states_file", help="A vcf file containing ancestral allele states. "
    )
    parser.add_argument("output_file", help="The tsinfer output file")
    parser.add_argument(
        "-m",
        "--metadata_file",
        default=None,
        help="The metadata file containing population and sample data",
    )
    parser.add_argument(
        "-n",
        "--max-variants",
        default=None,
        type=int,
        help="Keep only the first n variants",
    )
    parser.add_argument(
        "--target-samples",
        default=None,
        help="A target sampledata file, only variants present in this target file will \
            will be used in the resulting sampledata file",
    )
    parser.add_argument(
        "-p",
        "--progress",
        action="store_true",
        help="Show progress bars and output extra information when done",
    )
    parser.add_argument(
        "--ancestral-states-url",
        default=None,
        help="The source of ancestral state information for provenance.",
    )
    parser.add_argument(
        "--reference-name",
        default=None,
        help="The name of the reference for provenance.",
    )
    parser.add_argument(
        "--num-threads", type=int, default=1, help="Number of threads to use."
    )

    args = parser.parse_args()

    if args.num_threads > 1:
        run_multiprocessing(args, make_sampledata)
    else:
        report = make_sampledata(args)
        print(report)


if __name__ == "__main__":
    main()
