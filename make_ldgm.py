import argparse
import collections
import pandas as pd
import ld_graph
import numpy as np
import tskit
import networkx as nx
from tqdm import tqdm
import os


def get_mut_edges(ts):
    muts_to_brick = {}
    node_edge_dict = {}
    for tree, (_, edges_out, edges_in) in zip(ts.trees(), ts.edge_diffs()):
        for edge in edges_out:
            node_edge_dict.pop(edge.child)
        for edge in edges_in:
            node_edge_dict[edge.child] = edge.id
        for site in tree.sites():
            for mut in site.mutations:
                node = mut.node
                if node in node_edge_dict:
                    muts_to_brick[mut.id] = node_edge_dict[node]

    bricks_to_muts = collections.defaultdict(list)

    for mut, brick in muts_to_brick.items():
        bricks_to_muts[brick].append(mut)
    return bricks_to_muts


def simplify_ts(ts, superpop, start, stop, prune_snps):
    ts = ts.keep_intervals([[start, stop]])
    superpops = {
        "EAS": [0, 1, 2, 3, 4],
        "EUR": [5, 6, 7, 8, 9],
        "AFR": [10, 11, 12, 13, 14, 15, 16],
        "AMR": [17, 18, 19, 20],
        "SAS": [21, 22, 23, 24, 25],
    }
    if superpop in superpops:
        superpop_pops = superpops[superpop]

        # Simplify to the superpop samples and keep SNPs above given frequency
        superpop_samples = np.where(np.isin(ts.tables.nodes.population, superpop_pops))[
            0
        ]
        ts_superpop = ts.simplify(samples=superpop_samples)
        return_ts = ld_graph.prune_snps(ts_superpop, threshold=prune_snps)
    elif superpop == "ALL":
        # Keep all samples and SNPs above given frequency (in each superpopulation)
        print("Original ts has {} SNPs".format(ts.num_sites))
        sites_to_keep = []
        for superpop, superpop_pops in superpops.items():
            superpop_samples = np.where(
                np.isin(ts.tables.nodes.population, superpop_pops)
            )[0]
            ts_superpop = ts.simplify(samples=superpop_samples, filter_sites=False)
            a = ts_superpop.num_samples
            cur_pop_keep_sites = 0
            for tree in ts_superpop.trees():
                for site in tree.sites():
                    if len(site.mutations) == 1:
                        freq = tree.num_samples(site.mutations[0].node) / a
                        if freq >= prune_snps and freq <= 1 - prune_snps:
                            sites_to_keep.append(site.id)
                            cur_pop_keep_sites += 1
                    elif len(site.mutations) > 1:
                        raise ValueError("Site cannot have more than one mutation")
            print(
                "Superpop: {} has {} SNPs above 1% frequency".format(
                    superpop, cur_pop_keep_sites
                )
            )
        # tskit has a delete_sites function, so we need to only keep the sites that aren't in sites_to_keep
        all_sites = np.arange(0, ts.num_sites)
        sites_to_keep = np.sort(np.array(list(set(sites_to_keep))))
        sites_to_delete = all_sites[~np.in1d(all_sites, sites_to_keep)]
        print("Total number of SNPs to keep is {}".format(len(sites_to_keep)))
        print("Total number of SNPs being deleted is {}".format(len(sites_to_delete)))
        return_ts = ts.delete_sites(np.array(sites_to_delete))
    else:
        raise ValueError("Must specify a superpop or 'ALL'")

    return return_ts


def make_graphical_model(ts, recombination_threshold, path_threshold, softmin):
    if False:
        bts = ld_graph.brick_ts(
            ts, threshold=recombination_threshold, add_dummy_bricks=True
        )
        BG = ld_graph.brick_graph(bts, threshold=path_threshold)
        mutgraph, _ = ld_graph.reduce_graph(BG, bts, threshold=path_threshold)
    mutgraph, bts = ld_graph.reduce(
        ts,
        path_threshold=path_threshold,
        recombination_threshold=recombination_threshold,
        use_softmin=softmin,
    )

    #####
    # TODO: use code we've been using to create graphical model
    #####
    bricks_to_muts = get_mut_edges(bts)
    SNPs = list(bricks_to_muts.values())
    SNPs = [x[0] for x in SNPs]
    genos = ts.genotype_matrix()
    X = genos[
        np.array(
            [index for index, site in enumerate(ts.sites()) if len(site.mutations) > 0]
        )
    ]
    X = np.array(X, dtype=int)
    snplist, anc_alleles, deriv_alleles = ld_graph.return_snp_list(bts)
    return mutgraph, genos, (snplist, anc_alleles, deriv_alleles)


def iterate_intervals(
    chrom,
    superpop,
    bed_file,
    block,
    ts,
    prune_snps,
    recombination_threshold,
    path_threshold,
    softmin,
    path,
):
    dbsnp_bed_file = pd.read_csv(
        "/broad/oconnor/trees/nygc/bed_chr_" + chrom + ".bed.gz",
        delimiter="\t",
        skiprows=1,
        header=None,
    )
    # Subset to SNPs only in dbSNP
    dbsnp_bed_file = dbsnp_bed_file[
        (dbsnp_bed_file.iloc[:, 2] - dbsnp_bed_file.iloc[:, 1]) == 1
    ]
    # If multiple rsids refer to the same location, keep first
    dbsnp_bed_file = dbsnp_bed_file.drop_duplicates(subset=2, keep="first")
    rsids = dbsnp_bed_file.iloc[:, 3].to_numpy()
    dbsnp_bed_file = dbsnp_bed_file.iloc[:, 2].to_numpy()

    for row, interval in bed_file.iterrows():
        if row + 1 == block:
            print(row, interval)
            if not os.path.exists(
                path
                + "/1kg_"
                + superpop
                + "_chr"
                + chrom
                + "_"
                + str(interval["start"])
                + "_"
                + str(interval["end"])
                + "_MAF_"
                + str(prune_snps)
                + "_RF_"
                + str(recombination_threshold)
                + "_T_"
                + str(path_threshold)
                + ".genos"
            ):
                return_ts = simplify_ts(
                    ts, superpop, interval["start"], interval["end"], prune_snps
                )
                snp_pos = return_ts.tables.sites.position.astype(int)
                ts_ids = np.full(return_ts.num_sites, "NA", dtype=np.dtype("U100"))
                labeled_snps = np.isin(snp_pos, dbsnp_bed_file)
                index = np.argsort(dbsnp_bed_file)
                sorted_dbsnp_bed_file = dbsnp_bed_file[index]
                sorted_index = np.searchsorted(sorted_dbsnp_bed_file, snp_pos)

                yindex = np.take(index, sorted_index, mode="clip")
                mask = dbsnp_bed_file[yindex] != snp_pos

                result = np.ma.array(yindex, mask=mask)
                ts_ids[labeled_snps] = rsids[result.data[~result.mask]]

                (
                    mutgraph,
                    genotypes,
                    (_, anc_alleles, deriv_alleles),
                ) = make_graphical_model(
                    return_ts, recombination_threshold, path_threshold, softmin
                )
                assert (
                    genotypes.shape[0]
                    == len(ts_ids)
                    == len(anc_alleles)
                    == len(deriv_alleles)
                )
                np.savetxt(
                    path
                    + "/1kg_"
                    + superpop
                    + "_chr"
                    + chrom
                    + "_"
                    + str(interval["start"])
                    + "_"
                    + str(interval["end"])
                    + "_MAF_"
                    + str(prune_snps)
                    + "_RF_"
                    + str(recombination_threshold)
                    + "_T_"
                    + str(path_threshold)
                    + ".genos",
                    genotypes,
                )
                nx.write_edgelist(
                    mutgraph,
                    path
                    + "/1kg_"
                    + superpop
                    + "_chr"
                    + chrom
                    + "_"
                    + str(interval["start"])
                    + "_"
                    + str(interval["end"])
                    + "_MAF_"
                    + str(prune_snps)
                    + "_RF_"
                    + str(recombination_threshold)
                    + "_T_"
                    + str(path_threshold)
                    + ".adjlist",
                )
                pd.DataFrame(
                    {
                        "position": snp_pos,
                        "rsid": ts_ids,
                        "anc_allele": anc_alleles,
                        "deriv_allele": deriv_alleles,
                    }
                ).to_csv(
                    path
                    + "/1kg_"
                    + superpop
                    + "_chr"
                    + chrom
                    + "_"
                    + str(interval["start"])
                    + "_"
                    + str(interval["end"])
                    + "_MAF_"
                    + str(prune_snps)
                    + "_RF_"
                    + str(recombination_threshold)
                    + "_T_"
                    + str(path_threshold)
                    + ".snplist"
                )


def main():
    print("starting evaluation")
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "chr", type=int, help="Chromosome from which to build LDGM at intervals"
    )
    parser.add_argument(
        "superpop",
        type=str,
        help='Which superpop to build LDGM from (or ALL to use all populations). Must be "AFR", "AMR", "EAS", "EUR", "SAS" or "ALL"',
    )
    parser.add_argument("bed_input", type=str, help="Path to input bed file.")
    parser.add_argument("block", type=str, help="Which block to run.")
    parser.add_argument("ts_input", type=str, help="Path to input tree sequence.")
    parser.add_argument(
        "output",
        type=str,
        help="Path to output edge list, genotypes, and SNP list (TODO)",
    )
    parser.add_argument(
        "--recombination-threshold",
        "-r",
        type=float,
        default=0.01,
        help="Recombination threshold to use",
    )
    parser.add_argument(
        "--prune-snps",
        type=float,
        default=0.01,
        help="SNP frequency threshold to prune",
    )
    parser.add_argument(
        "--path-threshold", "-p", type=int, default=10, help="Path threshold"
    )
    parser.add_argument(
        "--softmin", "-s", action="store_true", help="If passed, uses softmin"
    )

    args = parser.parse_args()
    print("Starting interval {}".format(args.block))
    chrom = str(args.chr)
    assert chrom in args.bed_input
    assert chrom in args.ts_input
    bed_file = pd.read_csv(args.bed_input, delim_whitespace=True)
    ts = tskit.load(args.ts_input)
    print("Loaded tree sequence")
    iterate_intervals(
        chrom,
        args.superpop,
        bed_file,
        int(args.block),
        ts,
        args.prune_snps,
        args.recombination_threshold,
        args.path_threshold,
        args.softmin,
        args.output,
    )


if __name__ == "__main__":
    main()
