import argparse
import tskit
import tsinfer
import pandas as pd

# Script to run tsinfer
# Creates ancestors tree sequences with biallelic SNPs
# Then maps other variants onto the tree sequence during match_samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "chr", type=int, help="Chromosome from which to build LDGM at intervals"
    )
    parser.add_argument("num_threads", type=int, help="Number of threads to use")

    args = parser.parse_args()

    chr = str(args.chr)

    df = pd.read_csv(
        "1kg_phased_genotypes_chr" + chr + ".excluded_sites.csv", header=None
    )
    exclude = df.values.reshape(df.shape[0])
    samples = tsinfer.load("1kg_chr" + chr + ".samples")
    anc = tsinfer.generate_ancestors(
        samples,
        exclude_positions=exclude,
        progress_monitor=True,
        num_threads=args.num_threads,
    )
    copy = anc.copy("1kg_chr" + chr + ".ancestors")
    copy.finalise()
    anc = tsinfer.load("1kg_chr" + chr + ".ancestors")
    ats = tsinfer.match_ancestors(
        samples, anc, progress_monitor=True, num_threads=args.num_threads
    )
    ats.dump("1kg_chr" + chr + ".ats")
    ats = tskit.load("1kg_chr" + chr + ".ats")
    ts = tsinfer.match_samples(
        samples, ats, progress_monitor=True, num_threads=args.num_threads
    )
    ts.dump("1kg_chr" + chr + ".trees")


if __name__ == "__main__":
    main()
