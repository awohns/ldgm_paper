"""
Simple analysis script producing a simulated and inferred tree sequence.
LDGMs and LDGM precision matrices are subsequently created from these files.
"""

import argparse
import msprime
import stdpopsim
import tskit
import tsinfer
import ldgm
import numpy as np


def prune_sites(ts, prune_snps):
    sites_to_keep = []
    for tree in ts.trees():
        for site in tree.sites():
            if len(site.mutations) == 1 and len(ts.site(site.id).mutations) == 1:
                freq = tree.num_samples(site.mutations[0].node) / ts.num_samples
                if freq >= prune_snps and freq <= 1 - prune_snps:
                    sites_to_keep.append(site.id)
                elif len(site.mutations) > 1:
                    pass
    # tskit has a delete_sites function, so we need to only keep the sites that aren't in sites_to_keep
    all_sites = np.arange(0, ts.num_sites)
    sites_to_keep = np.sort(np.array(list(set(sites_to_keep))))
    sites_to_delete = all_sites[~np.in1d(all_sites, sites_to_keep)]
    return ts.delete_sites(np.array(sites_to_delete)).simplify()


def create_edgelist(
    ts, prefix, path_threshold, recombination_threshold, output_fn, seed
):
    mutgraph, bts = ldgm.make_ldgm(
        ts,
        path_weight_threshold=path_threshold,
        recombination_freq_threshold=recombination_threshold,
        progress=True,
        num_processes=10,
        chunksize=10,
    )
    edgelist = ldgm.return_edgelist(mutgraph)
    edgelist.to_csv(output_fn + ".ldgm.edgelist", index=None)
    population_dict = {"sim_pop": list(range(0, ts.num_samples))}
    snplist = ldgm.make_snplist(bts, population_dict=population_dict).round(decimals=4)
    snplist.to_csv(output_fn + ".snplist", index=None)


def main():
    parser = argparse.ArgumentParser(description="Script to run simulated analysis.")
    parser.add_argument(
        "-s",
        "--seed",
        default=1,
        type=int,
        help="Seed for simulation",
    )
    args = parser.parse_args()
    seed = args.seed

    ts = msprime.simulate(
        5008,
        mutation_rate=1.2e-8,
        recombination_rate=1e-8,
        Ne=10000,
        length=3e4,#3e6,
        random_seed=seed,
    )
    ts.dump("msprime_ts_seed_" + str(seed) + ".trees")

    path_threshold = 8
    recombination_threshold = 0.01
    mut_threshold = 0.01
    pruned_ts = prune_sites(ts, mut_threshold)
    genotypes = pruned_ts.genotype_matrix()

    output_fn = (
        "sim_"
        + str(seed)
        + "_MAF_"
        + str(mut_threshold)
        + "_RF_"
        + str(recombination_threshold)
        + "_T_"
        + str(path_threshold)
    )
    np.savetxt(output_fn + ".genos", genotypes.astype(int), fmt="%i", delimiter=",")
    print(
        "Simulated ts has {} sites, pruned ts has {} sites".format(
            ts.num_sites, pruned_ts.num_sites
        )
    )
    create_edgelist(
        pruned_ts,
        "msprime",
        path_threshold=path_threshold,
        recombination_threshold=recombination_threshold,
        output_fn=output_fn,
        seed=seed,
    )

    sampledata = tsinfer.SampleData.from_tree_sequence(
        ts, use_sites_time=False, use_individuals_time=False
    )
    inferred_ts = tsinfer.infer(sampledata, progress_monitor=True)
    inferred_ts.dump("tsinfer_ts_seed_" + str(seed) + ".trees")
    inferred_ts = tskit.load("tsinfer_ts_seed_" + str(seed) + ".trees")
    inferred_pruned_ts = prune_sites(inferred_ts, mut_threshold)
    genotypes = pruned_ts.genotype_matrix()
    output_fn = (
        "inferred_sim_"
        + str(seed)
        + "_MAF_"
        + str(mut_threshold)
        + "_RF_"
        + str(recombination_threshold)
        + "_T_"
        + str(path_threshold)
    )
    np.savetxt(output_fn + ".genos", genotypes.astype(int), fmt="%i", delimiter=",")

    print(
        "Inferred ts has {} sites, pruned inferred ts has {} sites".format(
            inferred_ts.num_sites, inferred_pruned_ts.num_sites
        )
    )
    create_edgelist(
        inferred_pruned_ts,
        "tsinfer",
        path_threshold=path_threshold,
        recombination_threshold=recombination_threshold,
        output_fn=output_fn,
        seed=seed,
    )


if __name__ == "__main__":
    main()
