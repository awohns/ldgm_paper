"""
Simple analysis script producing a simulated and inferred tree sequence.
LDGMs and LDGM precision matrices are subsequently created from these files.
"""


import msprime
import stdpopsim
import tskit
import tsinfer

for i in range(1, 21):
    ts = msprime.simulate(5008, mutation_rate=1.2e-8, recombination_rate=1e-8, Ne=10000, length=3e6, random_seed=i)
    ts.dump("simulated_data/msprime_ts_seed_" + str(i) + ".trees")

    sampledata = tsinfer.SampleData.from_tree_sequence(ts, use_sites_time=False, use_individuals_time=False)
    inferred_ts = tsinfer.infer(sampledata, progress_monitor=True)
    inferred_ts.dump("simulated_data/tsinfer_ts_seed_" + str(i) + ".trees")
