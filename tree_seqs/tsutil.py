"""
Various utilities for manipulating tree sequences and running tsinfer.
"""
import argparse
import json
import pandas as pd

import tskit
import tsinfer
import daiquiri
import numpy as np


def run_simplify(args):
    ts = tskit.load(args.input)
    ts = ts.simplify()
    ts.dump(args.output)


def remove_trios(args):
    sampledata = tsinfer.load(args.input)
    trio_info = pd.read_csv(args.trio_data, skiprows=23, delimiter="\t")
    trio_names = trio_info["SAMPLE_ID"].values
    remove_inds = []
    for index, ind in enumerate(sampledata.individuals()):
        if ind.metadata["individual_id"] in trio_names:
            remove_inds.append(index)
    retain_inds = np.where(
        ~np.isin(np.arange(0, sampledata.num_individuals), remove_inds)
    )[0]
    sampledata_notrios = sampledata.subset(individuals=retain_inds)
    output_notrios = sampledata_notrios.copy(args.output)
    output_notrios.finalise()


def remove_trios_ts(args):
    ts = tskit.load(args.input)
    trio_info = pd.read_csv(args.trio_data, skiprows=23, delimiter="\t")
    trio_names = trio_info["SAMPLE_ID"].values
    remove_nodes = []
    for index, ind in enumerate(ts.individuals()):
        if json.loads(ind.metadata)["individual_id"] in trio_names:
            remove_nodes.append(ind.nodes)
    remove_nodes = np.concatenate(remove_nodes)
    retain_inds = np.where(~np.isin(ts.samples(), remove_nodes))[0]
    ts_notrios = ts.simplify(retain_inds)
    ts_notrios.dump(args.output)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = "command"

    subparser = subparsers.add_parser("simplify")
    subparser.add_argument("input", type=str, help="Input tree sequence")
    subparser.add_argument("output", type=str, help="Output tree sequence")
    subparser.set_defaults(func=run_simplify)

    subparser = subparsers.add_parser("remove_trios")
    subparser.add_argument("input", type=str, help="Input sampledata")
    subparser.add_argument(
        "trio_data",
        type=str,
        help="Path to file containing information of individuals in trios",
    )
    subparser.add_argument("output", type=str, help="Output sampledata")
    subparser.set_defaults(func=remove_trios)

    subparser = subparsers.add_parser("remove_trios_ts")
    subparser.add_argument("input", type=str, help="Input ts")
    subparser.add_argument(
        "trio_data",
        type=str,
        help="Path to file containing information of individuals in trios",
    )
    subparser.add_argument("output", type=str, help="Output ts")
    subparser.set_defaults(func=remove_trios_ts)

    daiquiri.setup(level="INFO")

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
