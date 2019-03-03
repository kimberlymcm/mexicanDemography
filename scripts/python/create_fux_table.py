"""
   Creates the fux table based on:
     divergence: human and chimp divergence
     tri_freq: trinucleotide frequencies
     transition_mat: Trinuclotide transition rate matrix inferred along the human lineage by
     Hwang and Green, PNAS 101:13994 (2004)
  Fux_table: table of misidentification probabilities.
  More details at: https://bitbucket.org/gutenkunstlab/dadi
"""

import sys
import os.path
import argparse

import dadi


TRANS_MAT_FN = "/home/kfm/kfm_projects/NA/NA_data/getIntrons/dadi_other/Q.HwangGreen.human.dat"
TRI_FREQ_FN = "/home/kfm/kfm_projects/NA/NA_data/getIntrons/dadi_other/tri_freq.dat"


def main(args):

    outfn = args.out_prefix + str(args.divergence) + ".dat"
    if os.path.isfile(outfn):
        sys.exit("Output file already exists")

    Q = dadi.Numerics.array_from_file(TRANS_MAT_FN)
    tri_freq = dict((line.split()[0], float(line.split()[1]))
                    for line in file(TRI_FREQ_FN).readlines())

    outfn = args.out_prefix + str(args.divergence) + ".dat"
    dadi.Misc.make_fux_table(outfn, args.divergence, Q, tri_freq)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create fux table")
    parser.add_argument("--out_prefix",
        default="/home/kfm/kfm_projects/NA/NA_data/getIntrons/dadi_other/fux_table_")
    parser.add_argument("--divergence", default=0.0121, type=float)
    args = parser.parse_args()
    main(args)

