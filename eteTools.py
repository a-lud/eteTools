#!/usr/bin/env python3

# This script is designed to take as input a list of CodeML output files/directories
# and return a single table for comparisons.
import argparse
import logging
import os
from pathlib import Path

from lib import utility

# --------------------------------------------------------------------------- #
# Logging information
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s"
)

# getArgs Set up argument parser.
def getArgs():
    """Get user arguments and set up parser"""
    desc = """\
    # -------------------------------------------------------- #
    #                ETE3 Evol output-to-table                 #
    # -------------------------------------------------------- #

    This is a simple script that parses the output of ETE3 evol
    and returns a series of informative tables. This tool
    provides a little more flexibility than just using the std-
    out from the ETE3 evol tool.

    ------------------------------------------------------------
    """

    epi = """\
    Code written by Alastair J. Ludington
    University of Adelaide
    2022
    """

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=desc,
        epilog=epi,
    )

    # Required, positional input file arguments
    parser.add_argument(
        "input",
        help="Directory path to ETE3 evol results",
        metavar="/path/to/input",
    )
    parser.add_argument(
        "outdir", help="Pipeline output directory", metavar="/path/to/outdir"
    )

    # Optional argument - specify directories (models)
    parser.add_argument(
        "-m",
        "--models",
        help="Pass directories (models) you want to parse data for",
        nargs="*",
        default=None,
    )

    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = getArgs()
    outdirs = utility.listDirs(args.input)

    # Get summary information
    pd_lrt, dict_summary, pd_branches = utility.parseCodeMl(outdirs)

    # Write tables to file
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    utility.summaryDictToCsv(dict_summary, args.outdir)
    pd_lrt.to_csv(path_or_buf=os.path.join(args.outdir, "lrt.csv"), index=False)
    pd_branches.to_csv(
        path_or_buf=os.path.join(args.outdir, "branches.csv"), index=False
    )

    logging.info("Finished")
