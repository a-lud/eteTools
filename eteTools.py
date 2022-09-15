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

if __name__ == "__main__":

    args = utility.getArgs()
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
