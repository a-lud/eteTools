#!/usr/bin/env python3

# This script is designed to take as input a list of CodeML output files/directories
# and return a single table for comparisons.
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
    etedirs = utility.listDirs(
        args.input
    )  # [<DirEntry 'OGID'>, <DirEntry 'OGID'>, ...]

    # Summary objects for all orthologs
    pd_lrt, dict_summary, pd_branches, pd_beb = utility.parseCodeMl(etedirs)

    # Write tables to file
    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    utility.summaryDictToCsv(dict_summary, args.outdir)
    pd_lrt.to_csv(path_or_buf=os.path.join(args.outdir, "lrt.csv"), index=False)
    if not pd_branches.empty:
        pd_branches.to_csv(
            path_or_buf=os.path.join(args.outdir, "branches.csv"), index=False
        )
    else:
        logging.warning("Branches data frame empty. Will not write any file.")

    pd_beb.to_csv(path_or_buf=os.path.join(args.outdir, "beb.csv"), index=False)

    logging.info("Finished")
