#!/usr/bin/env python3

import logging
import argparse

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- #
# Logging information
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s"
)


def getArgs():
    """Get user arguments and set up parser"""
    desc = """\
    # -------------------------------------------------------- #
    #                 Drop-out LRT comparison                  #
    # -------------------------------------------------------- #

    Little script to take two LRT tables generated by
    eteTools and compare the branch-site LRT statistic to the
    site LRT statistic from a drop-out test. It assumes that the
    site models uses the same samples as the BS tests, but
    without the foreground species present.

    Significant signals of positive selection are reported on
    the foreground species when the BS models report a sig.
    LRT p-value and the Site models do not. This indicates
    that the foregorund branch experiences positive selection
    while no species on the background branch does.

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
        "branchsite",
        help="File path to Branch-Site models summary CSV",
        metavar="/path/to/model-branch-site.csv",
    )
    parser.add_argument(
        "site",
        help="File path to Site models summary CSV",
        metavar="/path/to/model-site.csv",
    )
    parser.add_argument(
        "outcsv", help="Pipeline output directory", metavar="/path/to/out.csv"
    )

    parser.add_argument(
        "-p",
        "--pvalue",
        help="P-value threshold to filter results (default: %(default)s). ",
        metavar=0.05,
        default=0.05,
    )

    args = parser.parse_args()
    return args


def compareLRT(bs, site, pval, outcsv):
    """Import LRT tables and compare p-values for Branch-Site and Site models. Only
    compares the 'bsA/bsA1' Branch-Site models to 'M1a/M2a' and 'M7/M8' models."""

    pd_bs = pd.read_csv(bs).drop(["df", "note"], axis=1)
    pd_site = pd.read_csv(site).drop(["df", "note"], axis=1)

    # Filter for Only M1-M2 and M7-M8 comparisons
    pd_site = pd_site.loc[
        ((pd_site.null == "M1") & (pd_site.alt == "M2"))
        | ((pd_site.null == "M7") & (pd_site.alt == "M8"))
    ]

    # Left-join dataframes on sample ID
    pd_merged = pd.merge(pd_bs, pd_site, on="file")
    print(pd_merged)

    # Outcomes
    conditions = [
        # Signal of PS
        (pd_merged["pval_x"] <= pval) & (pd_merged["pval_y"] > pval),
        # PS on foreground and background OR not significant
        (pd_merged["pval_x"] <= pval) & (pd_merged["pval_y"] <= pval),
        (pd_merged["pval_x"] > pval) & (pd_merged["pval_y"] <= pval),
        (pd_merged["pval_x"] > pval) & (pd_merged["pval_y"] > pval),
        # Poor model file
        (np.isnan(pd_merged["pval_x"])) | (np.isnan(pd_merged["pval_y"])),
    ]

    # Corresponding values for conditions
    values = ["PS_fg", "PS_fg_bg", "PS_bg", "no_PS", "poor_fit"]

    # Assign column values based on condition
    pd_merged["signal"] = np.select(conditions, values)

    # Write output file
    pd_merged.to_csv(outcsv)


if __name__ == "__main__":
    # 1. read in data csv files
    # 2. Filter null and alt columns for respective model comparisons
    # 3. compare LRT of BS to Sites
    #   - BS < alpha < Site
    #   - BS & Site < alpha - not what we want
    #   - Site < alpha < BS - not what we want
    #   - alpha < BS & Site

    args = getArgs()
    compareLRT(args.branchsite, args.site, args.pvalue, args.outcsv)
