import argparse
import logging
import re
from os import path, scandir
from pathlib import PurePath

import pandas as pd
from Bio.Phylo.PAML import codeml

from lib import EteResults


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

    args = parser.parse_args()
    return args


def listDirs(path):
    """Return the sub-directories in the user provided input path. These
    sub-directories should correspond to the MSA/Genes that were run, with
    each sub-directory containing the CodeML model outputs."""
    subdirs = []
    [subdirs.append(f) for f in scandir(path) if f.is_dir()]
    return subdirs


def getName(path):
    """Get the name of the current MSA"""
    path_split = PurePath(path).parts
    return path_split[-1]


def getModels(path):
    """Get the models that were run for each gene. ETE3 evol uses
    the model as the prefix of the output directory."""
    dirs = listDirs(path)

    # Iterate over each model and parse first section before
    # tilde
    models = []
    for d in dirs:
        tmp = d.name.split("~", 1)[0]

        if "." in tmp:
            m = tmp.split(".", 1)[0]
            models.append(m)
        else:
            models.append(tmp)

    return models


def readCodemlOut(path):
    """Read CodeML outputs into dictionary structures for current MSA"""
    dirs = listDirs(path)

    # Iterate over models
    cml_dict = {}
    for model in dirs:
        m = model.name.split("~", 1)[0]

        if "." in m:
            m = m.split(".", 1)[0]

        # Path to 'out' file
        p = f"{model.path}/out"

        # Codeml dict
        cml = codeml.read(p)

        # Get np values
        with open(p, "r") as f:
            lines = f.readlines()
            match = int([lines.index(i) for i in lines if i.startswith("lnL(")][0])
            match = lines[match].rstrip()
            np = parseNp(st=match)

        cml["np"] = np
        cml_dict[m] = cml

    return cml_dict


def parseNp(st):
    """Extract the NP values from CodeML output"""
    ex = re.search(".+np:(.*)\\):.+", st).group(1).lstrip()
    return ex


def getSiteClasses(input):
    """Convert the 'site classes' field into a pandas data frame for Site models."""
    lst_df = []
    for key, value in input.items():
        df = pd.DataFrame([value])
        df.columns = [str(col) + "_" + str(key) for col in df.columns]
        lst_df.append(df)
    df = pd.concat(lst_df, axis=1)
    return df


def getSiteClassesBranchSite(input):
    """Convert the 'site classes' field into a pandas data frame for Branch-Site models."""
    df_lst = []
    for key, value in input.items():
        for k, v in value.items():
            if isinstance(v, dict):
                df = pd.DataFrame([v])
                df.columns = [str(col) + "_" + str(key) for col in df.columns]
                df_lst.append(df)
            else:
                df = pd.DataFrame([{k: v}])
                df.columns = [str(col) + "_" + str(key) for col in df.columns]
                df_lst.append(df)
    df = pd.concat(df_lst, axis=1)
    return df


def getSiteClassesClade(input):
    """Convert the 'site classes' field into a pandas data frame for Clade models."""
    df_lst = []
    for key, value in input.items():
        for k, v in value.items():
            if isinstance(v, dict):
                df = pd.DataFrame([v])
                df.columns = [
                    k.replace(" ", "-") + "_" + str(key) + "_" + str(col)
                    for col in df.columns
                ]
                df_lst.append(df)
            else:
                df = pd.DataFrame([{k: v}])
                df.columns = [str(col) + "_" + str(key) for col in df.columns]
                df_lst.append(df)
    df = pd.concat(df_lst, axis=1)
    return df


def getBranchResults(input, file):
    """Build pandas dataframe from Branch information in CodeML output files for Null, Branch-Site and Site models"""
    branches_list = []
    for br, val in input.items():
        d = pd.DataFrame([val])
        d.insert(0, "file", file)
        d.insert(1, "branch", br)
        branches_list.append(d)
    return pd.concat(branches_list)


def buildSummaryTable(input):
    """Build a summary Pandas table for each model class. Empty fields will be removed"""

    # 1. clean dict of empty values
    ret = {}
    for key, value in input.items():
        if len(value) == 0:
            continue
        else:
            ret[key] = pd.concat(value)
    return ret


def mergeSummaryDicts(input):
    """Given an list of dictionaries of arbitary length, append the 'values' of each dict
    into a list for matching keys."""
    ret = {"null": [], "site": [], "branch-site": [], "clade": [], "branch": []}

    # Append each datatframe to list
    for d in input:
        for key, value in d.items():
            ret[key].append(value)

    # Concatenate all dataframes for each key
    for key, value in ret.items():
        if len(value) == 0:
            continue
        else:
            ret[key] = pd.concat(value)

    # Remove empty keys
    keys_to_pop = [k for k, v in ret.items() if len(v) == 0]
    for key in keys_to_pop:
        ret.pop(key, None)

    return ret


def parseCodeMl(input):
    """Wrapper function for the functions that do all the work."""

    logging.info("[parseCodeMl] Building LRT, summary and branch tables")

    # Aggregated output structures
    lrt = []
    branches = []
    summary = []

    # Iterate over each ortholog output directory
    for outdir in input:
        l = EteResults.EteResults(outdir.path).getLRT()
        s, b = EteResults.EteResults(outdir.path).getSummary()

        # Append to branches list object
        lrt.append(l)
        summary.append(s)
        branches.append(b)

    # Return LRT table, summary dicts by model type and concatenated branch dataframe
    return pd.concat(lrt), mergeSummaryDicts(summary), pd.concat(branches)


def summaryDictToCsv(input, outdir):
    """Write to file each table in the summary dictionary, using the key as the filename."""
    for key, table in input.items():
        table.to_csv(path_or_buf=path.join(outdir, f"model-{key}.csv"), index=False)
