from ctypes import util
import logging
from lib import utility
import pandas as pd
from scipy import stats


class EteResults:
    """EteResults class contains information relating to an ETE3 evol output directory"""

    def __init__(self, path=""):
        """Class to store the results from an ETE3 evol run. Expects that the
        parent directory is a unique name corresponding to the MSA that was
        used to generate the results - e.g. the ortholog identifier/gene name"""
        self.path = path
        self.name = utility.getName(path)
        self.models = utility.getModels(path)
        self.codeml = utility.readCodemlOut(path)

    def getLRT(self):
        """Build a table from the CodeML results for each model"""
        # Dictionary valid comparisons - {null_model: [alternative_models]}
        null_alt = {
            "M0": ['M1', 'M2', 'M8', 'M7', 'bsA', 'bsA1', 'bsC', 'bsD', 'b_free'],
            'M1': ['M2', 'M8', 'bsA', 'bsA1', 'bsC', 'bsD'],
            'M2': ['bsC', 'bsD'],
            'M8': ['bsC', 'bsD'],
            'M7': ['M2', 'M8', 'bsA', 'bsA1', 'bsC', 'bsD'],
            'bsA': ['bsC', 'bsD'],
            'bsA1': ['M2', 'M8', 'bsA', 'bsC', 'bsD'],
            'bsC': ['bsD'],
            'b_free': ['bsA', 'bsA1', 'bsC', 'bsD'],
            'b_neut': ['bsA', 'bsA1', 'bsC', 'bsD', 'b_free']
        }

        # 1. Build dictionary of information required for LRT statistic
        codeml = {}
        for model in self.codeml:
            cml = self.codeml[model]  # Codeml dict

            # Key value is dependent on number of sites for model
            k = list(cml['NSsites'].keys())[0]

            # Comparison information
            lnl = cml['NSsites'][k]['lnL']
            np = cml['np']

            # dictionary that can be turned into data frame
            data = {
                "np": np,
                "Lnl": lnl,
            }

            codeml[model] = data

        # 2. Run LRT statistic for each valid model comparison
        lrt_list = []
        for model in self.models:
            if model in null_alt.keys():
                # keep only comparisons that we have data for
                alt_models = {key: codeml[key] for key in null_alt[model]}

                # LRT statistic for null_model to all alt_models
                for key, values in alt_models.items():

                    df_alt = int(values['np'])
                    df_null = int(codeml[model]['np'])

                    lnl1 = values['Lnl']
                    lnl0 = codeml[model]['Lnl']

                    degFree = abs(df_alt - df_null)  # degrees of freedom

                    if lnl1 < lnl0:
                        lrt_list.append(pd.DataFrame([{
                            'file': self.name,
                            'null': model,
                            'alt': key,
                            'df': degFree,
                            'pval': "",
                            'note': 'lnl1 < lnl0'
                        }]))
                    else:
                        # Data we need to conduct the LRT statistic
                        diff = 2 * (lnl1 - lnl0)
                        # Taken from ETE3 source
                        pval = 1 - stats.chi2.cdf(diff, degFree)

                        # Append to dataframe
                        lrt_list.append(pd.DataFrame([{
                            'file': self.name,
                            'null': model,
                            'alt': key,
                            'df': degFree,
                            'pval': pval,
                            'note': ''
                        }]))

        # Concatenate all model comparisons into a single data frame (by alignment)
        lrt = pd.concat(lrt_list) 
        return lrt

    def getSummary(self):
        """Generate a summary table of all the models run for a MSA."""
        model_types = {
            'null': ['M0'],
            'site': ['M1', 'M2', 'M7', 'M8'],
            'branch-site': ['bsA', 'bsA1'],
            'clade': ['bsC', 'bsD'],
            'branch': ['b_free', 'b_neut']
        }

        
        # Output objects
        summary = {
            'null': [],
            'site': [],
            'branch-site': [],
            'clade': [],
            'branch': []
        }
        branch = []

        for model in self.codeml:

            # Get model type
            modelType = [k for k, v in model_types.items() if model in v][0]

            # Relevant CodeML dictionary
            cml = self.codeml[model]

            # Key value is dependent on number of sites for model
            k = list(cml['NSsites'].keys())[0]

            # Custom field
            np = cml['np']

            # Top level info
            mdl = cml['model']
            codonModel = cml['codon model']

            # NSsites info
            description = cml['NSsites'][k]['description']
            treeLen = cml['NSsites'][k]['tree length']
            lnl = cml['NSsites'][k]['lnL']

            # Parameters info
            kappa = cml['NSsites'][k]['parameters']['kappa']

            # Informative model information
            data = pd.DataFrame([{
                'file': self.name,
                'model': model,
                'model-long': mdl,
                'codon model': codonModel,
                'description': description,
                'tree-length': treeLen,
                "Lnl": lnl,
                "np": np,
                'kappa': kappa
            }])

            # Get model specific information
            match modelType:
                case 'null':
                    tmp = pd.DataFrame([{
                        'omega': cml['NSsites'][k]['parameters']['omega'],
                        'dN': cml['NSsites'][k]['parameters']['dN'],
                        'dS': cml['NSsites'][k]['parameters']['dS']
                    }])
                    tmp_branch = utility.getBranchResults(cml['NSsites'][k]['parameters']['branches'], self.name)
                    
                    # Join model specific data to objects and append to output object
                    branch.append(data.merge(tmp_branch, on='file', how='left'))
                    summary[modelType].append(pd.concat([data, tmp], axis=1))

                # Ugly as all heck but it gets the job done
                case 'site':
                    tmp = pd.DataFrame([{
                        'siteClassModel': cml['site-class model'],
                        'p0': cml['NSsites'][k]['parameters']['p0'] if model == 'M8' else '',
                        'w': cml['NSsites'][k]['parameters']['w'] if model == 'M8' else '',
                        'p': cml['NSsites'][k]['parameters']['p'] if model in ['M7', 'M8'] else '',
                        'q': cml['NSsites'][k]['parameters']['q'] if model in ['M7', 'M8'] else ''
                    }])

                    tmp_sc = utility.getSiteClasses(cml['NSsites'][k]['parameters']['site classes'])
                    tmp_branch = utility.getBranchResults(cml['NSsites'][k]['parameters']['branches'], self.name)
                    
                    # Join model specific data to objects and append to output object
                    branch.append(data.merge(tmp_branch, on='file', how='left'))

                    summary[modelType].append(pd.concat([data, tmp, tmp_sc], axis=1))
                    
                case 'branch-site':
                    tmp = pd.DataFrame([{'siteClassModel': cml['site-class model'] }])
                    tmp_sc = utility.getSiteClassesBranchSite(cml['NSsites'][k]['parameters']['site classes'])
                    
                    summary[modelType].append(pd.concat([data, tmp, tmp_sc], axis=1))

                case 'clade':
                    tmp = pd.DataFrame([{'siteClassModel': cml['site-class model'] }])
                    tmp_sc = utility.getSiteClassesClade(cml['NSsites'][k]['parameters']['site classes'])
                    summary[modelType].append(pd.concat([data, tmp, tmp_sc], axis=1))

                case 'branch':
                    tmp = pd.DataFrame([{
                    'omega1': cml['NSsites'][k]['parameters']['omega'][0],
                    'omega2': cml['NSsites'][k]['parameters']['omega'][1],
                    'dN': cml['NSsites'][k]['parameters']['dN'],
                    'dS': cml['NSsites'][k]['parameters']['dS']
                    }])
                    tmp_branch = utility.getBranchResults(cml['NSsites'][k]['parameters']['branches'], self.name)

                     # Join model specific data to objects and append to output object
                    branch.append(data.merge(tmp_branch, on='file', how='left'))
                    summary[modelType].append(pd.concat([data, tmp], axis=1))
        
        # Return summary tables dict and single branch df
        return utility.buildSummaryTable(summary), pd.concat(branch)
