# EteTools

This is a simple tool I've written to parse output from [`ETE3 evol`][ete] and write summary CSV files.
The tools uses [BioPython][bpy] to interface with the CodeML output files. Calculation of the LRT statistic
has been based off the [method used][lrt] in the ETE3 evol software.

# Required python modules

The tools requires the following non-standard python libraries.

- [Pandas][pd]
- [SciPy][spy]
- [BioPython][bpy]

## How to use

`EteTools` has been designed to be used in a custom Nextflow pipeline, but can work on any output directory
generated by `ete evol`. The tool parses the `CodeML` output file in each model directory to get key
information. It then uses this information to:

- Perform LRT statistics for all valid comparisons (comparisons the same as those from `ETE3 evol`)
- Create concatenated summary tables for each model type
- Report the sites under selection (BEB) - TODO

# Nextflow pipeline

I've implemented this tool in my nextflow pipeline that can be found [here][nf]


[ete]: https://github.com/etetoolkit/ete
[bpy]: https://biopython.org/wiki/PAML
[lrt]: https://github.com/etetoolkit/ete/blob/2b207357dc2a40ccad7bfd8f54964472c72e4726/ete3/evol/evoltree.py#L87
[pd]: https://pandas.pydata.org/
[spy]: https://scipy.org/
[nf]: https://github.com/a-lud/nf-pipelines