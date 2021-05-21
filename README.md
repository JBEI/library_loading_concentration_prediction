# DIVA DNA Seq Utilities

This package is a repository of tools used by the DIVA DNA Seq team at the Joint BioEnergy Institute (JBEI) to predict the library loading concentration for their MiSeq runs.  It employs a machine learning model to determine the relationship between library fragment sizes, loading conentration and cluster density.

## Installation 

The dependencies of this tool can be installed with:
```
pip install -r requirements.txt
```

## Use

This tool can then be used either in a Jupyter notebook or on the command line. A bioanalyzer ladder file and result file are needed. The tool generates a plot that relates library loading concentration for your library to the MiSeq flow cell cluster density.
