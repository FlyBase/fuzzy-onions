Fuzzy-Onions - Helper scripts for the FlyBase scRNAseq project
==============================================================

Fuzzy-Onions (*fzo*) is a set of Python modules and helper scripts to
facilitate the annotation and curation of single-cell RNA-seq (scRNAseq)
datasets for exploitation by [FlyBase](https://flybase.org/).


Configuration
-------------
Fuzzy-Onions expects a configuration file under the name
`$XDG_CONFIG_HOME/fuzzyonions/config` (on GNU/Linux) or
`~/Library/Application Support/fuzzyonions/config` (on Mac OS).

Here is a sample configuration file:

```
[store]
directory: <path to a directory containing the datasets>

[proformae]
directory: <path to a directory containing FlyBase's proformae files>
```


Basic usage
-----------
Call `fzo -h` to get a list of the available subcommands, and
`fzo <subcommand> -h` to get the help message for a given subcommand.

The main commands are:

* `download`, to download a dataset from the SCEA;

* `list`, to list the locally available datasets;

* `explorer`, to enter a mode allowing to interactively explore the
contents of a dataset;
  
* `extract`, `sumexpr`, and `proforma`, to extract useful data for
curation from a dataset.


Dataset description file
------------------------
Several subcommands expect a JSON-formatted file describing the dataset
to use.

Here is a sample description file:

```
{
    "Dataset ID": "<the SCEA dataset ID>",
    "Symbol": "<Base symbol to use for FlyBase FBlc entities>",
    "Cell types column": "<Name of the column containing cell types>",
    "Corrections": [
    	# List of "corrections" to apply to SCEA annotations
    	{
    		"Source": "<Source column to correct>",
    		"Destination:" "<Name of the column to add>",
    		"Values": [
    			# List of replacements to do
    			[ "Original value", "New value", "Comment" ]
    		]
    	}
    ],
    "Conditions": [
         # List of columns in the experiment design table
         # used to select cells belonging to different samples
    ],
    "Samples": [
    	# List of biosamples in the datasets
    	{
    		"Symbol": "<Suffix to identify the biosample>",
    		"Selectors": [
    			# List of values to select cells belonging to
    			# that biosample (in the same order as in the
    			# "Conditions" key above)
    		],
    		"Title": "<Free-text name of that biosample>",
    		"Stage": "<Developmental stage, as a FBdv term>"
    	}
    ]
}
```


Copying
-------
Fuzzy-Onions is distributed under the terms of the 1-clause BSD license.