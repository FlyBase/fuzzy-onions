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
production: <path to a directory containing the datasets>
staging: <same, but for datasets from the staging server>

[curation]
proformae: <path to a directory containing FlyBase's proformae files>
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
    "dataset_id": "<the SCEA dataset ID>",
    "symbol": "<Base symbol to use for FlyBase FBlc entities>",
    "cell_types_column": "<Name of the column containing cell types>",
    "corrections": [
    	# List of "corrections" to apply to SCEA annotations
    	{
    		"source": "<Source column to correct>",
    		"destination:" "<Name of the column to add>",
    		"values": [
    			# List of replacements to do
    			[ "Original value", "New value", "Comment" ]
    		]
    	}
    ],
    "conditions": [
         # List of columns in the experiment design table
         # used to select cells belonging to different samples
    ],
    "samples": [
    	# List of biosamples in the datasets
    	{
    		"symbol": "<Suffix to identify the biosample>",
    		"selectors": [
    			# List of values to select cells belonging to
    			# that biosample (in the same order as in the
    			# "Conditions" key above)
    		],
    		"title": "<Free-text name of that biosample>",
    		"stage": "<Developmental stage, as a FBdv term>"
    	}
    ]
}
```


Copying
-------
Fuzzy-Onions is distributed under the terms of the MIT license.