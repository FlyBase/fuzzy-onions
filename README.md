Fuzzy-Onions - Helper scripts for the FlyBase scRNAseq project
==============================================================

Fuzzy-Onions (*fzo*) is a set of Python modules and helper scripts to
facilitate the annotation and curation of single-cell RNA-seq (scRNAseq)
datasets for exploitation by [FlyBase](https://flybase.org/).


Setup
-----
Fuzzy-Onions is currently _not_ available on PyPI, and will likely never
be due to its highly specific nature.

Clone the repository locally, then run [uv sync](https://docs.astral.sh/uv/)
from within the local copy:

```sh
$ git clone https://github.com/FlyBase/fuzzy-onions.git
$ cd fuzzy-onions
$ uv sync
```

The tool can then be invoked through UV as follows:

```sh
$ uv run fzo --help
```

or, if you are not within the `fuzzy-onions` directory:

```sh
$ uv --project /path/to/fuzzy-onions run fzo --help
```

I recommend setting up an alias like:

```sh
alias fzo="uv --project /path/to/fuzzy-onions run fzo"
```

so that you can simply call `fzo` from anywhere:

```sh
$ fzo --help
```


Configuration
-------------
Fuzzy-Onions expects a configuration file under the name
`$XDG_CONFIG_HOME/fuzzyonions/config` (on GNU/Linux) or
`~/Library/Application Support/fuzzyonions/config` (on Mac OS).

Here is a minimal sample configuration file:

```
[store]
production: <path to a directory containing the datasets>
staging: <same, but for datasets from the staging server>

[curation]
proformae: <path to a directory containing FlyBase's proformae files>
```

See the [example configuration file](docs/example.conf) for more details
about the contents of that file.


Basic usage
-----------
Call `fzo -h` to get a list of the available subcommands, and
`fzo <subcommand> -h` to get the help message for a given subcommand.

The main commands are:

* `store`, to manage the local cache of datasets;

* `explorer`, to enter a mode allowing to interactively explore the
contents of a dataset;

* `curate`, to produce a proforma or a summarised expression table from
a dataset.

* `discover`, to automatically discover newly published scRNAseq papers.

* `tracker`, to keep track of the datasets and their status.


Copying
-------
Fuzzy-Onions is distributed under the terms of the MIT license.
