% Curation of scRNAseq datasets

# Prerequisites

* A MacOS computer installed according to the standard procedure in
  place in Cambridge. In particular:

  * a local checkout of the repository containing the proforma templates
    (hereafter assumed to be `~/SVN_folders/proformae`;

  * the ability to run Peeves;

  * write access to the staging folder on flybase-vm (hereafter assumed
    to be in `~/export-vm/staging/records`).

* Write access to the Harvard FTP server.


# Overview of the procedure

1. A paper reporting a scRNAseq experiment is published.

2. The paper is curated at FlyBase per the established procedures.

3. The scRNAseq curators check whether new terms will be needed in FBbt
   or FBcv to annotate the scRNAseq data, and if so request them from
   the ontology editors.

4. The scRNAseq data are processed by EBI curators and eventually made
   available on the Single-Cell Expression Atlas (SCEA).

5. The scRNAseq curators fetch and explore the data from the SCEA.

6. The scRNAseq curators fix the SCEA annotations if needed and send
   corrections back to the EBI curators.

7. The scRNAseq curators extract the data relevant for FlyBase from the
   raw data.

8. The scRNAseq curators produce a proforma record for the dataset, to
   be loaded using the standard record loading procedure.

9. The scRNAseq curators produce a “summarised expression table”, to be
   deposited on the Harvard FTP server.

10. Harvard developers load the summarised expression table into the
    database.

This document describes the role of the scRNAseq curators in steps 3 and
5–9.

Step 3 is isolated because it only requires the reading of the paper.
Steps 5–9 additionally require the actual scRNAseq data and therefore
can only be done once said data are made available by the EBI curators
in step 4.


# Adding new ontology terms

There are 3 cases where new ontology terms may be needed for curating a
scRNAseq paper/dataset.


## New methods

The paper is using a method for which we don’t already have a CV term in
FBcv. In this case, the scRNAseq curators should ask for the missing
term to be added at the appropriate place under the `dataset descriptor`
(FBcv:0003021) hierarchy.

We already have a selection of terms covering a large range of routinely
used methods, so the need for new terms here should arise infrequently.
We can expect, however, to need new terms for the precise scRNAseq
protocols. We currently (as of December 2021) have terms for three of
them (all under `single-cell library sequencing`, FBcv:0009006):
`DROP-Seq` (FBcv:0009006), `Smart-seq2` (FBcv:0009007), and `10x
sequencing` (FBcv:0009008). They correspond to protocols that have
already been used in fly scRNAseq papers, but there are several other
protocols that have been used in other species (at least 10 distinct
protocols have already been described) and we may need to add terms for
them if they start being used in fly studies as well.

If a method or experimental condition cannot be described by any
existing term in FBcv, the curators may have a look at what has been
done in some foreign ontologies to decide of new terms to propose to the
ontology editors. The following ontologies in particular may be of use:

* [Experimental condition ontology (XCO)](https://www.ebi.ac.uk/ols/ontologies/xco),

* [Evidence & Conclusion Ontology (ECO)](https://www.ebi.ac.uk/ols/ontologies/eco),

* [Zebrafish Experimental Conditions Ontology (ZECO)](https://www.ebi.ac.uk/ols/ontologies/zeco),

* [Experimental Factor Ontology (EFO)](https://www.ebi.ac.uk/ols/ontologies/efo).


## Missing terms for existing anatomical structures or cell types

This is the case where a scRNAseq paper refers to an anatomical
structure or cell that is already known to exist (it’s not a discovery
from the paper), but that for one reason or another is currently
insufficiently described (or possibly not described at all) in FBbt.

A recent example is a paper reporting a scRNAseq study on hemocytes
isolated from adult flies, which has required the addition of terms
representing adult hemocytes because at the time FBbt only had terms for
the embryonic/larval hemocytes.

This also encompasses the addition of “grouping terms” that might be
needed to collectively refer to a group of structures or cell types. An
example of that case was the addition of `endocrine cell of the ring
gland`, to regroup all the various types of endocrine cells found in the
ring gland.


## Terms for newly described cell types

A scRNAseq paper may contain a claim that the authors have discovered a
new cell type, which may then warrant the introduction of a new term for
that cell type.

However, the claim must be carefully examined before a corresponding
term is added. This is ultimately a judgement call from the scRNAseq
curators, but the following principles should be observed when deciding
whether a new term should be added:

* A cluster is not necessarily a cell type. Do _not_ create a new term for
  each cluster in a scRNAseq study.

* Be mindful of the distinction between a cell state and a cell type.
  Only cell types should be described in FBbt. For example, if a study
  reports a cluster of “latent hemocytes” and a cluster of “activated
  hemocytes”, both clusters should be annotated as “hemocytes” and no
  new term is needed.

* Only explicit claims should be considered. If a potentially new cell
  type only appears in the legend of a t-SNE plot without ever being
  mentioned in the text, that is not a claim. The authors have to
  explicitly state in the main text that they have found what they think
  is a new cell type.

* The claimed new cell type should be characterised by more than a mere
  transcriptional profile. Basically, there should be something to be
  said about that cell type, beyond the fact that the cells have a high
  level of expression of a given gene (or of a combination of genes).

* The claim ideally be backed by more evidence than just the scRNAseq
  experiment (for example an IF staining confirming the presence of
  cells expressing a given marker in tissues).

When in doubt, be conservative.


# Data from the EBI

This section describes what data are available from the EBI’s
Single-Cell Expression Atlas and what is needed by the scRNAseq
curators.

The main page for the EBI’s Single-Cell Expression Atlas is
<https://www.ebi.ac.uk/gxa/sc/experiments>. It lists all the currently
available datasets. You may filter by species so that you only get the
Drosophila melanogaster datasets:
<https://www.ebi.ac.uk/gxa/sc/experiments?species=%22drosophila%20melanogaster%22>.

The SCEA’s staging server is at
<https://wwwdev.ebi.ac.uk/gxa/sc/experiments>. It has exactly the same
structure and the same kind of contents as the production server, but
contains additional datasets that have not been published in production
yet. The FlyBase scRNAseq curators can use the staging server to know
which datasets are next in queue to reach the production server and to
start working on them before they are officially released.


## Informations from the SCEA website

Most of the data needed for curating the scRNAseq datasets are to be
found in downloadable files (described in the next section), but there
are a few pieces of information that need to be obtained from the web
site itself:

1. The **Dataset ID** for the dataset to be curated. This is available
   in the URL pointing to the dataset. For example, in
   <https://www.ebi.ac.uk/gxa/sc/experiments/E-MTAB-8698/>, the Dataset
   ID is `E-MTAB-8698`. This is also displayed at the bottom of the
   “Supplementary Information” page for the dataset.

2. The corresponding paper, whose reference is given as free text on the
   main page for the dataset. The scRNAseq curators should match that
   reference with a FlyBase FBrf.

3. The reference genome assembly used by the EBI curators to align the
   sequencing reads. This is available on the “Supplementary
   Information” page, under the name “Reference”.


## Data files

All data files provided by the EBI curators are available in the
“Downloads” section. The FlyBase scRNAseq curators will need the
Experiment Design table, the Raw Counts table, and the Normalised Counts
table.


### Experiment Design table

This table contains the cell-level annotations for the dataset. Each row
in that table represents a single cell.

The columns usually have relatively self-explanatory names. They may
vary slightly from one dataset to the next, but the following “core set”
of annotations is almost always present (though they may be empty):

* `Assay` – this is the unique identifier for the cell;

* `Sample Characteristics[organism]` – should obviously always be
  `Drosophila melanogaster` for the datasets that are of interest for
  FlyBase;

* `Sample Characteristics[strain]` – this one is sometimes misused by
  the EBI curators, who may use it to record the genotype of the sample
  (e.g. `w[1118]`) instead of the actual strain (e.g. `Oregon-R`);

* `Sample Characteristics[genotype]` – can either be an actual genotype,
  or a free-text description such as `wild type`;

* `Sample Characteristics[sex]`;

* `Sample Characteristics[developmental stage]`;

* `Sample Characteristics[age]` – may be used in combination with the
  previous field to give more precision;

* `Sample Characteristics[organism part]` – used to indicate which organ
  or tissue a cell has been dissected from.

Depending on the dataset and the actual experiments performed,
additional columns may be present to record further experimental
conditions. This may include, e.g.

* `Sample Characteristics[cell type]`, if cells of a specific type were
  selected by FACS or a similar method prior to being subjected to
  scRNAseq;

* `Sample Characteristics[stimulus]`, if the flies/tissues/cells were
  subjected to some kind of treatment as part of the experiments.

Inferred cell type annotations, if present, will typically be
represented using two distinct columns:

* `Factor Value[inferred cell type - authors labels]`, which contains
  the inferred cell types as provided by the authors;

* `Factor Value[inferred cell type - ontology labels]`, which represents
  the EBI curators’ attempt to match the original labels with known cell
  types.

Each column is normally accompanied by a sister column that is intended
to contain an identifier of some sort (typically the identifier for an
ontology term, though identifiers for other kinds of records may also be
used). For example, the `Sample Characteristics[organism]` column is
paired with a `Sample Characteristics Ontology Term[organism]` column
that will contain the value
`http://purl.obolibrary.org/obo/NCBITaxon_7227` (the PURL identifier for
Drosophila melanogaster in the NCBI Taxonomy Database).

The scRNAseq curators should explore thoroughly the contents of the
Experiment Design table, make note of the available annotations, and
check whether the contents match what they were expecting based on what
was described in the paper. They may decide to exclude the dataset from
curation if they can’t reconcile the Experiment Design table with the
paper. A real example of that involved a dataset containing only the
data from one experimental condition while the corresponding paper
described three distinct experimental conditions.

Another reason to temporarily exclude a dataset is if it does not
contain any inferred cell type annotations. This can happen either
because the authors did not even try to infer cell types from the
transcriptomic profiles (e.g. this was out of scope for their study), or
because they have not transmitted their inferred cell types to the EBI
curators.  When working on a dataset without cell type annotations, the
scRNAseq curators should check in the corresponding paper if cell types
have been inferred by the authors; if so, they should try requesting the
inferred cell types from them, on behalf of both FlyBase and the EBI,
and keep the dataset “on hold” until the cell type annotations have been
obtained. If cell types have not been inferred by the authors, or if it
has not been possible to obtain the annotations, then the dataset may be
curated as normal.

Exploring the Experiment Design table may be done with any tool capable
of manipulating/displaying TSV files. The [Fuzzy-onions program][Using
Fuzzy-onions], has some features to that effect.


### Counts tables

The Counts tables contain the number of reads for each gene and each
individual cell. They are basically huge matrixes where each column is a
single cell ID (which should match a row in the Experiment Design table)
and each row is a single gene ID (which should match a known FBgn entity
in FlyBase).

The Raw Counts table is to be used to extract the number of reads for
each sample. This is not expected to be much useful but we still record
that bit of data.

The Normalised Counts table is to be used to produce a “summarised
expression table” that associates cell types (instead of individual
cells) with gene expression levels.

The counts tables are in the MatrixMarket format and require some
specialised software to be manipulated (such as SciPy, in Python).
However the use of the aforementioned [Fuzzy-onions program][Using
Fuzzy-onions] is strongly recommended, as it can automatically perform
all the steps involving the counts files (indeed it has initially been
expressly designed for that purpose).


# Amending the SCEA annotations

During the curation process, there are two types of corrections that may
need to be applied to the Experiment Design table provided by the EBI
curators:

* adding FlyBase identifiers or CV terms wherever possible (e.g. for
  strains, developmental stages, etc.), when it has not already been
  done by the EBI curators;

* correcting the inferred cell types’ ontology labels, if those selected
  by the EBI curators are not appropriate.


## Adding FlyBase identifiers

This is done only for the benefit of SCEA users, as part of the
cooperation agreement between FlyBase and the EBI on the scRNAseq
project, and will have no consequence on what FlyBase does with the
scRNAseq data. The rationale is that when someone uses the SCEA website
to view a dataset, it could be useful for them to have direct links to
FlyBase report pages for the strains, transgenic markers, stocks, etc.
that have been used.

For example, if cells in the experiment design table are annotated as
coming from the Oregon-R strain, the scRNAseq curators may insert the
link to the corresponding report page
(http://flybase.org/reports/FBsn0000276) in a new column next to the
original `Sample Characteristics[strain]` column.

The scRNAseq curators should use the data from thin curation to get the
appropriate identifiers. If they find that a genetic object referenced
in the annotations is not represented in the database, they may ask the
thin curators to have another look at the associated paper and consider
adding the missing object.

Note that there is no need to add FBbt identifiers for cell types. Such
identifiers are automatically added by a specific process (Zooma) in the
EBI pipeline. For cell types, all the scRNAseq curators need to do is to
ensure that each free-form cell type name corresponds exactly to a FBbt
term label.

The experiment design table amended with such additional identifiers
should be sent back to the EBI curators (Nancy) by e-mail.


## Correcting cell types

There are several reasons why the cell types provided by the EBI
curators could need corrections.

* Plain error on the part of the EBI curators, who picked up an
  incorrect term from FBbt. The EBI curators are not drosophilists and
  may not be familiar with FBbt, so it can happen. This includes the
  case where they have chosen a term that is strictly speaking correct
  but without realising that a more precise term was also available.

* Use of an imprecise term because the appropriate more precise term has
  only been added recently to FBbt. The EBI curators only work from a
  released version of FBbt, so any new term added since the last release
  (e.g., during step 3 of the current procedure) will be temporarily
  invisible to them.

* Adult vs embryonic/larval distinction. In FBbt, we almost always have
  a triad of terms for any single cell type: one referring to the cell
  type independently of the developmental stage, one referring to the
  adult cell type, and one referring to the embryonic/larval cell type.
  In FlyBase, to annotate scRNAseq clusters, we always want to use the
  most precise term: if a dataset has been obtained from an adult
  sample, the corresponding clusters should be annotated using the terms
  for the adult cell types, if they exist.

A complication regarding the last point is that the SCEA is not
interested in the adult/larval distinction. They want to annotate cell
types using only the generic, stage-neutral terms. This means that a
dataset may sometimes require two distinct set of corrections: one set
intended for the SCEA (which should not include any correction for
adult-vs-larval cases), and one intended for FlyBase (which should
always include such corrections if needed).

Cell type corrections intended for the SCEA should be sent to the EBI
curators under the form of a separate TSV file containing one column for
the original term and one column for the suggested correction, with
optionally a third column containing a comment for the benefit of the
EBI curators.


# Extracting data for FlyBase

Before preparing the proforma to load into FlyBase, the following pieces
of data need to be extracted from [the files previously
described][Data files]:

* the number of individual cells in each sample;

* if cell types have been identified, the number of individual cells in
  each cluster of a given cell type;

* the number of sequencing reads in each sample.

This step is performed automatically by the [Fuzzy-onions program][Using
Fuzzy-onions] and this section is thus for information only.

As [described earlier][Experiment Design table], each row in the
Experiment Design table represents an individual cell, so extracting the
number of cells in a given sample can be done by filtering the table on
a discriminating value for that sample and counting the number of rows
after filtering.  For example, if a dataset contains two samples that
differ by a treatment applied to the larvae prior to dissection, with
one set of larvae being left untreated and the other one being infested
by a parasite, this might be represented by the presence of two
different values in the `Sample Characteristics[stimulus]` column:
`control` for the untreated sample, and `infeste`d for the infested
sample; in that case, the number of cells in the untreated sample will
simply be the number of rows with the value `control` in the `Sample
Characteristics[stimulus]` column.

Similarly, to get the number of cells in a cluster of a given cell type,
the Experiment Design table needs to be filtered to keep only the rows
containing the desired cell type in the `Factor Value[inferred cell type
- ontology labels]` column — after applying [corrections][Correcting
  cell types] if needed. Of course, if a dataset is divided in several
samples, the table must also be filtered to only get the cells from a
given sample. Following the example above, the number of cells in the
cluster of “lamellocytes” from the “infested” sample will be the number
of rows that have the value `infested` in the `Sample
Characteristics[stimulus]` column and the value `lamellocyte` in the
`Factor Value[inferred cell type - ontology labels]` column.

Extracting the number of reads from a given sample involves manipulating
both the Experiment Design table and the Raw Counts table. The procedure
is:

* Filter the Experiment Design table to get only the rows corresponding
  to the desired sample.

* Extract the corresponding cell identifiers (in the `Assay` column).

* Filter the Raw Counts table to only keep the columns corresponding to
  those cell identifiers.

* Sum all values in the resulting filtered raw count table. The sum is
  the number of reads across all genes and all cells in the sample.


# Producing the proforma

The dataset is to be described in a proforma assembled as follows:

* a publication proforma for the reference associated with the dataset;

* a **project**-type dataset proforma (that is, a dataset proforma with
  the `LC2a` field set to `project`) representing the overall dataset;

* for each sample in the dataset:

  * a **biosample**-type dataset proforma representing that particular
    sample;

  * a **assay**-type dataset proforma representing the scRNAseq
    experiment performed on that sample;

  * a **result**-type dataset proforma representing the subsequent
    analysis of the scRNAseq data;

  * for each identified cell cluster (if cell type annotations are available):

    * a **result**-type dataset proforma representing that cell cluster.

Because it is designed to be flexible enough to suit the needs of
different kinds of datasets, the dataset proforma may look complicated.
The following subsections will describe how to use the fields of that
proforma to describe scRNAseq datasets.

The final proforma should have a `.coll` extension, and be loaded by the
standard procedure used for all other proformae.

The [Fuzzy-onions program][Using Fuzzy-onions] will automatically assemble the
proforma and fill all the fields that can be automatically filled.


## Project-level proforma

### LC1f (Database ID for dataset)

Set to `new` for a newly curated dataset, or to an existing dataset ID
to amend an existing record.


### LC1a (Symbol) {#project-LC1a}

This is a curator-chosen symbol to uniquely refer to that dataset in the
entire database. The current convention is to use the following format:
`scRNAseq_YYYY_Authorname`, where `YYYY` is the publication year of the
associated reference and `Authorname` is the last name of the
reference’s first author.

Note: This convention makes colliding symbols possible, if two datasets are
published in the same year by the same author (or by two authors with the same
last name). The probability of such collisions is considered to be very low for
now, and it is deemed enough to rely on Peeves to avoid them – Peeves will warn
curators if the symbol they’re trying to use when creating a new dataset record
has already been used for another dataset, in which case the curators are free
to alter the symbol as they see fit to avoid the collision (e.g. by adding a
letter after the year, as in `scRNAseq_YYYYb_Authorname`). The convention may
be updated in the future if/when the number of published scRNAseq datasets
becomes high enough to make colliding symbols a recurring annoyance.


### LC6g (Dataset title)

A free text title for the dataset, should be less than 140 characters.
Current recommendation is something along the lines of `Single-cell
RNA-seq study of <tissue and stage information> [upon condition]`. Do
not mention Drosophila melanogaster (keep that for the `LC6a` field). Do
not include a final dot.

Examples from curated datasets:

* Single-cell RNA-seq study of stage 6 embryos

* Single-cell RNA-seq study of larval hemocytes upon wasp infestation

* Single-cell RNA-seq study of pupal olfactory projection neurons

* Single-cell RNA-seq study of larval lymph gland hemocytes over time
  and upon wasp infestation

* Single-cell RNA-seq study of third instar larval eye discs and the
  role of Rbf


### LC6a (Description)

A free text simple description of the dataset. This should give enough
information about what the dataset is about while at the same time being
kept relatively consistent across datasets. Current recommendation is to
always start with either `A characterization of the diverse populations
of cells in…` (if cell types have been identified) or` A single-cell
transcriptomic study of the diverse populations from…` (if there are no
cell types), and to complement with dataset-specific informations.
However this is a rough guideline only and deviating from it is fine if
it makes sense for a particular dataset.

Examples:

* A characterization of the diverse populations of hemocytes in the
  hemolymph of third instar Drosophila larvae infested with Leptopilina
  boulardi.

* A single-cell transcriptomic study of the diverse populations of
  olfactory projection neurons in the Drosophila pupal brain.

* A single-cell transcriptomic study of the diverse populations of cells
  making up the developing eye disc during the third instar larval
  stage, in wild-type and Rbf mutant context.

* A study using single-cell mRNA sequencing on cells from precisely
  staged embryos to produce a single-cell resolution map of the
  Drosophila embryo.


### LC2a (Type of dataset entity)

Should always be set to `project ; FBcv:0003023`.


### LC2b (Type of dataset data)

Most often will be set to `transcriptome ; FBcv:0003034`, though
depending on the dataset another CV term from the `project type`
(FBcv:0003029) hierarchy might be more suitable.


### LC13c (Key GO term(s) - Biological Process)

This is interpreted as the biological processes that are investigated in
the study.

Terms used so far in this field include:

* `embryonic morphogenesis ; GO:0048598`

* `brain development ; GO:0007420`

* `larval lymph gland hemocyte differentiation ; GO:0035168`

* `eye-antennal disc development ; GO:0035214`

* `wing disc pattern formation ; GO:0035222`

* `defense response to other organism ; GO:0098542`

Note: The other GO term fields (`LC13a` for Cellular Component, `LC13b`
for Molecular Function) are presumably unlikely to ever be needed given
the nature of scRNAseq data. However it is of course fine to use them if
it happens to be relevant for a particular dataset.

### LC4a (Species of derivation)

Should always be `Dmel`.


### LC4i (Other species of interest)

This field refers to any other species somehow involved in the
experiment. The value should come from the `organism` table of the
database.

So far, this has been used for scRNAseq experiments performed on larvae
infested with a parasite (value `Lbou` for Leptopilina boulardi).


### LC6d (Dataset/collection members stored in database)

At the project level, this should always say `N`.

Leave the subsequent fields `LC6e` and `LC6f` empty.


### LC11m (Experimental protocol – dataset)

At the project level, this field should be set to term(s) from the
`study design` (FBcv:0003130) hierarchy.

Examples:

* `biotic stimulus study ; FBcv:0003134`

* `developmental stage study ; FBcv:0003127`

* `tissue type study ; FBcv:0003130`


### LC11a/LC6b/LC11c/LC11e (Experimental protocol)

Those fields are intended to describe how the single cells were obtained
(`LC11a`), how the mRNA libraries were prepared from the cells (`LC6b`),
how the libraries were sequenced (`LC11c`), and how the resulting data
were analyzed (`LC11e`).

In principle, each of these fields should be used at a specific level in
the dataset model:

* `LC11a` at the biosample level;

* `LC6b` and `LC11c` at the assay level;

* `LC11e` at the result level.

However, it is possible to use any of them at the project level if the
same protocol has been used across all biosamples in a given dataset.
Their values will then be automatically propagated to all biosamples,
assays, and results in the dataset, avoiding the need for repetition.

Typically, it is expected that the protocols for preparing and
sequencing the mRNA libraries and for analyzing the data will be the
same across biosamples (since using different protocols to process
different samples would likely be a questionable study design), so the
fields `LC6b`, `LC11c`, and `LC11e` will most often be filled once at
the project level.

How to use these fields exactly will be described in the sections
covering the biosample-, assay-, and result-level proformae.


### LC7a (Types of additional data available)

This is a free text description of what data are available outside of
FlyBase.

Currently, we only describe what is available on the EBI’s Single Cell
Expression Atlas, using the following boilerplate: “The EMBL-EBI’s
Single Cell Expression Atlas provides cell-level annotations, clustering
data, raw and normalised read counts, and putative marker genes.”


### LC99a / LC99b (DataSet Accession Number / Database name)

These fields collectively point to an external database. `LC99b` is the
name of the external database (must correspond to a database name in
Chado’s `db` table), and `LC99a` is the accession number for the current
dataset in that database.

For now, we are only referencing the EBI’s Single Cell Expression Atlas,
so there should be only one instance of each of these fields. `LC99b`
should be set to `EMBL-EBI Single Cell Expression Atlas Datasets`, and
`LC99a` should be set to the SCEA’s dataset ID (obtained from the SCEA
website [as described here][Informations from the SCEA website].


### LC8c (Created by)

This is a Markdown-style link to the website of the research group that
produced the dataset (if the group has such a website), as in:

    [Name of the group that produced the dataset](URL of their lab).


## Biosample-level proforma

### LC1f (Database ID for dataset)

Set to `new` for a newly curated dataset, or to the existing dataset ID
to amend an existing record.


### LC1a (Symbol) {#biosample-LC1a}

This is a curator-chosen symbol to uniquely identify the biosample.
Current convention is to re-use the [symbol chosen at the project
level](#project-LC1a) and to append to it a short code to distinguish
between the different samples in the same dataset.

For example, in a dataset that contains one sample from non-infested
larvae and one sample from wasp-infested larvae, the corresponding
symbols could be `scRNAseq_YYYY_Authorname_NI` for the non-infested
sample and `scRNAseq_YYYY_Authorname_WI` for the wasp-infested sample.


### LC6g (Dataset title) {#biosample-LC6g}

A short title for the sample. Current recommendation is to follow the
template `Cells <from stage and tissue> [upon condition]`. A more
precise term than “Cells” should be used if possible, if it is known
that the sample is made of cells of a more precise type.

Examples:

* Cells from stage 6 embryos

* Circulating hemocytes from wasp-infested third instar larvae at 24h
  post-infestation

* GH146-labelled olfactory projection neurons from pupal brains

* Lymph gland cells from non-infested larvae at 72h after egg-laying

* Cells from wild-type eye discs


### LC2a (Type of dataset entity)

Should always be `biosample ; FBcv:0003024`.


### LC2b (Type of dataset data)

Should be a CV term from the `biosample type` (FBcv:0003044) hierarchy.
Will most often be `isolated cells ; FBcv:0003047`, or `isolated nuclei
; FBcv:0009004` for single-nucleus RNA-seq.


### LC3 (Dataset “belongs_to” this project)

Must be set to the same value as the [LC1a field at the project
level](#project-LC1a); this is what establishes the link between the
project-level object and the biosample-level object.


### LC4a (Species of derivation)

Should always be set to `Dmel`.


### LC4i (Other species of interest)

As at the project level, this should be the 4-letters abbreviation for
any non-Dmel species that is involved in the experiment, if any.


### LC4h (Strain used)

The FlyBase symbol for the fly strain used in this sample, if known.
Example: `Oregon-R`.


### LC4b (Strain used)

To be used instead of `LC4h` when the strain used (as reported in the
paper) is not present in the database and therefore has no associated
symbol. But check with thin curators whether the unknown strain
could/should be added to the database (which would then allow to use
`LC4h`) before using this field.


### LC4f (Genotype used)

The full genotype of the flies or larvae from which the sample was
obtained, if known. The genotype should be written as a single-line
using “standard” genotype notation, e.g. `y[1] v[1]; P{CaryP}attP2`.

Note: If parts of the genotype are significant to understand what the
sample is (for example, a mutant allele in a study that compares
scRNAseq profiles between wild-type and mutant flies, or a GAL4
transgene driving the expression of GFP that is used to collect cells of
a certain type), these parts must also be listed in `LC12a`.


### LC4g (Stage and tissue)

A [TAP](TAP.notes)-like statement describing the developmental stage and
tissue from which the cells of the sample come from. Should be of the
form

    <e><t>FBdv term<a>FBbt term<s><note>free text comment

The `<e>` field is not used for dataset objects and should be left
empty. The `<s>` field, intended for “cellular component” GO terms, is
not expected to ever be needed given the nature of scRNAseq data.

If the sample comes specifically from individuals of a given sex, this
should be specified as a `male` or `female` qualifier in the `<t>`
field.

If the individuals were collected at a specific time (e.g. “96h after
egg-laying”), this may be specified in the free text comment.

Examples:

* `<e><t>embryonic stage 6<a>embryo<s><note>`

* `<e><t>third instar larval stage<a>embryonic/larval hemolymph<s><note>`

* `<e><t>adult stage | female<a>ovary<s><note>`

* `<e><t>third instar larval stage<a>wing disc<s><note>`


### LC12a / LC12b (Experimental entity / Type of experimental entity)

Those two fields are used to list any FlyBase feature (gene, allele,
insertion, etc.) that can be considered an experimental factor in the
experiment from which the dataset is obtained. `LC12a` is the FlyBase
identifier for the feature and `LC12b` is a value taken from the
`experimental_design` CV. Repeat those two fields as many times as
necessary.

For transgenes, the identifier used here should be the one for the
corresponding allele object (`FBal...`), as identified by the thin
curators.

Of the allowed values for `LC12b`, the ones typically expected to be
used are `allele_used` and `transgene_used`, and the generic value
`experimental_design` for any non-gene entity.

Note that the LC12b field does not use the standard syntax for CV term
fields (`term name ; CV term ID`).


### LC6d (Dataset/collection members stored in database)

At the biosample level this should always say `N`.

Leave the subsequent fields `LC6e` and `LC6f` empty.


### 7.2.15 LC11m (Experimental protocol – dataset)

A set of CV terms describing the sample (or the method(s) used to
prepare the sample). The terms should be picked in the `biosample
attribute` (FBcv:0003137) hierarchy.

Examples:

* `tissue dissection ; FBcv:0003169`

* `cell isolation ; FBcv:0003170` (be careful that this term is only
  intended for isolation by means of affinity purification or FACS;
  physical dissociation of cells from a dissected tissue or organ does
  not count as `cell isolation`, even if the result is a suspension of
  isolated cells)

* `biotic treatment ; FBcv:0003166`

* `multi-individual sample ; FBcv:0003141` or its opposite `single
  individual sample ; FBcv0003140` (presumably one of those two terms
  will always be used)


### LC11a (Experimental protocol, source isolation and prep)

This field is intended to describe how the cells were prepared prior to
RNA extraction.

As for the other “protocol” fields (`LC6b`, `LC11c`, `LC11e`), this can
be described once at the project level if all samples in a dataset have
been obtained in exactly the same way. However, most of the time the
difference between the samples will precisely reside in the way they
have been obtained, whereas all the subsequent steps (RNA extraction,
sequencing, and analysis) will be the same across samples.

This field should only be used to capture critical informations about
the experimental protocol that are not captured elsewhere in the
biosample-level proforma. In particular, should *not* be recorded in
this field:

* the strain or genotype (captured in `LC4h`, `LC4b`, or `LC4f`, as
  appropriate; also in `LC12a` for significant alleles or transgenes);

* the stage at which the individuals were collected (captured in the
  `<t>` subfield of `LC4g`);

* the organism parts that were dissected (captured in the `<a>` subfield
  of `LC4g`);

* any part of the protocol that is already described by a FBcv term in
  `LC11m`.

If something is left to be described after excluding the above, care
should be taken to never include low-level technical details, such as
the composition of buffers, the centrifugation parameters, or the
providers’ names and catalogue numbers for any single reagent (reagents
themselves are most likely to be not critical enough to be described).


## Assay-level proforma

### LC1f (Database ID for dataset)

Set to `new` for a newly curated dataset, or to the existing dataset ID
to amend an existing record.


### LC1a (Symbol) {#assay-LC1a}

A curator-chosen symbol to uniquely identify the assay. The current
convention is to re-use the [symbol used for the corresponding
biosample](#biosample-LC1a) and to append to it the string `_seq`.
Following the examples of biosamples given above, the symbols for the
corresponding assays would be `scRNAseq_YYYY_Authorname_NI_seq` and
`scRNAseq_YYYY_Authorname_WI_seq`.


### LC6g (Dataset title)

The recommended format is `Single-cell RNA-seq of <title of
corresponding biosample>`.

Examples:

* Single-cell RNA-seq of cells from stage 6 embryos

* Single-cell RNA-seq of circulating hemocytes from wasp-infested third
  instar larvae at 24h post-infestation

* Single-cell RNA-seq of cells from wild-type eye discs


### LC2a (Type of dataset entity)

Should always be set to `assay ; FBcv:0003025`.


### LC2b (Type of dataset data)

Must be a CV term from the `assay type` (FBcv:0003050) hierarchy. The
two expected possible values are `single-cell RNA-Seq ; FBcv:0009000` or
`single-nucleus RNA-Seq ; FBcv:0009001`.


### LC3 (Dataset “belongs_to” this project)

Must be set to the same value as the [LC1a field in the project-level
proforma](#project-LC1a).


### LC14a (Dataset assay/collection is “assay_of” this biosample)

Must be set to the same value as the [LC1a field of the corresponding
biosample-level proforma](#biosample-LC1a).


### LC14d (Dataset “technical_reference_is”)

If the overall dataset contains an assay that serves as a technical
reference for this one, this field should be set to the same value as
the `LC1a` field of the technical reference assay object.


### LC14e (Dataset “biological_reference_is”)

If the overall dataset contains an assay that serves as a biological
reference for this one (e.g., an assay from an untreated control or from
a wild-type genotype), this field should be set to the same value as the
`LC1a` field of the biological reference assay object.


### LC4e (Species of derivation)

Should always be set to `Dmel`. The value of this field is not
propagated from the biosample-level object, so it should be explicitly
repeated at the assay level.


### LC4i (Other species of interest)

Should be set to the 4-letters abbreviation of any other species somehow
involved in the experiment. Like `LC4e`, this field is not propagated
from the biosample-level object, so its value (if it has one) should be
explicitly repeated at the assay level.


### LC6d (Dataset/collection members stored in database)

At the assay level this should always say `N`.


### LC6e (Number of entities in dataset/collection)

This field should be set to the number of sequencing reads in the
corresponding sample, extracted from the raw data [as described
here][Extracting data for FlyBase].


### LC6f (Comment on number of entities in dataset/collection)

This is a description of what the value stored in `LC6e` means. On the
dataset report, this will be displayed immediately next to the value of
`LC6e`, so this should be set to something like `reads in this sample`.


### LC11m (Experimental protocol – dataset)

A set of CV terms describing the assay. The terms should be picked up
from the `assay method` (FBcv:0003208) hierarchy. Expected values
include:

* `10x sequencing ; FBcv:0009008`

* `DROP-Seq ; FBcv:0009006`

* `Smart-seq2 ; FBcv:0009007`

* or the generic term `single-cell library sequencing ; FBcv:0009005` if
  another scRNAseq protocol was used (however in that case, consider
[requesting a FBcv term][Adding new ontology terms] representing that
method).

In addition to terms representing the scRNAseq protocol, terms
representing the sequencing technology may also be used, such as `454
sequencing ; FBcv:0003212` or `Illumina sequencing ; FBcv:0003213`.


### LC6b (Experimental protocol, dataset/collection preparation)

This field is intended to describe how the biological material obtained
from the sample was processed to obtain the mRNA library used for
downstream sequencing.

This field is better left empty, as the FBcv terms in the assay-level
`LC11m` field should normally be sufficient to describe the protocol
used for preparing the mRNA library. Only use that field if a critical
information about the protocol cannot be covered by a FBcv term (and in
that case, consider requesting a new FBcv term that would convey such an
information).

If it is used at all, this field may be used at the project-level, if
the same library preparation method has been used for all samples in a
given dataset.


### LC11c (Experimental protocol, mode of assay)

This field is intended to describe how the mRNA libraries (prepared as
described in `LC6b`) were sequenced. As for `LC6b`, this can be
described at the project level rather than at the assay level if the
sequencing method is the same across all samples.

This field is better left empty, as the FBcv terms in the assay-level
`LC11m` field should normally be sufficient to describe the sequencing
protocol. Only use that field if a critical information about the
protocol cannot be covered by a FBcv term (and in that case, consider
requesting a new FBcv term that would convey such an information).


## Result-level proforma

### LC1f (Database ID for dataset)

Set to `new` for a newly curated dataset, or to the existing dataset ID
to amend an existing record.


### LC1a (Symbol) {#result-LC1a}

A curator-chosen symbol to uniquely identify the result object. Current
convention is to reuse the [symbol of the corresponding
assay](#assay-LC1a) and to append to it the string `_clustering`, as in
`scRNAseq_YYYY_Authorname_NI_seq_clustering`.


### LC6g (Dataset title)

The recommended format is `Clustering analysis of <title of corresponding biosample>`.


### LC2a (Type of dataset entity)

Should always be set to `result ; FBcv:0003026`.


### LC2b (Type of dataset data)

Should always be set to `cell clustering analysis ; FBcv:0009002`.


### LC3 (Dataset “belongs_to” this project)

Must be set to the same value as the [LC1a field in the project-level
proforma](#project-LC1a).


### LC14b (Dataset result is “analysis_of”)

Must be set to the same value as the [LC1a field of the corresponding
assay-level dataset object](#assay-LC1a).


### LC14h (FlyBase reference gene model annotation set)

Should be set to a version of the FlyBase gene model annotation set, as
in `Dmel R6.28`. The version that has been used for a particular dataset
should be looked for [on the SCEA website][Informations from the SCEA
website].


### LC4a (Species of derivation)

Should always be set to `Dmel`. The value of this field is not
propagated from the biosample-level or assay-level object, so it should
be explicitly repeated at the result level.


### LC4i (Other species of interest)

Should be set to the 4-letters abbreviation of any other species somehow
involved in the experiment. Like `LC4e`, this field is not propagated
from the biosample-level or assay-level object, so its value (if it has
one) should be explicitly repeated at the result level.


### LC6d (Dataset/collection members stored in database)

At the result level this should always say `N`.


### LC6e (Number of entities in dataset/collection)

This field should be set to the number of isolated cells in the
corresponding sample, extracted from the raw data [as described
here][Extracting data for FlyBase].


### LC6f (Comment on number of entities in dataset/collection)

This is a description of what the value stored in `LC6e` means. On the
dataset report, this will be displayed immediately next to the value of
`LC6e`, so this should be set to something like `cells in this sample`.


### LC11e (Experimental protocol, data analysis)

This field is intended to describe how the sequencing raw data were
analysed.

Similarly to the `LC6b` and `LC11c` fields, this can be described once
at the project level if the data for all samples were analysed in the
same way (I have seen at least one dataset where it was not the case,
though).

This part of the protocol is typically described with many details in
scRNAseq papers. Not all details should be included in the dataset
report.

What should be included (if originally provided):

* the analysis programs and published protocols used (e.g., Drop-seq
  tools, STAR, Seurat, etc.) with an accompanying reference or URL;

* the reference genome used for alignment;

* the quality control and filtering parameters (which criteria were used
  to exclude cells);

* how dimensionality reduction was performed (e.g., how many variable
  genes were considered for principal component analysis);

* how clustering was performed (e.g., how many principal components were
  retained).

Fine details of software parameters (e.g., “-outFilterScoreMinOverLread
0.4 -outFilterMatchNMinOverLread 0.4…”) or names of software functions
called during the analysis should not be included.

For improved consistency across dataset reports, it is recommended not
to simply re-use the original text from the paper’s Methods section, but
to paraphrase it around the important pieces of information listed
above.

Examples:

* Reads were processed using the Picard suite and the Drop-seq tools
  v1.13 following the Drop-seq Core Computational Protocol v1.2. They
  were mapped to the reference genome BDGP 6.02. The quality filtered
  datasets were merged using the common genes into a single Seurat
  (v3.1) object and integrated using Harmony (Korsunsky et al., 2019).
  Dimensionality reduction and clustering was performed by computing 20
  principal components using the top 2000 most variable genes and a
  clustering resolution of 0.4.

* Reads were aligned to the Drosophila melanogaster reference genome
  using STAR v2.4.2a (Dobin et al., 2013). Uniquely mapped reads that
  overlap with genes were counted using HTSeq-count v0.7.1 (Anders et
  al., 2015), with cells having fewer than 300,000 uniquely mapped reads
  removed. Principal component analysis was performed on the top 500
  most variable genes; typically 7 to 12 principal components were used
  for further analysis.

* Reads were processed using the Drop-seq tools v1.13. They were mapped
  to the Drosophila melanogaster reference genome BDGP 6.02. Further
  analysis was performed using Seurat v3.0. For quality filtering, cells
  with more than 5,000 genes or less than 200 or 400 genes were removed,
  as well as cells where mitochondrial genes represented more than 10%
  of the total UMI count. Dimensionality reduction was performed by
  computing 52 principal components and with a clustering resolution of
  0.8.


## Cell cluster-level proforma

### LC1f (Database ID for dataset)

Set to `new` for a newly curated dataset, or to the existing dataset ID
to amend an existing record.


### LC1a (Symbol) {#cluster-LC1a}

A curator-chosen symbol to uniquely identify the cell cluster object.
Current convention is to reuse the [symbol of the corresponding
clustering analysis object](#result-LC1a) and to append to it a suffix
representing the cell type, as in
`scRNAseq_YYYY_Authorname_NI_seq_clustering_prohemocytes`.

The suffix representing the cell type (`prohemocytes` in the example
above) should come from a FBbt term, with the following transformations:

* any space in the cell type name is replaced by an underscore;

* any stage qualifier (adult, larval, embryonic) is removed;

* an underscore is added at the beginning;

* a final ‘s’ is appended.

For example, the suffix for a cluster of crystal cells from a larval
sample (FBbt term: `embryonic/larval crystal cell`) would be
`_crystal_cells`.


### LC6g (Dataset title)

The recommended format is `Clustering analysis of <title of
corresponding biosample>, <cell type>s cluster`.

Examples:

* Clustering analysis of hemocytes from wasp-infested third instar
  larvae, crystal cells cluster

* Clustering analysis of circulating hemocytes from wasp-infested third
  instar larvae at 24 h post-infestation, plasmatocytes cluster


### LC2a (Type of dataset entity)

Should always be set to `result ; FBcv:0003026`.


### LC2b (Type of dataset data)

Should always be set to `transcriptional cell cluster ; FBcv:0009003`.


### LC3 (Dataset “belongs_to” this project)

Should be set to the same value as the [LC1a field of the object
representing the clustering analysis](#result-LC1a) to which this
cluster belongs (_not_ the value from the project-level proforma!).


### LC4a (Species of derivation)

Should always be set to `Dmel`. The value of this field is not
propagated from the higher-level objects, so it should be explicitly
repeated at the cluster level.


### LC4i (Other species of interest)

Should be set to the 4-letters abbreviation of any other species somehow
involved in the experiment. Like `LC4e`, this field is not propagated
from the higher-level objects, so its value (if it has one) should be
explicitly repeated at the cluster level.


### LC4g (Stage and tissue)

A [TAP](TAP.notes)-like statement that is used to indicate the type of
the cells in this cluster. Should be of the form

    <e><t>FBdv term<a>FBbt term<s><note>free text comment

The `<a>` field must be set to a FBbt term representing the cell type
(note that Peeves only checks whether the term is a valid FBbt term; it
does not check whether it represents a cell type or any other kind of
anatomical object present in FBbt).

The `<t>` field should be set to the same value as in the `LC4g` field
at the corresponding biosample level, including any `male` or `female`
qualifier if relevant.

Examples:

* `<e><t>wandering third instar larval stage<a>embryonic/larval crystal cell<s><note>`

* `<e><t>third instar larval stage | male<a>prohemocyte<s><note>`


### LC6d (Dataset/collection members stored in database)

At the cluster level this should always say `N`.


### LC6e (Number of entities in dataset/collection)

This field should be set to the number of isolated cells in the
corresponding cluster, extracted from the raw data [as described
here][Extracting data for FlyBase].


### LC6f (Comment on number of entities in dataset/collection)

This is a description of what the value stored in `LC6e` means. On the
dataset report, this will be displayed immediately next to the value of
`LC6e`, so this should be set to something like `cells in this cluster`.


# Producing the summarised expression table

The “summarised expression table” is basically a condensed version of
the Normalised Counts table, in which expression data for individual
cells are conflated into expression data for clusters of cells.

This table can and should only be generated for datasets where inferred
cell type annotations are available.

This step is performed automatically by the [Fuzzy-onions program][Using
Fuzzy-onions], but for reference, the procedure to generate the table
fragment for one cluster is as follows:

1. Filter the Experiment Design table to keep only the rows
   corresponding to the desired cluster in the desired sample.

2. Extract the cell identifiers for those rows.

3. Filter the Normalised Counts table to keep only the columns
   corresponding to the retained cell identifiers.

4. For each gene (row) in the filtered normalised counts table:

   1. compute the sum of all values in the row, divided by the number of
      columns that have a non-null value: this is the _mean expression
      level_ of the gene across all cells of the cluster in which the
      gene is expressed;

   2. compute the count of columns that have a non-null value, divided
      by the count of all columns; this is the _spread_ or _extent of
      expression_, the proportion of cells in the cluster in which the
      gene is expressed.

   3. store in the summarised table:

      * the gene identifier, in a column named `genes`,

      * the mean expression level, in a column named `mean_expr`,

      * the extend of expression, in a column named `spread`,

      * the cell type associated to the cluster (as a FBbt term), in a
        column named `cell_type`,

      * the [symbol given to the cluster](#cluster-LC1a) in the proforma
        produced for the dataset, in a column named `sample`.

Note: The `cell_type` column is not, strictly speaking, necessary, as
the cell type information will be found in the `LC4g` field of the
dataset object representing the cluster, to which the row is linked
through the `sample` column. Its inclusion here is intended to
facilitate consistency checks on the produced table.

Repeat the procedure for all clusters in all samples in the dataset.
Write the complete table in a TSV-formatted file, with a header line
containing the column names and starting with a hash character (‘#’).
Name that file after the record that created (or will create, if it has
not been loaded yet) the dataset object in the database, with an added
`.scrnaseq.tsv` extension.

Once the proforma has been submitted for loading during an epicycle,
upload the TSV file to the Harvard FTP server, in the
`fromcam/scrnaseq/for_load` folder.  The table will be automatically
loaded at the end of the current release cycle.  The work of the
Cambridge scRNAseq curators stops here – until the next dataset, that
is.

# Using Fuzzy-onions

Fuzzy-onions (_fzo_) is a Python script designed to facilitate the
curation of scRNAseq data. It fully automatizes the [extraction of data
from the raw files][Extracting data for FlyBase] and the [production of
the summarised expression table][Producing the summarised expression
table], partially automatizes the [production of the proforma][Producing
the proforma], and provides tools to help with other steps of the
procedure.

When using Fuzzy-onions, the procedure to curate a dataset after it has
been made available on the SCEA becomes:

1. Fetch and explore the dataset using Fuzzy-onions.

2. Prepare a JSON file describing the dataset, including corrections if
   needed.

3. Have Fuzzy-onions generate

  * the correction files intended for the SCEA;

  * a partially pre-filled proforma;

  * the summarised expression table.

4. Complete the proforma automatically generated.

5. Dispatch the resulting files:

  * send the correction files to the EBI curators;

  * submit the proforma for loading;

  * upload the summarised expression table to the Harvard FTP server.


## Setup

Fuzzy-onions is maintained at <https://github.com/FlyBase/fuzzy-onions>.
To install it:

```sh
$ python -m pip install git+https://github.com/FlyBase/fuzzy-onions
```

Check the installation by running `fzo --version`. You should get a
version message.

Fuzzy-onions needs to know about two directories:

* A directory where the raw data from the EBI’s Single-Cell Expression
  Atlas will be downloaded; by default, this is `~/scRNAseq/raw`.
Fuzzy-onions will automatically create this directory if it does not
exist.

* A directory containing the proforma templates, by default
  `~/SVN_folders/proformae`.

If you need to change the default directories, use the `fzo conf`
command to edit Fuzzy-onions’ configuration file and set the correct
directories:

```ini
[store]
production: /Users/<your_username>/scRNAseq/raw

[proformae]
directory: /Users/<your_username>/SVN_folders/proformae
```

Replace the `production` key with a `staging` key to instruct
Fuzzy-onions to fetch the datasets from the EBI‘s staging server. You
may also use both the `production` and the `staging` key to fetch
datasets from both servers; in that case, Fuzzy-onions will look for a
dataset first on the production server, then on the staging server if
the dataset was not found on the production server.


## Fetching and exploring a dataset

Once a dataset is made available by the EBI curators on the SCEA
website, get the Dataset ID [as described here][Informations from the
SCEA website]. In this section, we’ll take dataset `E-MTAB-8698` as an
example.

Download the dataset with Fuzzy-onions:

```sh
$ fzo download E-MTAB-8698
```

[As described previously][Data files], you then need to explore the
dataset, and particularly the Experiment Design table provided by the
EBI curators. This table lists all the individual cells of the dataset,
with a number of cell-level annotations.

Fuzzy-onions provides a basic “dataset explorer” mode to get a quick
glimpse of the contents of the dataset. After a dataset has been
downloaded, start the explorer and open the dataset to explore:

```sh
$ fzo explorer
fzo-explorer> open E-MTAB-8698
```

The first thing to do is then to list the columns of the experiment
design table, to see which annotations are available:

```sh
fzo-explorer> columns
 0: Assay
 1: Sample Characteristic[organism]
 2: Sample Characteristic Ontology Term[organism]
 3: Sample Characteristic[strain]
 4: Sample Characteristic Ontology Term[strain]
 5: Sample Characteristic[individual]
 6: Sample Characteristic Ontology Term[individual]
 7: Sample Characteristic[age]
 8: Sample Characteristic Ontology Term[age]
 9: Sample Characteristic[developmental stage]
10: Sample Characteristic Ontology Term[developmental stage]
11: Sample Characteristic[sex]
12: Sample Characteristic Ontology Term[sex]
13: Sample Characteristic[genotype]
14: Sample Characteristic Ontology Term[genotype]
15: Sample Characteristic[stimulus]
16: Sample Characteristic Ontology Term[stimulus]
17: Sample Characteristic[organism part]
18: Sample Characteristic Ontology Term[organism part]
19: Sample Characteristic[cell type]
20: Sample Characteristic Ontology Term[cell type]
21: Factor Value[infect]
22: Factor Value Ontology Term[infect]
23: Factor Value[inferred cell type - ontology labels]
24: Factor Value Ontology Term[inferred cell type - ontology labels]
25: Factor Value[inferred cell type - authors labels]
26: Factor Value Ontology Term[inferred cell type - authors labels]
```

The “Sample Characteristics” columns collectively describe the sample a
given cell belongs to. You should check the values of these columns and
compare them to the paper. Use the `values -c X` command to list the
values of column `X` and the number of cells annotated with that value,
and compare with what you would expect from the paper.

For example, the Cattenoz et al. (2020) paper (FBrf0245988), which
reported the E-MTAB-8698 dataset, describes a scRNAseq experiment on
wild-type female third instar larvae split in two experimental
conditions: a control group, and a group of larvae infested with
Leptopilina boulardi. Let’s then check the values of the “Sample
Characteristics” columns for “genotype” (column 13), “sex” (column 11),
“developmental stage” (column 9), and “stimulus” (column 15):

```sh
fzo-explorer> values -c 13
wild-type genotype: 14180 cells

fzo-explorer> values -c 11
female: 14180 cells

fzo-explorer> values -c 9
wandering third instar larval stage: 14180 cells

fzo-explorer> values -c 15
none: 5966 cells
infested by the wasp Leptopilina boulardi: 8214 cells
```

All cells in the dataset are annotated as coming from wild-type, female,
third-instar larvae, and as having been subjected either to nothing
(“none”) or to wasp infestation. Therefore, the contents of the dataset
fit with the description from the paper, and everything is fine.

Another thing to check at this time is the presence of cell type
annotations. Here, there are two columns for cell type annotations:
`Factor Value[inferred cell type - authors labels]` (column 25) and
`Factor Value[inferred cell type - ontology labels]` (column 23). The
first one contains the cell types as provided by the authors, which
again should match what they report in the paper:

```sh
fzo-explorer> values -c 25
plasmatocyte-2: 2223 cells
plasmatocyte-1: 2062 cells
plasmatocyte-robo2: 838 cells
plasmatocyte-0: 2697 cells
plasmatocyte-3: 1413 cells
plasmatocyte-Rel: 936 cells
plasmatocyte-Lsp: 343 cells
plasmatocyte-Pcd: 370 cells
plasmatocyte-vir1: 746 cells
lamellocyte-2: 263 cells
plasmatocyte-Inos: 590 cells
plasmatocyte-prolif: 135 cells
plasmatocyte-AMP: 42 cells
crystal cell: 55 cells
lamellocyte-1: 1065 cells
plasmatocyte-Implasmatocyte2: 4 cells
```

The second column contains the cell types as inferred by the EBI
curators. They are based on the cell types provided by the authors, but
correspond to actual FBbt terms for cell types:

```sh
fzo-explorer> values -c 23
plasmatocyte: 12399 cells
lamellocyte: 1328 cells
crystal cell: 55 cells
```

Here, the EBI curators have regrouped all the different “plasmatocyte-X”
clusters into a single “plasmatocyte” cell type, the two “lamellocyte-1”
and “lamellocyte-2” clusters into a single “lamellocyte” cell type; the
original “crystal cell” annotation has been left as it is, since it
already corresponds to an actual cell type.

Use the `select` command to filter the experiment design table on
specific values. For example, to get only the cells from the “infested”
sample:

```sh
fzo-explorer> select “15:infested by the wasp Leptopilina boulardi”
infested by the wasp Leptopilina boulardi: 8214 cells
```

The filter remains in place until it is changed or removed by another
select command. While the filter is in place, you can then get the
distribution of cell types in the filtered set:

```sh
fzo-explorer> values -c 23
plasmatocyte: 6481 cells
lamellocyte: 1322 cells
crystal cell: 14 cells
```

## Preparing a dataset description file

In order to use Fuzzy-onions to accomplish the steps 5–9 of the curation
procedure, you need to prepare a JSON-formatted file containing some
information about the dataset. Fuzzy-onions will then use that file as
input.

### General structure

The file should contain a single JSON dictionary with the following
keys:

* `dataset_id`: This is the EBI’s Dataset ID for the dataset
  (E-MTAB-8698 in the example above).

* `symbol`: This will be the [symbol for the project-level dataset
  object](#project-LC1a); all other symbols will be automatically
derived from it.

* `cell_types_column`: The name of the column in the experiment design
  table containing the cell type annotations. This key may be absent if
  the dataset does not contain cell type annotations. For the
  E-MTAB-8698 dataset, this would be `Factor Value[inferred cell type -
  ontology labels]`.

* `conditions`: A list containing the names of the columns used to
  distribute individual cells into different samples. For the
  E-MTAB-8698 dataset, as we’ve seen above this would the column `Sample
  Characteristics[stimulus]`, which divides the cells between a control
  sample and an infested sample.

* `corrections`: A list of dictionaries describing the corrections to
  apply to the experiment design table. This will be described [in a
  following section][The corrections section].

* `samples`: A list of dictionaries describing the different samples.
  Each sample is in turn described by the following keys:

  * `symbol`: A suffix to append to the project-level symbol to make up
    the sample-level symbol.

  * `selectors`: A list of values used to select the cells that belong
    to this sample. There should be as many values as there are columns
    in the `Conditions` key, in the same order. For example, for the
    control sample of the E-MTAB-8698 dataset, the list should contain
    the single value `none`, whereas for the infested sample, it should
    contain the single value `infested by the wasp Leptopilina boulardi`.

  * `title`: A free-form title for the sample. It should be chosen
     following the [guidelines set forth above](#biosample-LC6g), but
     should start with a lowercase letter.

  * `stage`: A term from FBdv indicating the developmental stage of the
     sample.

Here is an example of a dataset description file for the E-MTAB-8698
dataset (minus the `corrections` section):

```json
{
    "dataset_id": "E-MTAB-8698",
    "symbol": "scRNAseq_2020_Cattenoz",
    "cell_types_column": "Factor Value[inferred cell type - ontology labels]",
    "conditions": [
        "Sample Characteristics[stimulus]"
    ],
    "samples": [
        {
            "symbol": "_NI",
            "selectors": [ "none" ],
            "title": "hemocytes from non-infested third instar larvae",
            "stage": "wandering third instar larval stage"
        },
        {
            "symbol": "_WI",
            "selectors": [ "infested by the wasp Leptopilina boulardi" ],
            "title": "hemocytes from wasp-infested third instar larvae",
            "stage": "wandering third instar larval stage"
        }
    ]
}
```

### The corrections section

The `corrections` section that may be specified in a dataset description
file is intended to describe the corrections to apply to the experiment
design table. See [the corresponding section][Amending the SCEA
annotations] for the types of corrections to do and the rationale behind
them.

The section contains an array of dictionaries, where each dictionary is
a correction set. The syntax for a single correction set is as follows:

```json
{
    "target": "<internal|scea>",
    "source": "<name of the column containing the values to correct>",
    "destination": "<name of the column where to store the corrected values>",
    "values": [
        [ "<value to correct>", "<corrected value>", "<comment>" ],
        [ "<value to correct>", "<corrected value>", "<comment>" ],
        ...
    ]
}
```

The `target` key indicates whether this correction set is intended for
FlyBase (`internal`) or for the SCEA. That key is optional, if it is
absent then the correction set is assumed to be intended for both
FlyBase and the SCEA.

The `source` key indicates which column contains the values that need to
be corrected. If several columns require corrections, each column should
have its own correction set.

The `destination` key indicates which column should contain the
corrected values. It can be the name of an existing column in the
experiment design table (in which case existing values in that column
will be overwritten), or the name of a new column to be added to the
table. If that key is absent, the `destination` column is assumed to be
the same as the `source` column.

Finally, the `values` key is an array of triples representing the
corrections themselves. The first element of each triple is the value to
search for in the `source` column, the second element is the corrected
value to write into the `destination` column, and the last element is an
optional comment.

#### Appending FlyBase identifiers to annotations

[As explained before][Adding FlyBase identifiers], whenever possible the
scRNAseq curators should propose FlyBase identifiers to complement the
annotations already made by the EBI curators.

For example, in the E-MTAB-8698 dataset, the column `Sample
Characteristics[strain]` contains the value `Oregon-R` (annotated by the
EBI curators, based on what they found in the paper). The Oregon-R
strain is represented in FlyBase, so we can add a supplementary column
containing its FlyBase identifier. This is done with the following
correction set:

```json
{
    "source": "Sample Characteristics[strain]",
    "destination": "Sample Characteristics[strain - flybase identifier]",
    "values": [
        [ "Oregon-R", "http://flybase.org/reports/FBsn0000276", "" ]
    ]
}
```

Basically, that example reads: “On every row containing `Oregon-R` in
the `Sample Characteristics[strain]` column, put
`http://flybase.org/reports/FBsn0000276` in the `Sample
Characteristics[strain - flybase identifier]` column (to be created).”

When adding FlyBase identifiers, they should always go into a column
with the same name as the column containing the original value, with the
suffix `- flybase identifier` added (per agreement with the EBI
curators). The exception is if/when we add FlyBase identifiers for
transgenes that have been somehow used to select cells (typically a GAL4
driver used to drive GFP expression for FACS) – then the identifiers
should always go into a column named `Sample Characteristics[FACS marker
- flybase identifier]`, regardless of how the original column was named.

#### Correcting cell types

Let’s have another look at the cell types provided by the EBI curators
in dataset E-MTAB-8698:

```sh
fzo-explorer> values -c 23
plasmatocyte: 12399 cells
lamellocyte: 1328 cells
crystal cell: 55 cells
```

This is all correct for the SCEA. However, because the samples of this
dataset have been obtained from third instar larvae, for FlyBase we want
to use `embryonic/larval plasmatocyte` and `embryonic/larval crystal
cell`, instead of `plasmatocyte` and `crystal cell` (`lamellocyte` is
fine; these cells only exist in larvae and therefore there is no
`embryonic/larval lamellocyte` term in FBbt).

This is done with the following correction set:

```json
{
    "target": "internal",
    "source": "Factor Value[inferred cell type - ontology labels]",
    "values": [
        [ "plasmatocyte", "embryonic/larval plasmatocyte", "" ],
        [ "crystal cell", "embryonic/larval crystal cell", "" ]
    ]
}
```

Basically, that example reads: “In the `Factor Value[inferred cell type
- ontology labels]` column, replace all instances of `plasmatocyte` by
`embryonic/larval plasmatocyte` and all instances of `crystal cell` by
`embryonic/larval crystal cell` – for FlyBase only.”

## Producing the output files

Once the dataset description file is complete (including the correction
sets, if needed), it can be given to Fuzzy-onions to automatically
generate

* the files to send back to the EBI curators with the corrected
  annotations;

* a pre-filled proforma for the dataset;

* the summarised expression table.

In the examples in the following sections, the dataset description file
will be assumed to be named `E-MTAB-8698.json`. The actual name of the
file is irrelevant and has no meaning for Fuzzy-onions.


### Correction files for EBI curators

Use the `fzo fixscea` command to generate the files that need to be sent
back to the EBI curators:

```sh
$ fzo fixscea E-MTAB-8698.json
```

This will generate two files:

* `experiment-design.with-fbids.tsv`, containing a version of the
  experiment design table enriched with FlyBase identifiers;

* `celltypes-fbterms.tsv`, containing a table of corrected cell type
  annotations in the format expected by EBI curators (if there were cell
  type annotations to correct in the dataset, which is actually not the
  case for the E-MTAB-8698 dataset).

Send those files to the EBI curators (Nancy) by e-mail.


### Pre-filled proforma

Use the `fzo proforma` command to generate a proforma for the dataset:

```sh
$ fzo proforma -o proforma.coll E-MTAB-8698.json
```

The resulting `proforma.coll` file will be pre-filled with all the
informations that can be automatically extracted from the data files and
the dataset description file. Fill the remaining fields according [to
the guidelines][Producing the proforma]. Rename it according to the
record naming guidelines in the Curation Manual (e.g. for Cambridge,
`dg99.coll` for the 99th record submitted by Damien Goutte-Gattat), have
it checked by Peeves, then submit it for loading.


### Summarised expression table

Use the `fzo sumexpr` command to generate the summarised expression table:

```sh
$ fzo sumexpr -o <recordname>.scrnaseq.tsv E-MTAB-8698
```

Make sure the output file is named according to the record submitted
above (e.g. `dg99.scrnaseq.tsv`) and upload it to the
`fromcam/scrnaseq/for_load` folder on the Harvard FTP server.
