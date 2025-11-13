# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2025 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import json
import logging
from datetime import datetime

import click
from click_shell.core import make_click_shell
from sssom.parsers import parse_sssom_table
from .export import FBCV_SSSOM_SET

@click.group(invoke_without_command=True)
@click.pass_context
def util(ctx):
    """Access the utility commands."""
    
    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-util> ")
        shell.cmdloop()

@util.command()
@click.argument('meta_in', type=click.Path(exists=True))
@click.argument('sumexpr_in', type=click.Path(exists=True))
@click.option(
    '--meta-out',
    type=click.File('w'),
    default='meta.json',
    help="Write the metadata to the specified file."
)
@click.option(
    '--sumexpr-out',
    type=click.File('w'),
    default='sumexpr.tsv',
    help="Write the summarised expression data to the specified file."
)
@click.pass_context
def bta(ctx, meta_in, sumexpr_in, meta_out, sumexpr_out):
    """Converts a data dump into a file suitable for the Alliance.
    
    Given the "flat files" for scRNAseq data that can be found downloaded
    from FlyBase (both the file containing the dataset metadata, and the
    file containing the summarised expression data), this command will
    produce files that _should_ be suitable for ingestion into the Alliance
    of Genome Resources.
    """
    
    # Metadata
    with open(meta_in, 'r') as f:
        dumped_datasets = json.load(f)['datasets'] 
    dc = DatasetConverter()
    json.dump(dc.convert(dumped_datasets), meta_out, indent=4)
    
    # Summarised expression data
    print("cluster_id\tcluster_symbol\tgene_id\tmean_expr\tspread", file=sumexpr_out)
    ignored = {}
    with open(sumexpr_in, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            values = line.strip().split('\t')
            cluster_id = values[7]
            cluster_symbol = values[8]
            if cluster_id not in dc._cluster_ids:
                if cluster_id not in ignored:
                    logging.warning(f"Ignoring cluster {cluster_id} ({cluster_symbol})")
                    ignored[cluster_id] = 1
            gene_id = values[11]
            mean = float(values[13])
            spread = float(values[14])
            print(f"FB:{cluster_id}\t{cluster_symbol}\t{gene_id}\t{mean:0.6f}\t{spread:0.6f}", file=sumexpr_out)
                
            
            
        
class DatasetConverter():
    """A helper class to convert a dumped dataset into an Alliance dataset."""
    
    def __init__(self):
        self._today = datetime.now().isoformat(timespec='seconds')
        self._fbcv_map = parse_sssom_table(FBCV_SSSOM_SET)
        self._cluster_ids = {}
    
    def convert(self, datasets):
        dataset_objects = []
        sample_objects = []
        for dataset in datasets:
            self._convert_project(dataset, dataset_objects, sample_objects)
        return {
            "high_throughput_expression_dataset_annotation_ingest_set": dataset_objects,
            "high_throughput_expression_dataset_sample_annotation_ingest_set": sample_objects
            }
        
    def _get_new_dict(self):
        return {
            "date_created": self._today,
            "internal": False,
            "obsolete": False
        }
        
    def _convert_project(self, dataset, dataset_objects, sample_objects):
        logging.info(f"Converting dataset {dataset['symbol']}")
        
        d = self._get_new_dict()
        d["primary_external_id"] = dataset["id"]
        # The dataset ID in "htp_expression_dataset_dto" is the same as the
        # primary_external_id (and is in fact the MOD ID), because apart from
        # from the top-level project, none of those objects correspond to an
        # external database entity.
        d["htp_expression_dataset_dto"] = {"curie": dataset["id"]}
        d["name"] = dataset["title"]
        # WARNING: The "symbol" slot is currently only defined on the
        # HTPExpressionDatasetAnnotation, _not_ on its *DTO counterpart, so
        # the following is not compliant
        d["symbol"] = dataset["symbol"]
        
        if "reference" in dataset:
            # This is a top-level project
            d["htp_expression_dataset_dto"]["secondary_identifiers"] = dataset["accessions"]
            d["references_curies"] = [dataset["reference"]] # FIXME: Should it be converted to a FB:FBrf?
            d["category_tag_names"] = self._fbcv_to_category_tags(dataset["study_cvterms"])
        
        dataset_objects.append(d)
        
        subseries = []
        for sample in dataset["samples"]:
            if "samples" in sample:
                # This is a "sub-project", to be rendered as a distinct
                # dataset object in the Alliance schema
                self._convert_project(sample, dataset_objects, sample_objects)
                subseries.append(sample["id"])
            else:
                self._convert_sample(sample, dataset, dataset_objects, sample_objects)
        if len(subseries) > 0:
            d["sub_series"] = subseries
            
        if "analysis" in dataset and len(dataset["analysis"]["clusters"]) > 0:
            d["analysis"] = self._convert_analysis(dataset["analysis"])
            
    def _convert_sample(self, sample, parent, dataset_objects, sample_objects):
        d = self._get_new_dict()
        
        d["primary_external_id"] = sample["id"]
        d["htp_expression_sample_dto"] = {"curie": sample["id"]}
        d["htp_expression_sample_title"] = sample["title"]
        # FIXME: The "symbol" slot is missing from the
        # HTPExpressionDatasetSampleAnnotationDTO class, so the following is
        # not compliant
        d["symbol"] = sample["symbol"]
        d["dataset_ids"] = [parent["id"]]
        # FIXME: Shouldn't that slot be called "taxon_curie" in the DTO class?
        d["taxon"] = "NCBITaxon:7227"
        # FIXME: Shouldn't that slot by "genetic_sex_name" in the DTO class?
        d["genetic_sex"] = sample.get("sex", "pooled sexes")
        # FIXME: The LinkML schema defines an enum value with only a small
        # subset of allowed OBI values, but that enum is not referenced from
        # anywhere within the schema and the "htp_expression_sample_type" is
        # defined with a range of "OBITerm", suggesting that _any_ OBI term
        # can be used -- unsure whether this is actually true or not.
        d["htp_expression_sample_type_curie"] = self._fbcv_to_obi(sample["sample_type"])
        # FIXME: MMO has only one term for all scRNAseq stuff; frustratingly,
        # OBI does have more precise terms for some specific subtypes of
        # scRNAseq, but we can't use OBI here...
        d["expression_assay_used"] = "MMO:0000862"
        d["htp_expression_sample_location_dtos"] = [{
            "anatomical_structure_curie": sample["tissue"]
            }]
        d["htp_expression_sample_age_dto"] = {
            "stage_dto": {
                "developmental_stage_start_curie": sample["stage"]
                }
            }
        if "genotype" in sample:
            # FIXME: Should ideally use AGM models instead of free text (much
            # more complex but would allow to represent the exact genetic
            # entities involved, when they are known).
            d["genomic_information_dto"] = {
                "biosample_text": sample["genotype"]
                }
            
        # FIXME: No way of storing the FBcv terms representing how the sample
        # was collected/analysed.
        # FIXME: No way of storing the collection protocol, though I guess
        # this could be done with a new Note Type since this is free text?
        # FIXME: No way of storing the biological reference, if any?
        
        if "analysis" in sample and len(sample["analysis"]["clusters"]) > 0:
            d["analysis"] = self._convert_analysis(sample["analysis"])
        
        sample_objects.append(d)
        
    def _convert_analysis(self, analysis):
        # FIXME: Shouldn't there be a AnalysisDTO class?
        # With something like a "analysis_type_curie" slot?
        d = self._get_new_dict()
        d["primary_external_id"] = analysis["id"]
        d["analysis_type"] = analysis["analysis_type"]
        d["cell_count"] = analysis["cell_count"]
        d["symbol"] = analysis["symbol"]
        if "analysis_protocol" in analysis:
            d["note_dtos"] = [{
                "free_text": analysis["analysis_protocol"],
                "note_type_name": "analysis_protocol"
                }]
        # FIXME: No way of storing the FBcv term(s) representing the analysis
        # pipeline (so no way of distinguishing between SCEA-analysed datasets
        # and other datasets).
        clusters = []
        for cluster in analysis["clusters"]:
            # FIXME: Shouldn't there be a SingleExpressionClusterDTO class?
            # With something like a "cell_type_curie" slot?
            c = self._get_new_dict()
            c["primary_external_id"] = cluster["id"]
            c["symbol"] = cluster["symbol"]
            c["title"] = cluster["title"]
            c["cell_count"] = cluster["cell_count"]
            c["cell_type"] = cluster["cell_type"]
            clusters.append(c)
            self._cluster_ids[cluster["id"][3:]] = 1
        d["cell_clusters"] = clusters
        return d
        
        
    def _fbcv_to_category_tags(self, fbcv_terms):
        tag_map = self._fbcv_map.df.loc[self._fbcv_map.df["object_category"] == "Data Set Category Tags"]
        tags = []
        for fbcv_term in fbcv_terms:
            tag = tag_map[tag_map["subject_id"] == fbcv_term]["object_label"]
            if len(tag) == 1:
                tags.append(tag.values[0])
            else:
                logging.warning(f"Missing (or ambiguous) FBcv-to-Tags mapping: {fbcv_term}")
        return tags
    
    def _fbcv_to_obi(self, fbcv_term):
        obi_map = self._fbcv_map.df.loc[self._fbcv_map.df["object_source"] == "obo:obi.owl"]
        obi_term = obi_map[obi_map["subject_id"] == fbcv_term]["object_id"]
        if len(obi_term) == 1:
            return obi_term.values[0]
        else:
            logging.warning(f"Missing FBcv-to-OBI mapping: {fbcv_term}")
            return None
