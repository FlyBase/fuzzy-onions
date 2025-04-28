# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2023 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import os.path
from shutil import copyfile

import click
from click_shell.core import make_click_shell
from pandas import DataFrame, read_csv
from pronto import Ontology

CELL_TYPES_COLUMN = 'Factor Value[inferred cell type - ontology labels]'


class AnnotationContext(object):
    """A helper object for all commands manipulating annotations."""

    def __init__(self, store, fbbt_path):
        self.store = store
        self._fbbt_path = fbbt_path
        self._fbbt = None
        self._fbbt_reverse_dict = None

    @property
    def fbbt(self):
        if self._fbbt is None:
            self._fbbt = Ontology(self._fbbt_path)
        return self._fbbt

    @property
    def fbbt_ids(self):
        if self._fbbt_reverse_dict is None:
            self._fbbt_reverse_dict = {}
            for term in self.fbbt.terms():
                self._fbbt_reverse_dict[term.name] = term.id
        return self._fbbt_reverse_dict

    def get_cell_types_for_dataset(self, dataset, invalid_only=False):
        if CELL_TYPES_COLUMN not in dataset.experiment_design:
            return []

        cell_types = []
        for cell_type in sorted(
            dataset.experiment_design[CELL_TYPES_COLUMN].dropna().unique()
        ):
            term_id = self.fbbt_ids.get(cell_type, None)
            if not invalid_only or term_id is None:
                cell_types.append((cell_type, term_id))

        return cell_types


@click.group(name="annots", invoke_without_command=True)
@click.option(
    '--fbbt',
    '-f',
    'fbbt_path',
    type=click.Path(exists=True),
    help="""Use the specified FBbt file. Default is to download the ontology
            from the PURL server.""",
)
@click.pass_context
def annots(ctx, fbbt_path):
    """Manipulate cell type annotations."""

    if fbbt_path is None:
        fbbt_path = 'http://purl.obolibrary.org/obo/fbbt.obo'
    ctx.obj = AnnotationContext(ctx.obj.raw_store, fbbt_path)

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-annots>")
        shell.cmdloop()


@annots.command()
@click.argument('dsids', nargs=-1)
@click.option(
    '--invalid-only',
    '-i',
    is_flag=True,
    default=False,
    help="List only invalid cell types.",
)
@click.pass_obj
def validate(ctx, dsids, invalid_only):
    """Validate cell type annotations in specified dataset(s).

    This command checks that cell type annotations in the specified
    dataset(s) are valid FBbt terms. It prints a list of all cell
    types used in a given dataset along with the corresponding FBbt
    term ID (or None if the annotation is not a valid FBbt term).
    """

    if 'all' in dsids:
        datasets = ctx.store.datasets
    else:
        datasets = []
        for dsid in dsids:
            datasets.append(ctx.store.get(dsid))

    print("dataset,cell type annotation,term id")
    for dataset in datasets:
        for cell_type, term_id in ctx.get_cell_types_for_dataset(dataset, invalid_only):
            print(f"{dataset.id},{cell_type},{term_id}")


@annots.command('invalid')
@click.argument('trackfile', type=click.Path())
@click.option(
    '--output',
    '-o',
    default=None,
    metavar='FILENAME',
    help="""Write the updated table to the specified file.
            The default is to write to the original file.""",
)
@click.pass_obj
def track_invalid_annotations(ctx, trackfile, output):
    """Track all invalid cell type annotations.

    This command checks all datasets for invalid cell type annotations
    and tracks them in a TSV file.
    """

    if os.path.exists(trackfile):
        trackdata = read_csv(trackfile, sep='\t')
        backup_needed = True
    else:
        trackdata = DataFrame(columns=['dataset', 'cell type annotation', 'new term'])
        backup_needed = False

    new_rows = []

    for dataset in ctx.store.datasets:
        # Get tracked annotations
        subset = trackdata[trackdata['dataset'] == dataset.id]

        # Get all current invalid annotations
        annots = ctx.get_cell_types_for_dataset(dataset, invalid_only=True)

        for cell_type, _ in annots:
            # Get replacement term, if the annotation was already tracked
            cell_subset = subset[subset['cell type annotation'] == cell_type]
            if len(cell_subset) == 1:
                new_term = cell_subset['new term'].iloc[0]
            else:
                new_term = ''

            # Add the annotation to the new tracked data
            new_rows.append(
                {
                    'dataset': dataset.id,
                    'cell type annotation': cell_type,
                    'new term': new_term,
                }
            )

    # Compile new tracked set
    trackdata = DataFrame(columns=trackdata.columns, data=new_rows)

    if output is None:
        output = trackfile
        if backup_needed:
            copyfile(output, f'{output}.bak')
    trackdata.to_csv(output, sep='\t', index=False)
