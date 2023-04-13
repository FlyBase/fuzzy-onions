# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2023 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import click
from click_shell.core import make_click_shell
from pronto import Ontology

CELL_TYPES_COLUMN = 'Factor Value[inferred cell type - ontology labels]'


def _get_fbbt_reverse_dict(path):
    fbbt = Ontology(path)
    reverse_dict = {}
    for term in fbbt.terms():
        reverse_dict[term.name] = term.id
    return reverse_dict


def _get_cell_types_for_dataset(dataset, fbbt):
    if CELL_TYPES_COLUMN not in dataset.experiment_design:
        return []

    cell_types = []
    for cell_type in sorted(
        dataset.experiment_design[CELL_TYPES_COLUMN].dropna().unique()
    ):
        term_id = fbbt.get(cell_type, None)
        cell_types.append((cell_type, term_id))
    return cell_types


@click.group(name="annots", invoke_without_command=True)
@click.pass_context
def annots(ctx):
    """Manipulate cell type annotations."""

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-annots>")
        shell.cmdloop()


@annots.command()
@click.argument('dsids', nargs=-1)
@click.option(
    '--fbbt',
    '-f',
    'fbbt_path',
    default='http://purl.obolibrary.org/obo/fbbt.obo',
    type=click.Path(exists=True),
    help="Validate against the specified ontology file.",
)
@click.option(
    '--invalid-only',
    '-i',
    is_flag=True,
    default=False,
    help="List only invalid cell types.",
)
@click.pass_obj
def validate(ctx, dsids, fbbt_path, invalid_only):
    """Validate cell type annotations in specified dataset(s).

    This command checks that cell type annotations in the specified
    dataset(s) are valid FBbt terms. It prints a list of all cell
    types used in a given dataset along with the corresponding FBbt
    term ID (or None if the annotation is not a valid FBbt term).
    """

    term_ids_by_name = _get_fbbt_reverse_dict(fbbt_path)

    if 'all' in dsids:
        datasets = ctx.raw_store.datasets
    else:
        datasets = []
        for dsid in dsids:
            datasets.append(ctx.raw_store.get(dsid))

    print("dataset,cell type annotation,term id")
    for dataset in datasets:
        for cell_type, term_id in _get_cell_types_for_dataset(
            dataset, term_ids_by_name
        ):
            if not invalid_only or term_id is None:
                print(f"{dataset.id},{cell_type},{term_id}")
