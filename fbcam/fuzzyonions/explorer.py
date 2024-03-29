# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import click
from click_shell import make_click_shell


class ExplorerContext(object):
    """Context holder for the dataset explorer.

    This object keeps track of which dataset is currently being
    explored and which subset of rows is selected.
    """

    def __init__(self, store):
        self._store = store
        self._dataset = None
        self._subset = None
        self._subset_filters = []

    @property
    def dataset(self):
        if self._dataset is None:
            self._dataset = self._store.datasets[0]
        return self._dataset

    @property
    def subset(self):
        if self._subset is None:
            return self.dataset.experiment_design
        else:
            return self._subset

    @subset.setter
    def subset(self, subset):
        self._subset = subset
        if subset is None:
            self._subset_filters.clear()

    def load_dataset(self, dsid):
        self._dataset = self._store.get(dsid)

    def filter_subset(self, column, value):
        self.subset = self.subset.loc[self.subset[column] == value]
        self._subset_filters.append([column, value])

    def get_filter_string(self):
        if len(self._subset_filters) == 0:
            return '(all)'
        else:
            return ' > '.join([b for _, b in self._subset_filters])


@click.group(invoke_without_command=True)
@click.pass_context
def explorer(ctx):
    """Explore a dataset."""

    ctx.obj = ExplorerContext(ctx.obj.raw_store)
    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-explorer> ")
        shell.cmdloop()


@explorer.command('open')
@click.argument('dsid')
@click.pass_obj
def open_dataset(ctx, dsid):
    """Open a dataset."""

    ctx.load_dataset(dsid)


@explorer.command()
@click.pass_obj
def columns(ctx):
    """List the columns from the experiment design."""

    for i, col in enumerate(ctx.dataset.experiment_design.columns):
        print(f"{i:2d}: {col}")


@explorer.command()
@click.argument('column')
@click.option(
    '--count-cells',
    '-c',
    is_flag=True,
    default=False,
    help="Print the number of cells for each value.",
)
@click.pass_obj
def values(ctx, column, count_cells):
    """Print unique values in a column."""

    expd = ctx.subset
    col = _get_column(ctx, column)
    if not col:
        print(f"Column '{column}' no found")
        return

    for value in expd[col].dropna().unique():
        if count_cells:
            count = len(expd.loc[expd[col] == value])
            print(f"{value}: {count} cells")
        else:
            print(value)


@explorer.command()
@click.argument('selectors', nargs=-1)
@click.option(
    '--clear', '-c', is_flag=True, default=False, help="Reset previous selection."
)
@click.pass_obj
def select(ctx, selectors, clear):
    """Select a subset of the experiment design."""

    if clear:
        ctx.subset = None

    for selector in selectors:
        column, value = selector.split(':', maxsplit=1)
        col = _get_column(ctx, column)
        if not col:
            print(f"Column '{column}' not found")
            return

        ctx.filter_subset(col, value)

    filter_string = ctx.get_filter_string()
    print(f"{filter_string}: {len(ctx.subset)} cells")


@explorer.command()
@click.pass_obj
def clear(ctx):
    """Clear any subset selection."""

    ctx.subset = None


def _get_column(ctx, spec):
    if spec in ctx.dataset.experiment_design.columns:
        return spec
    elif spec.isnumeric():
        spec = int(spec)
        if spec in range(len(ctx.dataset.experiment_design.columns)):
            return ctx.dataset.experiment_design.columns[spec]

    return None
