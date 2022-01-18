# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2021 Damien Goutte-Gattat
#
# Redistribution and use of this script, with or without modifications,
# is permitted provided that the following conditions are met:
#
# 1. Redistributions of this script must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import click
from click_shell import make_click_shell


@click.group(invoke_without_command=True)
@click.pass_context
def explorer(ctx):
    """Explore a dataset."""

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
@click.option('--count-cells', '-c', is_flag=True, default=False,
              help="Print the number of cells for each value.")
@click.pass_obj
def values(ctx, column, count_cells):
    """Print unique values in a column."""

    expd = ctx.subset
    col = _get_column(ctx, column)
    if not col:
        print(f"Column '{column}' no found")
        return

    for value in expd[col].unique():
        if count_cells:
            count = len(expd.loc[expd[col] == value])
            print(f"{value}: {count} cells")
        else:
            print(value)


@explorer.command()
@click.argument('selectors', nargs=-1)
@click.option('--clear', '-c', is_flag=True, default=False,
              help="Reset previous selection.")
@click.pass_obj
def select(ctx, selectors, clear):
    """Select a subset of the experiment design."""
    
    if clear:
        ctx.subset = None

    for selector in selectors:
        column, value = selector.split(':')
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
