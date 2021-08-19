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

    for col in ctx.dataset.experiment_design.columns:
        print(col)


@explorer.command()
@click.argument('column')
@click.option('--count-cells', '-c', is_flag=True, default=False,
              help="Print the number of cells for each value.")
@click.pass_obj
def values(ctx, column, count_cells):
    """Print unique values in a column."""

    expd = ctx.subset
    if column not in expd.columns:
        print(f"No column named '{column}'")

    for value in expd[column].unique():
        if count_cells:
            count = len(expd.loc[expd[column] == value])
            print(f"{value}: {count} cells")
        else:
            print(value)


@explorer.command()
@click.argument('column')
@click.argument('value')
@click.pass_obj
def select(ctx, column, value):
    """Select a subset of the experiment design."""

    expd = ctx.subset
    ctx.subset = expd.loc[expd[column] == value]
    print(f"{value}: {len(ctx.subset)} cells")


@explorer.command()
@click.pass_obj
def clear(ctx):
    """Clear any subset selection."""

    ctx.subset = None
