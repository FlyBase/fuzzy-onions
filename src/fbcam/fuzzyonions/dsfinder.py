# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import os.path
from shutil import copyfile
import re

import click
from click_shell.core import make_click_shell
from pandas import concat, read_csv, DataFrame
import numpy
import pymupdf
import llm
import logging
import openai


class DiscoverContext(object):
    """A helper object for all discovery operations."""

    def __init__(self, config, database):
        self._config = config
        self._database = database
        self._pattern = re.compile(self._config.get('textmining', 'pattern'), re.IGNORECASE)
        self._model = None

    @property
    def flybase_table_file(self):
        return self._config.get(
            'discovery', 'flybase_papers_list', fallback='flybase-papers.tsv'
        )

    @property
    def scea_dataset_file(self):
        return self._config.get('discovery', 'scea_dataset_list')

    @property
    def scea_staging_dataset_file(self):
        return self._config.get('discovery', 'scea_staging_list')

    @property
    def cursor(self):
        return self._database.cursor

    @property
    def model(self):
        if not self._model:
            self._model = llm.get_model('gpt-4o-mini')
            self._model.key = self._config.get('textmining', 'openai_key')
            logging.getLogger("openai").setLevel(logging.ERROR)
            logging.getLogger("httpx").setLevel(logging.ERROR)
        return self._model

    def filter_out_known_datasets(self, pmids):
        table = read_csv(self.scea_dataset_file, dtype=str)
        results = [p for p in pmids if p not in table['PubMed ID'].values]

        table = read_csv(self.scea_staging_dataset_file, dtype=str)
        return [p for p in results if p not in table['PMID'].values]

    def get_dataset_fbrfs(self):
        """Get all references that have been flagged with 'dataset'."""

        query = f'''SELECT pub.uniquename
                    FROM
                             pub
                        JOIN pubprop USING (pub_id)
                        JOIN cvterm  ON pubprop.type_id = cvterm.cvterm_id
                        JOIN cv      USING (cv_id)
                    WHERE
                             cv.name         LIKE 'pubprop type'
                         AND cvterm.name     LIKE 'harv_flag'
                         AND pubprop.value   LIKE 'dataset'
                         AND pub.is_obsolete = 'f';'''
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def get_internalnotes_for_fbrf(self, fbrf):
        """Get the internal notes for a given reference."""

        query = f'''SELECT pubprop.value
                    FROM
                             pub
                        JOIN pubprop USING (pub_id)
                        JOIN cvterm  ON pubprop.type_id = cvterm.cvterm_id
                        JOIN cv      USING (cv_id)
                    WHERE
                             cv.name        LIKE 'pubprop type'
                         AND cvterm.name    LIKE 'internalnotes'
                         AND pub.uniquename LIKE '{fbrf}';'''
        self.cursor.execute(query)
        return self.cursor.fetchall()

    def get_dataset_accessions_for_fbrf(self, fbrf):
        """Get dataset accession numbers associated with a reference."""

        notes = self.get_internalnotes_for_fbrf(fbrf)
        accessions = []
        for note in notes:
            for line in note[0].split('\n'):
                if not line.startswith('Dataset:'):
                    continue
                line = line[8:].strip().split('.')[0]
                for item in line.split(','):
                    item = item.strip()
                    if item.startswith('and '):
                        item = item[3:].strip()
                    if item == 'pheno':
                        continue
                    accessions.append(item)
        return ','.join(accessions)

    def get_pmid_for_fbrf(self, fbrf):
        """Get the PubMed ID for a given reference."""

        query = f'''SELECT dbxref.accession
                    FROM
                             pub
                        JOIN pub_dbxref USING (pub_id)
                        JOIN dbxref     USING (dbxref_id)
                        JOIN db         USING (db_id)
                    WHERE
                             db.name        LIKE 'pubmed'
                         AND pub.uniquename LIKE '{fbrf}';'''
        self.cursor.execute(query)
        return ','.join([a for a, in self.cursor.fetchall()])

    def get_citation_for_fbrf(self, fbrf):
        """Get a text citation for a given reference."""

        query = f'''SELECT pubauthor.surname, pubauthor.rank, pub.pyear
                    FROM
                             pubauthor
                        JOIN pub USING (pub_id)
                    WHERE
                             pub.uniquename LIKE '{fbrf}'
                    ORDER BY pubauthor.rank;'''
        self.cursor.execute(query)
        res = self.cursor.fetchall()
        year = res[0][2]
        n = len(res)
        if n == 1:
            return f'{res[0][0]}, {year}'
        elif n == 2:
            return f'{res[0][0]} and {res[1][0]}, {year}'
        else:
            return f'{res[0][0]} et al., {year}'

    def is_about_new_dataset(self, fbrf, pmid, use_llm=False):
        """Attempt to determine whether the given reference describes a
        new dataset.

        :param fbrf: FlyBase reference identifier
        :param pmid: PubMed identifier
        :return: a dictionary where the key represents the method used
                 to test the reference, and the value is the result;
                 empty if the full text of the reference is not
                 available.
        """
        fulltext = None
        res = {}
        pdf_file = self.get_fulltext_file(fbrf, pmid)
        if pdf_file:
            try:
                doc = pymupdf.open(pdf_file)
                fulltext = "".join([p.get_text() for p in doc])
            except:
                pass

        if fulltext:
            nmatch = len(self._pattern.findall(fulltext))
            res['grep'] = nmatch

            if use_llm and nmatch > 0:
                prompt = """
                    Can you tell whether, in the following text of a
                    scientific paper, the authors are describing
                    (possibly among other things) a new single-cell
                    RNA-sequencing dataset, or merely referencing a
                    dataset that has been published previously? If
                    they describe a new dataset, can you identify the
                    scRNAseq technology used (e.g. Chromium 10x)?
                    Answer with a simple, short sentence.
                    """
                try:
                    response = self.model.prompt(prompt, fragments=[fulltext])
                    res['llm'] = response.text()
                except openai.BadRequestError:
                    pass
        return res

    def get_fulltext_file(self, fbrf, pmid):
        """Get the full text of the given reference.

        :param fbrf: FlyBase reference identifier
        :param pmid: PubMed identifier
        :return: the full path to the PDF file, or None if the full text
                 is not available.
        """

        vm_base = self._config.get('textmining', 'vm_base')
        # First look up in the main PDF archive folder
        num = fbrf[4:7]
        fullpath = os.path.join(vm_base, f'pdf/{num}0000-{num}9999/{fbrf}.pdf')
        if os.path.exists(fullpath):
            return fullpath
        # Then look up in the staging area
        fullpath = os.path.join(vm_base, f'staging/pdf/{fbrf}.pdf')
        if os.path.exists(fullpath):
            return fullpath
        # Then look up in the PubGet folder
        fullpath = os.path.join(vm_base, f'staging/pubget/{pmid}.pdf')
        if os.path.exists(fullpath):
            return fullpath
        return None


class ResultsTable(object):
    """A class to encapsulate the results of a discovery operation."""

    def __init__(self, filename):
        self._filename = filename
        if os.path.exists(self._filename):
            self._table = read_csv(filename, sep='\t', dtype=str)
        else:
            self._table = DataFrame(
                columns=[
                    'FBrf',
                    'PMID',
                    'Citation',
                    'Accessions',
                    'Mentions',
                    'Confirmed',
                    'Organ/tissue',
                    'Comments'
                ]
            )

    def save(self, filename=None):
        if filename is None:
            filename = self._filename
            copyfile(filename, f'{filename}.bak')
        self._table.to_csv(filename, index=False, sep='\t')

    def add_new_references(self, fbrfs):
        new_rows = []
        for fbrf in fbrfs:
            if fbrf not in self._table['FBrf'].values:
                new_rows.append({'FBrf': fbrf})
        self._table = concat(self._table, DataFrame(data=new_rows))
        return len(new_rows)

    def fill_missing(self, ctx):
        incompletes = self._table[self._table['Citation'].isna()]['FBrf'].values
        with click.progressbar(incompletes) as bar:
            for fbrf in bar:
                pmid = ctx.get_pmid_for_fbrf(fbrf)
                accessions = ctx.get_dataset_accessions_for_fbrf(fbrf)
                citation = ctx.get_citation_for_fbrf(fbrf)

                mask = self._table['FBrf'] == fbrf
                self._table.loc[mask, 'PMID'] = pmid
                self._table.loc[mask, 'Citation'] = citation
                self._table.loc[mask, 'Accessions'] = accessions
        return len(incompletes)

    def fill_relevance(self, ctx, use_llm=False):
        candidates = self._table.loc[self._table['Mentions'].isna(), ('FBrf', 'PMID')].values
        pos = 0
        nofulltext = 0
        with click.progressbar(candidates) as bar:
            for fbrf, pmid in bar:
                res = ctx.is_about_new_dataset(fbrf, pmid, use_llm)
                nmatch = res.get('grep')
                if nmatch is None:
                    nofulltext += 1
                    nmatch = 'Full-text not available'
                elif nmatch > 0:
                    pos += 1
                self._table.loc[self._table['FBrf'] == fbrf, 'Mentions'] = nmatch
                if 'llm' in res:
                    self._table.loc[self._table['FBrf'] == fbrf, 'Comments'] = res.get('llm')
        return (pos, nofulltext)

    def clear_relevance_data(self):
        self._table.loc[:,'Mentions'] = numpy.nan


@click.group(invoke_without_command=True)
@click.pass_context
def discover(ctx):
    """Access the dataset discovery commands."""

    ctx.obj = DiscoverContext(ctx.obj.config, ctx.obj.database)

    if not ctx.invoked_subcommand:
        shell = make_click_shell(ctx, prompt="fzo-discover> ")
        shell.cmdloop()


@discover.command()
@click.option(
    '--filename',
    '-f',
    default=None,
    metavar='FILENAME',
    help="""Use the given table file instead of the one
                      specified in the configuration.""",
)
@click.option(
    '--output',
    '-o',
    default=None,
    metavar='FILENAME',
    help="""Write the updated table to the specified file.
            The default is to write to the original file.""",
)
@click.option('--use-llm',
    is_flag=True,
    default=False,
    help="""Use a LLM to try determinining whether a reference
            is about a new dataset""")
@click.pass_obj
def findnew(obj, filename, output, use_llm):
    """Find new FlyBase references that may be about scRNAseq."""

    if filename is None:
        filename = obj.flybase_table_file
    table = ResultsTable(filename)

    click.echo("Fetching dataset-flagged references...")
    n = table.add_new_references([f[0] for f in obj.get_dataset_fbrfs()])
    if n == 0:
        return

    click.echo(f"New dataset-flagged references: {n}")
    click.echo("Fetching additional data...")
    table.fill_missing(obj)

    click.echo("Checking for scRNAseq relevance...")
    n, nofulltext = table.fill_relevance(obj, use_llm=use_llm)
    click.echo(f"Number of potentially relevant references: {n}")
    click.echo(f"Number of references without full text: {nofulltext}")

    table.save(output)


@discover.command()
@click.argument('filename')
@click.option(
    '--output',
    '-o',
    type=click.File('w'),
    default='-',
    metavar='FILENAME',
    help="""Write the result to the specified file.
            The default is to write to standard output."""
)
@click.option('--force',
    is_flag=True,
    default=False,
    help="""Ignore previous relevance results.""")
@click.option('--use-llm',
    is_flag=True,
    default=False,
    help="""Use a LLM to try determinining whether a reference
            is about a new dataset""")
@click.pass_obj
def checknew(obj, filename, output, force, use_llm):
    """Check whether the provided references may be about scRNAseq."""

    table = ResultsTable(filename)

    click.echo("Fetching additional data..")
    n = table.fill_missing(obj)
    if n > 0:
        click.echo(f"Updated data for {n} references")

    if force:
        table.clear_relevance_data()
    click.echo("Checking for scRNAseq relevance...")
    n, nofulltext = table.fill_relevance(obj, use_llm)
    click.echo(f"Number of potentially relevant references: {n}")
    click.echo(f"Number of references without full text: {nofulltext}")

    table.save(output)


@discover.command()
@click.option(
    '--output',
    '-o',
    type=click.File('w'),
    default='-',
    help="Send output to the specified file.",
)
@click.pass_obj
def toscea(obj, output):
    """Make a list of datasets that are new to the SCEA."""

    table = read_csv(obj.flybase_table_file, sep='\t', dtype=str)
    subset = table[table['Confirmed'] == 'yes']
    unknown_pmids = obj.filter_out_known_datasets(subset['PMID'].values)

    subset = table[table['PMID'].isin(unknown_pmids)]
    subset = table.loc[
        table['PMID'].isin(unknown_pmids),
        ['FBrf', 'PMID', 'Accessions', 'Citation', 'Organ/tissue', 'Comments'],
    ]
    subset.to_csv(output, sep='\t', index=False)
