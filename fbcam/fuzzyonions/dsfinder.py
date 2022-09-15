# fuzzy-onions - FlyBase scRNAseq scripts
# Copyright © 2022 Damien Goutte-Gattat
#
# This file is part of the Fuzzy-onions project and distributed under
# the terms of the MIT license. See the LICENSE.md file in that project
# for the detailed conditions.

import os.path
import subprocess
from shutil import copyfile

import click
from click_shell.core import make_click_shell
from pandas import concat, read_csv, DataFrame


class DiscoverContext(object):
    """A helper object for all discovery operations."""

    def __init__(self, config, database):
        self._config = config
        self._database = database
        self._miner = None

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
    def textminer(self):
        if self._miner is None:
            hostname = self._config.get('textmining', 'host')
            directory = self._config.get('textmining', 'directory')
            script = self._config.get('textmining', 'grep_script')
            self._miner = TextMiner(hostname, directory, script)
        return self._miner

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

    def grep_fulltexts(self, references):
        self.textminer.open()
        regex = self._config.get('textmining', 'grep_regex')
        res = self.textminer.grep(regex, references)
        self.textminer.close()

        return res


class TextMiner(object):
    """A helper object to grep the fulltext archive."""

    def __init__(self, hostname, directory, script):
        self._hostname = hostname
        self._directory = directory
        self._script = script
        self._ready = False

    def open(self):
        if not self._ready:
            subprocess.run(['scp', self._script, f'{self._hostname}:grep-svm.sh'])
            self._ready = True

    def close(self):
        if self._ready:
            subprocess.run(['ssh', self._hostname, 'rm grep-svm.sh'])
            self._ready = False

    def grep(self, regex, references):
        lines = [f'{a} {b}' for a, b in references]
        command = f"bash ./grep-svm.sh {self._directory} '{regex}'"
        ssh = subprocess.Popen(
            ['ssh', self._hostname, command],
            shell=False,
            text=True,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        output, _ = ssh.communicate('\n'.join(lines) + '\n')

        res = {}
        for line in output.split('\n'):
            line = line.strip()
            if '\t' not in line:
                continue
            fbrf, nmatches = line.split('\t')
            res[fbrf] = nmatches
        return res


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
@click.pass_obj
def findnew(obj, filename, output):
    """Find new FlyBase references that may be about scRNAseq."""

    if filename is None:
        filename = obj.flybase_table_file

    if not os.path.exists(filename):
        table = DataFrame(
            columns=[
                'FBrf',
                'PMID',
                'Citation',
                'Accessions',
                'Mentions',
                'Confirmed',
                'Organ/tissue',
                'Comments',
            ],
            dtype=str,
        )
    else:
        table = read_csv(filename, sep='\t', dtype=str)

    click.echo("Fetching dataset-flagged references...")
    fbrfs = obj.get_dataset_fbrfs()

    click.echo("Fetching additional data...")
    newrows = []
    with click.progressbar(fbrfs) as bar:
        for (fbrf,) in bar:
            if fbrf in table['FBrf'].values:
                continue

            pmid = obj.get_pmid_for_fbrf(fbrf)
            accessions = obj.get_dataset_accessions_for_fbrf(fbrf)
            citation = obj.get_citation_for_fbrf(fbrf)
            newrows.append(
                {
                    'FBrf': fbrf,
                    'PMID': pmid,
                    'Citation': citation,
                    'Accessions': accessions,
                }
            )

    click.echo(f"New dataset-flagged references: {len(newrows)}")
    if len(newrows) == 0:
        return

    table = concat([table, DataFrame(data=newrows)])

    click.echo("Querying the fulltext archive...")
    queries = table.loc[table['Mentions'].isna(), ['FBrf', 'PMID']].values
    mentions = obj.grep_fulltexts(queries)
    for fbrf, nmatch in mentions.items():
        table.loc[table['FBrf'] == fbrf, 'Mentions'] = nmatch

    if output is None:
        output = filename
        copyfile(output, f'{output}.bak')
    table.to_csv(output, index=False, sep='\t')


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
