#!/usr/bin/env python
"""
Make the antismash output actually useful
"""

__author__ = "Matt Olm"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import glob
import os
import textwrap

import pandas as pd

from Bio import SeqIO

def main(args):
    """ Main entry point of the app """
    if not os.path.exists(args.out):
        os.makedirs(args.out)

    if not os.path.exists(args.smash):
        print('{0} is not an antismash directory- quitting'.format(args.smash))

    parse_antismash_folder(args.smash, args.out)

def parse_antismash_folder(folder, outdir, write=True):
    '''
    Take an antismash folder and parse it to make useful output
    '''
    gene_clusters_loc = os.path.join(folder, 'geneclusters.txt')
    ass_name = os.path.join(folder).split('/')[-2]

    # Parse geneclusters.txt
    if not os.path.isfile(gene_clusters_loc):
        print("{0} is empty!!!!! Returning nothing".format(\
            folder))
        return
    if os.stat(gene_clusters_loc).st_size == 0:
        print("{0} has an empty geneclusters- returning nothing".format(\
            folder))
        return
    gdb = load_geneclusters_txt(gene_clusters_loc, ass_name)

    # Load genbank file
    full_gb = glob.glob(os.path.join(folder, '*.final.gbk'))
    assert len(full_gb) == 1
    full_gb = full_gb[0]

    # Save geneclusters info
    gdb.to_csv(os.path.join(outdir, '{0}.geneclusters.txt'.format(ass_name)),\
              index=False, sep='\t')

    # Write fasta files
    Sdb = write_fasta_files(gdb, full_gb, outdir, write=True)
    
    return gdb, Sdb

def load_geneclusters_txt(file, ass_name):
    '''
    parse geneclusters file and return dataframe
    '''
    gdb = pd.read_table(file, header=None)
    cols = ['contig_full', 'scaffold_full', 'cluster_type', \
                   'cluster_genes', 'cluster_genes_again']
    gdb.columns = cols

    for c in ['cluster_genes', 'cluster_genes_again']:
        gdb[c] = [list(x.split(';')) for x in gdb[c]]

    gdb['contig'] = [x.split('_')[0] for x in gdb['contig_full']]
    gdb['scaffold'] = [x.split(' ')[0] for x in gdb['scaffold_full']]

    gdb.insert(0, 'cluster_number', range(1, len(gdb) + 1))

    # Make sure there are no double clusters
    assert len(gdb[gdb['cluster_number'].duplicated()]) == 0

    # Make sure cluster types are the same
    assert len(gdb[gdb['cluster_genes'] != gdb['cluster_genes_again']]) == 0

    # Make sure all clusters are continuous
    for l in gdb['cluster_genes'].tolist():
        it = (int(x.split('_')[1]) for x in l if x.startswith('ctg'))
        first = next(it)
        assert all(a == b for a, b in enumerate(it, first + 1))

    # Add assembly stuff
    gdb['assembly'] = ass_name
    gdb['cluster'] = ["cluster{1:03}".format(x, y) for x, y in zip(\
                        gdb['assembly'], gdb['cluster_number'])]

    del gdb['cluster_genes_again']
    cols.remove('cluster_genes_again')
    return gdb[['cluster', 'contig', 'scaffold', 'assembly', 'cluster_number'] \
               + cols]

def write_fasta_files(xdb, gb_file, outdir, write=True):
    Sdb = pd.DataFrame()
    gdb = xdb.copy()
    for seq_record in SeqIO.parse(gb_file, "genbank"):
        if seq_record.id in list(gdb['contig_full']):
            # For every cluster on this scaffold:
            for cluster, db in gdb[gdb['contig_full'] == seq_record.id].groupby('cluster'):

                # Set up .fasta files
                fasta_base = os.path.join(outdir, "{0}_{1}".format(\
                        cluster, db['assembly'].tolist()[0]))
                fna_handle = open(fasta_base + '.fna', 'w')
                faa_handle = open(fasta_base + '.faa', 'w')

                # Set up info table
                table = {'cluster':[], 'gene':[], 'scaffold':[], 'location':[], 'info':[]}

                # Figure out genes to find
                to_find = db['cluster_genes'].tolist()[0]
                assert(len(to_find) > 0)

                # Print necessary genes to file
                for feature in seq_record.features:
                    if feature.type != "CDS":
                        continue

                    # This gene needs to be printed
                    if feature.qualifiers['locus_tag'][0] in to_find:
                        gene = feature.qualifiers['locus_tag'][0]
                        scaffold = db['scaffold'].tolist()[0]

                        # mark gene as found
                        to_find.remove(gene)

                        # make gene header
                        header = ">{0}__{1}__{2}".format(cluster, gene, scaffold)

                        # write nucleotide
                        fna_handle.write("{0}\n{1}\n".format(header, \
                            textwrap.fill(str(feature.location.extract(seq_record).seq), 80)))

                        # write amino acid
                        faa_handle.write("{0}\n{1}\n".format(header, \
                            textwrap.fill(str(feature.qualifiers['translation'][0]), 80)))

                        # store information about gene
                        table['cluster'].append(cluster)
                        table['gene'].append(gene)
                        table['location'].append(str(feature.location))
                        table['scaffold'].append(scaffold)
                        if 'sec_met' in feature.qualifiers:
                            table['info'].append(feature.qualifiers['sec_met'])
                        else:
                            table['info'].append('')

                # Close handles
                faa_handle.close()
                fna_handle.close()

                # Make sure all genes were found
                assert(len(to_find) == 0)

                # Save datatable
                sdb = pd.DataFrame(table)
                sdb.to_csv(fasta_base + '.info.tsv', sep='\t', index=False)

                # Append datatable
                Sdb = pd.concat([Sdb, sdb])

    return Sdb

if __name__ == "__main__":
    """ From the output of antismash, this will make parse it out in a nice way"""
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("smash", help="location of antismash folder")
    parser.add_argument("out", help="location of output folder to store results")

    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)