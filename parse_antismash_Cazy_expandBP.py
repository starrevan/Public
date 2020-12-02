#!/usr/bin/env python
"""
Combine output from antismash_html_scrape.py and dbCan2 output
"""

__author__ = "Evan"
__version__ = "0.0.1"
__license__ = "MIT"

import csv
import pandas as pd
import argparse
import numpy as np



def read_csv(csvfile):
    print('read_csv(): type(csvfile)) = {}'.format(csvfile))
    print('')

    foo_df = pd.read_csv(csvfile)

    return foo_df

def main():
    parser = argparse.ArgumentParser(description = "Parse ouput of antismash_html_scrape.py and run_dbcan.py")
    parser.add_argument('-d', '--debug', help='Debugging output', action='store_true')
    parser.add_argument('Antismashcsv', type=argparse.FileType('r'), help='Input Antismash csv file')
    parser.add_argument('Cazycsv', type=argparse.FileType('r'), help='Input Cazy gff file')
    parser.add_argument("--version", action="version", version="%(prog)s (version {version})".format(version=__version__))
    parser.add_argument('-e', '--Expand', type=int, help='Optionally expand cazy search before and after anismash gene cluster (bp)', required=False)

    args = parser.parse_args()
    #read in files as dataframes
    anti_df = pd.read_csv(args.Antismashcsv)
    cazy_df = pd.read_csv(args.Cazycsv, sep='\t', header=None)

    #remove commas in 'start' and 'stop' from anti_df
    anti_df['start'] = anti_df['start'].str.replace(r'\D', '')
    anti_df['stop'] = anti_df['stop'].str.replace(r'\D', '')
    
    #rename cazy_df headers
    cazy_df = cazy_df.rename({0: 'sequence', 1: 'source', 2: 'feature', 3: 'start', 4: 'end', 5: 'score', 6: 'strand', 7: 'phase', 8: 'attributes'}, axis='columns')
    #keep only the CAZyme, TC, STP coding sequences
    cazy_df = cazy_df[cazy_df['feature'] != 'CDS'] 

    #Expand both + and - search by -e 
    expand = args.Expand 
    if expand is None: #set to 0 if not filled out
        expand = 0

    #combine antismash and cazy output
    #find cazy genes where 'start' is between the 'start' and 'stop' of the antismash output
    a = cazy_df.start.values
    bh = anti_df.stop.values.astype(np.int)+expand
    bl = anti_df.start.values.astype(np.int)-expand
    bl[bl < 0] = 0
    
    i, j = np.where((a[:, None] >= bl) & (a[:, None] <= bh))

    data = pd.DataFrame(
    np.column_stack([cazy_df.values[i], anti_df.values[j]]), 
    columns=cazy_df.columns.append(anti_df.columns))
    
    data = data[data.sequence == data.scaffold]
    data.to_csv('antiCazy_expand.csv', index=False)

#test output

    print(data)


if __name__ == '__main__':
    main()