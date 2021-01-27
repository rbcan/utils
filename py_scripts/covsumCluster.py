#!/usr/bin/env python3

import pandas as pd
import sys
import argparse
import datetime

def main(argv):

    # define user input
    parser = argparse.ArgumentParser(description='get cluster coverage sums of all genes in a cluster')
    parser.add_argument('-c', '--clustfile', help='file containing the gene clusters. assumes cluster ID in column a [default 2] and gene ID in column b [default 9]')
    parser.add_argument('-a', '--clustIDcolumn', type=int, default=2, help='column with cluster ID in file containing the gene clusters [default 2]')
    parser.add_argument('-b', '--geneIDcolumn', type=int, default=9, help='column with gene ID in file containing the gene clusters [default 9]')
    parser.add_argument('-g', '--genecov', help='file containing each gene (column d [default 1]) and its read coverage normalized by gene length (column e [default 2])')
    parser.add_argument('-d', '--genecolumn', type=int, default=1, help='column with gene ID in file containing coverage per gene [default 1]')
    parser.add_argument('-e', '--genecovcolumn', type=int, default=2, help='column with read coverage normalized by gene length [default 2]')   

    args = parser.parse_args()

    A = args.clustIDcolumn - 1
    B = args.geneIDcolumn - 1
    D = args.genecolumn - 1 
    E = args.genecovcolumn - 1 

    now = datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")
    print(now," Summing read coverage for genes within a gene cluster")

    if A < B:
        columnnames=['clusterID','geneID']
    else:
        columnnames=['geneID','clusterID']

    # import cluster file
    input = pd.read_csv(args.clustfile, 
        sep='\s|t', 
        usecols=[A,B],
        index_col=False, 
        header=None, 
        names=columnnames,
        engine='python',
        comment='#').drop_duplicates()

    # import gene coverage file
    cov = pd.read_csv(args.genecov, 
        sep='\t', 
        usecols=[D,E], 
        index_col=False, 
        header=None, 
        names=['geneID','coverage'])

    # combine cluster file and coverage table
    combi = pd.merge(input, cov, how='outer', on='geneID')

    # sum coverage within each gene cluster
    cov_df = combi.groupby(['clusterID']).agg(lambda x: x.tolist())
    cov_df['sum_coverage'] = combi.groupby(['clusterID']).sum()

    # subset to retriev summed coverage only
    cov_df_sumonly = cov_df['sum_coverage']

    # export files to csv
    cov_df.to_csv('cluster_coverage.csv', index=True, sep="\t")
    cov_df_sumonly.to_csv('cluster_only-cov.csv', index=True, sep="\t")

    # bye bye
    print('\nAll done - bye bye!')

if __name__ == "__main__": main(sys.argv[1:])