#!/usr/bin/env python3

import pandas as pd
import sys
import argparse
import datetime
import os

def main(argv):
    # define user input
    parser = argparse.ArgumentParser(description='get cluster coverage sums of all genes in a category of EggNOGmapper output')
    parser.add_argument('-c', '--clustfile', help='eggnog mapper output annotation file')
    parser.add_argument('-a', '--clustIDcolumn', type=int, required=True, nargs='+', help='column with cluster ID in file containing the gene clusters. Valid columns are 7: Gene Ontology terms, 8: EC number, 9: KEGG_ko, 10: KEGG_Pathway, 11: KEGG_Module, 12: KEGG_Reaction, 13: KEGG_rclass, 14: BRITE, 15: KEGG_TC, 16.: CAZy, 17: BiGG Reaction, 19.: eggNOG OGs, 21: COG Functional Category,')
    parser.add_argument('-b', '--geneIDcolumn', type=int, default=1, help='column with gene ID in eggNOG annotation file [default 1]')
    parser.add_argument('-g', '--genecov', help='file containing gene ID and its read coverage normalized by gene length')
    parser.add_argument('-d', '--genecolumn', type=int, default=1, help='column with gene ID in file containing coverage per gene [default 1]')
    parser.add_argument('-e', '--genecovcolumn', type=int, default=2, help='column with read coverage normalized by gene length [default 2]')   

    args = parser.parse_args()

    # overlap of two lists
    def intersection(lst1, lst2): 
        lst3 = [value for value in lst1 if value in lst2] 
        return lst3 
  
    list1 = [2,3,4,5,6,18,20] # exempt columns from emapper output
    intsect = intersection(list1, args.clustIDcolumn)
    if len(intsect) > 0:
        print(" ")
        print("Column(s)",intsect,"Not a valid emapper category to sum coverage - please choose only valid columns")
        print(" == EXITING == ")
        exit()

    for inlist in args.clustIDcolumn:

        A = inlist - 1
        B = args.geneIDcolumn - 1
        D = args.genecolumn - 1 
        E = args.genecovcolumn - 1 

        columnnames=["geneID","clusterID"] 
        # import cluster (emapper) file
        input = pd.read_csv(args.clustfile,
            sep='\t', 
            usecols=[B,A],
            index_col=False, 
            header=None, 
            names=columnnames,
            engine='python',
            comment='#').dropna().replace(to_replace =',map.*', value = '', regex = True) 

        # import headline from cluster (emapper) file to name category
        headline = pd.read_csv(args.clustfile, 
            sep='\t', 
            index_col=False, 
            header=None, 
            engine='python',
            skiprows=3,
            nrows=1)
        
        now = datetime.datetime.now().strftime("[%Y-%m-%d %H:%M:%S]")

        cat = headline[A].iloc[0]
        print(now," Summing read coverage for genes belonging to categories of",cat)

        # split categories into separate rows
        if args.clustIDcolumn == 21 : # for COG categories split after every 'letter'
            splitted = pd.concat([pd.Series(row['geneID'], list(row['clusterID']))              
                for _, row in input.iterrows()]).reset_index()
        else: # for all others split at the comma
            splitted = pd.concat([pd.Series(row['geneID'], row['clusterID'].split(','))              
                for _, row in input.iterrows()]).reset_index()

        # rename columns
        splitted = splitted.rename(columns={"index":"clusterID", 0:"geneID"})

        # import gene coverage table
        cov = pd.read_csv(args.genecov, 
            sep='\t', 
            usecols=[D,E], 
            index_col=False, 
            header=None, 
            names=['geneID','coverage'])

        # combine cluster file and coverage table
        combi = pd.merge(splitted, cov, how='outer', on='geneID')

        # sum coverage per gene cluster or category
        cov_df = combi.groupby(['clusterID']).agg(lambda x: x.tolist())
        cov_df['sum_coverage'] = combi.groupby(['clusterID']).sum()

        # subset to retriev summed coverage only
        cov_df_sumonly = cov_df['sum_coverage']

        # extort csv files
        name = os.path.basename(args.genecov)
        name = os.path.splitext(name)[0]
        cov_df.to_csv(name + '_' + cat + '_cov.csv', index=True, sep="\t")
        cov_df_sumonly.to_csv(name + '_' + cat + '_only-cov.csv', index=True, sep="\t")

    # bye bye
    print('\nAll done - bye bye!')

if __name__ == "__main__": main(sys.argv[1:])
