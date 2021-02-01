#!/usr/bin/python3

import pandas as pd
import os
import sys
import argparse

def main(argv):
 
 # define user input
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('-i', '--inputdir', default='.', help='set the input directory containing emapper output *.annotations')
    parser.add_argument('-s', '--single', action='store_true', help='produce additional output tsv for KO and COG counts of each single dataset')

    args = parser.parse_args()

    # set wd
    workingdir = os.chdir(args.inputdir)
    # list all files in wd and save as entries
    entries = os.listdir(workingdir)

    #######

    # create combi dataframe with column name used for merging in loop
    combi = pd.DataFrame(columns = ['KEGG_ko'])

    # loop over all files in wd that end with annotations
    for emapfile in entries:
        if emapfile.endswith('.annotations'):
            print('processing dataset ==> ',emapfile)
            # open emapper.annotations file (skip first 3 rows)
            sample = pd.read_csv(emapfile, sep='\t', skiprows=3)
            # select only query_name and KEGG_ko column
            select = sample.loc[:,['#query_name','KEGG_ko']]
            # drop all lines where KEGG_ko is NaN
            select = select[select['KEGG_ko'].notna()]
            select = select.replace(to_replace ='ko:', value = '', regex = True) 

            # iterate over each row to split KEGG_ko column by ',' and concatenate into new file 
            splitted = pd.concat([pd.Series(row['#query_name'], row['KEGG_ko'].split(','))              
                for _, row in select.iterrows()]).reset_index()
            # rename columns
            splitted.columns = ['KEGG_ko', 'query_name']
            
            # calculate counts for every KO 
            ko_counts = pd.DataFrame(splitted['KEGG_ko'].value_counts()).rename_axis('KEGG_ko', axis=0)
            ko_counts.columns = [str(emapfile).replace('.emapper.annotations','')]

            # loop over merging the KO counts of all input files into the dataframe combi
            combi = pd.merge(combi, ko_counts, how='outer', on='KEGG_ko')
            
            # if set save KOs and counts to separate files
            if args.single is True:
                print('saving KO counts for single dataset: ', str(emapfile).replace('.emapper.annotations',''))
                ko_counts.to_csv(emapfile + '.KEGG-ko-counts.csv')

        else:
            continue

    # rename column to comply with microbiomeanalyst input
    combi.rename(columns={'KEGG_ko': '#NAME'}, inplace=True)
    # replace missing values with 0
    combi = combi.fillna(0)
    # save KO counts
    print('\nsaving merged KO counts to KEGG_ko_counts.csv\n')
    combi.to_csv('KEGG_ko_counts.csv', index=False)

    ########

    # create combi dataframe with column name used for merging in loop
    combicog = pd.DataFrame(columns = ['COG'])

    # loop over all files in wd that end with annotations
    for emapfile in entries:
        if emapfile.endswith('.annotations'):
            print('processing dataset ==> ',emapfile)
            # open emapper.annotations file (skip first 3 rows)
            sample = pd.read_csv(emapfile, sep='\t', skiprows=3)
            #print(sample.head())
            # select only query_name and COG Functional cat. column
            select = sample.loc[:,['#query_name','COG Functional cat.']]
            # drop all lines where COG Functional cat. is NaN
            select = select[select['COG Functional cat.'].notna()]

            # iterate over each row to split COG Functional cat. column by '' and concatenate into new file 
            splitted = pd.concat([pd.Series(row['#query_name'], list(row['COG Functional cat.']))              
                for _, row in select.iterrows()]).reset_index()
            # rename columns
            splitted.columns = ['COG Functional cat.', 'query_name']
            
            # calculate counts for every COG 
            cog_counts = pd.DataFrame(splitted['COG Functional cat.'].value_counts()).rename_axis('COG', axis=0)
            cog_counts.columns = [str(emapfile).replace('.emapper.annotations','')]

            # loop over merging the KO counts of all input files into the dataframe combi
            combicog = pd.merge(combicog, cog_counts, how='outer', on='COG')
      
            # if set save COGs and counts to separate files
            if args.single is True:
                print('saving COG counts for single dataset: ', str(emapfile).replace('.emapper.annotations',''))
                cog_counts.to_csv(emapfile + '.COG-counts.csv')

        else:
            continue

    # rename column to comply with microbiomeanalyst input
    combicog.rename(columns={'COG': '#NAME'}, inplace=True)
    # replace missing values with 0
    combicog = combicog.fillna(0)
    # save COG counts
    print('\nsaving merged COG counts to COG_counts.csv\n')
    combicog.to_csv('COG_counts.csv', index=False)

    #####

    # create combi dataframe with column name used for merging in loop
    combi_pw = pd.DataFrame(columns = ['KEGG_Pathway'])

    # loop over all files in wd that end with annotations
    for emapfile in entries:
        if emapfile.endswith('.annotations'):
            print('processing dataset ==> ',emapfile)
            # open emapper.annotations file (skip first 3 rows)
            sample = pd.read_csv(emapfile, sep='\t', skiprows=3)
            # select only query_name and KEGG_ko column
            select = sample.loc[:,['#query_name','KEGG_Pathway']]
            # drop all lines where KEGG_ko is NaN
            select = select[select['KEGG_Pathway'].notna()]
            select = select.replace(to_replace =',map.*', value = '', regex = True) 

            # iterate over each row to split KEGG_ko column by ',' and concatenate into new file 
            splitted = pd.concat([pd.Series(row['#query_name'], row['KEGG_Pathway'].split(','))              
                for _, row in select.iterrows()]).reset_index()
            # rename columns
            splitted.columns = ['KEGG_Pathway', 'query_name']
            
            # calculate counts for every KO 
            pw_counts = pd.DataFrame(splitted['KEGG_Pathway'].value_counts()).rename_axis('KEGG_Pathway', axis=0)
            pw_counts.columns = [str(emapfile).replace('.emapper.annotations','')]

            # loop over merging the KO counts of all input files into the dataframe combi
            combi_pw = pd.merge(combi_pw, pw_counts, how='outer', on='KEGG_Pathway')
            
            # if set save KOs and counts to separate files
            if args.single is True:
                print('saving KEGG Pathway counts for single dataset: ', str(emapfile).replace('.emapper.annotations',''))
                pw_counts.to_csv(emapfile + '.KEGG-Pathway-counts.csv')

        else:
            continue

    # rename column to comply with microbiomeanalyst input
    combi_pw.rename(columns={'KEGG_Pathway': '#NAME'}, inplace=True)
    # replace missing values with 0
    combi_pw = combi_pw.fillna(0)
    # save KO counts
    print('\nsaving merged KEGG Pathway counts to KEGG_Pathway_counts.csv\n')
    combi_pw.to_csv('KEGG_Pathway_counts.csv', index=False)

    #####

    # create combi dataframe with column name used for merging in loop
    combi_cazy = pd.DataFrame(columns = ['CAZy'])

    # loop over all files in wd that end with annotations
    for emapfile in entries:
        if emapfile.endswith('.annotations'):
            print('processing dataset ==> ',emapfile)
            # open emapper.annotations file (skip first 3 rows)
            sample = pd.read_csv(emapfile, sep='\t', skiprows=3)
            # select only query_name and KEGG_ko column
            select = sample.loc[:,['#query_name','CAZy']]
            # drop all lines where KEGG_ko is NaN
            select = select[select['CAZy'].notna()]
            #select = select.replace(to_replace ='', value = '', regex = True) 

            # iterate over each row to split KEGG_ko column by ',' and concatenate into new file 
            splitted = pd.concat([pd.Series(row['#query_name'], row['CAZy'].split(','))              
                for _, row in select.iterrows()]).reset_index()
            # rename columns
            splitted.columns = ['CAZy', 'query_name']
            
            # calculate counts for every KO 
            cazy_counts = pd.DataFrame(splitted['CAZy'].value_counts()).rename_axis('CAZy', axis=0)
            cazy_counts.columns = [str(emapfile).replace('.emapper.annotations','')]

            # loop over merging the KO counts of all input files into the dataframe combi
            combi_cazy = pd.merge(combi_cazy, cazy_counts, how='outer', on='CAZy')
            
            # if set save KOs and counts to separate files
            if args.single is True:
                print('saving CAZy counts for single dataset: ', str(emapfile).replace('.emapper.annotations',''))
                cazy_counts.to_csv(emapfile + '.CAZy-counts.csv')

        else:
            continue

    # rename column to comply with microbiomeanalyst input
    combi_cazy.rename(columns={'CAZy': '#NAME'}, inplace=True)
    # replace missing values with 0
    combi_cazy = combi_cazy.fillna(0)
    # save KO counts
    print('\nsaving merged CAZy counts to CAZy_counts.csv\n')
    combi_cazy.to_csv('CAZy_counts.csv', index=False)

    # bye bye
    print('\nAll done - bye bye!')


if __name__ == "__main__": main(sys.argv[1:])
