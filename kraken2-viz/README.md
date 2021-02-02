# From kraken2 and bracken reports to visualisation

## 1. Merge report files

In order to visualise multiple taxonomy profiles obtained with kraken2 (and bracken) we need to merge the per-sample report files. 
Check out the different options with

    ./combine_kreports_modified.py -h

To merge kraken2 report files (using all *.report files):

    ./combine_kreports_modified.py -r *.report -o kraken2-merged.report

To merge bracken report files (using all *.report files):

    ./combine_kreports_modified.py -r *.report -o bracken-merged.report --bracken

To merge kraken2 report files (using all *.report files) and in addition to merged report also obtain a separate report for each taxonomic level ('kingdom':'K', 'superkingdom':'D','phylum':'P','class':'C','order':'O','family':'F','genus':'G','species':'S'):

    ./combine_kreports_modified.py -r *.report -o kraken2-merged.report --single-tax-level

To merge bracken report files (using all *.report files) and in addition to merged report also obtain a separate report for each taxonomic level ('kingdom':'K', 'superkingdom':'D','phylum':'P','class':'C','order':'O','family':'F','genus':'G','species':'S'):

    ./combine_kreports_modified.py -r *.report -o bracken-merged.report --bracken --single-tax-level


## 2. Visualise profiles

In order to visualise the abundance profiles we use kraken2 input to report the number of unclassified reads and the bracken report to report taxonomic abundances.

Input files are *kraken2-merged.report* and *bracken-merged.report* from **step 1**. Alternatively you can provide your own abundance tables in the following tab-separated formats (any lines starting with '#' will be ignored):

*kraken2-merged.report*

    Sample1 Sample2 Sample3 lvl_type	taxid	name
    1873876	1096341	1851559 U	0	unclassified
    2133580	524258	1627158 R	1	root
    [...]

*bracken-merged.report*

    Sample1 Sample2 Sample3 lvl_type	taxid	name
    2133580	524258	1627158 R	1	root
    [...]

Run visualisation R script replace *bracken-merged.report* and *kraken2-merged.report* with your own inoput files and run the following command

    Rscript render.R --inbracken bracken-merged.report --inkraken kraken2-merged.report -o ./
