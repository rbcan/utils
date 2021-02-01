#!/usr/bin/env python
################################################################
# combine_kreports_modified.py takes multiple kraken-style kraken2 
# or braken reports and combines them into a single report file
#
# MODIFIED AFTER combine_kreports.py from KrakenTools by
# Jennifer Lu under this link:
# https://github.com/jenniferlu717/KrakenTools/blob/master/combine_kreports.py
#
# modifications were done by by Rebecca Ansorge resulting in this script
# https://github.com/rbcan/utils/tree/master/kraken2-viz/combined_kreports_modified.py

#################################################################
# Rebecca Ansorge, rebecca.ansorge@quadram.ac.uk
# Modified: 01/02/2021
#
# This program reads in multiple Kraken2 or Braken report files and generates
# a combined Kraken report with columns for read counts 
#
# Parameters:
#   -h, --help................show help message.
#   -r X, --report-file X.....all input kraken reports (separated by spaces)
#   -o X, --output X..........output kraken report filename
#   --sample-names............sample names for each kraken report (separated by spaces)
#                             [if none are given, each sample named from input report file name]
#   --braken..................Set this for braken input - omits representation of unclassified fraction
#   --single-tax-level........Set this to create single tables for each major taxonomic level 
#                             in addition to overall merged report
# Each Input report file format (tab-delimited)
#   - percentage of total reads
#   - number of reads (including reads within subtree)
#   - number of reads (only at this level)
#   - taxonomic classification level (U, D, P, C, O, F, G, S,...etc)
#   - NCBI taxonomic ID
#   - name of level
# Output file format (tab-delimited)
#   - S1_all_reads, S2_all_reads, ...etc.
#   - taxonomic classification level (U, D, P, C, O, F, G, S,...etc)
#   - NCBI taxonomic ID
#   - name of level
# Methods 
#   - main
#   - process_kraken_report
####################################################################
import os, sys, argparse
import operator
from time import gmtime 
from time import strftime 

#Tree Class 
#usage: tree node used in constructing a taxonomy tree
#   including only the taxonomy levels and genomes identified in the Kraken report
class Tree(object):
    'Tree node.'
    def __init__(self, name, taxid, level_num, level_id, all_reads, lvl_reads, children=None, parent=None):
        self.name = name
        self.taxid = taxid
        self.level_num = level_num
        self.level_id = level_id
        self.tot_all = all_reads
        self.tot_lvl = lvl_reads
        self.all_reads = {}
        self.lvl_reads = {}
        self.children = []
        self.parent = parent
        if children is not None:
            for child in children:
                self.add_child(child)
    def add_child(self,node):
        assert isinstance(node,Tree)
        self.children.append(node)
    def add_reads(self, sample, all_reads, lvl_reads):
        self.all_reads[sample] = all_reads
        self.lvl_reads[sample] = lvl_reads
        self.tot_all += all_reads
        self.tot_lvl += lvl_reads
    def __lt__(self,other):
        return self.tot_all < other.tot_all
         
####################################################################
# process_kraken_report
# usage: parses a single line in the kraken report and extracts relevant information
# input: kraken report file with the following tab delimited lines 
#   - percent of total reads   
#   - number of reads (including at lower levels)
#   - number of reads (only at this level)
#   - taxonomy classification of level 
#       (U, - (root), - (cellular org), D, P, C, O, F, G, S) 
#   - taxonomy ID (0 = unclassified, 1 = root, 2 = Bacteria...etc)
#   - spaces + name 
# returns:
#   - classification/genome name
#   - taxonomy ID for this classification
#   - level for this classification (number)
#   - level name (U, -, D, P, C, O, F, G, S)
#   - all reads classified at this level and below in the tree
#   - reads classified only at this level
def process_kraken_report(curr_str):
    split_str = curr_str.strip().split('\t')
    if len(split_str) < 5:
        return []
    try:
        int(split_str[1])
    except ValueError:
        return []
    #Extract relevant information
    all_reads =  int(split_str[1])
    level_reads = int(split_str[2])
    level_type = split_str[3]
    taxid = split_str[4] 
    #Get name and spaces
    spaces = 0
    name = split_str[-1]
    for char in name:
        if char == ' ':
            name = name[1:]
            spaces += 1 
        else:
            break 
    #Determine which level based on number of spaces
    level_num = int(spaces/2)
    return [name, taxid, level_num, level_type, all_reads, level_reads]
    
####################################################################
# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--report-file','--report-files',
        '--report','--reports', required=True,dest='r_files',nargs='+',
        help='Input kraken report files to combine (separate by spaces)') 
    parser.add_argument('-o','--output', required=True,dest='output',
        help='Output kraken report file with combined information')
    parser.add_argument('--sample-names',required=False,nargs='+',
        dest='s_names',default=[],help='Sample names to use as headers in the new report')
    parser.add_argument('--braken',required=False,dest='braken',
        action='store_true', default=False,
        help='Set this for braken input - omits representation of unclassified fraction')
    parser.add_argument('--single-tax-level',required=False,dest='singletaxa',
        action='store_true', default=False,
        help='Set this to create single tables for each taxonomic level')
    args=parser.parse_args()
    

    # Initialize combined values 
    main_lvls = ['U','R','D','K','P','C','O','F','G','S']
    map_lvls = {'kingdom':'K', 'superkingdom':'D','phylum':'P','class':'C','order':'O','family':'F','genus':'G','species':'S'}
    count_samples = 0
    num_samples = len(args.r_files)
    if args.s_names:
        sample_names = args.s_names
    else:
        sample_names = []
        for r_file in args.r_files:
            name = os.path.basename(r_file)
            name = os.path.splitext(name)[0]
            sample_names.append(name)
    root_node = -1 
    prev_node = -1
    curr_node = -1
    u_reads = {0:0} 
    taxid2node = {}

    # Check input values 
    if len(sample_names) > 0 and len(sample_names) != num_samples: 
        sys.stderr.write("Number of sample names provided does not match number of reports\n")
        sys.exit(1)
    # Map names
    id2names = {} 
    id2files = {} 
    if len(sample_names) == 0:
        for i in range(num_samples):
            id2names[i+1] = "S" + str(i+1)
            id2files[i+1] = ""
    else:
        for i in range(num_samples):
            id2names[i+1] = sample_names[i] 
            id2files[i+1] = ""
    
    #################################################
    # STEP 1: READ IN REPORTS
    # Iterate through reports and make combined tree! 
    sys.stdout.write(">>STEP 1: READING REPORTS\n")
    sys.stdout.write("\t%i/%i samples processed" % (count_samples, num_samples))
    sys.stdout.flush()
    for r_file in args.r_files:
        count_samples += 1 
        sys.stdout.write("\r\t%i/%i samples processed" % (count_samples, num_samples))
        sys.stdout.flush()
        id2files[count_samples] = r_file
        # Open File 
        curr_file = open(r_file,'r')
        for line in curr_file: 
            report_vals = process_kraken_report(line)
            if len(report_vals) < 5:
                continue
            [name, taxid, level_num, level_id, all_reads, level_reads] = report_vals
            if level_id in map_lvls:
                level_id = map_lvls[level_id]
            # Unclassified 
            if level_id == 'U' or taxid == '0':
                u_reads[0] += all_reads
                u_reads[count_samples] = all_reads 
                continue
            # Tree Root 
            if taxid == '1': 
                if count_samples == 1:
                    root_node = Tree(name, taxid, level_num, 'R', 0,0)
                    taxid2node[taxid] = root_node 
                root_node.add_reads(count_samples, all_reads, level_reads) 
                prev_node = root_node
                continue 
            # Move to correct parent
            while level_num != (prev_node.level_num + 1):
                prev_node = prev_node.parent
            # IF NODE EXISTS 
            if taxid in taxid2node: 
                taxid2node[taxid].add_reads(count_samples, all_reads, level_reads) 
                prev_node = taxid2node[taxid]
                continue 
            # OTHERWISE
            # Determine correct level ID
            if level_id == '-' or len(level_id)> 1:
                if prev_node.level_id in main_lvls:
                    level_id = prev_node.level_id + '1'
                else:
                    num = int(prev_node.level_id[-1]) + 1
                    level_id = prev_node.level_id[:-1] + str(num)
            # Add node to tree
            curr_node = Tree(name, taxid, level_num, level_id, 0, 0, None, prev_node)
            curr_node.add_reads(count_samples, all_reads, level_reads)
            taxid2node[taxid] = curr_node
            prev_node.add_child(curr_node)
            prev_node = curr_node 
        curr_file.close()

    sys.stdout.write("\r\t%i/%i samples processed\n" % (count_samples, num_samples))
    sys.stdout.flush()

    #################################################
    # STEP 2: SETUP OUTPUT FILE
    sys.stdout.write(">>STEP 2: WRITING NEW REPORT HEADERS\n")
    o_file = open(args.output,'w') 
    # Lines mapping sample ids to filenames
    o_file.write("#Number of Samples: %i\n" % num_samples) 
    for i in id2names:
        o_file.write("#")
        o_file.write("%s\t" % id2names[i])
        o_file.write("%s\n" % id2files[i])
        #Report columns
    for i in sample_names:
        o_file.write("%s\t" % i)
    o_file.write("lvl_type\ttaxid\tname\n")
    #################################################
    # STEP 3a: PRINT TREE
    if args.braken:
        print(">>NOTE: braken input (--braken) - ignore unclassified")
    sys.stdout.write(">>STEP 3a: PRINTING OVERALL REPORT\n")
    # Print line for unclassified reads
    if args.braken:
        pass
    else:
        for i in u_reads:
            if (i > 0): #i == 0 or 
                o_file.write("%i\t" % u_reads[i])
        o_file.write("U\t0\tunclassified\n")
    # Print for all remaining reads 
    all_nodes = [root_node]
    curr_node = -1
    prev_node = -1
    while len(all_nodes) > 0:
        # Remove node and insert children
        curr_node = all_nodes.pop()
        if len(curr_node.children) > 0:
            curr_node.children.sort()
            for node in curr_node.children:
                all_nodes.append(node)
        # Print information for this node 
        for i in range(num_samples):
            if (i+1) not in curr_node.all_reads: 
                o_file.write("0\t")
            else:
                o_file.write("%i\t" % curr_node.all_reads[i+1])
        o_file.write("%s\t" % curr_node.level_id)
        o_file.write("%s\t" % curr_node.taxid)
        o_file.write(" "*curr_node.level_num*2)
        o_file.write("%s\n" % curr_node.name)

    o_file.close() 

#########################################################################################
    # STEP 3b: PRINT SEPARATE REPORTS FOR EVERY TAXONOMIC LEVEL
    if args.singletaxa:
        for osplit in main_lvls:
            if osplit in ["U","R"]:
                pass
            else:
                osplit_name = args.output + "-" + osplit
                osplit_file = open(osplit_name,'w') 
                # Lines mapping sample ids to filenames

                for i in sample_names:
                    osplit_file.write("%s\t" % i)
                osplit_file.write("lvl_type\ttaxid\tname\n")

                #################################################
                # STEP 3: PRINT TREE
                sys.stdout.write(">>STEP 3b: PRINTING SEPARATE REPORTS PER TAXONMIC LEVEL \n")
                # Print line for unclassified reads
                if args.braken:
                    pass
                else:
                    for i in u_reads:
                        if (i > 0): #i == 0 or 
                            osplit_file.write("%i\t" % u_reads[i])
                    osplit_file.write("U\t0\tunclassified\n")
                # Print for all remaining reads 
                all_nodes = [root_node]
                curr_node = -1
                prev_node = -1
                while len(all_nodes) > 0:
                    # Remove node and insert children
                    curr_node = all_nodes.pop()
                    if len(curr_node.children) > 0:
                        curr_node.children.sort()
                        for node in curr_node.children:
                            all_nodes.append(node)
                    # Print information for this node 
                        if curr_node.level_id == osplit:
                            for i in range(num_samples):
                                if (i+1) not in curr_node.all_reads: 
                                    osplit_file.write("0\t")
                                else:
                                    osplit_file.write("%i\t" % curr_node.all_reads[i+1])
                            osplit_file.write("%s\t" % curr_node.level_id)
                            osplit_file.write("%s\t" % curr_node.taxid)
                            osplit_file.write(" "*curr_node.level_num*2)
                            osplit_file.write("%s\n" % curr_node.name)

                osplit_file.close() 
        
####################################################################
if __name__ == "__main__":
    main()
