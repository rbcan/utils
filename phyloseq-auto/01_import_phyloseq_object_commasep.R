### CREATE PHYLOSEQ OBJECT ###

# FIRST THINGS FIRST #

# load required packages
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}

list.of.bioc <- c("phyloseq", "ape", "gtools", "plyr", "dplyr","tibble")
new.packages <- list.of.bioc[!(list.of.bioc %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

library("phyloseq")
library("gtools")
library("plyr")
library("dplyr")
library("argparser")
library("tibble")
library("ape")

# Create a parser for commandline arguments
par <- arg_parser("relative abundances features/OTUs")

# Add command line arguments
par <- add_argument(par, "--outdir", help="output directory", default="./")
par <- add_argument(par, "--feature", help="comma-separated feature/OTU table")
par <- add_argument(par, "--tax", help="comma-separated taxonomy file including header")
par <- add_argument(par, "--tree", help="tree in newick format (optional)", default="none")
par <- add_argument(par, "--meta", help="comma-separated metadata file (optional)")
par <- add_argument(par, "--sampleID", help="Column number (e.g. 1) with the sample ID", default=1)
par <- add_argument(par, "--mcat", help="comma-separated (e.g. Time,Treatment) column headers of metadata to include in plots", default="NoMetadata")

# Parse the command line arguments
argv <- parse_args(par)
MCAT = unlist(strsplit(argv$mcat,","))

#######################

##### DATA IMPORT #####

# import OTU data
in_otu = as.matrix(read.table(argv$feature, header=TRUE, sep=",", row.names = 1, comment.char=""))
head(in_otu)
class(in_otu)
ncol(in_otu)

# import taxonomy
in_tax = as.matrix(read.table(argv$tax, header=TRUE, sep=",", row.names = 1, fill=TRUE, na.strings=c("","NA","d__","p__","c__","o__","f__","g__","s__"), comment.char="")) %>%
  as.matrix()
head(in_tax)
class(in_tax)
ncol(in_tax)

# import tree
if (argv$tree == "none") {
  print("no tree file provided - proceeding without tree")
} else {
in_tree <- read.tree(file = argv$tree)
class(in_tree)
}

# import metadata
SAMPLEID = argv$sampleID # assign sample name column

if (is.na(argv$meta)){        
  print("no metadata provided - proceeding without")
  SampleName<-colnames(in_otu)
  metaIn<-data.frame(SampleName)
  SAMPLEID = 1 # overwrite sample name column if no metadata are provided
} else {
  metaIn <- argv$meta
  metaIn = read.table(metaIn, header=TRUE, sep=",", comment.char="") # import metadata table
} # check if metadata file was provided 

metaIn_tibble = tibble(metaIn) # create metadata tibble
metaIn_tibble = metaIn_tibble %>% filter( !grepl("^#",metaIn_tibble[[1]]))

if (MCAT == "NoMetadata") { 
  print("No metadata provided - proceed without")
  metaIn_tibble <- metaIn_tibble %>%
    mutate(NoMetadata="0")
} # if no metadata 'categories' are specified, then metadata only contain sample names and NoMetadata column for downstream plotting purposes

metaIn_tibble$rName = metaIn_tibble[[SAMPLEID]] # get sample name column 
metaIn_tibble = metaIn_tibble %>% column_to_rownames('rName') # copy sample names into rownames
metaIn_tibble = metaIn_tibble %>% rename(sampleIDs = all_of(SAMPLEID)) # rename sample name column to 'sampleIDs'
metaIn_tibble = metaIn_tibble %>% select(sampleIDs, all_of(MCAT)) # select metadata categories to be used in phyloseq object
in_metad = sample_data(metaIn_tibble)
class(in_metad)
head(in_metad)

my_OTU = otu_table(in_otu, taxa_are_rows = TRUE)
my_TAX = tax_table(in_tax)
head(my_OTU)
head(my_TAX)

# combine to create phyloseq object
if (argv$tree == "none") {
  my_physeq = phyloseq(my_OTU, my_TAX, in_metad)
} else {
  my_physeq = phyloseq(my_OTU, my_TAX, in_metad, in_tree)
}
my_physeq

# set working directory
dir.create(file.path(argv$outdir), showWarnings = FALSE)
setwd(file.path(argv$outdir))

# save phyloseq object as rds
saveRDS(my_physeq, file = "my_physeq.rds")
