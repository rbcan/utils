library("rmarkdown")
library("argparser")
library("tools")

# Parser
par <- arg_parser("render kraken2 visualisation - report")
# Arguments
par <- add_argument(par, "--inbracken", help="merged bracken input file")
par <- add_argument(par, "--inkraken", help="merged kraken input file")
par <- add_argument(par, "--topspecies", help="the number of top species to be displayes", default=10)
par <- add_argument(par, "--outdir", help="output directory", default="./")

# Parse the command line arguments
argv <- parse_args(par)

if (pandoc_available())
  cat("pandoc", as.character(pandoc_version()), "is available!\n")
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

inbracken_file<-file_path_as_absolute(argv$inbracken)
inbracken_file
inkraken_file<-file_path_as_absolute(argv$inkraken)
inkraken_file
xxx<-file_path_as_absolute(argv$outdir)
xxx

rmarkdown::render("../kraken2-viz_modified.R", params=list(inbracken=inbracken_file, inkraken=inkraken_file, topspecies=argv$topspecies, outdir=xxx), output_dir=argv$outdir)
