library("rmarkdown")
library("argparser")
library("tools")
# tidyverse ggdendro

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

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
inkraken_file<-file_path_as_absolute(argv$inkraken)
cat("Kraken merged report:  ", inkraken_file, "\n")
cat("Bracken merged report: ", inbracken_file, "\n")

# Create if not existing (warn if exists)
dir.create(file.path(argv$outdir))
outputPath<-file_path_as_absolute(argv$outdir)
cat("Output path:           ", outputPath, "\n")


vizScript  <- file.path(script.basename, "kraken2-viz.R")
rmarkdown::render(vizScript, params=list(inbracken=inbracken_file, inkraken=inkraken_file, topspecies=argv$topspecies, outdir=outputPath), output_dir=argv$outdir)
