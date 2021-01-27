library("rmarkdown")
library("argparser")

# Parser
par <- arg_parser("render kraken2 visualisation - report")
# Arguments
par <- add_argument(par, "--inbraken", help="merged braken input file")
par <- add_argument(par, "--inkraken", help="merged kraken input file")
par <- add_argument(par, "--outdir", help="output directory")
par <- add_argument(par, "--topspecies", help="the number of top species to be displayes", default=10)

# Parse the command line arguments
argv <- parse_args(par)

if (pandoc_available())
  cat("pandoc", as.character(pandoc_version()), "is available!\n")
Sys.setenv(RSTUDIO_PANDOC="/Applications/RStudio.app/Contents/MacOS/pandoc")

rmarkdown::render("kraken2-viz.R", params=list(inbraken=argv$inbraken, inkraken=argv$inkraken, outdir=argv$outdir, topspecies=argv$topspecies))
