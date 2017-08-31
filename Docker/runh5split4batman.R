#!/usr/bin/Rscript

## import optparse library
suppressPackageStartupMessages(library("optparse"))

## specify our desired options in a list
## by default OptionParser will add an help option equivalent to
## make_option(c("-h", "--help"), action="store_true", default=FALSE,
## help="Show this help message and exit")
option_list <- list(
   make_option(c("-i", "--inputData"), 
               help=".h5 file as NMR raw data"),
   make_option(c("-t", "--interval"), 
               help="An integer as interval spectra range")
)

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

if(!("inputData" %in% names(opt))) {
  print("no input argument given!")
  print_help(parser)
  q(status = 1,save = "no")
}

## Run .h5 splitter for BATMAN
library(h5)
library(h5split4batman)
cat("start splitting....\n")
h5split4batman(BrukerZippedFile=opt$inputData, spec_interval=opt$interval)

