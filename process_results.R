#!/usr/bin/env Rscript
library("optparse")


# SPECIFY PATH -------------------------------------------------------------
path_to_repo = "/Users/d/git_for_cv/phase_analysis_public"

# Prepare Option List -----------------------------------------------------

option_list = list(
  make_option(c("-d","--working_dir"), type="character", default=NULL, 
              help="Absolute Path of Working Directory", metavar="character"),
  make_option(c("-n", "--nd_folder_path"), type="character", default=NULL, 
              help="Absolute Path to Folder that contains .nd file", metavar="character"),
  make_option(c("-p", "--pdef_path"), type="character", default=NULL, 
              help="Absolute Path to Plate Definition File", metavar="character"),
  make_option(c("-r", "--results_txt_path"), type="character", default=NULL, 
              help="Absolute Path to folder that contains result .txt files", metavar="character"),
  make_option(c("-w", "--wells"), type="integer", default=384, 
              help="Amount of wells in the imaging plate", metavar="number"),
  make_option(c("-c", "--column_to_measure"), type="character", default="_int_b5", 
              help="Which bin to use for Fluorescence Intensity. Default is Mean Intensity (_int_b5)", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=NULL, 
              help="Absolute Path to Output Folder", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#Set the working directory, if not specified use the current working directory
if(is.null(opt$working_dir)) {
  opt$working_dir <- getwd()
}


cat(
  "\n",
  "#######", "\n\n",
  "Looking for required files in.........................", "\t", opt$working_dir, "\n",
  "For custom input files use options. To see them run...",  "\t", "./process_results.R --help", "\n\n",
  "#######",
  "\n\n",
  sep = ""
)


cat("Continue? [y/n] ")
ANSWER <- readLines("stdin",n=1)
if (ANSWER != "y" & ANSWER != "Y") {stop("User interrupt")}


# Source Scripts ----------------------------------------------------------
#(Only after user confirms to continue)

source(paste(path_to_repo,"private_scripts/Multiwell.R",sep = "/"))
source(paste(path_to_repo,"private_scripts/general.R",sep = "/"))
source(paste(path_to_repo,"private_scripts/MicroscopeToolBox.R",sep = "/"))
source(paste(path_to_repo,"private_scripts/raw_data_tools.R",sep = "/"))



# Populate Option List ----------------------------------------------------
# Check whether the user specified any elements of the option list via
# command line options. If not, search for it with list.files()

#search for .nd file. save to opt$nd_folder_path
if (is.null(opt$nd_folder_path)) {
  opt$nd_folder_path <-
    list.files(
      path = opt$working_dir,
      pattern = ".nd$",
      full.names = TRUE,
      recursive = TRUE
    )
  if (length(opt$nd_folder_path) > 1)
    stop(
      "More than 1 .nd file detected. Please keep only 1 .nd file in the directory or specify it manually"
    )
  if (length(opt$nd_folder_path) < 1)
    stop(
      "No .nd file detected. Please keep exactly 1 file in the directory. The name has to contain the word .nd"
    )
  opt$nd_folder_path <- dirname(opt$nd_folder_path)
}

#Get the plate definition file
if (is.null(opt$pdef_path)) {
  #get all .csv files that have 'pdef' in their name
  opt$pdef_path <-
    list.files(
      path = opt$working_dir,
      pattern = ".*pdef.*\\.csv",
      full.names = TRUE,
      recursive = TRUE
    )
  if (length(opt$pdef_path) > 1)
    stop(
      "More than 1 pdef file detected. Please keep only 1 pdef file in the directory or specify it manually"
    )
  if (length(opt$pdef_path) < 1)
    stop(
      "No pdef file detected. Please keep exactly 1 file in the directory. The name has to contain the word 'pdef' and it must be a .csv file."
    )
}


#Get the folder that contains the raw data (txt files)
if (is.null(opt$results_txt_path)) {
  opt$results_txt_path <-
    list.files(
      path = opt$working_dir,
      pattern = "images_output",
      full.names = TRUE,
      recursive = TRUE,
      include.dirs = TRUE,
      ignore.case = TRUE
    )
  if (length(opt$results_txt_path) > 1)
    stop(
      "More than 1 results folder detected. Please keep only 1 folder in the directory or specify it manually"
    )
  if (length(opt$results_txt_path) < 1)
    stop(
      "No results folder detected. Please keep exactly 1 folder in the directory. The name has to be images_output (not case sensitive)"
    )
}

#Check whether folder "results" exist, create if not and assign to opt$out_dir
if (is.null(opt$out_dir)) {
  OUT_PATH = paste(opt$working_dir, "results", sep = "/")
  if (!file.exists (OUT_PATH)) {
    dir.create(OUT_PATH, showWarnings = TRUE)
  }
  opt$out_dir = OUT_PATH
}



cat(
  "\n",
  "*********************************************", "\n\n",
  "working directory..", "\t", opt$working_dir, "\n",
  ".nd folder path....", "\t", opt$nd_folder_path, "\n",
  "pdef file path.....", "\t", opt$pdef_path, "\n",
  "result files path..", "\t", opt$results_txt_path, "\n\n",
  "*********************************************",
  "\n\n",
  sep = ""
)


pdef_head <- colnames(read.csv(opt$pdef_path,sep=",")) #get column names of pdef file
CHANNELS_INPUT <- unlist(lapply(pdef_head,function(x) if(x=="GFP" | x=="RFP" | x=="BFP") x))


# Create Experiment Info -------------------------------------------------------

design = microscope.get.design(
  F = opt$nd_folder_path,
  D = c("comp"),
  PDEF = opt$pdef_path,
  FORMAT = opt$wells,
  OUT = opt$results_txt_path,
  CHANELS = CHANNELS_INPUT,
  MEASURE.COL = opt$column_to_measure,
  DIR.res = opt$out_dir
)


# Process Raw Data ------------------------------------------------------

df.for.analysis <-
  f.process.raw.data(
    design.file = design,
    min.cell.size = 1000,
    max.cell.size = 2500,
    brightfield.cutoff = 0.8,
    normalize.concentration = TRUE
  )

# Save Processed Data -----------------------------------------------------

save(list=c("df.for.analysis","design"),file = paste(opt$working_dir,"data_ready_to_analyze.RData",sep="/"))
