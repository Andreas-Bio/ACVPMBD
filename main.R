#This is the main script for the automated ITS reference database curation
# Version 1.00, Author: Andreas Kolter

if (!"this.path" %in% installed.packages()) install.packages("this.path")

library("this.path")
setwd(this.path::here())

source("functions.R")

#load required packages, this installs packages from CRAN if needed
load_package("rstudioapi")
load_package("magrittr")
load_package("rentrez")
load_package("xml2")
load_package("utils")
load_package("dplyr")
load_package("BiocManager")
load_package("rgbif")
load_package("foreach")
load_package("parallel")
load_package("doParallel")
load_package("Biostrings", style="BiocManager")
#package countrycode used to create countrycodes.csv
#see: https://en.wikipedia.org/wiki/United_Nations_geoscheme for more information
load_package("ggplot2") #only needed for visualization
load_package("data.tree") #only needed for visualization

source("config.R")

#check internet connection and functionality of packages
connection_check()


genbank_rentrez <- entrez_search(db="nuccore", term=genbank_searchstr, use_history = T, retmax=10^8)
paste0("Query resulted in ", genbank_rentrez$count, " GenBank hits.") %>% print

#downloads raw data from GenBank, in case of any connection interrupt re-download from start is recommended, although partial download is supported (see functions.R)
#creates a temp directory which contains the raw data in xml format, old temp directory is overwritten automatically
#this step can take multiple hours due to the bandwidth limit imposed by NCBI

download_rentrez_search(genbank_rentrez) #this step should take ~1 hour

#process the xml files, extract metadata and sequence
raw_input <- list.files(path="./temp", pattern = "its_from_[0-9]*_raw.xml", full.names = T) 
extract_information_raw(raw_input)

#read all .fasta files and concatenate them then run ITSx to determine ITS regions, also filters for special characters in species names
#depending on your processor speed ITSx takes between 2 and 4 hours
#please make sure ITSx has multiple GB of hard disk space available, it creates rather huge temporary files
#please make sure to run the ITSx test example to check functionality
step1_combine_fasta()

file.copy(from="./temp/step1.fasta", to= paste0(itsx_directory,"step1.fasta"))
system(command = paste0(itsx_directory,"ITSx", args=paste0(" -i ",itsx_directory,"step1.fasta -o step1 --reset T",
                                                           " --multi_thread T --cpu ",number_of_cores," -t ",itsx_profiles," --complement F --save_regions none",
                                                           " --preserve T --fasta F --graphical F --summary F")))
if (file.exists(paste0(itsx_directory,"step1.fasta"))) file.remove(paste0(itsx_directory,"step1.fasta"))

#creates a step1 folder and copies the fasta, ITSx and metadata csv file
migitate_to_step1_folder()

#add ITSx information to metadata, filter for chimeric and short ITS1/2 sequences
#create step2 folder and delete it if it previously existed
filter_by_ITSx(input_folder = "./step1/", output_folder = "./step2/", its_minlen=its_minlen,
               meta_in = "step1_meta.csv",meta_out = "step2_meta.csv",
               itsx_file = "step1.positions.txt", dna_file = "step1.fasta",
               min_58S_length=min_58S_length, max_58S_length=max_58S_length)

#optional: rescue sequences and consolidate meta-files #move to step3 folder
rescue_consolidate(rescue_me=T) #this step should take less than 10 minutes

#consolidate taxonomy with the GBIF backbone taxonomy, report errors to GBIF (they will be fixed eventually)
taxonomy_assignment_step4() #this step should take less than 20 minutes

#only keep a maximum number of sequences per species to avoid cluttering the reference database
seq_dereplicate(max_seqspec=10) #this step should take less than 10 minutes

#compares downloaded sequences to fungi database built with selected UNITE sequences (one per fungal genus)
remove_fungi(whichITS=0) #0 = complete ITS
remove_fungi(whichITS=1)
remove_fungi(whichITS=2)


#adds spike-in sequences to increase fidelity in case of fungal amplification during metabarcoding
add_spike_ins()

#sanity checks, this only affects Tracheophyta, not the spike-ins as they are expected to have larger gaps
dist_clean(whichITS=0, family_forced=T, max_intraspecific=10, max_global=30, max_IUPAC=0.025) #each dist_clean step should finish in approx. 30 minutes
dist_clean(whichITS=1, family_forced=T, max_intraspecific=10, max_global=30, max_IUPAC=0.025)
dist_clean(whichITS=2, family_forced=T, max_intraspecific=10, max_global=30, max_IUPAC=0.025)
dist_clean(whichITS=5, family_forced=F, max_intraspecific=5, max_global=10, max_IUPAC=0.01)

migrate_to_step9()

get_distribution()

export_seq(include_GBIF_taxID=F, whichITS=0, style=export_tax_style, out_dir_name=out_dir_name) #0 = complete ITS
export_seq(include_GBIF_taxID=F, whichITS=1, style=export_tax_style, out_dir_name=out_dir_name) 
export_seq(include_GBIF_taxID=F, whichITS=2, style=export_tax_style, out_dir_name=out_dir_name) 
export_seq(include_GBIF_taxID=F, whichITS=0, style="RDP", out_dir_name=out_dir_name) #0 = complete ITS
export_seq(include_GBIF_taxID=F, whichITS=1, style="RDP", out_dir_name=out_dir_name) 
export_seq(include_GBIF_taxID=F, whichITS=2, style="RDP", out_dir_name=out_dir_name) 

viz_length(minseq_per_family=50, width_graphic_scale=0.156,
           height_graphic_absolute=14.4,y_axis_upper_limit=2000,
           y_axis_number_ticks=20, x_axis_font_size=10,
           y_axis_font_size=14, legend_offset_left=0.1)


krona_plot()


# #examples how to subset the database, requires output from export_seq ITS0 (full ITS)
# meta <- paste0("./", out_dir_name, "/metadata_full.csv") %>% 
#   read.table(., header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
# refdat <- paste0("./", out_dir_name, "/ITS0_SINTAX.fasta") %>% readDNAStringSet()
# 
# #region names that can be used for subsetting
# meta[,"region"] %>% gsub("\\(.*?\\)","",.) %>% strsplit(.," , ") %>% unlist %>% gsub(" $","",.) %>% table %>% sort %>% rev %>% names
# refdat[refdat %>% names %>% gsub("_.*|;.*","",.) %in% ( grep("N-Africa", meta[,"region"]) %>% meta[.,"acc"] ) ] %>% 
#   writeXStringSet(filepath = paste0("./", out_dir_name, "/ITS0_SINTAX_N-Africa.fasta"))
# 







