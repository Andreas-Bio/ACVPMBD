#disable scientific notation
options(scipen = 999)

#name of the run
out_dir_name <- "output"

#possible options SINTAX RDP
#any other value will just result in the metadata to be written (no sequence files)
export_tax_style <- "SINTAX"

#set directory to vsearch and ITSx executables
#don't forget the final /
vsearch_directory <- "C:/R/"

itsx_directory <- "/home/bio/Software/ITSx_1.1.3/" 

#define number of cores
number_of_cores <- 10

#This requires at least 20GB of free hard drive space.
delete_temporary_files <- FALSE

#define search string to download data from GenBank
genbank_searchstr <- "(internal transcribed spacer[Title] OR ITS2[Title] OR ITS1[Title]) AND Tracheophyta[Organism] NOT UNVERIFIED[All Fields] AND 200[SLEN] : 50000[SLEN] AND is_nuccore[filter] NOT plastid[filter] NOT patent[Title]"
                         
#If you have an entrez API key set it here to speed up the download. Make sure it is active and correctly entered or the download will fail.
#for more details see: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
set_entrez_key("")

#minimum and maximum accepted lengths, these values have been tuned for Tracheopyhta
min_58S_length <- 150
max_58S_length <- 170
its_minlen <- 100  #both ITS1 and ITS2

itsx_profiles <- "T,L" #Tracheophyta and Liverworts (needed for ferns). Change to anything else if other sequences are included (e.g. algae). See ITSx manual for codes.

