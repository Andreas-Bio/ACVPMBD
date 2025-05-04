# === PATHS ===
# Path to the ITSx executable (version ≥1.1.3 recommended). Ensure the final slash is included.
# Change this if ITSx is installed elsewhere on your system.
# Use an absolute path if running from a different location or environment.
itsx_directory <- "/home/bio/Software/ITSx_1.1.3/" 
vsearch_directory <- "./"

# === OUTPUT SETTINGS ===
# Name of the output folder. All final files will be stored here.
# Choose a descriptive name for each run to avoid overwriting previous outputs.
out_dir_name <- "output"

# Taxonomy annotation format for exported FASTA headers:
# Options:
#   "SINTAX" – compatible with VSEARCH and USEARCH-style classifiers.
#   "RDP"    – compatible with RDP classifier format.
# Any other value disables sequence export and only generates metadata.
export_tax_style <- "SINTAX"

# === PARALLELIZATION ===
# Number of CPU cores used for parallel processing.
# Use one less than total available to avoid freezing the system.
# Set manually if needed (e.g., number_of_cores <- 8).
number_of_cores <- parallel::detectCores()-1

# If TRUE, deletes temporary intermediate folders (step1–step9) after processing.
# Requires at least 20 GB of free disk space if set to FALSE.
# Recommended: TRUE for stable servers or routine runs; FALSE for debugging or audits.
delete_temporary_files <- FALSE

# === GENBANK SEARCH QUERY ===
# Query string passed to GenBank's `entrez_search()`. Adjust for different target groups.
genbank_searchstr <- "(internal transcribed spacer[Title] OR ITS2[Title] OR ITS1[Title]) AND Tracheophyta[Organism] NOT UNVERIFIED[All Fields] AND 200[SLEN] : 50000[SLEN] AND is_nuccore[filter] NOT plastid[filter] NOT patent[Title]"
                         
#If you have an entrez API key set it here to speed up the download. Make sure it is active and correctly entered or the download will fail.
#for more details see: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html
set_entrez_key("")

# Length thresholds for the 5.8S region (X58S).
# These are based on known 5.8S lengths in Tracheophyta (typically 155–167 bp).
# Relax only for algae or highly divergent groups.
min_58S_length <- 150
max_58S_length <- 170
its_minlen <- 100  #both ITS1 and ITS2

# ITSx profiles to use for detecting ITS regions.
# "T" = Tracheophyta, "L" = Liverworts (includes non-angiosperm lineages like ferns).
# Use "A" for algae, "F" for fungi, etc. See ITSx manual for full list.
# Adjust only if working with non-vascular plants or additional kingdoms.
itsx_profiles <- "T,L"

# Disable scientific notation in output (e.g. 1e+05 becomes 100000).
# Leave unchanged for compatibility with output formatting and plotting.
options(scipen = 999)
