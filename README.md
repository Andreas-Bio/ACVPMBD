under construction

# ACVPMBD

Automated Curation of Vascular Plant (Meta-)Barcoding Databases

This R-based workflow will create reference databases for ITS plant
barcoding (ITS1 & ITS2 & full ITS) in an unattended fashion. Sequence
files from GenBank in .xml format will be locally downloaded and
processed by R and VSEARCH. The GBIF Application Programming Interface
(API) is used for taxonomic checks.

<u>Dependencies & Requirements</u>

ITSx \>= 1.1.3, https://microbiology.se/software/itsx/

> HMMER \>= 3.4, http://hmmer.org/

VSEARCH \>= 2.26, https://github.com/torognes/vsearch

R \>=4.4.0, https://cran.r-project.org/, if Ubuntu: https://cran.r-project.org/bin/linux/ubuntu/
It might be necessary to install additional libraries: `sudo apt-get install libfreetype6-dev libpng-dev cmake libtiff5-dev libjpeg-dev libsodium-dev libharfbuzz-dev libfribidi-dev`

20GB of free disk space

8GB of RAM

Unix-like platform (WSL2 requires additional configuration, not supported)

<u>Installation</u>

1.  The GitHub repository needs to be copied in a user-created folder
    with read/write access.

The script will automatically determine its location and set the working
directory in R accordingly.

2.  Paths to dependent software (VSEARCH and ITSx) must be set in the
    config.R file.

There are checks at the start of the script to check if the paths have
been set correctly. However, it is recommended to run these tools
independently and check if they are working properly. It is not
necessary to include these dependencies in the PATH, downloading and
dropping the files & folders anywhere on the system or in the resources
folder is sufficient.

3.  Additional parameters within the config.R file include the
    specification of the number of threads to be utilized by the script.
    This can be manually defined and is set to a dynamic detection mode
    by default: parallel::detectCores()-1. Please note that this default
    setting is not adequate for HPC clusters.

4.  Required packages will be downloaded automatically during the first
    execution of the script. Updating existing packages first is
    recommended. In case the package installation results in an error it
    is necessary to re-start the script. If the error persists the renv
    information included in the resources directory may be used to
    install the exact package versions which were used during the
    creation of this script.

<u>Usage</u>

The script itself can run on Unix-like systems and Windows, however ITSx
requires a UNIX-like environment. This also restricts the script to
UNIX-like systems for all practical reasons. Should parts of the script
be executed on a Windows machine, it is recommended to install
additional software first:
https://cran.r-project.org/bin/windows/Rtools/

The script can be run by the terminal command: Rscript
pathtoscript/main.R and should be able to complete its run without any
user intervention. Optional: If RStudio is installed the script can be
run with Alt+Ctrl+R from within RStudio.

Internet connection is required throughout the script runtime. Although
all downloads are split into batches and each batch download is
attempted twice, a loss of internet connection longer than 5 minutes
will cause the script to throw an error and terminate.

Please note that each execution of the script automatically overrides
temporary files eventually left over from the last run. Backing those
files up is not necessary but eventually can provide some means of
troubleshooting if unexpected results are found in the final database. A
unique name for each run should be defined in the config.R file to avoid
overwriting the script output of the final database. This name will be
used to create a unique output folder. It is highly recommended to
assign a DOI (<https://www.doi.org>) to any database used in
publications (e.g., <https://zenodo.org>).

The metadata and final database are stored in a tab-delimited file using
the .csv file ending encoded with UTF-8. Note that to open this file in
Excel it is required to use Data->From Text/CSV instead of double
clicking. Any errors can be safely ignored, these are conversion errors
if the date is set to NA and can not be converted by Excel.

<u>Step-by-step explanation</u>

1.  Startup checks

The script performs start-up and connection checks. In detail, (A) the
rentrez package is used to query GenBank for Ephedra glauca\[Organism\].
The first hit is subsequently downloaded, but not saved to the hard
disk. This check only passes if the rentrez package has been installed
correctly and a connection to GenBank is successfully established. (B)
The MD5 checkfile of the NCBI taxonomy is downloaded to the hard drive.
This test only succeeds if the user has writing access to the script
folder and NCBI did not change their ftp-server address (has been stable
since 5+ years). The taxonomy file is being referenced twice in the
functions.R file:
<https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip.md5> during the
startup tests and <https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip>
in the function taxonomy_assignment_step4. (C) The GBIF API is queried
by checking the backbone taxonomy for *Ephedra glauca*. This checks if
the rgbif package has been installed correctly and if a connection to
GBIF could be established. (D) The VSERACH test first checks if a file
called vsearch can be found in the directory specified in the config.R
file. Second, vsearch is run with the parameter --usearch_global using
the ./resources/Marchantiophyta_final_addins_ITSfull.fasta file as a
query and database. This verifies that VSEARCH can be executed and the
correct version (depending on the operating system) has been installed.
This test also fails if the script directory is not writeable. (E) The
ITSx test checks if a file called ITSx can be found in the directory
specified in the config.R file. (F) Finally, the script checks for the
permission to create folders and to delete them.

2.  Search sequences in GenBank

The search term used as a query to find sequences in GenBank is
specified in the config.R file. Note that although the theoretical
minimum length assumed for ITS1/2 is 100bp, the lower size limit is set
to 200bp. This excludes sequences which are incomplete and/or have been
trimmed. Initial tests showed that decreasing this size limit to 100bp
results in an increase of just 1% more sequences. ITS sequences which
have been trimmed are extremely difficult to verify, as without the
conserved flanking regions it is often impossible to verify and detect
pseudogenes. Note that the keyword UNVERIFIED is used by GenBank to tag
sequences which have incorrect annotations (e.g., a plant sequence
annotated as fungi), if the incorrect annotation has been reported. The
upper size limit of 50000 bp has been implemented to exclude large
genomic contig sequences. This size limit currently only excludes one
sequence, however in the near future that number is expected to grow.
The upper size filter should not be smaller than 15000 bp, as this would
exclude nrDNA repeat contigs, often featuring complete 26S sequences. It
is important to consider that ribosomal DNA can also be found in the
mitochondrion and the chloroplast. Although these sequences have a very
low chance of being picked up by ITSx, as it has only been trained on
nuclear ribosomal DNA, it is safer to exclude them from the start. It is
crucial for as many parameters as possible to define exactly where they
are expected to be found. Simply searching for its2 without specifying
\[Title\] will result in \~20000 more matches, most of which are not ITS
sequences. This is because this search also includes paper which have
the word ITS in their title, e.g.
<https://www.ncbi.nlm.nih.gov/nuccore/AJ512245>. For the same reason the
term environmental is not excluded in the search, but during data
curation and only if it is part of the taxonomic descriptor. The default
search string therefore is: (internal transcribed spacer\[Title\] OR
ITS2\[Title\] OR ITS1\[Title\]) AND Tracheophyta\[Organism\] NOT
UNVERIFIED\[All Fields\] AND 200\[SLEN\] : 50000\[SLEN\] AND
is_nuccore\[filter\] NOT plastid\[filter\] NOT patent\[Title\]

3.  Download sequences from GenBank

Downloading large chunks of data, such as in this workflow, from GenBank
can be challenging for multiple reasons. First, there is a limit on the
number of parallel connections and second, the download speed is
limited. Third, the number of accessions downloaded via the package
rentrez, and ultimately by the NCBI API, is limited. The maximum number
of parallel connections is 3 for unregistered users. This limit can be
increased to 10 parallel connections by providing an API key in the
config.R file. The API key can be obtained for free by following these
instructions:
<https://support.nlm.nih.gov/knowledgebase/article/KA-05317/en-us>. This
potentially, depending on other factors, can cut the download time by 15
minutes. If one of the batch download chunks fails (each download chunk
is 500 accession numbers), it is re-tried once after 5 minutes have
elapsed. If it fails again the script will terminate with an error. This
will increase the download time by multiple hours if the internet
connection is instable and fails frequently (i.e., instable WLAN). In
extreme cases, should 5 minutes not prove to be sufficient, this time
limit can be increased by searching for error = function(e){ throughout
the script and change Sys.sleep(300) to the desired value.

4.  Parsing downloaded xml files

The xml files are processed in parallel and filtered. The GBSeq_organism
descriptor, which usually holds the binomial name is filtered for the
following keywords: environmental unverified spec cultivated spp. sp.
cf. and unexpected special characters (e.g., ! or &) or numbers. Some of
these filters only apply if the keyword was found at the end of the
string to avoid removing legitimate species, such as *Cheilocostus
speciosus*. The only special character which should be allowed is a
dash, as found, for example, in *Saxifraga federici-augusti*. The next
filter removes sequences if the GBSeq_taxonomy descriptor, usually
containing the whole taxonomic lineage does not contain a plant family
name, identified by ending on ceae. Although sequences containing more
than 1% of ambiguous nucleotides (=IUPAC nucleotide ambiguities) can
potentially be considered noisy, this script excludes sequences if they
have more than 2%. The reason behind accepting up to 2% IUPAC
nucleotides is that the intragenomic variations of ITS often make
ambiguous base calls of Sanger sequencing trace files necessary and do
not reflect noisy sequences but biological meaningful information.
Metadata extracted by the script, for example the sample location,
currently only contains information for \~44% of the sequences, however,
this number is expected to increase. For more information, see:
<https://ncbiinsights.ncbi.nlm.nih.gov/2023/05/01/sequences-genbank-sra/>.

5.  ITSx detection

ITSx has been published in 2013 it is continuously updated by the
author. ITSx outputs multiple files, however, the most comprehensive
evaluation requires parsing of the positions file. Although it may look
like ITSx has poor detection rates in some taxonomic groups, such as
Bromeliaceae (>90% of sequences with no detections), manual inspection
of these taxa reveals that they usually contain a large number of
untrustworthy sequences which are either too short, contain a lot of
unexpected A/T substitutions (strong indicator of pseudogenes), are
incorrectly annotated or contain inconsistent gaps in conserved regions.
This is due to the fact that most plant barcoding protocols are unfit to
deal with high GC sequences, leading to systematic taxonomic bias. For
details, see: DMSO and betaine significantly enhance the PCR
amplification of ITS2 DNA barcodes from plants
(<https://cdnsciencepub.com/doi/10.1139/gen-2019-0221>). ITSx is able to
detect chimeras by checking the order and count of ITS regions (18S,
ITS1, 5.8S, ITS2, 26S), which, if detected, are removed by this script.
Although sequences are not kicked from the database if the detection of
the 5.8S region failed, it is recommended to exercise caution. These
ITS1/2 sequences potentially still have very short 5.8S snippets
attached to them (too short for ITSx to detect) or have been generated
by internal ITS primers, resulting in partial sequences. The relaxed
upper length threshold of 500bp is only applied to ITS2, as Gymnosperm
ITS1 sequences can be up to 1500bp long (e.g., Pinaceae). The lower
length threshold can be defined in the config.R file and is also set to
a relaxed 100bp. This potentially allows partial ITS1/2 sequences to
pass this filter. A more aggressive minimum filter setting would be
150bp, however even partial ITS1/2 sequences can still hold valuable
information which can be used to infer at least a genus-level in most
applications. As the 5.8S region is highly conserved, a more narrowly
defined filter can be applied which only allows sequences between 150bp
and 170bp to pass. 99% of 5.8S sequences are in the range of 158bp to
162 bp, so a stricter filter should be applied if the downstream
application depends on clearly defined 5.8S sequences, which however is
not the case for barcoding.

6.  Rescue sequences

As ITSx depends on flanking regions (18S, 5.8S, 26S) to detect ITS
sequences, but some sequences on GenBank have been trimmed prior to
their upload, these sequences would be automatically removed. In an
attempt to keep these sequences in the database, already verified
sequences are used as a reference to compare sequences to which were not
successfully annotated so far. The threshold of 85% (minimum required
similarity including terminal gaps), is expected to only match sequences
which already are represented on genus level in the database. To exclude
partial sequences, a maximum terminal gap length of 10 bp (each) in the
alignment of query (potentially truncated unannotated sequence) and
match (annotated sequences from previous step) is enforced. In addition,
an upper size limit of 300 bp filters sequences which are unannotated
despite containing a flanking region.

7.  Taxonomic harmonization

The NCBI taxonomy (<https://www.ncbi.nlm.nih.gov/taxonomy>), as
mentioned in their disclaimer, is not an authoritative source of
taxonomic names. It contains entries, such as Magnoliophyta sp. MP 525,
which are given the rank of species. Recent (2022) taxonomic changes,
such as in the genus *Psoralea*, are not reflected in the database as of
Dec 2023. Another example is the genus *Isotrema*, which has been
revised in 2019 but not updated in the NCBI taxonomy. The scientificName
is copied from GBIF, meaning that synonyms are also represented in this
data descriptor. All other information taken from GBIF (species, genus,
family, order, phylum and taxid_gbif) point to the current accepted
taxonomic name. NCBI taxonomy identification is retained for downstream
verification (taxid_ncbi, taxonomy_ncbi and authority_ncbi).

8.  Dereplication

The script reduces sequences to a default of 10 per species. It
prioritizes complete ITS sequences, making the database suitable for
projects targeting both ITS1 and 2. While considering the country of
origin would ideally represent maximum intraspecific diversity, GenBank
often lacks this data. Thus, the script uses the first 3 characters of
the accession numbers to consider different studies, if possible.

9.  Removing off-target sequences

Despite ITSx using HMMER profiles trained on Tracheophyta and Marchantiophyta, a small number of fungal sequences may still pass due to misannotations, contamination, or non-specific primers. To remove these, the script aligns sequences against a curated fungal ITS reference using vsearch, retaining only hits longer than 100 bp and flagging those exceeding an adaptive similarity threshold (median + 4×IQR). Users may replace the fungal reference files in the resources folder if needed, but caution is advised, as public fungal databases like UNITE can contain mislabeled plant-origin sequences. This step helps ensure the final output includes only verified vascular plant and Marchantiophyte ITS regions.
