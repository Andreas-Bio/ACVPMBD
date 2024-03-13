#this contains the functions of the refdat script to increase readability
load_package <- function(my_package, style="CRAN")
{
  suppressWarnings({
    
    if (style=="BiocManager")
    {
      if ( library(package = my_package, logical.return = TRUE, character.only = TRUE) )
      {
        library(my_package, quiet = TRUE, character.only = TRUE)
      } else {
        BiocManager::install(my_package, ask=FALSE); library(my_package, quiet = TRUE, character.only = TRUE)
      }
    }
    
    
    if ( library(package = my_package, logical.return = TRUE, character.only = TRUE) )
    {
      library(my_package, quiet = TRUE, character.only = TRUE)
    } else {
      install.packages(my_package, quiet = TRUE, character.only = TRUE); library(my_package, quiet = TRUE, character.only = TRUE)
    }
    
    if (library(package = my_package, logical.return = TRUE, character.only = TRUE)) return(paste("package", my_package, "successfully loaded"))
    
  })
}

connection_check <- function(checkme = c("genbank", "ncbi_tax", "gbif", "create_folder", "vsearch", "itsx"))
{
  if ("genbank" %>% grepl(., checkme) %>% sum > 0)
  {
    genbank_rentrez <- entrez_search(db="nuccore", term="Ephedra glauca[Organism]", use_history = T, retmax=1)
    temp <- entrez_fetch(db="nuccore", web_history = genbank_rentrez$web_history,
                         retstart = 1, rettype = "xml", retmode="xml", retmax = 1, parse=F)
    print("NCBI GenBank connection test successful")
  }
  
  if ("ncbi_tax" %>% grepl(., checkme) %>% sum > 0)
  {
  curl::curl_download("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip.md5", destfile="./ncbi_taxonomy.zip.md5")
  if (file.exists("./ncbi_taxonomy.zip.md5")) temp_silent <- file.remove("./ncbi_taxonomy.zip.md5")
  print("NCBI taxonomy database connection test successful")
  }
  
  if ("gbif" %>% grepl(., checkme) %>% sum > 0)
  {
    temp <- name_backbone_checklist(name_data = "Ephedra glauca", rank = "SPECIES", kingdom = "PLANTAE", verbose=TRUE)
    if (temp[1,"matchType"]=="EXACT") { print("GBIF taxonomy database connection test successful") } else {
      stop("GBIF taxonomy database returned unexpected result")
    }
  }
  
  if ("vsearch" %>% grepl(., checkme) %>% sum > 0)
  {
    if ( list.files(path=vsearch_directory, pattern="^vsearch") %>% length > 0 ) 
    { print("vsearch executable found") 
    } else { stop("vsearch executable not found") }
    
    if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
    print("running vsearch executable test")
    temp <- system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./resources/Marchantiophyta_final_addins_ITSfull.fasta --db ./resources/Marchantiophyta_final_addins_ITSfull.fasta",
    " --id 0.85 --threads ",number_of_cores," --userout ./resources/vsearch_testrun.out --iddef 1 --userfields query+target+id+alnlen+mism+opens+qilo+qihi+tilo+tihi+ql+tl+qrow"))
    
    if (file.exists("./resources/vsearch_testrun.out") & file.size("./resources/vsearch_testrun.out") > 1000)
    { print("running vsearch executable test successfully completed")
    } else { stop("vsearch executable test failed")  }
    
    if (file.exists("./resources/vsearch_testrun.out")) temp <- file.remove("./resources/vsearch_testrun.out")
  }
  

  if ("itsx" %>% grepl(., checkme) %>% sum > 0)
  {
    if ( list.files(path=itsx_directory, pattern="^ITSx") %>% length > 0 ) { print("ITSx executable found") } else { stop("ITSx executable not found") }
  }
  
  if ("create_folder" %>% grepl(., checkme) %>% sum > 0)
  {
    if ("testdir" %in% dir()) unlink("testdir", force = TRUE, recursive = TRUE)
    dir.create("testdir")
    cat("ABC", file="./testdir/test.txt")
    unlink("testdir", force = TRUE, recursive = TRUE)
    print("creating a test folder and writing a file successful")
  }
  
}

#use file number as download start and stop, not sequence count
download_rentrez_search <- function(genbank_rentrez, download_start = 1, download_stop= ceiling(genbank_rentrez$count / 500) )
{
  closeAllConnections()
  
  if ( Sys.getenv("ENTREZ_KEY") != "")
    cl <- parallel::makeCluster(9, type = "PSOCK")
  if ( Sys.getenv("ENTREZ_KEY") == "")
    cl <- parallel::makeCluster(3, type = "PSOCK")
  
  doParallel::registerDoParallel(cl)
  if (foreach::getDoParRegistered() == FALSE) stop("CLuster could not be registered.")
  
  if (download_start == 1 && download_stop == ceiling(genbank_rentrez$count / 500) )
  {
    if ("temp" %in% dir()) unlink("temp", force = TRUE, recursive = TRUE)
    dir.create("temp")
  }
  
  foreach (j = download_start:download_stop, .inorder = FALSE, .packages=c("rentrez","magrittr","xml2")) %dopar%
    {
      options(scipen = 999)
      
      sample(11:33/100, size=1) %>% Sys.sleep()
      if (( j %% 20 == 0) | (j==download_start))
        paste(j * 500 - 500,"out of", download_stop*500 ,"downloaded") %>% print
      
      tryCatch(
        expr = {
          temp <- entrez_fetch(db="nuccore", web_history = genbank_rentrez$web_history,
                               retstart = j * 500 - 500, rettype = "xml", retmode="xml", retmax = 500, parse=F)
        },
        error = function(e){ 
          Sys.sleep(300)
          temp <- entrez_fetch(db="nuccore", web_history = genbank_rentrez$web_history,
                               retstart = j * 500 - 500, rettype = "xml", retmode="xml", retmax = 500, parse=F)
        }
      )

      temp_fileout <- paste0("./temp/its_from_",j * 500 - 500,"_raw.xml",collapse="")
      if (file.exists(temp_fileout)) file.remove(temp_fileout)
      xml2::write_xml(temp %>% xml2::as_xml_document(x = .), file = temp_fileout)
    }
  
  stopCluster(cl)
  closeAllConnections()
}

step1_combine_fasta <- function()
{
  fasta_input <- list.files(path="./temp", pattern = "its_filtered_startfrom_.*.fasta", full.names = T) 
  x <- readDNAStringSet(fasta_input)
  x <- x[x %>% names %>% gsub("^.*_|[A-Za-z_:-]*","",.) == ""]
  x <- x[x %>% names %>% grepl(":x_",.) %>% not]
  if (file.exists("./temp/step1.fasta")) file.remove("./temp/step1.fasta")
  writeXStringSet(x = x, filepath = "./temp/step1.fasta")
}

#use file number for start_from, not sequence number
extract_information_raw <- function(raw_files=raw_input, start_from = 1, disable_IUPAC_filter=F, disable_keyword_and_noSpecies_removal=F,
                                    disable_check_for_family_CEAE=F)
{
  cl <- parallel::makeCluster(number_of_cores, type = "PSOCK")
  doParallel::registerDoParallel(cl)
  if (foreach::getDoParRegistered() == FALSE) stop("CLuster could not be registered.")
  
  raw_files %<>% sub("^\\.","",.)
  temp_input <- raw_files
  
  foreach (i = start_from:length(temp_input), .inorder = FALSE, .packages=c("Biostrings","magrittr","xml2")) %dopar% 
    {
      result <- NULL
      metadata <- NULL
      options(scipen = 999)
      
      if (i %% 100 == 0) print(paste(i,"of",length(temp_input),"xml files processed"))
      temp_in <- read_xml(x = paste0( getwd(), temp_input[i]) ) %>% xml2::as_list(x = .) %>% .[[1]] 
      
      for (j in 1:length(temp_in))
      {
        temp_subset <- temp_fam <- temp_org <- temp_acc <- NA
        temp_date <- temp_auth <- temp_title <- temp_orig_title <- temp_orig_taxonomy <- temp_journal <- NA
        temp_country <- temp_voucher <- temp_collection <- temp_xref <-  temp_isolation <-  temp_collectedby <- NA
        
        temp <- temp_in[[j]]
        if (length(temp)<5) next
        
        temp_id <- temp$`GBSeq_primary-accession`[[1]]
        if (nchar(temp_id)<6) next
        
        temp_org <- temp$GBSeq_organism[[1]] %>% gsub(" var\\..*| subsp\\..*","",.)
        if (disable_keyword_and_noSpecies_removal==FALSE)
          if ( grepl('environmental|unverified|spec$|sp$|cultivated|spp.$| sp. | cf\\. | x |[0-9]|"|§|\\$|\\%|\\&|\\/|\\(|\\)|\\=|\\?|\\!', temp_org) ) next
        
        temp_spec <- temp_org %>% gsub(".* ","",.)
        temp_gen <- temp_org %>% gsub(" .*","",.)
        
        temp_fam <- temp$GBSeq_taxonomy[[1]] %>% gsub("^.*; (.*aceae);.*$","\\1",.)
        
        if(disable_check_for_family_CEAE==FALSE)
        {
        if (nchar(temp_fam)<3) next; if (temp_fam %>% is.character %>% not) next
        if (grepl(" |\\;|Eukaryota",temp_fam)) next
        }
        
        temp_seq <- temp$GBSeq_sequence[[1]] %>% DNAStringSet
        if (temp_seq %>% length == 0) next
        
        names(temp_seq) <- paste0(temp_id,"_f:",temp_fam,"_g:",temp_gen,"_s:",temp_spec)
        
        if (disable_IUPAC_filter==FALSE)
          if(temp_seq %>% alphabetFrequency(.,as.prob=T) %>% .[5:16] %>% sum > 0.02) next
        
        if (!is.null(temp$`GBSeq_primary-accession`[[1]]))
          temp_acc <- temp$`GBSeq_primary-accession`[[1]]
        
        if (!is.null(temp$`GBSeq_create-date`[[1]]))
          temp_date <- temp$`GBSeq_create-date`[[1]]
        
        if (!is.null(temp$GBSeq_references$GBReference$GBReference_authors))
          temp_auth <- temp$GBSeq_references$GBReference$GBReference_authors %>% unlist %>% paste0(.,collapse="; ")
        
        if (!is.null(temp$GBSeq_references$GBReference$GBReference_title[[1]]))
          temp_title <- temp$GBSeq_references$GBReference$GBReference_title[[1]]
        
        if (!is.null(temp$GBSeq_definition[[1]]))
          temp_orig_title <- temp$GBSeq_definition[[1]]

        if (!is.null(temp$GBSeq_taxonomy[[1]]))
          temp_orig_taxonomy <- temp$GBSeq_taxonomy[[1]]
        
        if (!is.null(temp$GBSeq_references$GBReference$GBReference_journal[[1]]))
          temp_journal <- temp$GBSeq_references$GBReference$GBReference_journal[[1]]
        
        temp_subset <- temp$`GBSeq_feature-table`$GBFeature$GBFeature_quals
        
        if (grepl('country"', temp_subset) %>% sum == 1)
          temp_country <- grep('country"', temp_subset) %>% temp_subset[[.]] %>% .[[2]] %>% .[[1]]
        
        if (grepl('specimen_voucher"', temp_subset) %>% sum == 1)
          temp_voucher <- grep('specimen_voucher"', temp_subset) %>% temp_subset[[.]] %>% .[[2]] %>% .[[1]]
        
        if (grepl('collection_date"', temp_subset) %>% sum == 1)
          temp_collection <- grep('collection_date"', temp_subset) %>% temp_subset[[.]] %>% .[[2]] %>% .[[1]]
        
        if (grepl('"taxon:', temp_subset) %>% sum == 1)
          temp_xref <- grep('"taxon:', temp_subset) %>% temp_subset[[.]] %>% .[[2]] %>% .[[1]]
        
        #if (grepl('isolation_source(?=")', temp_subset) %>% sum == 1)
        #  temp_isolation <- grep('isolation_source', temp_subset) %>% temp_subset[[.]] %>% .[[2]] %>% .[[1]]
        
        if (grepl('collected_by"', temp_subset) %>% sum == 1)
          temp_collectedby <- grep('collected_by"', temp_subset) %>% temp_subset[[.]] %>% .[[2]] %>% .[[1]]
        
        metadata <-  cbind(temp_acc, temp_date, temp_auth, temp_title, temp_orig_title, temp_orig_taxonomy,
                           temp_journal, temp_country, temp_voucher, temp_collection, temp_xref, temp_isolation, temp_collectedby ) %>%
          gsub("\n"," ",.) %>% cbind.data.frame %>% rbind.data.frame(metadata, .)
        
        result <- result %>% append(.,temp_seq)
        
      }
      
      if(length(result)<1) next
      
      temp_filename <- paste0("./temp/its_filtered_startfrom_",(i-1)*500,".fasta",collapse="")
      if(file.exists(temp_filename)) file.remove(temp_filename)
      writeXStringSet(x=result, filepath=temp_filename)
      
      temp_filename <- paste0("./temp/its_metadata_startfrom_",(i-1)*500,".csv",collapse="")
      if(file.exists(temp_filename)) file.remove(temp_filename)
      write.table(x=metadata , temp_filename, quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")
      
    }
  
  
  stopCluster(cl)
  closeAllConnections()
}

migitate_to_step1_folder <- function()
{
  if ("step1" %in% dir()) unlink("step1", force = TRUE, recursive = TRUE)
  dir.create("step1")
  
  file.copy(from="./temp/step1.fasta", to= paste0("./step1/step1.fasta"))
  temp <- list.files(path="./temp",pattern = "its_metadata_startfrom.*\\.csv" ,full.names=TRUE) 
  temp_csv <- NULL
  for (i in 1:length(temp))
    temp_csv <- read.csv(file=temp[i], sep='\t') %>% rbind.data.frame(.,temp_csv)
  temp_csv <- temp_csv
  colnames(temp_csv) %<>% gsub("temp_","",.)
  colnames(temp_csv) %<>% gsub("^title","publication",.)
  if(file.exists("./step1/step1_meta.csv")) file.remove("./step1/step1_meta.csv")
  write.table(x=temp_csv , "./step1/step1_meta.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")   
  file.copy(from ="step1.positions.txt", to="./step1/step1.positions.txt")
  file.copy(from= "step1.positions.txt", to="./temp/step1.positions.txt")
  file.copy(from= "./temp/step1.positions.txt", to="./step1/step1.positions.txt")
  if(file.exists("./step1/step1.positions.txt") & file.exists("step1.positions.txt")) file.remove("step1.positions.txt")
  list.files(path=itsx_directory, pattern="step1", full.names = T) %>% file.remove(.)
}

filter_by_ITSx <- function(input_folder = "./step1/", output_folder = "./step2/", its_minlen, meta_in = "step1_meta.csv",
                           meta_out = "step2_meta.csv", itsx_file = "step1.positions.txt", dna_file = "step1.fasta", min_58S_length, max_58S_length)
{
  itsx <- read.csv2(file = paste0(input_folder,itsx_file), sep="\t", header = FALSE)
  dna <- readDNAStringSet(filepath = paste0(input_folder,dna_file))
  meta <- read.table(file=paste0(input_folder,meta_in), sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE)   
  
  print(paste0(nrow(itsx)," ITSx detections in ", length(dna), " sequences" ))
  
  #remove chimeric detections
  itsx <- itsx[!grepl("Chimeric",itsx[,"V8"]),]
  
  #reduce all datasets to the common acc
  temp_common <- intersect(itsx[,1] %>% gsub("_.*","",.), dna %>% names %>% gsub("_.*","",.) ) %>% intersect(., meta[,"acc"] )
  
  itsx <- itsx[itsx[,1] %>% gsub("_.*","",.) %in% temp_common, ]
  dna <- dna[dna %>% names %>% gsub("_.*","",.) %in% temp_common]
  meta <- meta[,"acc"] %in% temp_common %>% meta[.,]
  if (nrow(meta) != nrow(itsx)) print("ITSx table mismatch")
  
  #add informtion to metadata
  meta <- meta[match(itsx[,1] %>% gsub("_.*","",.), meta[,"acc"]),]
  dna <- dna[match(itsx[,1] %>% gsub("_.*","",.), dna %>% names %>% gsub("_.*","",.))]
  
  ##fix 5.8S missing detection due to a deletion at the start of ITS2
  #select all hits which have ITS1 and ITS2 detected
  temp_subset <- ( itsx[,"V4"] %>% gsub("[A-Za-z: ]","",.) %>% nchar > 1 ) &
    ( itsx[,"V6"] %>% gsub("[A-Za-z: ]","",.) %>% nchar > 1 )
  
  #select all hits which have ITS1 and ITS2 detected but no 5.8S
  temp_subset_no58S <- temp_subset & ( itsx[,"V5"] %>% gsub("[A-Za-z: ]|5.8S:","",.) %>% nchar == 0 )

  #apply the fix only to genera which show this in at least 33% of their sequences to avoid over-fixing random pseudogenes
  #get a data frame of genera which are in temp_subset_no58S
  temp_no58S <- itsx[temp_subset_no58S,"V1"] %>% gsub("^.*_g:(.*)_s:.*","\\1",.) %>% table %>% cbind.data.frame
  if (temp_no58S %>% nrow > 0) colnames(temp_no58S) <- c("gen", "count")
  
  #get a data frame of genera which are in temp_subset
  temp_its12 <- itsx[temp_subset,"V1"] %>% gsub("^.*_g:(.*)_s:.*","\\1",.) %>% table %>% cbind.data.frame
  if (temp_its12 %>% nrow > 0) colnames(temp_its12) <- c("gen", "count")
  
  if ( (temp_no58S %>% nrow > 0) & (temp_its12 %>% nrow > 0) )
  {
    #check overlap and merge by names
    temp <- merge(temp_no58S, temp_its12, by.x="gen", by.y="gen")
    
    #get a list of genera
    temp_genera_select <- temp[round(temp[,2]/temp[,3]*100)>33,1] %>% as.character
    
    #subset previously made selection
    temp_subset_no58S  <- temp_subset_no58S & ( itsx[,"V1"] %>% gsub("^.*_g:(.*)_s:.*","\\1",.) %in% temp_genera_select )
    
    temp <- itsx[temp_subset_no58S,]  
    
    temp_no_start <- temp[,"V5"]=="5.8S: No start"
    temp_no_end <- temp[,"V5"]=="5.8S: No end"
    temp_proposed_ITS1_end <- temp[temp_no_start,"V6"] %>% gsub("[A-Za-z: ]|ITS2","",.) %>% strsplit(.,"-") %>% lapply(.,`[`,1) %>% as.integer - 161
    temp_max_seq_length <- temp[temp_no_start,"V2"] %>% gsub(" bp.","",.) %>% as.integer
    temp_no_start[temp_proposed_ITS1_end + 260 > temp_max_seq_length] <- FALSE
    temp_proposed_ITS1_end <- temp_proposed_ITS1_end[temp_proposed_ITS1_end < temp_max_seq_length + 260]
    
    temp[temp_no_start,"V5"] <- paste0("5.8S: ",temp_proposed_ITS1_end+1,"-",temp_proposed_ITS1_end+160)
    temp[temp_no_start,"V4"] <- temp[temp_no_start,"V4"] %>% gsub("(ITS1: .*)-.*$","\\1",.) %>% paste0(.,"-",temp_proposed_ITS1_end)
    
    temp_proposed_ITS2_start <- temp[temp_no_end,"V4"] %>% gsub("[A-Za-z: ]|ITS1","",.) %>% strsplit(.,"-") %>% lapply(.,`[`,2) %>% as.integer + 161
    #make sure to not set a start of ITS2 which would be > end of ITS2 or a start which is > length of sequence
    temp_current_ITS2_end <- temp[temp_no_end,"V6"] %>% gsub("[A-Za-z: ]|ITS2","",.) %>% strsplit(.,"-") %>% lapply(.,`[`,2) %>% as.integer
    temp_max_seq_length <- temp[temp_no_end,"V2"] %>% gsub(" bp.","",.) %>% as.integer
    temp_kickme <- ( temp_current_ITS2_end < temp_proposed_ITS2_start + 100 ) & ( temp_proposed_ITS2_start + 100 < temp_max_seq_length )
    #subset the selection by all positive values and then set the kickme values as negative
    temp_no_end[temp_no_end][temp_kickme] <- FALSE
    temp_proposed_ITS2_start <- temp[temp_no_end,"V4"] %>% gsub("[A-Za-z: ]|ITS1","",.) %>% strsplit(.,"-") %>% lapply(.,`[`,2) %>% as.integer + 161
    
    temp[temp_no_end,"V5"] <- paste0("5.8S: ",temp_proposed_ITS2_start-160,"-",temp_proposed_ITS2_start-1)
    temp_constructor <- temp[temp_no_end,"V6"] %>% gsub("(ITS2: ).*(-.*)$","\\1\\2",.) %>% strsplit(.,"-")
    temp[temp_no_end,"V6"] <- paste0(lapply(temp_constructor,`[`,1),temp_proposed_ITS2_start, "-",lapply(temp_constructor,`[`,2))
    
    itsx[temp_subset_no58S,] <- temp
  }
  
  meta <- grepl("Not found", itsx[,"V3"]) %>% not %>% cbind.data.frame(meta, LSU_found=.)
  paste0(sum(meta[,"LSU_found"], na.rm=T), " 26S (LSU) regions detected by ITSx") %>% print
  temp_26S <- itsx[,"V3"] %>% gsub("[A-Za-z: ]","",.) %>% strsplit(.,"-")
  temp_26_start <- temp_26S %>% lapply(.,`[`,1) %>% as.integer
  temp_26_stop <- temp_26S %>% lapply(.,`[`,2) %>% as.integer
  meta <- ( temp_26_stop - temp_26_start + 1 )  %>% cbind.data.frame(meta, LSU_length=.)
  meta <- subseq(dna, start=temp_26_start, end=temp_26_stop) %>% cbind.data.frame(meta, LSU_seq=.)
  meta[meta[,"LSU_found"]==FALSE,"LSU_seq"] <- NA
  
  meta <- grepl("Not found", itsx[,"V7"]) %>% not %>% cbind.data.frame(meta, SSU_found=.)
  paste0(sum(meta[,"SSU_found"], na.rm=T), " 18S (SSU) regions detected by ITSx") %>% print
  temp_18S <- itsx[,"V7"] %>% gsub("[A-Za-z: ]","",.) %>% strsplit(.,"-")
  temp_18_start <- temp_18S %>% lapply(.,`[`,1) %>% as.integer
  temp_18_stop <- temp_18S %>% lapply(.,`[`,2) %>% as.integer
  meta <- ( temp_18_stop - temp_18_start + 1 )  %>% cbind.data.frame(meta, SSU_length=.)
  meta <- subseq(dna, start=temp_18_start, end=temp_18_stop) %>% cbind.data.frame(meta, SSU_seq=.)
  meta[meta[,"SSU_found"]==FALSE,"SSU_seq"] <- NA
  
  meta <- grepl("Not found", itsx[,"V4"]) %>% not %>% cbind.data.frame(meta, ITS1_found=.)
  paste0(sum(meta[,"ITS1_found"], na.rm=T), " ITS1 regions detected by ITSx") %>% print
  temp_ITS1 <- itsx[,"V4"] %>% gsub("[A-Za-z: ]|ITS1","",.) %>% strsplit(.,"-")
  temp_ITS1_start <- temp_ITS1 %>% lapply(.,`[`,1) %>% as.integer
  temp_ITS1_stop <- temp_ITS1 %>% lapply(.,`[`,2) %>% as.integer
  meta <- ( temp_ITS1_stop - temp_ITS1_start + 1 )  %>% cbind.data.frame(meta, ITS1_length=.)
  meta <- subseq(dna, start=temp_ITS1_start, end=temp_ITS1_stop) %>% cbind.data.frame(meta, ITS1_seq=.)
  meta[meta[,"ITS1_found"]==FALSE,"ITS1_seq"] <- NA
  
  meta <- grepl("Not found| No ", itsx[,"V5"]) %>% not %>% cbind.data.frame(meta, "58S_found"=.)
  paste0(sum(meta[,"58S_found"], na.rm=T), " 5.8S regions detected by ITSx") %>% print
  temp_58S <- itsx[,"V5"] %>% gsub("[A-Za-z: ]|5.8S:","",.) %>% strsplit(.,"-")
  temp_58S_start <- temp_58S %>% lapply(.,`[`,1) %>% as.integer
  temp_58S_stop <- temp_58S %>% lapply(.,`[`,2) %>% as.integer
  meta <- ( temp_58S_stop - temp_58S_start + 1 )  %>% cbind.data.frame(meta, "58S_length"=.)
  meta <- subseq(dna, start=temp_58S_start, end=temp_58S_stop) %>% cbind.data.frame(meta, "58S_seq"=.)
  meta[meta[,"58S_found"]==FALSE,"58S_seq"] <- NA
  
  meta <- grepl("Not found", itsx[,"V6"]) %>% not %>% cbind.data.frame(meta, "ITS2_found"=.)
  paste0(sum(meta[,"ITS2_found"], na.rm=T), " ITS2 regions detected by ITSx") %>% print
  temp_ITS2 <- itsx[,"V6"] %>% gsub("[A-Za-z: ]|ITS2","",.) %>% strsplit(.,"-")
  temp_ITS2_start <- temp_ITS2 %>% lapply(.,`[`,1) %>% as.integer
  temp_ITS2_stop <- temp_ITS2 %>% lapply(.,`[`,2) %>% as.integer
  meta <- ( temp_ITS2_stop - temp_ITS2_start + 1 )  %>% cbind.data.frame(meta, "ITS2_length"=.)
  meta <- subseq(dna, start=temp_ITS2_start, end=temp_ITS2_stop) %>% cbind.data.frame(meta, "ITS2_seq"=.)
  meta[meta[,"ITS2_found"]==FALSE,"ITS2_seq"] <- NA
  
  meta <- ( meta[,"ITS2_found"] & meta[,"ITS1_found"] ) %>% cbind.data.frame(meta, "58S_complete"=.)
  meta <- ( meta[,"SSU_found"] & meta[,"58S_found"] ) %>% cbind.data.frame(meta, "ITS1_complete"=.)
  meta <- ( meta[,"LSU_found"] & meta[,"58S_found"] ) %>% cbind.data.frame(meta, "ITS2_complete"=.)
  
  paste0("removed ", which(meta[,"ITS2_length"] > 500 | meta[,"ITS2_length"] < its_minlen) %>% length, " ITS2 sequences longer than 500bp or shorter than ",its_minlen,"bp") %>% print
  meta[which(meta[,"ITS2_length"] > 500 | meta[,"ITS2_length"] < its_minlen),"ITS2_found"] <- FALSE
  meta[which(meta[,"ITS2_length"] > 500 | meta[,"ITS2_length"] < its_minlen),c("ITS2_length","ITS2_seq")] <- NA
  
  paste0("removed ", which(meta[,"ITS1_length"] < its_minlen) %>% length, " ITS1 sequences shorter than ",its_minlen,"bp") %>% print
  meta[which(meta[,"ITS1_length"] < its_minlen),"ITS1_found"] <- FALSE
  meta[which(meta[,"ITS1_length"] < its_minlen),c("ITS1_length","ITS1_seq")] <- NA
  
  paste0("removed ", which(!(meta[,"ITS2_found"] | meta[,"ITS1_found"])) %>% length, " metadata entries where neither ITS1 or ITS2 were detected") %>% print
  meta <- meta[meta[,"ITS2_found"] | meta[,"ITS1_found"],] 
  
  paste0(nrow(meta), " metadata entries where either ITS1 or ITS2 were detected") %>% print
  
  output_clean <- output_folder %>% gsub("[^a-z0-9]","",.)
  
  closeAllConnections()
  if (output_clean %in% dir()) unlink(output_clean, force = TRUE, recursive = TRUE)
  if (output_clean %in% dir() %>% not) dir.create(output_clean)
  
  meta <- meta[!(meta[,"ITS2_found"]==TRUE & meta[,"ITS1_found"]==TRUE & meta[,"58S_length"] %>% is.na),]
  meta <- meta[(meta[,"58S_length"] > min_58S_length & meta[,"58S_length"] < max_58S_length) | meta[,"58S_length"] %>% is.na,]

  temp_dna_select <- dna[dna %>% names %>% gsub("_.*","",.) %in% meta[,"acc"]]
  dna <- dna[match(meta[,"acc"], dna %>% names %>% gsub("_.*","",.)) ] 
  
  writeXStringSet(dna, filepath = paste0(output_folder,"step2.fasta"))
  write.table(x=meta , paste0(output_folder,meta_out), quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")  
  
  #before 316524
}

rescue_sequences <- function(whichITS=1)
{
  
  if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
  #Run vsearch to find matches with more than 85% identity between the stray and reference database.
  if(!file.exists(paste0(vsearch_directory, vsearch_exec))) stop("vsearch not found")
  
  #rescue sequences which have been discarded by ITSx due to having no flanking regions
  meta2 <- read.table(file="./step2/step2_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE)   
  dna2 <- readDNAStringSet(filepath = "./step2/step2.fasta")
  if ( sum (( dna2 %>% names %>% gsub("_.*","",.) ) != meta2[,"acc"]) == nrow(meta2)) print("Sequence order in fasta and metadata not identical, repeat previous step.")
  rm(dna2)
  
  itsx <- read.csv2(file = "./step1/step1.positions.txt", sep="\t", header = FALSE)
  meta1 <- read.table(file="./step1/step1_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE) 
  dna1 <- readDNAStringSet(filepath = "./step1/step1.fasta")
  
  #which accession numbers can be found in the GenBank download but not in the ITSx file (without ITS detection)
  temp_lost <- meta1[,"acc"] %in% ( itsx[,1] %>% gsub("_.*","",.) ) %>% not
  temp_lost_meta1 <- meta1[temp_lost,]
  #Which of those GenBank sequence descriptions indicate that it is an ITS1/2 sequence?
  temp <- paste0("internal transcribed spacer ",whichITS,"|ITS",whichITS)
  temp_lost_its <- temp_lost_meta1[grepl(temp, temp_lost_meta1[,"orig_title"]),"acc"]
  
  dna1_its_lost <- dna1[dna1 %>% names %>% gsub("_.*","",.) %in% temp_lost_its] 
  dna1_its_lost <- dna1_its_lost[dna1_its_lost %>% width < 350]
  if (file.exists("./step2/its_potential_recover.fasta")) file.remove("./step2/its_potential_recover.fasta")
  writeXStringSet(x = dna1_its_lost, filepath = "./step2/its_potential_recover.fasta")
  
  #note: It is expected that ITSx was not able to detect those sequences because the flanking regions were missing.
  #     This is not intended to recover sequences with flanking regions which did not match in ITSx because those will have a high proportion of pseudogenes and chimeras.
  #Which of the ITS1 sequences is shoter than a fixed threshold (see above)? Write the results to a file.
  temp_its_ref <- paste0("ITS",whichITS,"_found") %>% meta2[,.] %>% meta2[.,paste0("ITS",whichITS,"_seq")] %>% DNAStringSet
  
  names(temp_its_ref) <- paste0("ITS",whichITS,"_found") %>% meta2[,.] %>% meta2[.,"acc"]
  if (file.exists("./step2/temp_its_extracted.fasta")) temp_silent <- file.remove("./step2/temp_its_extracted.fasta")
  writeXStringSet(x = temp_its_ref, filepath = "./step2/temp_its_extracted.fasta")
  
  if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
  #Run vsearch to find matches with more than 85% identity between the stray and reference database.
  system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./step2/its_potential_recover.fasta --db ./step2/temp_its_extracted.fasta --id 0.85 --threads ",number_of_cores,
                          " --userout ./step2/its_potential_recover.out --iddef 1 --userfields query+target+id+alnlen+mism+opens+qilo+qihi+tilo+tihi+ql+tl+qrow"))
  
  
  meta1[setdiff(names(meta2), names(meta1))] <- NA
  
  #vsearch output:
  #V1:query V2:target V3:id V4:alnlen V5:mism V6:opens V7:first nuc of query aligned V8:last nuc of query aligned 
  #V9: first nuc of target aliged V10: last nuc of target aligned ##all alignment positions ignore terminal gaps
  #V11: query sequence length V12:target sequence length V13:alignment of query sequence to target
  
  temp_its_rescue <- read.csv("./step2/its_potential_recover.out", sep="\t", header = F)
  
  temp <- temp_its_rescue[,"V7"] <= 10 & ( temp_its_rescue[,"V11"] - temp_its_rescue[,"V8"] <= 10 ) &
    temp_its_rescue[,"V11"] < 350 & nchar(temp_its_rescue[,"V13"])  < 300
  
  temp_its_rescue <- temp_its_rescue[temp,]
  
  temp_its_rescue_acc <- temp_its_rescue[,"V1"] %>% gsub("_.*","",.)
  temp_its_rescue[,"V13"] %<>% gsub("-","",.)
  
  meta_r <- meta1[meta1[,"acc"] %in% temp_its_rescue_acc,]
  meta_r <- meta_r[match(temp_its_rescue_acc, meta_r[,"acc"]),] 
  meta_r[,paste0("ITS",whichITS,"_found")] <- TRUE
  meta_r[,paste0("ITS",whichITS,"_length")] <- nchar(temp_its_rescue[,"V13"])
  meta_r[,paste0("ITS",whichITS,"_seq")] <- temp_its_rescue[,"V13"]
  
  meta_r
  
}

rescue_consolidate <- function(rescue_me=T)
{
  
  if ("step3" %in% dir()) unlink("step3", force = TRUE, recursive = TRUE)
  dir.create("step3")
  
  if(rescue_me==T)
  {
    x <- rbind.data.frame(rescue_sequences(whichITS=2), rescue_sequences(whichITS=1))
    x <- x[x[,"acc"] %>% duplicated %>% not,]
    meta2 <- read.table(file="./step2/step2_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE) 
    meta2 <- rbind.data.frame(meta2, x)
    meta2[meta2[,"LSU_found"] %>% is.na,"LSU_found"] <- FALSE
    meta2[meta2[,"SSU_found"] %>% is.na,"SSU_found"] <- FALSE
    meta2[meta2[,"ITS1_found"] %>% is.na,"ITS1_found"] <- FALSE
    meta2[meta2[,"X58S_found"] %>% is.na,"X58S_found"] <- FALSE
    meta2[meta2[,"ITS2_found"] %>% is.na,"ITS2_found"] <- FALSE
    meta2[meta2[,"ITS1_complete"] %>% is.na,"ITS1_complete"] <- FALSE
    meta2[meta2[,"ITS2_complete"] %>% is.na,"ITS2_complete"] <- FALSE
    meta2[meta2[,"X58S_complete"] %>% is.na,"X58S_complete"] <- FALSE
    meta2 <- meta2[!(meta2[,"X58S_complete"]==TRUE & meta2[,"X58S_seq"] %>% is.na),]
    
    write.table(x=meta2 , "./step3/step3_meta.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8") 
  } else {
    file.copy(from="./step2/step2_meta.csv", to= paste0("./step3/step3_meta.csv"))
  }
}

taxonomy_assignment_step4 <- function()
{
  delete_temporary_files <<- delete_temporary_files
  
  if (delete_temporary_files) 
  {
    if ("temp" %in% dir()) unlink("temp", force = TRUE, recursive = TRUE)
    if ("step1" %in% dir()) unlink("step1", force = TRUE, recursive = TRUE)
    if ("step2" %in% dir()) unlink("step2", force = TRUE, recursive = TRUE)
  }
  
  curl::curl_download("https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip", destfile="./step3/ncbi_taxonomy.zip")
  utils::unzip(zipfile = "./step3/ncbi_taxonomy.zip", overwrite = TRUE, exdir = "./step3")
  
  #load data for subspecies and variety linkage
  temp_nodes <- read.table(file = "./step3/nodes.dmp", sep="\t")
  temp_nodes <- temp_nodes[,temp_nodes[1,]!="|"]
  
  meta3 <- read.table(file="./step3/step3_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE) 
  temp_taxref <- meta3[,"xref"] %>% gsub("taxon:","",.) %>% as.integer
  
  xref <- cbind.data.frame("step0"=temp_taxref %>% unique)
  xref <- cbind.data.frame(xref, "step0_tax"=temp_nodes[match(xref[,"step0"],temp_nodes[,"V1"]),"V5"])
  
  for (i in 1:15)
  {
    this_step <- i
    last_step <- i-1
    last_xref <- xref %>% colnames %>% grepl(last_step,.) %>% xref[,.] %>% .[,1]
    this_xref <- match(last_xref, temp_nodes[,"V1"]) %>% temp_nodes[.,"V3"]
    this_tax <- match(this_xref, temp_nodes[,"V1"]) %>% temp_nodes[.,"V5"]
    xref <- cbind.data.frame(xref, this_xref, this_tax)
    colnames(xref)[(ncol(xref)-1):ncol(xref)] <- c(paste0("step",i),paste0("step",i,"_tax"))
  }
  
  tax_levels <- c("species","genus","family","order")
  
  xtax <- cbind.data.frame("species"=rep(NA, nrow(xref)))
  
  for (i in 1:length(tax_levels))
  {
    this_tax <- tax_levels[i]
    tax_hits <- which(xref==this_tax, arr.ind=T)
    tax_hits[,2] <- tax_hits[,2]-1
    xtax[tax_hits[,1],i] <- xref[tax_hits] %>% as.numeric
  }
  
  colnames(xtax) <- tax_levels
  xtax <- cbind.data.frame("orig_xref"=xref[,"step0"], xtax)
  xtax <- xtax[complete.cases(xtax),]
  xtax <- cbind.data.frame(xtax, xtax) 
  xtax <- xtax[,-1]
  
  tax_names <- read.table("./step3/names.dmp", sep="\t", comment.char = "",
                          quote = "", dec = ".", stringsAsFactors = F, blank.lines.skip = T)
  
  tax_sci <- tax_names[tax_names[,"V7"]=="scientific name",]
  
  for (i in 1:4)
    xtax[,i] <- match(xtax[,i], tax_sci[,"V1"]) %>% tax_sci[.,"V3"]
  
  tax_names <- tax_names[tax_names[,"V1"] %in% xtax[,"species.1"],]
  tax_names <- tax_names[!grepl("common", tax_names[,"V7"]),]
  tax_names[,"V3"] %<>% gsub("\\[|\\]","",.)
  
  tax_auth <- tax_names[tax_names[,"V7"]=="authority",]
  xtax <- match(xtax[,"species"], tax_auth[,"V3"] %>% gsub("([A-Z]{1,}[a-z-]{1,} [,a-z-]{2,}) .*|([A-Z]{1,}[a-z-]{1,} [x]{1,1} [,a-z-]{2,}) .*","\\1\\2",.)) %>%
    tax_auth[.,"V3"] %>% cbind.data.frame(xtax,"scientificName"=.)
  
  xtax <- xtax[match(meta3[,"xref"] %>% gsub("taxon:","",.) %>% as.integer, xtax[,"orig_xref"]),]  
  
  meta3[,"authority_ncbi"] <- xtax[,"scientificName"]
  meta3[meta3[,"authority_ncbi"] %>% is.na,"authority_ncbi"] <- xtax[meta3[,"authority_ncbi"] %>% is.na,"species"]

  #fill the gaps with incorrectly annotated species names in NCBI
  temp_target <- meta3[,"authority_ncbi"] %>% is.na
  temp_replacement <- meta3[temp_target,c("orig_title")] %>% gsub("^([A-Z]{1,1}[a-z]{1,} [a-z-]{2,}) .*","\\1",.)
  
  temp_replacement[grepl(" cv| aff|hybrid| sp | sp\\.", temp_replacement)] <- ""
  temp_replacement[temp_replacement %>% gsub("[^ ]","",.) %>% nchar > 1] <- ""
  
  #remove all entries with no proper name detected 
  meta3[temp_target,"authority_ncbi"] <- temp_replacement
  xtax <- xtax[!(meta3[,"authority_ncbi"] %>% is.na | meta3[,"authority_ncbi"] ==""),]
  meta3 <- meta3[!(meta3[,"authority_ncbi"] %>% is.na | meta3[,"authority_ncbi"] ==""),]

  #WARNING for some reason some species names contain backslashes, remove those together with single quotation marks and special characters. 
  meta3[,"authority_ncbi"] %<>% gsub("[^-.()&, 0-9A-Za-z]|  ","",.) 
  
  write.table(x=meta3 , "./step3/step3_meta_edit.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")   
  
  temp_query <- cbind.data.frame("name"=meta3[,"authority_ncbi"], rank="SPECIES",
                                 "family"=xtax[,"family"], "order"=xtax[,"order"],
                                 "genus"= meta3[,"authority_ncbi"] %>% gsub(" .*","",.),
                                 "kingdom" = "PLANTAE")
  temp_query <- temp_query[temp_query %>% duplicated %>% not,]
  
  closeAllConnections()
  
  cl <- parallel::makeCluster(9, type = "PSOCK")
  
  doParallel::registerDoParallel(cl)
  if (foreach::getDoParRegistered() == FALSE) stop("CLuster could not be registered.")
  
  download_stop <- ceiling(temp_query %>% nrow / 1000)
  
  temp_tax <- foreach (j = 1:download_stop, .inorder = FALSE, .combine = rbind, .packages=c("rgbif")) %dopar%
    {
      options(scipen = 999)
      Sys.sleep(sample(111:333/100, size=1))
      
      tryCatch(
        expr = {
          temp <- name_backbone_checklist(name_data = temp_query[(j * 1000 - 1000 + 1):(j * 1000),], verbose=TRUE)
        },
        error = function(e){ 
          Sys.sleep(300)
          temp <- name_backbone_checklist(name_data = temp_query[(j * 1000 - 1000 + 1):(j * 1000),], verbose=TRUE)
        }
      )

      temp
    }
  
  stopCluster(cl)
  closeAllConnections()
  
  write.table(x=temp_tax %>% as.data.frame, file = "./step3/tax_results.csv",
              quote=FALSE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")   
  
  #temp_tax <- read.table(file = "./step3/tax_results.csv", sep='\t', dec=".",
  #           fileEncoding = "UTF-8", header=TRUE, quote = "")
  
  temp_tax <- temp_tax[temp_tax[,"kingdom"] != "Animalia",]
  temp_tax <- temp_tax[temp_tax[,"confidence"] >= 95,]
  temp_tax <- temp_tax[temp_tax[,"matchType"] != "NONE",]
  temp_tax <- temp_tax[temp_tax[,"matchType"] != "HIGHERRANK",]
  temp_tax <- temp_tax[temp_tax[,"verbatim_name"] %>% duplicated %>% not,] %>% data.frame
  
  meta3 <- read.table("./step3/step3_meta_edit.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE)   
  
  temp_select <- match(meta3[,"authority_ncbi"], temp_tax[,"verbatim_name"])
  meta3[,c("phylum","order","family","genus","species","speciesKey","scientificName")] <- temp_tax[temp_select,c("phylum","order","family","genus","species","speciesKey","scientificName")]
  
  meta3 <- meta3[!meta3[,"species"] %>% is.na,]
  meta3 <- meta3[meta3[,"phylum"] == "Tracheophyta",]
  
  #remove accidental hybrid characters in species names
  meta3[,"species"] <- meta3[,"species"] %>% sub(" ",".",.) %>% gsub(" .*","",.) %>% gsub("\\."," ",.)
  
  write.table(x=meta3 , file = "./step3/step3_meta_tax.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")   
  
  if ("step4" %in% dir()) unlink("step4", force = TRUE, recursive = TRUE)
  dir.create("step4")
  
  file.copy(from="./step3/step3_meta_tax.csv", to= paste0("./step4/step4_meta.csv"))
  if (file.exists("./step3/step3_meta_tax.csv")) temp_silent <- file.remove("./step3/step3_meta_tax.csv")
  
}

seq_dereplicate <- function(max_seqspec=10)
{
  meta4 <- read.table("./step4/step4_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  
  temp_abund <- table(meta4[,"scientificName"])
  temp_over <- temp_abund[temp_abund > 10]
  temp_dropme <- NULL
  
  for (i in 1:length(temp_over))
  {
    temp <- meta4[meta4[,"scientificName"] %in% ( temp_over %>% names %>% .[i]),]
    temp <- temp[order(sample(temp[,"acc"])),] #randomize
    temp <- temp[order(
      temp[,c("ITS1_complete","ITS2_complete")] %>% rowSums,
      temp[,"ITS1_complete"],
      temp[,"ITS2_complete"],
      temp[,"acc"] %>% gsub("^([A-Z]{1,2}[0-9]{1,1}).*","\\1",.) %>% duplicated %>% not,
      temp[,"ITS1_found"],
      temp[,"ITS2_found"],
      temp[,"X58S_complete"],
      decreasing = T),]
    
    temp_ITS1_found <- temp[,"ITS1_found"] %>% cumsum
    temp_ITS2_found <- temp[,"ITS2_found"] %>% cumsum
    
    temp_keep <- ( temp_ITS2_found <= max_seqspec & temp_ITS2_found %>% duplicated %>% not & temp_ITS2_found > 0 ) |
      ( temp_ITS1_found <= max_seqspec & temp_ITS1_found %>% duplicated %>% not & temp_ITS1_found > 0 )
    
    temp_dropme <- c(temp_dropme, temp[!temp_keep,"acc"])
  }
  
  paste0(temp_dropme %>% length, " sequences dropped due to >",max_seqspec," ITS1/2 sequences per species" ) %>% print
  
  meta4 <- meta4[meta4[,"acc"] %in% temp_dropme %>% not,]
  
  if ("step5" %in% dir()) unlink("step5", force = TRUE, recursive = TRUE)
  dir.create("step5")
  
  write.table(x=meta4, file = "./step5/step5_meta.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")   
  
}

remove_fungi <- function(whichITS=1)
{
  if(sum(whichITS!=c(0,1,2))>2) return("bad arg")
  
  if(file.exists("./step5/step5_meta_edit.csv"))
  {
    meta5 <- read.table("./step5/step5_meta_edit.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  } else { 
    meta5 <- read.table("./step5/step5_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  }
  
  if(whichITS==1)
    temp_select <- meta5[,"ITS1_found"]
  if(whichITS==2)
    temp_select <- meta5[,"ITS2_found"]
  if(whichITS==0)
    temp_select <- meta5[,"X58S_complete"] & meta5[,"ITS1_found"] & meta5[,"ITS2_found"]
  
  temp_names <- paste0(meta5[temp_select,"acc"],";tax=p:",meta5[temp_select,"phylum",],"_o:",
                       meta5[temp_select,"order",], "_f:", meta5[temp_select,"family",],
                       "_g:", meta5[temp_select,"genus",], "_s:", meta5[temp_select,"species",] %>% gsub(" ",".",.))
  
  temp_meta <- meta5[temp_select,]
  
  if(whichITS==1)
  {
    temp_out <- temp_meta[,"ITS1_seq"] %>% DNAStringSet
    temp_fungi_dat <- "./resources/Fungi_final_addins_ITS1_extended.fasta"
  }
  if(whichITS==2)
  {
    temp_out <- temp_meta[,"ITS2_seq"] %>% DNAStringSet
    temp_fungi_dat <- "./resources/Fungi_final_addins_ITS2_extended.fasta"
  }
  if(whichITS==0)
  {
    temp_out <- paste0(temp_meta[,"ITS1_seq"], temp_meta[,"X58S_seq"], temp_meta[,"ITS2_seq"]) %>% DNAStringSet
    temp_fungi_dat <- "./resources/Fungi_final_addins_ITSfull_extended.fasta"
  }
  
  names(temp_out) <- temp_names
  
  writeXStringSet(temp_out, filepath = "./step5/temp.fasta")
  
  if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
  temp <- system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./step5/temp.fasta --db ",temp_fungi_dat," --id 0.1 --threads ",number_of_cores,
                                  " --userout ./step5/temp_fungi_score.csv --iddef 2 --userfields query+target+id+alnlen+mism+opens+qilo+qihi+tilo+tihi+ql+tl+qrow"), 
                 invisible = TRUE, intern = TRUE)
  
  temp_its_check <- read.csv("./step5/temp_fungi_score.csv", sep="\t", header = F)
  
  temp_its_check <- temp_its_check[temp_its_check[,"V4"]>100,]
  
  #temp_its_check[,"V3"] %>% sort %>% plot
  
  temp_thresh <- temp_its_check[,"V3"] > ( temp_its_check[,"V3"] %>% median ) + (temp_its_check[,"V3"] %>% IQR * 4)
  #temp_its_check[temp_thresh,] %>% .[,1] %>% gsub("\\;.*","",.) %>% cbind %>% print(quote=F)
  
  temp_kickme <- temp_its_check[temp_thresh,] %>% .[,1] %>% gsub("\\;.*","",.) 
  print(paste0(length(temp_kickme)," fungal sequences removed"))
  
  meta5 <- meta5[meta5[,"acc"] %in% temp_kickme %>% not,] 
  
  write.table(x=meta5 , "./step5/step5_meta_edit.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")
  
}

add_spike_ins <- function()
{
  if ("step6" %in% dir()) unlink("step6", force = TRUE, recursive = TRUE)
  dir.create("step6")

  delete_temporary_files <<- delete_temporary_files
  
  if (delete_temporary_files) 
  {
    if ("step2" %in% dir()) unlink("step2", force = TRUE, recursive = TRUE)
    if ("step3" %in% dir()) unlink("step3", force = TRUE, recursive = TRUE)
    if ("step4" %in% dir()) unlink("step4", force = TRUE, recursive = TRUE)
  }
  
  
  temp_readme <- list.files(path="./resources", pattern="ITSfull", full.names = T)
  temp_input <- NULL
  
  for (i in 1:length(temp_readme))
  temp_input <- readDNAStringSet(filepath = temp_readme[i]) %>% append(temp_input, .)
  
  if (file.exists("./step6/spike_ins.fasta")) temp_silent <- file.remove("./step6/spike_ins.fasta")
  writeXStringSet(temp_input, filepath = "./step6/spike_ins.fasta", compress = F, append = F)
  
  # file.copy(from="./step6/spike_ins.fasta", to= paste0(itsx_directory,"spike_ins.fasta"))
  # system(command = paste0(itsx_directory,"ITSx", args=paste0(" -i ",itsx_directory,"spike_ins.fasta -o spike_ins --reset T",
  #                                                            " --multi_thread T --cpu ",number_of_cores," --complement F --save_regions none",
  #                                                            " --preserve T --fasta F --graphical F --summary F")))
  # if (file.exists(paste0(itsx_directory,"spike_ins.fasta"))) file.remove(paste0(itsx_directory,"spike_ins.fasta"))
  # file.copy(from= paste0(itsx_directory,"spike_ins.positions.txt"), to="./step6/spike_ins.positions.txt")
  # 
  # list.files(path=itsx_directory, pattern="spike_ins", full.names = T) %>% file.remove(.)
  
  metas <- cbind("acc"=temp_input %>% names %>% gsub("_.*","",.), "orig_title"= temp_input %>% names %>% as.character)
  write.table(x=metas , "./step6/spikein_meta.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8") 
  
  file.copy(from="./step6/spike_ins.fasta", to= paste0(itsx_directory,"spike_ins.fasta"))
  system(command = paste0(itsx_directory,"ITSx", args=paste0(" -i ",itsx_directory,"spike_ins.fasta -o spike_ins --reset T",
                                                             " --multi_thread T --cpu ",number_of_cores," --complement F --save_regions none",
                                                             " --preserve T --fasta F --graphical F --summary F")))
  
  if (file.exists(paste0(itsx_directory,"spike_ins.fasta.fasta"))) 
    temp_silent <- file.remove(paste0(itsx_directory,"spike_ins.fasta.fasta"))
  
  if (file.exists("./step6/spike_ins.positions.txt")) temp_silent <- file.remove("./step6/spike_ins.positions.txt")
  file.copy(from="./spike_ins.positions.txt", to= "./step6/spike_ins.positions.txt")
  list.files(path=".", pattern="spike_ins", full.names = T) %>% file.remove()
  
  filter_by_ITSx(input_folder = "./step6/", output_folder = "./step7/", its_minlen, meta_in = "spikein_meta.csv",
                 meta_out = "spikein2_meta.csv", itsx_file = "spike_ins.positions.txt", dna_file = "spike_ins.fasta",
                 min_58S_length=min_58S_length, max_58S_length=max_58S_length)

  # start here for manual update without ITSx!
  file.copy(from = "./step5/step5_meta_edit.csv", to = "./step7/step7_meta.csv", overwrite = T)
  
  metas <- read.table(file="./step7/spikein2_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE) 
  
  temp_spec <- metas[,"orig_title"] %>% gsub(".*_g:","",.) %>% gsub("\\."," ",.)
  temp_select <- temp_spec %>% grepl(" ",.) %>% not
  temp_spec[temp_select] <- temp_spec[temp_select] %>% gsub("_s:"," ",.)
  temp_gen <- metas[,"orig_title"] %>% gsub("^.*_g:(.*)_s:.*","\\1",.) %>% gsub("^(.*sedis$)|^NA","",.)
  temp_fam <- metas[,"orig_title"] %>% gsub("^.*_f:(.*)_g:.*","\\1",.) %>% gsub("^(.*sedis$)|^NA","",.)
  temp_ord <- metas[,"orig_title"] %>% gsub("^.*_o:(.*)_f:.*","\\1",.) %>% gsub("^(.*sedis$)|^NA","",.)
  temp_class <- metas[,"orig_title"] %>% gsub("^.*_c:(.*)_o:.*","\\1",.) %>% gsub("^(.*sedis$)|^NA","",.)
  temp_spec <- temp_spec %>% gsub(".*s:","",.) %>% gsub(" sp$","",.) %>% gsub("^(.*sedis$)|^NA","",.)
  temp_query <- cbind.data.frame("name"=temp_spec, "genus"=temp_gen, "family"=temp_fam, "order"=temp_ord )

  
  closeAllConnections()
  
  cl <- parallel::makeCluster(6, type = "PSOCK")
  
  doParallel::registerDoParallel(cl)
  if (foreach::getDoParRegistered() == FALSE) stop("CLuster could not be registered.")
  
  download_stop <- ceiling(temp_query %>% nrow / 1000)
  
  temp_tax <- foreach (j = 1:download_stop, .inorder = FALSE, .combine = rbind, .packages=c("rgbif")) %dopar%
    {
      options(scipen = 999)
      Sys.sleep(sample(111:333/100, size=1))
      temp <- name_backbone_checklist(name_data = temp_query[(j * 1000 - 1000 + 1):(j * 1000),], verbose=TRUE)
      temp
    }
  
  stopCluster(cl)
  closeAllConnections()
  
  write.table(x=temp_tax , file = "./step7/tax_results.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8")   
  temp_tax <- read.table(file = "./step7/tax_results.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header=T)  
  
  temp_tax <- temp_tax[temp_tax[,"verbatim_name"] %>% duplicated %>% not,] %>% data.frame
  
  metas <- read.table(file="./step7/spikein2_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE) 
  
  temp_select <- match(temp_query[,"name"], temp_tax[,"verbatim_name"])
  if (temp_select %>% length != metas %>% nrow) print("tax assignment error")
  metas[,c("kingdom","phylum","order","family","genus","species","speciesKey","scientificName")] <- temp_tax[temp_select,c("kingdom", "phylum","order","family","genus","species","usageKey","scientificName")]
  
  metas <- metas[metas[,"kingdom"] %>% is.na %>% not,]
  metas <- metas[metas[,"family"] %>% is.na %>% not,]
  
  metas[metas[,"genus"] %>% is.na,"genus"] <- paste0("indet",metas[metas[,"genus"] %>% is.na,"family"])
  metas[metas[,"species"] %>% is.na,"species"] <- paste0("indet", " ",metas[metas[,"species"] %>% is.na,"genus"] %>% gsub("indet","",.) ) 

  meta7 <- read.table(file="./step7/step7_meta.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE) 
  meta7 <- cbind.data.frame(meta7, "kingdom"="Plantae")
  
  metas[,colnames(meta7)[colnames(meta7) %in% colnames(metas) %>% not]] <- NA
  
  meta7 <- rbind.data.frame(meta7, metas)

  if ("step8" %in% dir()) unlink("step8", force = TRUE, recursive = TRUE)
  dir.create("step8")
  
  write.table(x=meta7 , file = "./step8/step8_meta.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8") 
  
}

dist_clean <- function(whichITS=0, family_forced, max_intraspecific, max_global, max_IUPAC)
{
  if(sum(whichITS != c(0,1,2,5)) > 3 ) return("bad arg")
  
  if(file.exists("./step8/step8_meta_edit.csv"))
  {
    meta8 <- read.table("./step8/step8_meta_edit.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  } else { 
    meta8 <- read.table("./step8/step8_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  }

  if(whichITS==1)
    temp_select <- meta8[,"ITS1_found"]
  if(whichITS==2)
    temp_select <- meta8[,"ITS2_found"]
  if(whichITS==0)
    temp_select <- meta8[,"X58S_complete"] & meta8[,"ITS1_found"] & meta8[,"ITS2_found"]
  if(whichITS==5)
    temp_select <- meta8[,"X58S_complete"] 
  
  
  temp_names <- paste0(meta8[temp_select,"acc"],";tax=k:",meta8[temp_select,"kingdom",] ,"_p:",meta8[temp_select,"phylum",],"_o:",
                       meta8[temp_select,"order",], "_f:", meta8[temp_select,"family",],
                       "_g:", meta8[temp_select,"genus",], "_s:", meta8[temp_select,"species",] %>% gsub(" ",".",.))
  
  temp_meta <- meta8[temp_select,]
  
  if(whichITS==1)
    temp_out <- temp_meta[,"ITS1_seq"] %>% DNAStringSet
  if(whichITS==2)
    temp_out <- temp_meta[,"ITS2_seq"] %>% DNAStringSet
  if(whichITS==0)
    temp_out <- paste0(temp_meta[,"ITS1_seq"], temp_meta[,"X58S_seq"], temp_meta[,"ITS2_seq"]) %>% DNAStringSet
  if(whichITS==5)
    temp_out <- temp_meta[,"X58S_seq"] %>% DNAStringSet
  
  names(temp_out) <- temp_names
  
  writeXStringSet(temp_out, filepath = "./step8/temp.fasta")
  
  if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
  temp_silent <- system(command = paste0(vsearch_directory, vsearch_exec, " --makeudb_usearch ./step8/temp.fasta --output ./step8/temp.fasta.udb"))
  
  temp <- system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./step8/temp.fasta --self -db ./step8/temp.fasta.udb --id 0.1 --threads ",number_of_cores,
                                  " --userout ./step8/temp_score.csv --iddef 3 --userfields query+target+id+alnlen+mism+opens+qilo+qihi+tilo+tihi+ql+tl+qrow"))
  
  temp_its_check <- read.csv("./step8/temp_score.csv", sep="\t", header = F)
  
  #don't check spike-ins from fungi and other organisms
  temp_its_check <- temp_its_check[grepl("p:Tracheophyta",temp_its_check[,"V1"] ),]
  
  #global distance treshold
  temp_select <-  temp_its_check[,"V3"] < ( 100 - max_global ) 
  
  kickme3 <- NULL
  kickme3 <- c(temp_its_check[temp_select,"V1"] %>% gsub(";.*","",.), kickme3)
  
  paste0(kickme3 %>% length, " sequences tagged due to a max. global distance threshold of ", max_global, "%") %>% print
  
  
  kickme1 <- NULL
  
  if (family_forced == TRUE )
  {
    #should match own family if not only sequence in family or only sequence in genus
    
    temp_select <- ( temp_its_check[,"V1"] %>% gsub("^.*_f:(.*)_g:.*","\\1",.) ) != ( temp_its_check[,"V2"] %>% gsub("^.*_f:(.*)_g:.*","\\1",.) )
    temp_select <- (temp_out %>% names) %in% c(temp_its_check[temp_select,"V1"], temp_its_check[temp_select,"V2"]) 
    
    writeXStringSet(temp_out[temp_select], filepath = "./step8/temp_extracted.fasta", append = F)
    writeXStringSet(temp_out[!temp_select], filepath = "./step8/temp_suspicious_removed.fasta", append = F)
    
    if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
    temp <- system(command = paste0(vsearch_directory, vsearch_exec , " --usearch_global ./step8/temp_extracted.fasta --self -db ./step8/temp_suspicious_removed.fasta --id 0.1 --threads ",number_of_cores,
                                    " --userout ./step8/temp_score_fam.csv --iddef 3 --userfields query+target+id+alnlen+mism+opens+qilo+qihi+tilo+tihi+ql+tl+qrow"))
    
    temp_its_check_fam <- read.csv("./step8/temp_score_fam.csv", sep="\t", header = F)
    temp_its_check_fam <- temp_its_check_fam[grepl("p:Tracheophyta",temp_its_check_fam[,"V1"] ),]
    
    temp_select <- ( temp_its_check_fam[,"V1"] %>% gsub("^.*_f:(.*)_g:.*","\\1",.) ) != ( temp_its_check_fam[,"V2"] %>% gsub("^.*_f:(.*)_g:.*","\\1",.) )
    
    temp_its_check_fam <- temp_its_check_fam[temp_select,]
    
    #family still left in reference database?
    temp_reference <- readDNAStringSet(filepath = "./step8/temp_suspicious_removed.fasta")
    
    temp_no_ref_family <- ( temp_its_check_fam[,"V1"] %>% gsub("^.*_f:(.*)_g:.*","\\1",.) ) %in% ( temp_reference %>% names %>% gsub("^.*_f:(.*)_g:.*","\\1",.) )
    
    temp_no_ref_family <- temp_its_check_fam[temp_no_ref_family,]
    
    kickme1 <- temp_its_check_fam[,"V1"] %>% gsub(";.*","",.) %>% unique
    
    paste0(kickme1 %>% length, " sequences tagged due to a family mismatch") %>% print
    
  }
  

  if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
  temp_silent <- system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./step8/temp.fasta --self -db ./step8/temp.fasta.udb --id 0.65 --threads ",
                                         number_of_cores, " --userout ./step8/temp_score.csv --iddef 3 --maxaccepts 50 --userfields query+target+id+alnlen"))
  
  temp_its_check <- read.csv("./step8/temp_score.csv", sep="\t", header = F)

  #don't check spike-ins from fungi and other organisms
  temp_its_check <- temp_its_check[grepl("p:Tracheophyta",temp_its_check[,"V1"] ),]
  
  kickme2 <- NULL

  #create a helper data.frame of top-hits only
  temp_tophits <- temp_its_check[temp_its_check[,1] %>% duplicated %>% not,]
  
  #limit to only hits to same species
  temp_reduced <- temp_its_check[(temp_its_check[,"V1"] %>% gsub("^.*_s:","",.) ) == ( temp_its_check[,"V2"] %>% gsub("^.*_s:","",.) ),]
  #which of the sequences hit with an median of below intraspecific filter
  temp_median <- temp_reduced %>% group_by(V1) %>% summarize(median(V3)) %>% data.frame

  
  
  #select only hits which have an median hit percentage similarity below the threshold
  temp_tophits_select <- temp_tophits[match(temp_median[temp_median[,2] < 100 - max_intraspecific,1], temp_tophits[,1]),]

  #select from those only sequences where the top hit is not the same species
  kickme2 <- temp_tophits_select[
     ( temp_tophits_select[,"V1"] %>% gsub("^.*_s:","",.) ) != ( temp_tophits_select[,"V2"] %>% gsub("^.*_s:","",.) )
     ,"V1"]
  
  ##remove: the sequence is more than the threshold away from the same species and the tophit is not the same species
  #keep only sequences where the top hit is not to the same species
  
  select_mismatches <- temp_tophits[( temp_tophits[,"V1"] %>% gsub("^.*_s:","",.) ) != ( temp_tophits[,"V2"] %>% gsub("^.*_s:","",.) ),"V1"]
  temp_mism <- temp_reduced[temp_reduced[,"V1"] %in% select_mismatches,] 
  temp_mism <- temp_mism[temp_mism[,"V1"] %>% duplicated %>% not,]
  temp_mism <- temp_mism[temp_mism[,"V3"] < 100 - max_intraspecific,]
  #temp_mism[,"V1"] %>% gsub(".*f:(.*)_s:.*","\\1",.) %>% table %>% sort %>% rev %>% head(.,n=50)
   
  kickme2 <- c(kickme2, temp_mism[,"V1"] %>% gsub("^.*_s:","",.)) 
  kickme2 %<>% unique 
   
 #meta8[meta8[,"acc"] %in% kickme,"species"] %>% table %>% sort %>% rev %>% head(.,n=20)
  
  paste0(kickme2 %>% length, " sequences tagged due to a intraspecies distance threshold of ", max_intraspecific, "%") %>% print
  
  temp_reduced <- temp_out[( temp_out %>% names %>% gsub(";.*","",.) )  %in% ( temp_its_check[,"V1"] %>% gsub(";.*","",.) )]
  temp_freq <- temp_reduced %>% alphabetFrequency(., as.prob=T)
  
  kickme4 <- temp_reduced[temp_freq[,5:15] %>% rowSums > max_IUPAC] %>% names %>% gsub(";.*","",.)
  paste0(kickme4 %>% length, " sequences tagged due to a max. allowed IUPAC content of ", max_IUPAC*100, "%") %>% print
  
  
  
  kickme <- c(kickme1, kickme2, kickme3, kickme4) %>% unique
  #meta8[meta8[,"acc"] %in% kickme,"species"] %>% table %>% sort %>% rev %>% head(.,n=20)
  
  paste0(kickme %>% length, " sequences removed, total") %>% print
  
  
  meta8 <- meta8[meta8[,"acc"] %in% kickme %>% not,] 
  
  write.table(x=meta8 , file = "./step8/step8_meta_edit.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8") 
  
  
}

migrate_to_step9 <- function()
{
  if ("step9" %in% dir()) unlink("step9", force = TRUE, recursive = TRUE)
  if ("step9" %in% dir() %>% not) dir.create("step9")
  
  if(file.exists("./step8/step8_meta_edit.csv"))
  {
    file.copy(from="./step8/step8_meta_edit.csv", to= paste0("./step9/step9_meta.csv"))
  } else { 
    file.copy(from="./step8/step8_meta.csv", to= paste0("./step9/step9_meta.csv"))
  }
  
  delete_temporary_files <<- delete_temporary_files
  
  if (delete_temporary_files) 
  {
    if ("step5" %in% dir()) unlink("step5", force = TRUE, recursive = TRUE)
    if ("step6" %in% dir()) unlink("step6", force = TRUE, recursive = TRUE)
    if ("step7" %in% dir()) unlink("step7", force = TRUE, recursive = TRUE)
  }
  
  
}

get_distribution <- function()
{

  #getting the taxonomic information takes approximately 1 hour and 20 minutes
  
  meta9 <- read.table("./step9/step9_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  
  temp_query <- meta9[,"speciesKey"] %>% unique
  
  closeAllConnections()
  
  cl <- parallel::makeCluster(9, type = "PSOCK")
  
  doParallel::registerDoParallel(cl)
  if (foreach::getDoParRegistered() == FALSE) stop("CLuster could not be registered.")
  
  occ_count <- NULL
  
  occ_count  <- foreach (j = 1:length(temp_query), .inorder = TRUE, .combine = append, .packages=c("rgbif")) %dopar%
    {
      options(scipen = 999)
      Sys.sleep(1/10)

      tryCatch(
        expr = {
          temp <- occ_count(taxonKey= temp_query[j], facet='country', facetLimit=999) 
        },
        error = function(e){ 
          Sys.sleep(300)
          temp <- occ_count(taxonKey= temp_query[j], facet='country', facetLimit=999) 
        }
      )

      list(temp)
    }
  
  stopCluster(cl)
  closeAllConnections()
  
  save <- occ_count 
  names(occ_count) <- temp_query

  occ_df <- occ_count %>% lapply(., as.data.frame) %>% do.call(rbind,.) 
  occ_df <- cbind.data.frame(occ_df, taxid = occ_df %>% rownames %>% gsub("\\..*","",.))
  
  #eliminate ZZ (not designated) and counts below 3, summarize on continent level
  occ_df <- occ_df[occ_df[,"country"]!="ZZ",]
  occ_df <- occ_df[occ_df[,"count"]>2,]
  
  #replace ISO2 with ISO3 codes
  temp_iso <- enumeration_country()
  occ_df[,"iso3"] <- temp_iso[match(occ_df[,"country"], temp_iso[,"iso2"]),"iso3"]
  
  # library(countrycode)
  # temp <- cbind.data.frame(isocodes,
  #                          UNsubregion=countrycode(sourcevar = isocodes[, "name"], origin = "country.name", destination = "un.regionsub.name"),
  #                          UNinterregion=countrycode(sourcevar = isocodes[, "name"], origin = "country.name", destination = "un.regionintermediate.name")
  #                          )
  # write.table(x=temp , file = "./resources/countrycodes.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8") 
  #fixed Kosovo and Curacao spelling and added region for Antartica
  
  temp_country <- read.table("./resources/countrycodes.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  temp_country[temp_country[,"code"] %>% is.na,"code"] <- "NA" #fix Namibia
  temp_country[,"iso3"] <- temp_iso[match(temp_country[,"code"], temp_iso[,"iso2"]),"iso3"]
  
  temp_select <- (temp_country[,"UNsubregion"]=="Latin America and the Caribbean" | 
                    temp_country[,"UNsubregion"]=="Sub-Saharan Africa" ) & 
                    temp_country[,"UNinterregion"] %>% is.na %>% not
  
  temp_country[temp_select,"UNsubregion"] <- temp_country[temp_select,"UNinterregion"]
  
  temp_select <- match(occ_df[,"country"], temp_country[,"code"])
  temp_region_name <- temp_country[temp_select,"UNsubregion"]
  
  temp_region_name  %<>% gsub("Eastern ","E-",.) %>% gsub("Northern ","N-",.) %>% gsub("South-eastern ","SE-",.) %>%
    gsub("South ","S-",.) %>% gsub("Southern ","S-",.) %>% gsub("Western ","W-",.) %>% 
    gsub("Central ","Cent. ",.) %>% gsub("Middle ","Cent. ",.) %>% gsub("Australia and New Zealand","AUS/NZ",.) 
 
  occ_df <- cbind.data.frame(occ_df, region = temp_region_name)
  rownames(occ_df) <- NULL
  
  occ_df <- cbind.data.frame(occ_df, region_naked = occ_df[,"region"] %>% gsub(".* |.*-","",.) )
  
  region_agglomerate <- occ_df %>% group_by(taxid, region, region_naked) %>%
    arrange(region_naked, region, iso3, taxid) %>% 
    reframe(joined=paste(iso3,count, sep=" ", collapse = "")) %>%
    group_by(taxid) %>% arrange(region_naked) %>%
    reframe(joined=paste(region,joined, sep="(", collapse = ")")) %>% data.frame
  
  region_agglomerate[,2] %<>% gsub("(?<=[0-9])([A-Z]{1,})",", \\1",.,perl=T) %>%
    gsub(" ([0-9])",": \\1",.) %>% gsub("\\("," \\(",.) %>% gsub("\\)","\\), ",.) %>% 
    gsub("([0-9]{1,})$","\\1)",.)
  
  colnames(meta9)[meta9 %>% colnames %in% "speciesKey"] <- "taxid_gbif"
  colnames(meta9)[meta9 %>% colnames %in% "country"] <- "country_ncbi"
  colnames(meta9)[meta9 %>% colnames %in% "xref"] <- "taxid_ncbi"
  colnames(meta9)[meta9 %>% colnames %in% "orig_taxonomy"] <- "taxonomy_ncbi"
  colnames(meta9)[meta9 %>% colnames %in% "orig_title"] <- "title_ncbi"
    
  meta9 <- cbind.data.frame(meta9, "region"=match(meta9[,"taxid_gbif"], region_agglomerate[,"taxid"]) %>% region_agglomerate[.,2]) 
  
  meta9 <- meta9[,match(
         c("acc", "taxid_gbif", "voucher", "kingdom", "phylum", "order", "family", "genus", "species", "scientificName", "region", "ITS1_found", "ITS1_length",
           "X58S_found", "X58S_length", "ITS2_found", "ITS2_length" , "ITS1_complete", "ITS1_seq", "X58S_complete", "X58S_seq", "ITS2_complete", "ITS2_seq",
           "LSU_found", "LSU_length", "LSU_seq", "SSU_found", "SSU_length", "SSU_seq", "date" ,"collection", "auth", "publication", "title_ncbi", "taxid_ncbi",
           "taxonomy_ncbi", "authority_ncbi", "journal", "country_ncbi", "collectedby"), colnames(meta9)  
  )]                
  
  if ("step10" %in% dir()) unlink("step10", force = TRUE, recursive = TRUE)
  if ("step10" %in% dir() %>% not) dir.create("step10")
  
  write.table(x=meta9 , file = "./step10/step10_meta.csv", quote = TRUE, row.names = FALSE, sep='\t', dec=".", append=FALSE, fileEncoding = "UTF-8") 
  
  delete_temporary_files <<- delete_temporary_files
  
  if (delete_temporary_files) 
  {
    if ("step8" %in% dir()) unlink("step8", force = TRUE, recursive = TRUE)
    if ("step9" %in% dir()) unlink("step9", force = TRUE, recursive = TRUE)
  }
  
  if(file.exists("./step8/temp.fasta.udb")) temp_silent <- file.remove("./step8/temp.fasta.udb")
  
}


export_seq <- function(include_GBIF_taxID=F, whichITS=0, style="SINTAX", out_dir_name)
{
 
  if (out_dir_name %in% dir() %>% not) dir.create(path=out_dir_name)
  
  meta10 <- read.table("./step10/step10_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  
  if(whichITS==1)
    temp_select <- meta10[,"ITS1_found"]
  if(whichITS==2)
    temp_select <- meta10[,"ITS2_found"]
  if(whichITS==0)
    temp_select <- meta10[,"X58S_complete"] & meta10[,"ITS1_found"] & meta10[,"ITS2_found"]
  
 if (style=="SINTAX" & include_GBIF_taxID==F)
  temp_names <- paste0(meta10[temp_select,"acc"],";tax=p:",meta10[temp_select,"phylum",],",o:",
                       meta10[temp_select,"order",], ",f:", meta10[temp_select,"family",],
                       ",g:", meta10[temp_select,"genus",],
                       ",s:", meta10[temp_select,"species",] %>% gsub(" ","_",.)) %>%
    gsub(",[a-z]{1,1}_NA.*","",.)
  
  if (style=="SINTAX" & include_GBIF_taxID==T)
    temp_names <- paste0(meta10[temp_select,"acc"],"_",meta10[temp_select,"taxid_gbif"],";tax=p:",meta10[temp_select,"phylum",],",o:",
                         meta10[temp_select,"order",], ",f:", meta10[temp_select,"family",],
                         ",g:", meta10[temp_select,"genus",],
                         ",s:", meta10[temp_select,"species",] %>% gsub(" ","_",.))

  if (style=="RDP" & include_GBIF_taxID==F)
  {
    temp_names <- paste0(meta10[temp_select,"acc"]," p_",meta10[temp_select,"phylum",],";o_",
                         meta10[temp_select,"order",], ";f_", meta10[temp_select,"family",],
                         ";g_", meta10[temp_select,"genus",], ";s_",
                         meta10[temp_select,"species",] %>% gsub(" ","_",.))
    
     temp_tax_names <- paste0(meta10[temp_select,"acc"],"\tp_",meta10[temp_select,"phylum",],";o_",
                         meta10[temp_select,"order",], ";f_", meta10[temp_select,"family",],
                         ";g_", meta10[temp_select,"genus",], ";s_",
                         meta10[temp_select,"species",] %>% gsub(" ","_",.))
     
     temp_tax_names %<>% gsub(";[a-z]{1,1}_NA.*","",.)

     writeLines(text=temp_tax_names,
                con = paste0("./",out_dir_name,"./ITS",whichITS,"_",style,".fasta"))
   }
  
  
  if (style=="RDP" & include_GBIF_taxID==T)
  {
    temp_names <- paste0(meta10[temp_select,"acc"], "_", meta10[temp_select,"taxid_gbif"],
                         " p_",meta10[temp_select,"phylum",],";o_",
                         meta10[temp_select,"order",], ";f_", meta10[temp_select,"family",],
                         ";g_", meta10[temp_select,"genus",], ";s_",
                         meta10[temp_select,"species",] %>% gsub(" ","_",.))
    
    temp_tax_names <- paste0(meta10[temp_select,"acc"],"_", meta10[temp_select,"taxid_gbif"],
                             "\tp_",meta10[temp_select,"phylum",],";o_",
                             meta10[temp_select,"order",], ";f_", meta10[temp_select,"family",],
                             ";g_", meta10[temp_select,"genus",], ";s_",
                             meta10[temp_select,"species",] %>% gsub(" ","_",.))
    
    temp_tax_names %<>% gsub(";[a-z]{1,1}_NA.*","",.)
    
    writeLines(text=temp_tax_names,
               con = paste0("./",out_dir_name,"./ITS",which,"ITS_",style,".fasta"))
  }
  
  temp_meta <- meta10[temp_select,]
  
  if(whichITS==1)
    temp_out <- temp_meta[,"ITS1_seq"] %>% DNAStringSet
  
  if(whichITS==2)
    temp_out <- temp_meta[,"ITS2_seq"] %>% DNAStringSet
  
  if(whichITS==0)
    temp_out <- paste0(temp_meta[,"ITS1_seq"], temp_meta[,"X58S_seq"], temp_meta[,"ITS2_seq"]) %>% DNAStringSet
  
  names(temp_out) <- temp_names
  
  if(whichITS==1)
    writeXStringSet(temp_out, filepath = paste0("./",out_dir_name,"./ITS1_",style,".fasta"))
  if(whichITS==2)
    writeXStringSet(temp_out, filepath = paste0("./",out_dir_name,"./ITS2_",style,".fasta"))
  if(whichITS==0)
    writeXStringSet(temp_out, filepath = paste0("./",out_dir_name,"./ITS0_",style,".fasta"))
  
  temp_metafilename <- paste0("./",out_dir_name,"./metadata_full.csv")
  
  if (file.exists(temp_metafilename)) file.remove(temp_metafilename)
  temp_silent <- file.copy(from="./step10/step10_meta.csv", to=temp_metafilename )

  if (temp_metafilename %>% gsub("full","reduced",.) %>% file.exists)
    temp_metafilename %>% gsub("full","reduced",.) %>% file.remove
  
  meta10 %>% colnames %in% c("ITS1_found","X58S_found","ITS2_found",
                             "ITS1_complete", "ITS1_seq",
                             "X58S_complete","X58S_seq",
                             "ITS2_complete","ITS2_seq",
                             "LSU_found","LSU_length",
                             "LSU_seq","SSU_found",
                             "SSU_length","SSU_seq") %>%
    not %>% meta10[,.] %>% write.table(x=.,
                                       file = temp_metafilename %>% gsub("full","reduced",.),
                                       quote = TRUE, row.names = FALSE,
                                       sep='\t', dec=".", append=FALSE,
                                       fileEncoding = "UTF-8") 
  
}


viz_length <- function(minseq_per_family=50, width_graphic_scale=0.156,
                       height_graphic_absolute=14.4,y_axis_upper_limit=2000,
                       y_axis_number_ticks=20, x_axis_font_size=10,
                       y_axis_font_size=14, legend_offset_left=0.1)
{  
  
  meta10 <- read.table("./step10/step10_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 

  meta <- meta10[meta10[,"X58S_found"] & meta10[,"phylum"]=="Tracheophyta",]
  temp_fam_keep <- meta[,"family"] %>% table %>% names %>% .[meta[,"family"] %>% table > minseq_per_family]
  meta <- meta[meta[,"family"] %in% temp_fam_keep,]
  
  meta_fam <- meta %>% group_by(family) %>% 
    select_if(is.numeric) %>%
    summarise_all(median, na.rm = T) %>% 
    as.data.frame
  
  meta_fam <- meta_fam[,colnames(meta_fam) %in% c("family", "ITS1_length", "ITS2_length", "X58S_length")]
  
  meta_fam <- meta_fam[meta_fam[,"ITS1_length"] %>% order,]
  meta_fam <- meta_fam[complete.cases(meta_fam),]
  
  meta_dat <- meta[,colnames(meta) %in% c("family", "order", "ITS1_length", "ITS2_length", "X58S_length")]
  
  meta_dat[,"ITS0_length"] <- meta_dat[,"ITS1_length"] + meta_dat[,"ITS2_length"] + meta_dat[,"X58S_length"]
  meta_dat <- meta_dat[meta_dat[,"family"] %in% meta_fam[,"family"],]
  
  meta_dat[,"family"] <- factor(meta_dat[,"family"] %>% as.character, levels=meta_fam[,"family"], ordered=TRUE) 
  
  pdf(file = paste0("./",out_dir_name,"/its_length_minseq",minseq_per_family,".pdf"), width = width_graphic_scale*nrow(meta_fam), height = height_graphic_absolute)
  #to get legend comment out fill outside aes
  print( ggplot(dat=meta_dat,aes(x=family,group=family)
  ) + geom_violin(
    width=1,
    trim=T,
    draw_quantiles=F, #c(0.25,0.5,0.75)
    na.rm=T,
    position = "dodge",
    scale="width",
    alpha=1,
    linewidth=.5,
    color=alpha("black",.25),
    show.legend=T,
    adjust=.5,
    #fill="#00A0B0",
    aes(y=ITS1_length,fill="ITS1")
  ) + geom_violin(
    width=1,
    trim=T,
    draw_quantiles=F, #c(0.25,0.5,0.75)
    na.rm=T,
    position = "dodge",
    scale="width",
    linewidth=.5,
    alpha=1,
    color=alpha("black",.25),
    show.legend=T,
    adjust=.5,
    #fill="#00A0B0",
    aes(y=ITS0_length,fill="ITS full")
  ) + geom_violin(
    width=1,
    trim=F,
    draw_quantiles=F,
    na.rm=T,
    linewidth=.5,
    scale="width",
    alpha=.7,
    adjust=.5,
    position = "dodge",
    show.legend=T,
    color=alpha("black",.25),
    #fill="#EDC951",
    aes(y=ITS2_length,fill="ITS2")
  ) + theme(
    axis.title = element_text(face="bold",size=16),	
    axis.text.y = element_text(colour="black",size=y_axis_font_size),
    axis.text.x = element_text(angle= 90, hjust=1, vjust=.25, size=x_axis_font_size,color="black"),
    strip.text.x = element_text(angle= 90, hjust=.5, vjust=.5, size=16,color="black"),
    axis.line.y = element_line(colour="black", linetype = "solid"),
    axis.line.x = element_line(colour="black", linetype = "solid"),
    strip.background = element_rect(colour = "black", fill = NA),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(colour = alpha("black",.25)),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position = c(legend_offset_left,1),
    legend.text= element_text(colour="black",size=14),
    legend.title= element_text(face="bold",size=16),	
    legend.justification = c(1,1),
    legend.background = element_rect(colour=NA,fill=alpha("grey90", .85)),
    legend.key=element_blank(),
    panel.spacing=unit(0, "lines"),
    strip.placement="outside"
  ) + scale_y_continuous(
    breaks=scales::pretty_breaks(n=y_axis_number_ticks),
    expand = c(0,0),
    limit= c(0,y_axis_upper_limit),
    name="ITS length [bp]"
  ) + scale_x_discrete(
    expand = c(0,.55),
    name=NULL,
    labels=waiver()
  # ) + geom_vline(
  #   xintercept=seq(0.5+1,100,1),
  #   colour=alpha("black",.1),
  #   linetype=5
  ) + scale_fill_manual(
    values=c("#CC333F","#EDC951","#00A0B0"),
    name="ITS region"
  ) ) # ) from print
  dev.off()
}



krona_plot <- function()
{
  meta10 <- read.table("./step10/step10_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 
  meta10 <- meta10[meta10[,"phylum"]=="Tracheophyta",]

  plotme <- cbind.data.frame(freq=1,
                             ord.n = meta10[,"order"],
                           fam.n = meta10[,"family"],
                           gen.n= meta10[,"genus"],
                           spec.n= meta10[,"species"])
  
  plotme <- plotme %>% group_by(spec.n, gen.n, fam.n, ord.n) %>%
    summarise("freq"=sum(freq)) %>% as.data.frame

  plotme$pathString <- paste(" ",plotme$ord.n,plotme$fam.n,plotme$gen.n,plotme$spec.n,sep = "/")
  plotme <- as.Node(plotme)
  plotme$Do(function(x) x$freq <- ifelse(is.null(x$freq), 0, x$freq) + sum(Get(x$children, "freq")), traversal = "post-order")
  
  #print(plotme,"freq")
  repval <- plotme$Get("freq")
  
  plotme <- plotme %>% as.list(.,mode="simple",unname=T)
  plotme <- as_xml_document(list(root=plotme))%>%as.character

  temp <- strsplit(plotme,"\n")[[1]][2] %>% strsplit(.,"(?=<)",perl=T) %>% .[[1]]
  temp <- paste0(temp[c(T,F)],temp[c(F,T)]) %>% gsub("[0-9]","",.)
  temp <- gsub("</(.*)?>","</node>",temp) %>% gsub(">.*$",">",.)
  temp <- gsub("<(?=[^/])(.*)?>",
               '<node name="\\1"><magnitude><val>__\\1</val></magnitude><score><val>__\\1</val></score>',temp,perl=T)
  temp %<>% gsub("root"," ",.)
  repos <- grep("__",temp)
  for (i in repos)
    temp[i] <- gsub("<val>__.*?<",paste0("<val>",repval[which(repos==i)],"<"),temp[i],perl=T)
  
  col_attr <- paste0('<color attribute="score" valueStart="1" valueEnd="10" hueStart="0" hueEnd="120"></color>')
  temp <- c(col_attr, temp,"</krona></div>","</body>","</html>")  %>% unlist
  tout <- paste(temp,collapse="\n") #%>% writeLines(.,con="C:/R/test.txt")
  
  if (paste0("./",out_dir_name,"/Krona.html") %>% file.exists) paste0("./",out_dir_name,"/Krona.html") %>% file.remove
  file.copy(from="./resources/Krona1.html", to=paste0("./",out_dir_name,"/Krona.html"))

  write(tout,file=paste0("./",out_dir_name,"/Krona.html"),append=TRUE) 
  
}
















