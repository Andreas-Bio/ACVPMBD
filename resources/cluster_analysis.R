





meta8 <- read.table("./step9/step9_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8") 

whichITS <- 0

if(whichITS==1)
  temp_select <- meta8[,"ITS1_found"]
if(whichITS==2)
  temp_select <- meta8[,"ITS2_found"]
if(whichITS==0)
  temp_select <- meta8[,"X58S_complete"] & meta8[,"ITS1_found"] & meta8[,"ITS2_found"]

temp_names <- paste0(meta8[temp_select,"acc"],";tax=p:",meta8[temp_select,"phylum",],"_o:",
                     meta8[temp_select,"order",], "_f:", meta8[temp_select,"family",],
                     "_g:", meta8[temp_select,"genus",], "_s:", meta8[temp_select,"species",] %>% gsub(" ",".",.))

temp_meta <- meta8[temp_select,]

if(whichITS==1)
  temp_out <- temp_meta[,"ITS1_seq"] %>% DNAStringSet

if(whichITS==2)
  temp_out <- temp_meta[,"ITS2_seq"] %>% DNAStringSet

if(whichITS==0)
  temp_out <- paste0(temp_meta[,"ITS1_seq"], temp_meta[,"X58S_seq"], temp_meta[,"ITS2_seq"]) %>% DNAStringSet


names(temp_out) <- temp_names
temp_out <- temp_out[grep("Tracheophyta", temp_out %>% names)]

if(whichITS==1)
  writeXStringSet(temp_out, filepath = "./step10/ITS1.fasta")
if(whichITS==2)
  writeXStringSet(temp_out, filepath = "./step10/ITS2.fasta")
if(whichITS==0)
  writeXStringSet(temp_out, filepath = "./step10/ITS0.fasta")


for (j in 1:length(whichITS)) 
  system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./step10/ITS",whichITS[j],".fasta --db ./step10/ITS",whichITS[j],".fasta",
                          " --id 0.80 --threads ",number_of_cores," --maxaccepts 50 --userout ./step10/ITS",whichITS[j],"_LCA_cluster_0_hits.csv --iddef 2", 
                          " --userfields query+target+id+alnlen"))


cluster_at <- c(1, 0.995, 0.99, 0.98, 0.97, 0.96, 0.95)
whichITS <- c(0,1,2)
vsearch_info <- vector(mode="list", length=length(whichITS)*length(cluster_at))

for (i in 1:length(cluster_at))
{
  cluster_actual <- cluster_at[i]
  for (j in 1:length(whichITS)) 
  {
    ITS_actual <- whichITS[j]
    if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
    vsearch_info[[i*j]] <- system(command = paste0(vsearch_directory, vsearch_exec, " --cluster_size ./step10/ITS",ITS_actual,".fasta --id ",cluster_actual," --threads ",number_of_cores,
                                                   " --consout ./step10/ITS",ITS_actual,"_consensus_cluster_",cluster_actual,".fasta --uc ./step10/ITS",ITS_actual,
                                                   "_table_cluster_",cluster_actual,".csv --iddef 2 --clusterout_id"))
  }
}

cluster_file <- list.files(path = "./step10", pattern = "ITS._table_cluster_", full.names = T)

for (i in 1:length(cluster_file))
{
  clust_level <- cluster_file[i] %>% gsub(".*cluster_(.*).csv","\\1",.) %>% as.numeric
  if (clust_level==0) next
  temp <- cluster_file[i] %>% read.table(.,sep="\t", header=F)
  temp <- temp[grep("Tracheophyta", temp[,"V9"]),]
  temp_c <- temp[temp[,"V1"]=="C",]
  temp <- temp[temp[,"V1"]!="C",]
  clust_seeds <- temp[temp[,"V1"]=="S","V9"]
  clust_seeds <- clust_seeds [ temp_c[match(clust_seeds , temp_c[,"V9"]),"V3"] > 1] #only take sequences from clusters with 2 or more sequences in them
  
  ITS_actual <- cluster_file[i] %>% gsub(".*ITS(.*)_table.*","\\1",.)
  clust_seq <- paste0("./step10/ITS",ITS_actual,".fasta") %>% readDNAStringSet
  clust_seq <- clust_seq[names(clust_seq) %in% clust_seeds %>% not]
  writeXStringSet(clust_seq, filepath = "./step10/temp_reduced.fasta")
  
  if (.Platform$OS.type=="windows") {vsearch_exec <- "vsearch.exe"} else { vsearch_exec <- "vsearch" }
  temp_silent <- system(command = paste0(vsearch_directory, vsearch_exec, " --cluster_size ./step10/temp_reduced.fasta --id ",clust_level," --threads ",number_of_cores,
                                         " --consout ./step10/temp_cons.fasta --uc ./step10/temp_clust.csv --iddef 2 --clusterout_id"))
  
  temp2 <- read.table("./step10/temp_clust.csv",sep="\t", header=F)
  temp2 <- temp2[temp2[,"V1"]!="C",]
  temp2 <- temp2[grep("Tracheophyta", temp2[,"V9"]),]
  
  temp2[,"V2"] <- temp2[,"V2"] + ( temp[,"V2"] %>% max + 1)
  temp2[,"V9"] <- temp2[,"V9"] %>% gsub("^(.*);","\\1a:",.)
  
  seq2 <- readDNAStringSet("./step10/temp_cons.fasta")
  seq_clusterid <- names(seq2) %>% gsub(".*clusterid=","",.) %>% as.integer
  seq_clusterid <- seq_clusterid + ( temp[,"V2"] %>% max + 1)
  names(seq2) <- names(seq2) %>% gsub("clusterid=(.*)","clusterid=",.) %>% paste0(.,seq_clusterid)
  names(seq2) <- names(seq2) %>% gsub("centroid=(.*);tax","centroid=\\1a;tax",.)
  paste0("./step10/ITS",ITS_actual,"_consensusa_cluster_",clust_level,".fasta") %>% writeXStringSet(x = seq2,filepath = .)
  
  res <- rbind.data.frame(temp, temp2)
  
  cluster_file[i] %>% gsub("_table_","_tablea_",.) %>% write.csv(x=res,file=.,row.names=F)
}  


cluster_file <- list.files(path = "./step10", pattern = "ITS._tablea_cluster_", full.names = T)

for (i in 1:length(cluster_file))
{
  temp <- cluster_file[i] %>% read.table(.,sep=",", header=T)
  temp <- temp[grep("Tracheophyta", temp[,"V9"]),]
  temp <- temp[temp[,"V1"]!="C",]
  
  temp <- temp %>% group_by(V2) %>% reframe("seq"=paste0(V9, sep=" ", collapse="")) %>% data.frame
  
  cluster <- cbind.data.frame(temp, "centroid_id"=temp[,2] %>% sub(";.*","",.), "centroid"= temp[,2] %>% sub(" .*","",.))
  
  for (j in 1:nrow(cluster))
  {
    temp <- cluster[j,2] %>% as.character %>% strsplit(x=., split=" ") %>% unlist %>% gsub("^.*_f:","",.) %>% strsplit(.,"_")
    cluster[j,"same_fam"] <- temp %>% lapply(.,`[[`,1) %>% unlist %>% unique %>% length == 1
    cluster[j,"same_gen"] <- temp %>% lapply(.,`[[`,2) %>% unlist %>% unique %>% length == 1
    cluster[j,"same_spec"] <- temp %>% lapply(.,`[[`,3) %>% unlist %>% unique %>% length == 1
  }
  
  cluster[,"centroid_org"] <- cluster[,"centroid"]
  
  for (j in 1:nrow(cluster))
  {
    if (cluster[j,"same_fam"]==F)  
    {
      cluster[j,"centroid"] <- cluster[j,"centroid"] %>%
        sub(".*_o:","o:",.) %>%
        sub("_f:.*","",.) %>%
        paste0(cluster[j,"V2"],";tax=",.)
      next
    }
    
    if (cluster[j,"same_gen"]==F)  
    {
      cluster[j,"centroid"] <- cluster[j,"centroid"] %>%
        sub(".*_o:","o:",.) %>%
        sub("_g:.*","",.) %>%
        paste0(cluster[j,"V2"],";tax=",.)
      next
    }
    
    if (cluster[j,"same_spec"]==F)  
    {
      cluster[j,"centroid"] <- cluster[j,"centroid"] %>%
        sub(".*_o:","o:",.) %>%
        sub("_s:.*","",.) %>%
        paste0(cluster[j,"V2"],";tax=",.)
      next
    }
    
    cluster[j,"centroid"] <- cluster[j,"centroid"] %>%
      sub(".*_o:","o:",.) %>%
      paste0(cluster[j,"V2"],";tax=",.)
    
  }
  
  write.csv(cluster,file = cluster_file[i] %>% gsub("_tablea_","_processed_",.))  
  
  seq_file <- cluster_file[i] %>% gsub("tablea","consensus",.) %>% gsub(".csv",".fasta",.)
  seq_file2 <- cluster_file[i] %>% gsub("tablea","consensusa",.) %>% gsub(".csv",".fasta",.)
  
  seq <- seq_file %>% readDNAStringSet 
  seq2 <- seq_file2 %>% readDNAStringSet 
  
  seq <- append(seq, seq2)
  seq <- seq[match( seq %>% names %>% gsub(".*clusterid=","",.) %>% as.numeric, cluster[,"V2"] )]
  
  for (j in 1:length(seq))
    names(seq)[j] <- names(seq)[j] %>% gsub("^centroid=(.*);tax=.*(seqs=.*);.*","\\2_cent:\\1",.) %>% paste0(cluster[j,"centroid"],";",.)
  
  writeXStringSet(x = seq, filepath= seq_file %>% gsub("_consensus_","_LCA_",.) )
  
}


cluster_file <- list.files(path = "./step10", pattern = "ITS._LCA_cluster_", full.names = T)
whichITS <- c(0,1,2)
cluster_at <- c(1, 0.995, 0.99, 0.98, 0.97, 0.96, 0.95)

for (i in 1:length(cluster_at))
  for (j in 1:length(whichITS)) 
  {
    actual_file <- paste0("ITS",whichITS[j],"_LCA_cluster_",cluster_at[i],".fasta", collapse="") %>% grep(.,cluster_file,value=T)
    
    actual_out <- actual_file %>% gsub(".fasta","_hits.csv",.)
    
    system(command = paste0(vsearch_directory, vsearch_exec, " --usearch_global ./step10/ITS",whichITS[j],".fasta --db ",actual_file,
                            " --id 0.70 --threads ",number_of_cores," --maxaccepts 50 --userout ",actual_out," --iddef 2", 
                            " --userfields query+target+id+alnlen"))
  }




#build a master table
whichITS <- c(0,1,2)
cluster_at <- c(1, 0.995, 0.99, 0.98, 0.97, 0.96, 0.95)
cluster_file <- list.files(path = "./step10", pattern = "ITS._table_cluster_", full.names = T)


for (w in 1:length(whichITS))
  for (j in 1:length(cluster_at))
  {
    #read raw reads after vsearch matching
    raw <- paste0("./step10/ITS",whichITS[w],"_LCA_cluster_",cluster_at[j],"_hits.csv",collapse="") %>% read.table
    
    #get the start of each query results
    get_start <- raw[,"V1"] %>% duplicated %>% not %>% which
    
    #create a list with sequences to keep (top sequences) and filter alternatives if not needed
    keep_me <- NULL
    
    #for each of the query sequences do
    for (i in 1:length(get_start))
    {
      is_self_singleton_cluster <- NULL
      #get the start and stop of each query sequence and write the content to temp
      if(i == length(get_start)) { temp_stop <- nrow(raw) 
      } else {
        temp_stop <- get_start[i+1]-1
      }
      temp <- raw[get_start[i]:temp_stop,]
      temp_orig <- temp
      #identify if the top hit is a self-hit or not (in the latter case we need to remove the alternatives)
      temp_no_alternatives <-  !grepl("[0-9]a$", temp[,"V2"])
      identify_top_hits <- temp[temp_no_alternatives,"V3"] == ( temp[temp_no_alternatives,"V3"] %>% max )
      is_top_self <- ( temp[temp_no_alternatives,][identify_top_hits,"V1"] %>% gsub(";.*","",.) ) ==( temp[temp_no_alternatives,][identify_top_hits,"V2"] %>% gsub(".*_cent:","",.) )
      #if (is_top_self %>% sum > 0) print(i)
      
      if (is_top_self %>% sum > 0) is_self_singleton_cluster <- grep("seqs=1", temp[temp_no_alternatives,][identify_top_hits,"V2"][is_top_self], value=T)
      
      if (is_self_singleton_cluster %>% length > 0) 
      { 
        temp <- temp[! ( ( temp[,"V2"] %in% is_self_singleton_cluster) | grepl("[0-9]a$", temp[,"V2"]) ),] #simple LOOCV possible, discard alternatives
      } else {
        if (is_top_self %>% sum > 0) temp <- temp[!temp[,"V2"] %in% temp[temp_no_alternatives,][identify_top_hits,"V2"],] #keep alternatives
      }
      
      if (is_top_self %>% sum == 0) temp <- temp[!grepl("[0-9]a$", temp[,"V2"]),] #no LOOCV needed, no self-hit. discard alternatives
      
      if ( nrow(temp)==0 ) next
      
      temp <- temp[( temp[,"V1"] %>% gsub(";.*","",.)  ) != ( temp[,"V2"] %>% gsub(".*_cent:","",.) ),]
      
      if ( nrow(temp)==0 ) next
      
      #get the new top hit after the cleaning and keep it
      identify_top_hits <- temp[,"V3"] == ( temp[,"V3"] %>% max )
      
      keep_me <- c(keep_me,temp[identify_top_hits,] %>% rownames %>% as.numeric)
      
      #i 87 cluster 1 size self-hit #i 51 self-hit multiple seqs in cluster 
      
    }
    
    actual <- raw[keep_me,]
    
    #create a new cluster tax description based on the remval of the query sequence if in cluster but not centroid
    
    clust_inf <- paste0("./step10/ITS",whichITS[w],"_processed_cluster_",cluster_at[j],".csv",collapse="") %>% read.csv %>% .[,-1]
    clust_inf <- clust_inf[match(actual[,"V2"] %>% sub(";.*","",.) %>% as.numeric,clust_inf[,"V2"]),]
    
    query_id <- actual[,"V1"] %>% sub(";.*","",.)
    
    for(i in 1:nrow(clust_inf))
      clust_inf[i,"revised_seq"] <- query_id[i] %>% paste0(" ",.,".*?( |$)") %>% gsub(., " ",clust_inf[i,"seq"]) %>% gsub(" $","",.)
    
    nchar(clust_inf[,"revised_seq"]) / nchar(clust_inf[,"seq"])
    
    
    clust_inf[nchar(clust_inf[,"revised_seq"]) / nchar(clust_inf[,"seq"]) < 0.5,c("seq","revised_seq")]
    
    # clust_inf[which(rownames(clust_inf)==112654),]
    # i<-2600 
    # 
    cluster <- clust_inf
    
    cluster[,"centroid"] <- cluster[,"centroid_org"]
    
    for (i in 1:nrow(cluster))
    {
      temp <- cluster[i,"revised_seq"] %>% as.character %>% strsplit(x=., split=" ") %>% unlist %>% gsub("^.*_f:","",.) %>% strsplit(.,"_")
      cluster[i,"same_fam"] <- temp %>% lapply(.,`[[`,1) %>% unlist %>% unique %>% length == 1
      cluster[i,"same_gen"] <- temp %>% lapply(.,`[[`,2) %>% unlist %>% unique %>% length == 1
      cluster[i,"same_spec"] <- temp %>% lapply(.,`[[`,3) %>% unlist %>% unique %>% length == 1
    }
    
    
    for (i in 1:nrow(cluster))
    {
      if (cluster[i,"same_fam"]==F)  
      {
        cluster[i,"centroid"] <- cluster[i,"centroid"] %>%
          sub(".*_o:","o:",.) %>%
          sub("_f:.*","",.) %>%
          paste0(cluster[i,"V2"],";tax=",.)
        next
      }
      
      if (cluster[i,"same_gen"]==F)  
      {
        cluster[i,"centroid"] <- cluster[i,"centroid"] %>%
          sub(".*_o:","o:",.) %>%
          sub("_g:.*","",.) %>%
          paste0(cluster[i,"V2"],";tax=",.)
        next
      }
      
      if (cluster[i,"same_spec"]==F)  
      {
        cluster[i,"centroid"] <- cluster[i,"centroid"] %>%
          sub(".*_o:","o:",.) %>%
          sub("_s:.*","",.) %>%
          paste0(cluster[i,"V2"],";tax=",.)
        next
      }
      
      cluster[i,"centroid"] <- cluster[i,"centroid"] %>%
        sub(".*_o:","o:",.) %>%
        paste0(cluster[i,"V2"],";tax=",.)
      
    }
    
    cluster[,"centroid"] <- sapply(cluster[,"centroid"],gsub,pattern="^.*tax=o:",replacement="tax=o:")
    
    for (i in 1:nrow(actual))
      cluster[i,"tax_revised"] <- actual[i,"V2"] %>% gsub(";tax=(.*);seqs=",paste0(";",cluster[i,"centroid"],";seqs="),.) 
    
    actual[,"V2"] <- cluster[,"tax_revised"]
    
    actual <- actual[,"V1"] %>% gsub(";.*","",.) %>% cbind(actual, "q_acc"=.)
    actual <- actual[,"V1"] %>% gsub(".*_o:(.*)_f:.*","\\1",.) %>% cbind(actual, "q_ord"=.)
    actual <- actual[,"V1"] %>% gsub(".*_f:(.*)_g:.*","\\1",.)  %>% cbind(actual, "q_fam"=.)
    actual <- actual[,"V1"] %>% gsub(".*_g:(.*)_s:.*","\\1",.)  %>% cbind(actual, "q_gen"=.)
    actual <- actual[,"V1"] %>% gsub(".*_s:","",.) %>% cbind(actual, "q_spec"=.)
    
    if (cluster_at[j]==0)
    {
      actual <- actual[,"V2"] %>% gsub(";.*","",.) %>% cbind(actual, "h_acc"=.)
      actual <- actual[,"V2"] %>% gsub(".*_o:(.*)_f:.*","\\1",.) %>% cbind(actual, "h_ord"=.)
      actual <- actual[,"V2"] %>% gsub(".*_f:(.*)_g:.*","\\1",.)  %>% cbind(actual, "h_fam"=.)
      actual <- actual[,"V2"] %>% gsub(".*_g:(.*)_s:.*","\\1",.)  %>% cbind(actual, "h_gen"=.)
      actual <- actual[,"V2"] %>% gsub(".*_s:","",.) %>% cbind(actual, "h_spec"=.)
    } else {
      actual <- actual[,"V2"] %>% gsub(".*cent:","",.) %>% cbind(actual, "h_acc"=.)
      actual <- actual[,"V2"] %>% gsub(".*o:(.*?)_f:.*","\\1",.) %>% gsub(".*o:(.*);seq.*","\\1",.) %>% gsub("^[0-9].*","",.) %>% cbind(actual, "h_ord"=.)
      actual <- actual[,"V2"] %>% gsub(".*f:(.*?)_g:.*","\\1",.) %>% gsub(".*f:(.*);seq.*","\\1",.) %>% gsub("^[0-9].*","",.) %>% cbind(actual, "h_fam"=.)
      actual <- actual[,"V2"] %>% gsub(".*g:(.*?)_s:.*","\\1",.) %>% gsub(".*g:(.*);seq.*","\\1",.) %>% gsub("^[0-9].*","",.) %>% cbind(actual, "h_gen"=.)
      actual <- actual[,"V2"] %>% gsub(".*s:(.*?);seq.*","\\1",.) %>% gsub("^[0-9].*","",.) %>% cbind(actual, "h_spec"=.)
    }
    
    get_start <- actual[,"V1"] %>% duplicated %>% not %>% which
    how_many_max_values <- tapply(actual[,"V3"], actual[,"V1"] %>% factor, function(x){which(x==max(x))}) %>% lapply(.,max)
    how_many_max_values <- how_many_max_values[match(actual[get_start,"V1"], how_many_max_values %>% names)]
    
    actual[,"fam_match"] <- NA; actual[,"gen_match"] <- NA; actual[,"spec_match"] <- NA
    actual[,"ambig"] <- NA
    
    for (i in 1:length(get_start))
    {
      temp <- actual[get_start[i] : (get_start[i] + how_many_max_values[[i]] - 1) ,]
      
      if ( identical(temp[,"h_spec"], temp[,"q_spec"]) )
      {
        actual[get_start[i],"spec_match"] <- TRUE 
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- TRUE 
        actual[get_start[i],"ambig"] <- FALSE
        next
      } 
      
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( temp[1,"q_spec"] %in% temp[,"h_spec"] ) & 
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) ) ) 
      {
        actual[get_start[i],"spec_match"] <- NA
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- TRUE 
        actual[get_start[i],"ambig"] <- TRUE
        next
      } 
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( "" %in% temp[,"h_spec"] %>% not ) & 
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) ) ) 
      {
        actual[get_start[i],"spec_match"] <- FALSE
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- TRUE 
        actual[get_start[i],"ambig"] <- FALSE
        next
      } 
      
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( "" %in% temp[,"h_spec"] ) & 
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) ) ) 
      {
        actual[get_start[i],"spec_match"] <- NA
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- TRUE 
        actual[get_start[i],"ambig"] <- TRUE
        next
      } 
      
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) %>% not) & 
           ( temp[1,"q_gen"] %in% temp[,"h_gen"] ) & 
           ( identical(temp[,"h_fam"], temp[,"q_fam"])) ) 
      {
        actual[get_start[i],"spec_match"] <- NA
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- NA
        actual[get_start[i],"ambig"] <- TRUE
        next
      } 
      
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) %>% not) & 
           ( "" %in% temp[,"h_gen"] %>% not ) & 
           ( identical(temp[,"h_fam"], temp[,"q_fam"])) ) 
      {
        actual[get_start[i],"spec_match"] <- FALSE
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- FALSE
        actual[get_start[i],"ambig"] <- FALSE
        next
      } 
      
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) %>% not) & 
           ( "" %in% temp[,"h_gen"]) & 
           ( identical(temp[,"h_fam"], temp[,"q_fam"])) ) 
      {
        actual[get_start[i],"spec_match"] <- NA
        actual[get_start[i],"fam_match"] <- TRUE 
        actual[get_start[i],"gen_match"] <- NA
        actual[get_start[i],"ambig"] <- TRUE
        next
      } 
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) %>% not) & 
           ( identical(temp[,"h_fam"], temp[,"q_fam"]) %>% not)  & 
           ( "" %in% temp[,"h_fam"] %>% not ) )
      {
        actual[get_start[i],"spec_match"] <- FALSE
        actual[get_start[i],"fam_match"] <- FALSE 
        actual[get_start[i],"gen_match"] <- FALSE
        actual[get_start[i],"ambig"] <- FALSE
        next
      } 
      
      
      if ( ( identical(temp[,"h_spec"], temp[,"q_spec"]) %>% not) &
           ( identical(temp[,"h_gen"], temp[,"q_gen"]) %>% not) & 
           ( identical(temp[,"h_fam"], temp[,"q_fam"]) %>% not)  & 
           ( "" %in% temp[,"h_fam"]  ) )
      {
        actual[get_start[i],"spec_match"] <- NA
        actual[get_start[i],"fam_match"] <- NA
        actual[get_start[i],"gen_match"] <- NA
        actual[get_start[i],"ambig"] <- TRUE
        next
      } 
      
    }
    
    actual <- actual[get_start,]
    
    write.csv(x = actual, file = paste0("./step10/ITS",whichITS[w],"_LCA_cluster_",cluster_at[j],"_evaluated.csv",collapse=""))
  }

whichITS <- c(0,1,2)
files <- paste0("ITS.*_LCA_cluster_.*_evaluated.csv",collapse="") %>% list.files(path="./step10",., full.names = T)

ggal <- NULL

for (i in 1:length(files))
{
  actual <- files[i] %>% read.table(.,sep=",", header=T, row.names = NULL) %>% .[,-(1:3)]
  temp_add_ITS <- files[i] %>% gsub(".*ITS(.*)_LCA.*","\\1",.) %>% as.integer
  temp_add_clus <- files[i] %>% gsub(".*cluster_(.*)_eval.*","\\1",.) %>% as.numeric
  ggal <- cbind.data.frame(actual, "ITS"=temp_add_ITS, "thresh"=temp_add_clus) %>% rbind.data.frame(.,ggal)
}

ggal <- ggal[ggal[,"thresh"]!="0",]

#one word classification

ggal[,"dec"] <- "mism"
ggal[which(ggal[,"ambig"]==TRUE),"dec"] <- "ambig"
ggal[which(ggal[,"spec_match"]==TRUE),"dec"] <- "hit"

ggal[,"dec"] %<>% factor(., levels=c("hit","ambig","mism"), ordered=T)
ggal[,"thresh"] %<>% factor(., levels= c( 1, 0.995, 0.99, 0.98, 0.97, 0.96, 0.95), ordered=T)

for (i in 1:length(whichITS))
{
  temp <- ggal[ggal[,"ITS"] %in% whichITS[i],]
  duplicated_species <- ( temp[temp[,"thresh"]=="1","q_spec"] %>% duplicated | temp[temp[,"thresh"]=="1","q_spec"] %>% duplicated(.,fromLast = T) ) %>% temp[temp[,"thresh"]=="1",][.,"q_acc"]
  temp_delete_me <- temp[temp[,"q_acc"] %in% duplicated_species %>% not,"q_acc"] %>% unique
  ggal <- ggal[!(ggal[,"ITS"] %in% whichITS[i] & ggal[,"q_acc"] %in% temp_delete_me),]
}

ggal[ggal[,"ITS"]==2,"ITS"] <- "ITS2"
ggal[ggal[,"ITS"]==1,"ITS"] <- "ITS1"
ggal[ggal[,"ITS"]==0,"ITS"] <- "ITS full"

ggplot(data=ggal, aes(x=factor(thresh), alluvium=q_acc, stratum=dec, fill=dec)
) + geom_stratum(
) + geom_flow(
  na.rm= TRUE,
  aes.flow= "backward",
  color="black",
  linewidth=.5,
  alpha=0.7
) + scale_x_discrete(
  expand = c(0,0),
  #breaks=cluster_at,
  #labels=c(1:5) %>% as.character,
  name="cluster thresholds"#,
  #breaks=scales::pretty_breaks(n=5)
) + theme(
  axis.title = element_text(face="bold",size=16),
  axis.text.y = element_text(colour="black",size=14),
  axis.text.x = element_text(angle= 90, hjust=1, vjust=.25, size=16,color="black"),
  strip.text.x = element_text(angle= 0, hjust=.5, vjust=.5, size=16,color="black"),
  axis.line.y = element_line(colour="black", linetype = "solid"),
  axis.line.x = element_line(colour="black", linetype = "solid"),
  strip.background = element_rect(colour = "black", fill = NA),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = alpha("black",.25)),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.position = c(0,1),#c(1,1),
  legend.text= element_text(colour="black",size=14),
  legend.title= element_text(face="bold",size=16),
  legend.justification = c(0,1),#c(1,1),
  legend.background = element_rect(colour=NA,fill=alpha("grey90", .95)),
  legend.key=element_blank(),
  panel.spacing=unit(0, "lines"),
  strip.placement="outside"
) + scale_y_continuous(
  breaks=function(x) { seq(0,max(x),max(x)/10) },
  #trans="pseudo_log",
  expand = c(0,0),
  labels=paste0(seq(0,100,10),"%"),       #function(x) paste0(x/nrow(ggal)*5*100, "%"),
  #limits=c(0,0.5),
  name="total sequences [%]"
) + facet_wrap (
  ~ ITS,  
  scales="free_y",
  nrow=1
) + scale_fill_manual(
  values=c("#00A0B0","#EDC951","#CC333F"), #blue#yellow#red
  labels=c("correct","ambig","false"),
  name=NULL
)



#accept ident thresholds

whichITS <- c(0,1,2)

load("ggal.list") %>% print

nrow(ggal)
#head(ggal)

ggal <- ggal[ggal[,"thresh"]==1,]

ggal[,"reject"] <- NA
ggal[,"thresh"] <- NA

rej_thresh <- c(1,0.995,0.99,0.98,0.97,0.96,0.95,0)
add_me <- NULL

for (i in 1:length(rej_thresh))
{
  temp <- ggal
  temp[,"thresh"] <- rej_thresh[i]
  add_me <- rbind.data.frame(add_me, temp)
}

ggal <- add_me

ggal[,"reject"] <- ggal[,"V3"] < (ggal[,"thresh"]*100)

ggal[,"dec"] <- "mism"
ggal[which(ggal[,"ambig"]==TRUE),"dec"] <- "ambig"
ggal[which(ggal[,"spec_match"]==TRUE),"dec"] <- "hit"
ggal[which(ggal[,"reject"]==TRUE),"dec"] <- "reject"

ggal[,"dec"] %<>% factor(., levels=c("hit","ambig","mism","reject"), ordered=T)
ggal[,"thresh"] %<>% factor(., levels= c(0, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 1), ordered=T)

for (i in 1:length(whichITS))
{
  temp <- ggal[ggal[,"ITS"] %in% whichITS[i],]
  duplicated_species <- ( temp[temp[,"thresh"]=="1","q_spec"] %>% duplicated | temp[temp[,"thresh"]=="1","q_spec"] %>% duplicated(.,fromLast = T) ) %>% temp[temp[,"thresh"]=="1",][.,"q_acc"]
  temp_delete_me <- temp[temp[,"q_acc"] %in% duplicated_species %>% not,"q_acc"] %>% unique
  ggal <- ggal[!(ggal[,"ITS"] %in% whichITS[i] & ggal[,"q_acc"] %in% temp_delete_me),]
}

ggal[ggal[,"ITS"]==2,"ITS"] <- "ITS2"
ggal[ggal[,"ITS"]==1,"ITS"] <- "ITS1"
ggal[ggal[,"ITS"]==0,"ITS"] <- "ITS full"


outplot <- ggplot(data=ggal, aes(x=factor(thresh), alluvium=q_acc, stratum=dec, fill=dec)
) + geom_stratum(
) + geom_flow(
  na.rm= TRUE,
  aes.flow= "backward",
  color="black",
  linewidth=.5,
  alpha=0.7
) + scale_x_discrete(
  expand = c(0,0),
  #breaks=cluster_at,
  #labels=c(1:5) %>% as.character,
  name="reject thresholds"#,
  #breaks=scales::pretty_breaks(n=5)
) + theme(
  axis.title = element_text(face="bold",size=16),
  axis.text.y = element_text(colour="black",size=14),
  axis.text.x = element_text(angle= 90, hjust=1, vjust=.25, size=16,color="black"),
  strip.text.x = element_text(angle= 0, hjust=.5, vjust=.5, size=16,color="black"),
  axis.line.y = element_line(colour="black", linetype = "solid"),
  axis.line.x = element_line(colour="black", linetype = "solid"),
  strip.background = element_rect(colour = "black", fill = NA),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = alpha("black",.25)),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.position = c(0,1),#c(1,1),
  legend.text= element_text(colour="black",size=14),
  legend.title= element_text(face="bold",size=16),
  legend.justification = c(0,1),#c(1,1),
  legend.background = element_rect(colour=NA,fill=alpha("grey90", .95)),
  legend.key=element_blank(),
  panel.spacing=unit(0, "lines"),
  strip.placement="outside"
) + scale_y_continuous(
  breaks=function(x) { seq(0,max(x),max(x)/20) },
  #trans="pseudo_log",
  expand = c(0,0),
  labels=paste0(seq(0,100,5),"%"),       #function(x) paste0(x/nrow(ggal)*5*100, "%"),
  #limits=c(0,0.5),
  name="total sequences [%]"
) + facet_wrap (
  ~ ITS,  
  scales="free_y",
  nrow=1
) + scale_fill_manual(
  values=c("#00A0B0","#EDC951","#CC333F","#93896b"), #blue#yellow#red
  labels=c("correct","ambig","false","reject"),
  name=NULL
)

plot(outplot)





####on genus level


whichITS <- c(0,1,2)

load("ggal.list") %>% print

nrow(ggal)
#head(ggal)

ggal <- ggal[ggal[,"thresh"]==1,]

ggal[( ggal[,"spec_match"] %>% is.na ) & ( ggal[,"gen_match"] == TRUE ) & ( ggal[,"gen_match"] %>% is.na %>% not ),"ambig"] <- FALSE

ggal[,"reject"] <- NA
ggal[,"thresh"] <- NA

rej_thresh <- c(0.99,0.97,0.95,0.90,0.85,0.8,0)
add_me <- NULL

for (i in 1:length(rej_thresh))
{
  temp <- ggal
  temp[,"thresh"] <- rej_thresh[i]
  add_me <- rbind.data.frame(add_me, temp)
}

ggal <- add_me

ggal[,"reject"] <- ggal[,"V3"] < (ggal[,"thresh"]*100)

ggal[,"dec"] <- "mism"
ggal[which(ggal[,"ambig"]==TRUE),"dec"] <- "ambig"
ggal[which(ggal[,"gen_match"]==TRUE),"dec"] <- "hit"
ggal[which(ggal[,"reject"]==TRUE),"dec"] <- "reject"

ggal[,"dec"] %<>% factor(., levels=c("hit","ambig","mism","reject"), ordered=T)
ggal[,"thresh"] %<>% factor(., levels= rej_thresh %>% rev, ordered=T)

# for (i in 1:length(whichITS))
# {
#   temp <- ggal[ggal[,"ITS"] %in% whichITS[i],]
#   duplicated_species <- ( temp[temp[,"thresh"]=="0","q_spec"] %>% duplicated | temp[temp[,"thresh"]=="0","q_spec"] %>% duplicated(.,fromLast = T) ) %>% temp[temp[,"thresh"]=="0",][.,"q_acc"]
#   temp_delete_me <- temp[temp[,"q_acc"] %in% duplicated_species %>% not,"q_acc"] %>% unique
#   ggal <- ggal[!(ggal[,"ITS"] %in% whichITS[i] & ggal[,"q_acc"] %in% temp_delete_me),]
# }

ggal[ggal[,"ITS"]==2,"ITS"] <- "ITS2"
ggal[ggal[,"ITS"]==1,"ITS"] <- "ITS1"
ggal[ggal[,"ITS"]==0,"ITS"] <- "ITS full"


outplot <- ggplot(data=ggal, aes(x=factor(thresh), alluvium=q_acc, stratum=dec, fill=dec)
) + geom_stratum(
) + geom_flow(
  na.rm= TRUE,
  aes.flow= "backward",
  color="black",
  linewidth=.5,
  alpha=0.7
) + scale_x_discrete(
  expand = c(0,0),
  #breaks=cluster_at,
  #labels=c(1:5) %>% as.character,
  name="reject thresholds"#,
  #breaks=scales::pretty_breaks(n=5)
) + theme(
  axis.title = element_text(face="bold",size=16),
  axis.text.y = element_text(colour="black",size=14),
  axis.text.x = element_text(angle= 90, hjust=1, vjust=.25, size=16,color="black"),
  strip.text.x = element_text(angle= 0, hjust=.5, vjust=.5, size=16,color="black"),
  axis.line.y = element_line(colour="black", linetype = "solid"),
  axis.line.x = element_line(colour="black", linetype = "solid"),
  strip.background = element_rect(colour = "black", fill = NA),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = alpha("black",.25)),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.position = c(0,1),#c(1,1),
  legend.text= element_text(colour="black",size=14),
  legend.title= element_text(face="bold",size=16),
  legend.justification = c(0,1),#c(1,1),
  legend.background = element_rect(colour=NA,fill=alpha("grey90", .95)),
  legend.key=element_blank(),
  panel.spacing=unit(0, "lines"),
  strip.placement="outside"
) + scale_y_continuous(
  breaks=function(x) { seq(0,max(x),max(x)/20) },
  #trans="pseudo_log",
  expand = c(0,0),
  labels=paste0(seq(0,100,5),"%"),       #function(x) paste0(x/nrow(ggal)*5*100, "%"),
  #limits=c(0,0.5),
  name="total sequences [%]"
) + facet_wrap (
  ~ ITS,  
  scales="free_y",
  nrow=1
) + scale_fill_manual(
  values=c("#00A0B0","#EDC951","#CC333F","#93896b"), #blue#yellow#red
  labels=c("correct","ambig","false","reject"),
  name=NULL
)

plot(outplot)




##regular clustering on genus level


whichITS <- c(0,1,2)

load("ggal.list") %>% print

ggal[( ggal[,"spec_match"] %>% is.na ) & ( ggal[,"gen_match"] == TRUE ) & ( ggal[,"gen_match"] %>% is.na %>% not ),"ambig"] <- FALSE
ggal <- ggal[ggal[,"thresh"]!="0",]

#one word classification

ggal[,"dec"] <- "mism"
ggal[which(ggal[,"ambig"]==TRUE),"dec"] <- "ambig"
ggal[which(ggal[,"gen_match"]==TRUE),"dec"] <- "hit"

ggal[,"dec"] %<>% factor(., levels=c("hit","ambig","mism"), ordered=T)
ggal[,"thresh"] %<>% factor(., levels= c( 1, 0.995, 0.99, 0.98, 0.97, 0.96, 0.95), ordered=T)

for (i in 1:length(whichITS))
{
  temp <- ggal[ggal[,"ITS"] %in% whichITS[i],]
  duplicated_species <- ( temp[temp[,"thresh"]=="1","q_spec"] %>% duplicated | temp[temp[,"thresh"]=="1","q_spec"] %>% duplicated(.,fromLast = T) ) %>% temp[temp[,"thresh"]=="1",][.,"q_acc"]
  temp_delete_me <- temp[temp[,"q_acc"] %in% duplicated_species %>% not,"q_acc"] %>% unique
  ggal <- ggal[!(ggal[,"ITS"] %in% whichITS[i] & ggal[,"q_acc"] %in% temp_delete_me),]
}

ggal[ggal[,"ITS"]==2,"ITS"] <- "ITS2"
ggal[ggal[,"ITS"]==1,"ITS"] <- "ITS1"
ggal[ggal[,"ITS"]==0,"ITS"] <- "ITS full"

ggplot(data=ggal, aes(x=factor(thresh), alluvium=q_acc, stratum=dec, fill=dec)
) + geom_stratum(
) + geom_flow(
  na.rm= TRUE,
  aes.flow= "backward",
  color="black",
  linewidth=.5,
  alpha=0.7
) + scale_x_discrete(
  expand = c(0,0),
  #breaks=cluster_at,
  #labels=c(1:5) %>% as.character,
  name="cluster thresholds"#,
  #breaks=scales::pretty_breaks(n=5)
) + theme(
  axis.title = element_text(face="bold",size=16),
  axis.text.y = element_text(colour="black",size=14),
  axis.text.x = element_text(angle= 90, hjust=1, vjust=.25, size=16,color="black"),
  strip.text.x = element_text(angle= 0, hjust=.5, vjust=.5, size=16,color="black"),
  axis.line.y = element_line(colour="black", linetype = "solid"),
  axis.line.x = element_line(colour="black", linetype = "solid"),
  strip.background = element_rect(colour = "black", fill = NA),
  panel.background = element_rect(fill = "white"),
  panel.grid.major.y = element_line(colour = alpha("black",.25)),
  panel.grid.minor.y = element_blank(),
  panel.grid.major.x = element_blank(),
  legend.position = c(0,1),#c(1,1),
  legend.text= element_text(colour="black",size=14),
  legend.title= element_text(face="bold",size=16),
  legend.justification = c(0,1),#c(1,1),
  legend.background = element_rect(colour=NA,fill=alpha("grey90", .95)),
  legend.key=element_blank(),
  panel.spacing=unit(0, "lines"),
  strip.placement="outside"
) + scale_y_continuous(
  breaks=function(x) { seq(0,max(x),max(x)/10) },
  #trans="pseudo_log",
  expand = c(0,0),
  labels=paste0(seq(0,100,10),"%"),       #function(x) paste0(x/nrow(ggal)*5*100, "%"),
  #limits=c(0,0.5),
  name="total sequences [%]"
) + facet_wrap (
  ~ ITS,  
  scales="free_y",
  nrow=1
) + scale_fill_manual(
  values=c("#00A0B0","#EDC951","#CC333F"), #blue#yellow#red
  labels=c("correct","ambig","false"),
  name=NULL
)




