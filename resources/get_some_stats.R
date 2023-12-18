library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)



x <- read.table(file = "step10_meta.csv", sep = "\t", header=TRUE)

y <- x[,"date"] %>% as.Date(.,format="%d-%B-%Y") %>% table %>% cbind.data.frame
colnames(y) <- c("date","count")
y[,1] %<>% as.Date

y %<>% arrange(date) %>% mutate("cumulative"=cumsum(count))

#number of ITS sequences by year
ggplot(data=y, aes(date, cumulative)) +
  geom_line()


#number of total species sequenced with ITS
w <- x[,c("species","date")]

w <- w[w[,"species"] %>% duplicated %>% not,]

w %<>% arrange(date) %>% group_by(date) %>% summarise("count"=length((species)))
w %<>% data.frame
w[,"date"] %<>% as.Date(.,format="%d-%B-%Y")

w %<>% arrange(date) %>% mutate("cumulative"=cumsum(count))

ggplot(data=w, aes(date, cumulative)) +
  geom_line()

  
#number of sequences per step


meta <- read.table("./step1/step1_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8")
#meta6 does not exist, setp6 is only spike in folder

temp_spec <- meta[,"orig_title"] %>% gsub("^([A-Z]{1,1}[a-z]{1,} [a-z-]{1,}) .*|.*","\\1",.)
temp_gen <- temp_spec %>% gsub(" .*","",.)
temp_fam <- meta[,"orig_taxonomy"] %>% gsub("^.*([A-Z]{1,1}[a-z]{1,}ceae); .*","\\1",.)
temp_ord <- meta[,"orig_taxonomy"] %>% gsub("^.*([A-Z]{1,1}[a-z]{1,}ales); .*|.*","\\1",.)

temp <- cbind.data.frame("acc"=meta[,"acc"],"fam"=temp_fam, "gen"=temp_gen, "spec"=temp_spec, "ord"=temp_ord)
temp <- temp[temp[,"gen"]!= "" & temp[,"spec"]!= "" & temp[,"fam"] != "",]

paste0("number of species: ", temp[,"spec"] %>% unique %>% length) #meta1 112545
paste0("number of genera: ", temp[,"gen"] %>% unique %>% length) #meta1 11339 
paste0("number of families: ", temp[,"fam"] %>% unique %>% length) #meta1 449 
paste0("number of orders: ", temp[,"ord"] %>% unique %>% length) #meta1 449 
paste0("number of sequences: ", nrow(meta)) #meta1 409034 


meta <- read.table("./step3/step3_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8")
#meta6 does not exist, setp6 is only spike in folder

temp_spec <- meta[,"orig_title"] %>% gsub("^([A-Z]{1,1}[a-z]{1,} [a-z-]{1,}) .*|.*","\\1",.)
temp_gen <- temp_spec %>% gsub(" .*","",.)
temp_fam <- meta[,"orig_taxonomy"] %>% gsub("^.*([A-Z]{1,1}[a-z]{1,}ceae); .*","\\1",.)
temp_ord <- meta[,"orig_taxonomy"] %>% gsub("^.*([A-Z]{1,1}[a-z]{1,}ales); .*|.*","\\1",.)

temp <- cbind.data.frame("acc"=meta[,"acc"],"fam"=temp_fam,
                         "gen"=temp_gen,
                         "spec"=temp_spec,
                         "ord"=temp_ord,
                         "ITS1"=meta[,"ITS1_found"],
                         "ITS2"=meta[,"ITS2_found"],
                         "ITS0"=(meta[,"ITS1_found"] & meta[,"ITS2_found"]))

temp <- temp[temp[,"gen"]!= "" & temp[,"spec"]!= "" & temp[,"fam"] != "",]

paste0("number of species ITS1/ITS2/ITS0/all: ", temp[temp[,"ITS1"]==TRUE,"spec"] %>% unique %>% length,
       " / ", temp[temp[,"ITS2"]==TRUE,"spec"] %>% unique %>% length,
       " / ", temp[temp[,"ITS0"]==TRUE,"spec"] %>% unique %>% length,
       " / ", temp[,"spec"] %>% unique %>% length) #meta2 112545

paste0("number of genera ITS1/ITS2/ITS0/all: ", temp[temp[,"ITS1"]==TRUE,"gen"] %>% unique %>% length,
       " / ", temp[temp[,"ITS2"]==TRUE,"gen"] %>% unique %>% length,
       " / ", temp[temp[,"ITS0"]==TRUE,"gen"] %>% unique %>% length,
       " / ", temp[,"gen"] %>% unique %>% length) #meta2 112545

paste0("number of families ITS1/ITS2/ITS0/all: ", temp[temp[,"ITS1"]==TRUE,"fam"] %>% unique %>% length,
       " / ", temp[temp[,"ITS2"]==TRUE,"fam"] %>% unique %>% length,
       " / ", temp[temp[,"ITS0"]==TRUE,"fam"] %>% unique %>% length,
       " / ", temp[,"fam"] %>% unique %>% length) #meta2 112545

paste0("number of order ITS1/ITS2/ITS0/all: ", temp[temp[,"ITS1"]==TRUE,"ord"] %>% unique %>% length,
       " / ", temp[temp[,"ITS2"]==TRUE,"ord"] %>% unique %>% length,
       " / ", temp[temp[,"ITS0"]==TRUE,"ord"] %>% unique %>% length,
       " / ", temp[,"ord"] %>% unique %>% length) #meta2 112545

paste0("number of sequences ITS1/ITS2/ITS0/all: ",
       meta[meta[,"ITS1_found"]==TRUE,] %>% nrow,
       " / ", meta[meta[,"ITS2_found"]==TRUE,] %>% nrow,
       " / ", meta[meta[,"ITS1_found"] & meta[,"ITS2_found"],] %>% nrow,
       " / ", meta[,] %>% nrow) #meta2 112545


meta <- read.table("./step9/step9_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8")
#meta_o <- read.table("./step8/step8_meta.csv", header = TRUE, sep='\t', dec=".", fileEncoding = "UTF-8")
#meta6 does not exist, setp6 is only spike in folder
#fungal step use: "./step5/step5_meta_edit.csv"
#for distance based analysis jump from 5_edit to meta9


meta <- meta[meta[,"phylum"]=="Tracheophyta",]

paste0("number of sequences ITS1/ITS2/ITS0/all: ",
       meta[meta[,"ITS1_found"]==TRUE,"species"] %>% length,
       " / ", meta[meta[,"ITS2_found"]==TRUE,"species"] %>% length,
       " / ", meta[meta[,"ITS1_found"] & meta[,"ITS2_found"],"species"] %>% length,
       " / ", meta[,"species"] %>% length) #meta2 112545

paste0("number of species ITS1/ITS2/ITS0/all: ",
       meta[meta[,"ITS1_found"]==TRUE,"species"] %>% unique %>% length,
       " / ", meta[meta[,"ITS2_found"]==TRUE,"species"] %>% unique %>% length,
       " / ", meta[meta[,"ITS1_found"] & meta[,"ITS2_found"],"species"] %>% unique %>% length,
       " / ", meta[,"species"] %>% unique %>% length) #meta2 112545

paste0("number of genera ITS1/ITS2/ITS0/all: ",
       meta[meta[,"ITS1_found"]==TRUE,"genus"] %>% unique %>% length,
       " / ", meta[meta[,"ITS2_found"]==TRUE,"genus"] %>% unique %>% length,
       " / ", meta[meta[,"ITS1_found"] & meta[,"ITS2_found"],"genus"] %>% unique %>% length,
       " / ", meta[,"genus"] %>% unique %>% length) #meta2 112545

paste0("number of families ITS1/ITS2/ITS0/all: ",
       meta[meta[,"ITS1_found"]==TRUE,"family"] %>% unique %>% length,
       " / ", meta[meta[,"ITS2_found"]==TRUE,"family"] %>% unique %>% length,
       " / ", meta[meta[,"ITS1_found"] & meta[,"ITS2_found"],"family"] %>% unique %>% length,
       " / ", meta[,"family"] %>% unique %>% length) #meta2 112545

paste0("number of order ITS1/ITS2/ITS0/all: ",
       meta[meta[,"ITS1_found"]==TRUE,"order"] %>% unique %>% length,
       " / ", meta[meta[,"ITS2_found"]==TRUE,"order"] %>% unique %>% length,
       " / ", meta[meta[,"ITS1_found"] & meta[,"ITS2_found"],"order"] %>% unique %>% length,
       " / ", meta[,"order"] %>% unique %>% length) #meta2 112545


meta4[meta4[,"order"] %in% meta5[,"order"] %>% not,]




paste0("number of species: ", meta[,"species"] %>% unique %>% length) #meta4 102273 #meta5  #meta7 
paste0("number of genera: ", meta[,"genus"] %>% unique %>% length) #meta4 10218 #meta5  #meta7 
paste0("number of families: ", meta[,"family"] %>% unique %>% length) #meta4 421 #meta5  #meta7 


temp_tax <- read.table(file = "./step3/tax_results.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE)   

temp_tax <- temp_tax[temp_tax[,"kingdom"] != "Animalia",]
temp_tax <- temp_tax[temp_tax[,"matchType"] != "NONE",]
temp_tax <- temp_tax[temp_tax[,"matchType"] != "HIGHERRANK",]
temp_tax <- temp_tax[temp_tax[,"verbatim_name"] %>% duplicated %>% not,] %>% data.frame

temp_tax <- temp_tax[,"verbatim_name"] %>%
  gsub("^([A-Z]{1,1}[a-z]{1,} [a-z-]{2,}) .*","\\1",.) %>%
  cbind.data.frame(temp_tax, "verbatim_species"=.)

#howManySpecMoreThanOne <- temp_tax %>% group_by(verbatim_species) %>% summarise(length(species)) %>% as.data.frame


  
temp_tax[,"status"] %>% table
temp_tax[,"confidence"] %>% table

meta3 <- read.table("./step3/step3_meta_edit.csv", sep='\t', dec=".", fileEncoding = "UTF-8", header = TRUE)  


















  
  