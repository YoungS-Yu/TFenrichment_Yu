library(stringr)
library(tidyr)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ChIPseeker)
library(dplyr)
library(tibble)


### Load variables ###
TF <- read.delim("TF.txt",sep = "\t", header =FALSE, quote = "") %>% pull()
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
annotations <- read.table("tx2gene_human_EnsDBl105.txt", sep = "\t", quote = "", header = T)

#total gene number
annotations  %>% filter(symbol != "<NA>") %>% dplyr::select(symbol) %>% filter(!duplicated(symbol)) %>% nrow
#result:39473


### Filtering TF ###
samples <- list.files(path = "./sample_bed",pattern="bed", full.names = T)
files <- file.path(samples)
names(files) <- str_replace(samples, "./sample_bed/", "") %>% str_replace(".sorted.bed", "") 

peakAnnoNames <- c()
for (var in names(files)){
  newname <- as.character(str_split(var,"_", simplify = T))[1]
  ifelse(newname %in% TF, peakAnnoNames <- append(peakAnnoNames,var), NA)  
}
filtered_files <- files[names(files) %in% peakAnnoNames]


### Peak annotation using ChIPseeker ###
# Promoter region [-1000,1000]
peakAnnoList <- lapply(filtered_files, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-1000, 1000), verbose=FALSE)

### Add annotations for each peaks ###
AnnoList <- list()
for (var in peakAnnoNames){
  anno <- data.frame(peakAnnoList[[var]]@anno)
  anno$transcriptId <- sub('\\.[0-9]*$', '', anno$transcriptId)
  newlist <- list(left_join(anno, annotations, by= c("transcriptId"="tx_id")))
  names(newlist) <- var
  AnnoList <- append(AnnoList,newlist)
}

### IDR filtering ###
AnnoList2=list()
for (var in peakAnnoNames){
  newlist <- AnnoList[[var]] %>% dplyr::filter(V5 >= 830.482) %>% arrange(-V7) %>% list() #%>% dplyr::filter(AnnoList_t[[var]]$width <= summary(AnnoList_t[[var]]$width)[5]) -> newdf
  names(newlist) <- var
  AnnoList2 <- append(AnnoList2,newlist)
}

### save RDS files ###
saveRDS(AnnoList2, file="AnnoList2.RData")
saveRDS(peakAnnoNames, file="peakAnnoNames.RData")





### ChIP-seq based TF enrichment ###
Anno_table_maker <- function ( AnnoList = AnnoList2, genelist = ALgenes, peakAnnoNames = peakAnnoNames ) {
  Anno_df <- data.frame()
  for (var in peakAnnoNames){
    newdf = data.frame(
      target.genes = AnnoList[[var]] %>% filter(annotation=="Promoter") %>% filter(symbol != "<NA>") %>% filter(symbol %in% genelist) %>% dplyr::select(symbol) %>% filter(!duplicated(symbol)) %>% nrow,
      total.genes = AnnoList[[var]] %>% filter(annotation=="Promoter") %>% filter(symbol != "<NA>") %>% dplyr::select(symbol) %>% filter(!duplicated(symbol)) %>% nrow
    )
    rownames(newdf) <- var
    Anno_df <- rbind(Anno_df,newdf)
  }
  
  Anno_fish_table.ls <- list()
  for (var in peakAnnoNames){
    Anno_fish_table = data.frame("TF"=c(Anno_df[var,]$target.genes , Anno_df[var,]$total.genes - Anno_df[var,]$target.genes),
                                 "nonTF"=c(length(genelist) - Anno_df[var,]$target.genes  , 39473 - Anno_df[var,]$total.genes -length(genelist) + Anno_df[var,]$target.genes ), 
                                 row.names= c("Target","nonTarget"))
    Anno_fish_table.ls[[length(Anno_fish_table.ls) + 1]] <- Anno_fish_table
  }
  names(Anno_fish_table.ls) <- peakAnnoNames
  
  Anno_fish_p <- data.frame()
  for (var in peakAnnoNames){
    fishTest<-fisher.test(Anno_fish_table.ls[[var]])
    newdf <- data.frame( "logPvalue" = -log(fishTest$p.value, base=10), "OddsRatio" = fishTest$estimate, "ci_1" = fishTest$conf.int[1], "ci_2" = fishTest$conf.int[2] )
    rownames(newdf) <- var
    Anno_fish_p <- rbind(Anno_fish_p,newdf)
  }
  
  Anno_df_p <- cbind(Anno_df,Anno_fish_p)
  Anno_df_order <- Anno_df_p %>% arrange(-logPvalue) %>% mutate(generatio = target.genes / total.genes) %>% mutate(rank = row_number(-logPvalue)) %>% tibble::rownames_to_column(var = "Names")
  
  return(Anno_df_order)
}



### main ###
### Read RDS files ###
annoList_f <- readRDS("AnnoList2.RData")
peakAnnoNames <- readRDS("peakAnnoNames.RData")

### Input gene lists ###
Lys_bio <- read.delim(file = "Lysosomal_biogenesis_Bordi.txt", sep = "\t", header = F, quote = "")
DE_Lys <- read.table("Differentially_expressed_lysosomal_genes.txt", sep = "\t", header = F) 

Lys_bio <- Lys_bio %>% pull()
DE_Lys <- DE_Lys %>% pull()

### results ###
Lys_bio_table <- Anno_table_maker(annoList_f,Lys_bio,peakAnnoNames)
DE_Lys_table <- Anno_table_maker(annoList_f,DE_Lys,peakAnnoNames)




sessionInfo()




