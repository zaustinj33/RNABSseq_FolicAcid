library(dplyr)

# Find files to merge
# change suffix to find appropriate files
suffix <- "_Genome10xCall_3_Cutoff_0ML.txt"  #depth sites 0ML
#suffix <- "_3strucFilter.txt"  # sites after secondary structure filter
path_fdr <- "/call_sites/"

mRNAfiles <- list.files(path=paste0(getwd(),path_fdr), pattern = suffix,
                        full.names = T, recursive = T)

mRNA_dfs <- lapply(mRNAfiles, function(x) {
  # read in raw files and grab only relevant columns to merge together
  y <- read.csv(x, comment.char = '',header = T,sep='\t')[,c(1,2,5,6,7,3,18)]
  name <- gsub(suffix,"",basename(x))
  name <- gsub("^LF", "tLF", name)
  name <- gsub("^MF", "tMF", name)
  name <- gsub("^HFold", "tHF", name)

  # add name from column for easy merging
  names(y)[-6:-8] <- paste0(names(y)[-6:-8],"_",name)
  
  
  #20x, 3m5C filter for replicates
  # comment out if preparing depth file 
  #y <- y %>%
   # filter_if(grepl("C_count",colnames(y)), any_vars(. >= 6)) %>%
    #filter_if(grepl("cov",colnames(y)), any_vars(. >= 20)) %>%
    #filter_if(grepl("methRate",colnames(y)), any_vars(. >= 0.15))


  cat(name,": ",nrow(y),"\n")
  y$group <- paste(y$X.SeqID, y$refPos)
  y[,-1:-2]
})


# Combined sites dataframe
merge.total <- Reduce(function(...) merge(..., by = c('group','seqContext','refStrand'), all = T), mRNA_dfs)
merge.total <- merge.total[!duplicated(merge.total$group),]
merge.total$group <- gsub("chr","",merge.total$group)
merge.total$group <- gsub("M","MT",merge.total$group)
merge.total$chrom <- gsub(" .*","",merge.total$group)
merge.total$refPos <- as.integer(gsub(".* ","",merge.total$group))

#merge.total.final <- merge.total
#write.csv(merge.total, "allm5C_libraries.csv",row.names = F)
write.csv(merge.total, "allm5C_libraries_6Cfiltered.csv",row.names = F)

## Write bed file for calling sites at depth
all_bed <- merge.total %>% select(chrom, refPos)
all_bed$chrom <- paste0("chr", all_bed$chrom)
all_bed$chrom <- gsub("chrMT", "chrM", all_bed$chrom)
all_bed$end <- all_bed$refPos+1
all_bed$refPos <- all_bed$refPos-1
write.table(all_bed, "all_m5C.bed", quote = F, sep='\t', row.names = F, col.names = F)

## Filter unfiltered df by filtered df to get depth at all sites in the
unfiltered <- read.csv("allm5C_libraries_0MLfilteredDepthAnno.csv", comment.char = '',header = T)

filtered <- read.csv("allm5C_libraries_6Cfiltered.csv", comment.char = '',header = T)
Anno <- unfiltered %>% select(group, gene, type, start, end, trans,  gene_name, length, binLoc)

filtered_depth <- subset(merge.total,group %in% filtered$group)
write.csv(filtered_depth, "allm5C_libraries_filteredDepth.csv",row.names = F)

filtered_Anno <- left_join(filtered, Anno, by="group")
write.csv(filtered_depthAnno, "allm5C_libraries_filteredAnno.csv",row.names = F)
