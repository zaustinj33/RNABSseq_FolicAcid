library(ggpubr)
library(dplyr)
library(tidyr)
library(ensembldb)
library(AnnotationHub)

all_m5C <- read.csv("allm5C_libraries_0MLfilteredDepthAnno.csv")
overlap_files <- list.files(pattern = "_allOverlap.csv$")

#txdb <- makeTxDbFromGFF('gencode.v28.annotation.gtf.gz', tx_type=TRUE)
#saveDb(txdb, 'txdb.gencode24.sqlite')
#txdb <- loadDb(file = 'txdb.gencode24.sqlite')
#transcripts(txdb)

# annotation track
ah <- AnnotationHub()
edb <- ah[[names(query(ah, "EnsDb.Mmusculus.v104"))]]
gene_ids <- all_m5C$gene[!is.na(all_m5C$gene)]
anno <- as.data.frame(genes(edb, filter = ~ gene_id == gene_ids, columns = c("gene_biotype",'gene_name'),
              return.type = "DataFrame"))

biotype_order <- c("protein_coding","processed_pseudogene", "transcribed_processed_pseudogene",
                   "transcribed_unprocessed_pseudogene","lncRNA", "Mt_rRNA", "Mt_tRNA")
sample_order <- c("tLF","tMF","tHF", "pLF","pMF","pHF")
mypal <- pal_npg("nrc", alpha = 0.9)(9)[1:7]


## Individual files
plot_biotype <- function(input_file){
  overlap <- read.csv(input_file)  #all_m5C
  condition <- gsub("_allOverlap.csv","",input_file)
  overlap$gene_name <- overlap$gene
  overlap <- left_join(overlap, all_m5C, by='group')
  overlap <- left_join(overlap, anno, by='gene_name')


  biotypes <- as.data.frame(table(overlap$gene_biotype))
  biotypes <- biotypes[order(-biotypes$Freq),]

  biotypes$Var1 <- factor(biotypes$Var1, levels = rev(biotype_order))
  biotypes$perc <- biotypes$Freq / sum(biotypes$Freq)

  #p <- ggbarplot(biotypes, x="Var1", y="Freq", sort.val = "desc", x.text.angle = 90, fill = "Var1", sort.by.groups=FALSE)
  #p
  write.csv(biotypes, paste0(condition,"_biotypeTable.csv"), row.names = F)
  p<-ggplot(biotypes, aes(x="", y=Freq, fill=Var1)) +
    geom_bar(stat='identity', width=1,color='black') +
    coord_polar("y", start=0) +
    theme_void() +
    scale_fill_manual(breaks = (biotype_order),
                      values=mypal) +
    theme(legend.position = "none")
  p
  ggsave(filename = paste0('plot_biotype_',condition,'_pie.png'),height = 15,width = 15,units = 'cm',dpi = 200)
}
plot_biotype(overlap_files[[4]])
lapply(overlap_files, plot_biotype)

## Composite stacked bar plot
gather_biotype <- function(input_file){
  overlap <- read.csv(input_file)
  overlap$gene_name <- overlap$gene
  overlap <- left_join(overlap, all_m5C, by='group')
  overlap <- left_join(overlap, anno, by='gene_name')


  biotypes <- t(as.data.frame(table(overlap$gene_biotype)))
  colnames(biotypes) <- biotypes[1,]
  biotypes <- biotypes[-1,]
  return(biotypes)
}

test <- lapply(overlap_files, gather_biotype)
test_df <- bind_rows(test)

file_names <- lapply(overlap_files, basename)
file_names <- gsub("_alloverlap.csv","",file_names)
test_df$sample <- file_names
test

long_plot <- test_df  %>%
  pivot_longer(!sample, names_to = "biotype", values_to = "count") %>%
  mutate_at(3,as.numeric) %>%
  replace(is.na(.), 0) %>%
  group_by(sample) %>%
  mutate(per = prop.table(count)*100)

#long_plot$biotype <- gsub("lincRNA","lncRNA", long_plot$biotype)
#long_plot$biotype <- gsub("lincRNA","lncRNA", long_plot$biotype)

long_plot$biotype <- factor(long_plot$biotype, levels = rev(biotype_order))
long_plot$sample <- factor(long_plot$sample, levels = (sample_order))


ggplot(long_plot, aes(x=sample, y=per, fill=biotype)) + geom_bar(stat='identity', position='stack') +
  theme_bw() + scale_fill_manual(breaks = (biotype_order),
                                 values=mypal)

ggsave(filename = 'plot_biotype_legend.png',height = 15,width = 15,units = 'cm',dpi = 200)
