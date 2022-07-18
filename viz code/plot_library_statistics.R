## m5C, mRNA counts plots, to be combined with fig01a_boxplot

library(ggplot2)
library(ggsci)
library(dplyr)
library(ggpubr)
library(rstatix)
library(scales)

# set colors
#my_pal <- pal_npg('nrc')(10)
my_cols_mRNA <- c("#00FF00","#0000FF","#FF0000")#my_pal[c(2,8,5)] # mRNA all, MT,nonMT colors
# Import overlap files

overlap_files <- list.files(pattern = "Overlap.csv$", recursive = T)
anno_file <- read.csv('allm5C_libraries_0MLfilteredDepthAnno.csv')
# only grab methylLevel columns
counts_DF <- function(name){
  overlap <- read.csv(name)
  print(nrow(overlap))
  overlap <- anno_file[anno_file$group %in% overlap$group,]
  DF <- data.frame(group=overlap$group,
                   mRNA_count = overlap$trans)
  return(DF)

}
ML_list <- lapply(overlap_files, counts_DF)
stats_DF <- bind_rows(ML_list[[4]],ML_list[[5]],ML_list[[6]], .id = "sample")

stats_DF$condition <- ifelse(grepl("1",stats_DF$sample),"tHF",
                             ifelse(grepl("2",stats_DF$sample),"tLF",
                                    "tMF"))
stats_DF$condition <- factor(stats_DF$condition, levels = c('tLF','tMF','tHF'))

m5C_sum <- stats_DF %>%
  group_by(condition) %>%
  summarise(m5C_count = n_distinct(group),
            mRNA = n_distinct(mRNA_count))

m5C <- ggplot(m5C_sum, aes(x = condition, fill = condition, y=m5C_count)) +
  geom_bar(stat='identity',color='black',size = 1,width = 0.6)+
  theme_classic() +
  scale_fill_manual(values=c(my_cols_mRNA)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0,50,100),
                     labels = c(0,50,100),
                     limits = c(0,100)) +
  theme(legend.position = 'none') +
  ylab("Unique m5C sites")

m5C

mRNA <- ggplot(m5C_sum, aes(x = condition, fill = condition, y=mRNA)) +
  geom_bar(stat='identity',color = 'black',size=1,width = 0.6)+
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c(my_cols_mRNA)) +
  scale_y_continuous(expand = c(0,0),
                     breaks = c(0,50,100),
                     labels = c(0,50,100),
                     limits = c(0,100))+
  ylab("Unique Transcripts")
mRNA

library(gridExtra)
plot <- grid.arrange(m5C, mRNA, nrow=1, widths = c(2,2))
#library(gtable)
#plot <- rbind(ggplotGrob(p),ggplotGrob(m5C),ggplotGrob(mRNA))

ggsave(filename = 'm5C_stats_total.png',plot = plot, height = 3,width = 5,dpi = 200)
write.csv(m5C_sum, "m5C_stats_total.csv", row.names = F)
