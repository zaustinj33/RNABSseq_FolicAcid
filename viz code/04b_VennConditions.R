library(grid)
library(VennDiagram)
library(scales)
library(gplots)

## Venn diagram to compare m5C replicates and overlap

# set directories #
home_dir <- paste0(getwd())
setwd(home_dir)

#union_files <- list.files(pattern = "allUnion.csv", full.names = T)  # Union
OLdepth_files <- list.files(pattern = "allOverlap.csv$", full.names = T)  # Overlap with depth

# Write high-confidence sites only
library(dplyr)
df <- OLdepth_files %>%
  lapply(read.csv) %>%
  bind_rows

all_m5C <- read.csv("allm5C_libraries_0MLfilteredDepthAnno.csv")[c(1,40:49)]
df_anno <- left_join(df, all_m5C, by='group')
write.csv(df_anno, "allm5C_UnionConditions.csv", row.names = F)

generate_gene_list <- function(file){
  overlap <- read.csv(file)
  name <-  gsub("all.*$","",basename(file))

  site_list <- overlap$group
  return(site_list)
}
####### CHANGE THIS TO MATCH COMPARISON ###############
names = OLdepth_files
input <- lapply(names, generate_gene_list)
names(input) <- lapply(names, basename)

png(file=paste0("MF_6Coverlap_venn.png"))
venn_table <- venn(input, intersections = TRUE)
dev.off()
#########################################################
#Union of true overlap dataframe
merge_input <- unlist(input)
all_m5C <- read.csv("allm5C_libraries_0MLfilteredDepthAnno.csv")
union_df <- all_m5C[all_m5C$group %in% merge_input,]
write.csv(union_df, "allm5C_UnionConditions.csv", row.names = F)

#########################################################
library('ggplot2')
library('dplyr')
library('ggpubr')
## Intersection lists
all_m5C <- read.csv("allm5C_libraries_filteredDepth.csv")
intersections <- unlist(attr(venn_table,"intersections"))

# Compare methylation levels
plot_methylation_level_DMS <- function(input_compare, name){
  input_compare <- input_compare[input_compare$group %in% intersections,]
  info <- input_compare[, grepl(name, colnames(input_compare))]

  m5C_info <- (data.frame('total_ML' = (info[,2]+info[,5])/(info[,1]+info[,4]),
                          'poly_ML' = (info[,8]+info[,11])/(info[,7]+info[,10])))

  m5C_info <- m5C_info[complete.cases(m5C_info),]
  print(wilcox.test(m5C_info$total_ML, m5C_info$poly_ML, paired=T))

  boxplot_input <- stack(m5C_info)
  print(paste("Depth overlap: ", nrow(m5C_info)))
  print(paste("non-NA positions: ", nrow(m5C_info)))

  color <- ifelse(name == 'HF', 'red',
                 ifelse(name == 'MF', 'blue',
                       ifelse(name == 'LF', 'green', 'black')))
  boxplot_input$color <- "fill"

  g <- ggplot(boxplot_input, aes(x=ind, y=values, fill=color)) + geom_boxplot() + theme_classic() +
    scale_fill_manual(values= color) +
    scale_y_continuous(limits = c(0,1,1), breaks = c(0,0.5,1)) +
    theme(axis.ticks.length=unit(.25, "cm")) +
    scale_x_discrete(labels=c("total_ML" = paste('total',name), "poly_ML" = paste('poly',name))) +
    xlab("Fraction") +
    ylab("Methylation level") +
    stat_compare_means(comparisons = list(c('total_ML', 'poly_ML')), label.y = 0.9)
  g

  ggsave(filename = paste0(name,'all_depthOL_boxplot.png'),width = 3,height = 3,dpi = 200)

}
plot_methylation_level_DMS(all_m5C, 'LF')

