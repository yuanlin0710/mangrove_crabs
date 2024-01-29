# DESeq2
library(DESeq2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(apeglm)
library(vegan)
library(pheatmap)
library(dplyr)


setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/DA/deseq2")
metadata <- read.table("metadata.txt", header = TRUE, stringsAsFactors=FALSE)
all4DA<-read.csv("rescale_readcount_deseq2_forboxplot.csv", row.names = "X")



# how comes the all4DA:
# the re-scale value is got by dds$sizeFactor in deseq2
# rescaled count in each entry = its raw count in merged_table_DA / rescale value
# so that the overall read count is normalized

metadata$species = as.factor(metadata$species)
metadata$species = relevel(metadata$species, ref = "TD")



rawcount<-all4DA

rawcount <- rawcount[rowMeans(rawcount)>5,]


boxplotNoColorWdodgedJitter <- function(STRINGgenename){
  library(dplyr)
  library(ggplot2)
  # use STRINGgenename to pass the STRINGgenename
  # remember to pass the STRINGgenename as a STRING (use quotation marks "")!!!!
  
  sum_df<-rawcount[rownames(rawcount)==STRINGgenename,]
  sum_df<-t(sum_df)
  sum_df<-data.frame(sum_df)
  
  metadata <- read.table("metadata.txt", header = TRUE, stringsAsFactors=FALSE)
  merged_df<-merge(metadata, sum_df, by.x=0, by.y=0)
  colnames(merged_df)[4] <- "abundance"
  saved_filename<-paste0(STRINGgenename,"_abundBoxplot.png")
  p<-ggplot(merged_df, aes(x=species, y=abundance)) +
    geom_boxplot(width=0.4) 
  #fill = c("#4E84C4","#CE5143","#C3D7A4")
  p <- p + ggtitle(STRINGgenename)
  p <- p + ylab("abundance")
  p <- p + xlab(NULL)
  p <- p +
    theme(
      panel.background = element_rect(fill='transparent'), #transparent panel bg
      #plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
      panel.grid.major = element_blank(), #remove major gridlines
      panel.grid.minor = element_blank(), #remove minor gridlines
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      #legend.background = element_rect(fill='transparent'), #transparent legend bg
      #legend.box.background = element_rect(fill='transparent') #transparent legend panel
      plot.title = element_text(size=70, face = "italic"), #change plot title font size
      axis.text = element_text(size = 35),
      axis.title.y = element_text(size=50), # change y title font size
      axis.title.x = element_text(size=50), # change x title font size
    )
  p<- p + geom_point(aes(color=diet),size=8)+scale_color_manual(values=c("#8b0000","#008000","#964B00"))
  p<- p + theme (legend.position = "none")
  #ggsave(filename=saved_filename,plot=p, width = 10, height = 10)
  #fill = c("#4E84C4","#D16103","#C3D7A4")
  p
}

#boxplotNoColorWdodgedJitter("Proteobacteria")
boxplotNoColorWdodgedJitter("Acidobacteria")
#boxplotNoColorWdodgedJitter("Bacteroidetes")
boxplotNoColorWdodgedJitter("Chloroflexi")
boxplotNoColorWdodgedJitter("Actinobacteria")
boxplotNoColorWdodgedJitter("Candidatus Staskawiczbacteria")
boxplotNoColorWdodgedJitter("Planctomycetes")
