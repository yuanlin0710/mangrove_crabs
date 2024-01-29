setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/merge_DEDAtable_from_SQM")

# merge tables together and filter out "unclassified" stuff

EV<-read.csv("EV_subbac_phylum_table.csv")
MF<-read.csv("MF_subbac_phylum_table.csv")
PB<-read.csv("PB_subbac_phylum_table.csv")
TD<-read.csv("TD_subbac_phylum_table.csv")

all4DA<-merge(EV, MF, by="X", all=TRUE)
all4DA<-merge(all4DA, PB, by="X", all=TRUE)
all4DA<-merge(all4DA, TD, by="X", all=TRUE)


all4DA[is.na(all4DA)]<-0

write.csv(all4DA, "merged_table_DA.csv")

#manually remove rows with "unclassified" & change names if the first char is special char or number (may cause DESeq2 running error)]
#also, correct the name of TD

all4DA<-read.csv("merged_table_DA.csv", row.names = "X.1")




# DESeq2
library(DESeq2)
library(data.table)
library(ggplot2)
library(ggrepel)
library(apeglm)
library(vegan)
library(pheatmap)




setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/DA/deseq2")
metadata <- read.table("metadata.txt", header = TRUE, stringsAsFactors=FALSE)

metadata$species = as.factor(metadata$species)
metadata$species = relevel(metadata$species, ref = "TD")


rawcount<-all4DA
row.names(rawcount)<-rawcount$X
rawcount<-subset(rawcount, select=-c(X))

rawcount <- rawcount[rowMeans(rawcount)>5,]

dds <- DESeqDataSetFromMatrix(countData = rawcount,
                              colData = metadata,
                              design = ~ species)


dds = DESeq(dds)
cbind(resultsNames(dds))
dds

res1 <- results(dds, name = "species_EV_vs_TD", alpha = 0.05)
summary(res1)
resLFC1 <- lfcShrink(dds, coef="species_EV_vs_TD", type="apeglm")
View(as.data.frame(resLFC1))

resOrdered1<- res1[order(res1$padj),]
summary(resOrdered1)


res2 <- results(dds, name = "species_MF_vs_TD", alpha = 0.05)
summary(res2)
resLFC2 <- lfcShrink(dds, coef="species_MF_vs_TD", type="apeglm")
View(as.data.frame(resLFC2))

res3 <- results(dds, name = "species_PB_vs_TD", alpha = 0.05)
summary(res3)
resLFC3 <- lfcShrink(dds, coef="species_PB_vs_TD", type="apeglm")
View(as.data.frame(resLFC3))


res_signif <- results(dds, alpha=0.05)
summary(res_signif)
View(as.data.frame(res_signif))



#plotMA(res1, ylim=c(-10,10))
#plotMA(resLFC1, ylim=c(-10,10))

#write.csv(resOrdered, file = "DESeq2_taxa_all.csv")
#write.csv(res_signif, file = "DESeq2_taxa_alpha_sig.csv")
#write.csv(resLFC, file = "DESeq2_taxa_LFC.csv")



write.csv(resLFC1,"DA_EVvsTD_deseq2.csv")
write.csv(resLFC2,"DA_MFvsTD_deseq2.csv")
write.csv(resLFC3,"DA_PBvsTD_deseq2.csv")

EV_DA<-as.data.frame(resLFC1)
MF_DA<-as.data.frame(resLFC2)
PB_DA<-as.data.frame(resLFC3)

EV_DA<-EV_DA[EV_DA$padj<0.05,]
EV_DA<-EV_DA[!is.na(EV_DA$padj),]

MF_DA<-MF_DA[MF_DA$padj<0.05,]
MF_DA<-MF_DA[!is.na(MF_DA$padj),]

PB_DA<-PB_DA[PB_DA$padj<0.05,]
PB_DA<-PB_DA[!is.na(PB_DA$padj),]


write.csv(EV_DA,"DA_EVvsTD_sigOnly_deseq2.csv")
write.csv(MF_DA,"DA_MFvsTD_sigOnly_deseq2.csv")
write.csv(PB_DA,"DA_PBvsTD_sigOnly_deseq2.csv")


# Venn diagram
library("ggvenn")

A <-list('EV'=rownames(EV_DA),'MF'=rownames(MF_DA), 'PB' = rownames(PB_DA))
ggvenn(A)
#ggsave("DA_venn_3mangrove_deseq2.pdf")


# PCA
rld <- rlog(dds)
plotPCA(rld, intgroup ="species")->PCA_w_label
PCA_w_label<-PCA_w_label+geom_text(aes(label=colnames(assay(rld))), size=2)
PCA_w_label

#sample distances & heatmap
rld.t = t(assay(rld))
sampleDists <-  vegdist(rld.t,method="euclidean")
geneDists <- dist(t(rld.t))
clustergene <- hclust(sampleDists, method="average")
clustersample <- hclust(geneDists, method="average")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rownames(metadata)
colnames(sampleDistMatrix) <- NULL
deseq.sample2geneheatmap <- pheatmap(as.data.frame(rld.t),cluster_rows = TRUE,cluster_cols = TRUE,main = "Clustered heatmap of samples.DESeq2",border_color = NA,annotation_legend = T, show_colnames = F)


rm(list=ls())







#aldex2
setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/merge_DEDAtable_from_SQM")
all4DA<-read.csv("merged_table_DA.csv", row.names = "X.1")

setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/DA/aldex2")



library(ALDEx2)
library(CoDaSeq)


# sample input format: selex
#data(selex)
abun_table<-all4DA

row.names(abun_table)<-abun_table$X
abun_table$X<-NULL


conds <- c(rep("mangrove",9),rep("non-mangrove",3))
mm<-model.matrix(~conds)
x.all <- aldex(abun_table, conds, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")


filtered_abun_table<- codaSeq.filter(abun_table, min.reads=1000, min.occurrence=0.99, samples.by.row=FALSE)
f.x <- aldex.clr(filtered_abun_table, mm, denom="all", verbose=F)

glm.test <- aldex.glm(f.x,mm)
featurelist <- rownames(glm.test)


clr <- aldex.effect(f.x, conds,include.sample.summary=TRUE,CI = T,useMC = F,verbose = T)
#log abundance
E.E.clr <- t(clr[,grep("rab.sample", colnames(clr))])
rownames(E.E.clr) <- gsub("rab.sample.", "", rownames(E.E.clr))
exp <- apply(E.E.clr, 1, function(x) 2^x)
E.clr <- t(apply(exp, 2, function(x) log2(x) - mean(log2(x)) )) %>% as.matrix()
#filter genes for plotting graph

#PCA
pcx <- prcomp(E.clr)
pc1.aldex <- round(pcx$sdev[1]^2/sum(pcx$sdev^2),2)
pc2.aldex <- round(pcx$sdev[2]^2/sum(pcx$sdev^2),2)
xlab.aldex <- paste("PC1: ", pc1.aldex, sep="")
ylab.aldex <- paste("PC2: ", pc2.aldex, sep="")
en.aldex = envfit(pcx,E.clr, permutations = 999, na.rm = TRUE)
scores.aldex <- as.data.frame(scores(pcx$x))
loadings.aldex <- as.data.frame(pcx$rotation) 
en_coord_cont.aldex = as.data.frame(scores(en.aldex, "vectors")) * ordiArrowMul(en.aldex,fill = scores.aldex$PC1)
speciescol=row.names(scores.aldex)
gg.aldex = ggplot(data = scores.aldex, aes(x = PC1, y = PC2, label= speciescol)) + 
  geom_point(data = scores.aldex, aes(colour = speciescol), size = 3, alpha = 0.5) +
  geom_text(hjust=1, vjust=0, size=2.5)+
  #  stat_ellipse(aes(fill= speciescol), alpha=.2,type='t',linewidth=1, geom="polygon")+ 
  theme_minimal()+
  #  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
  #               data = en_coord_cont.aldex, linewidth =1, alpha = 0.5, colour = "grey30") +
  #  geom_text(data = en_coord_cont.aldex, aes(x = PC1, y = PC2), colour = "grey30", 
  #            fontface = "bold", label = row.names(en_coord_cont.aldex),size = 2) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"),    
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"))

color_condition<-c("EV1"="#7CAE00","EV2"="#7CAE00","EV3"="#7CAE00",
                   "PB1"="#7CAE00","PB2"="#7CAE00","PB3"="#7CAE00",
                   "MF1"="#964B00","MF2"="#964B00","MF3"="#964B00",
                   "TD1"="#F8766D","TD2"="#F8766D","TD3"="#F8766D")

gg.aldex+scale_color_manual(values=color_condition)+scale_fill_manual(values=color_condition)+xlab(xlab.aldex)+ylab(ylab.aldex)+ ggtitle("PCA plot.ALDEx2") -> gg.aldex
gg.aldex



#heatmap for abundance
aldex.sample2geneheatmap <- pheatmap(E.clr,cluster_rows = TRUE,cluster_cols = TRUE,border_color = NA,annotation_legend = T)


# DE table
DErow<-clr[,"effect"]
DErow<-data.frame(DErow)
row.names(DErow)<-row.names(clr)
colnames(DErow)[1]<-"aldex2_effect"
DErow$significant_effect<-FALSE


# for taxa
for (row in 1:nrow(DErow)){
  value_tmp<-DErow[row,1]
  DErow[row, ]$significant_effect<-(-1>value_tmp | value_tmp>1)
}

write.csv(DErow,"DA_sigOnly_aldex2.csv")
write.csv(glm.test,"DA_glmtest_all_aldex2.csv")




rm(list=ls())









# ancombc
setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/merge_DEDAtable_from_SQM")
all4DA<-read.csv("merged_table_DA.csv", row.names = "X.1")

setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/DA/ancombc")


