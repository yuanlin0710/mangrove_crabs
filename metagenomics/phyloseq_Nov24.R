setwd("C:/Users/All users.DESKTOP-B09KNMV/Documents/FYP_lynn/metagenome analysis/phyloseq")
library("dplyr")
library("ggplot2")
phyloseq_table<-function(SPECIES_NAME){
  taxfilename<-paste0(SPECIES_NAME,"_tax_family.csv")
  otufilename<-paste0(SPECIES_NAME,"_OTU_family.csv")
  colname1<-paste0(SPECIES_NAME,"_tax_long.family")
  CRtax<-read.csv(taxfilename)
  CROTU<-read.csv(otufilename)
  
  CR_taxtable<-merge(CROTU, CRtax, by.x="X", by.y="X")
  CR_taxtable<-subset(CR_taxtable, select = c("X", colname1))
  colnames(CR_taxtable)[2]<-"tax_long"
  CR_taxtable$k<-gsub("\\;p_.*","",CR_taxtable$tax_long)
  CR_taxtable$p<-gsub("\\;c_.*","",CR_taxtable$tax_long)
  CR_taxtable$c<-gsub("\\;o_.*","",CR_taxtable$tax_long)
  CR_taxtable$o<-gsub("\\;f_.*","",CR_taxtable$tax_long)
  CR_taxtable$f<-CR_taxtable$tax_long
  CR_taxtable$p<-gsub(".*\\;p_","p_",CR_taxtable$p)
  CR_taxtable$c<-gsub(".*\\;c_","c_",CR_taxtable$c)
  CR_taxtable$o<-gsub(".*\\;o_","o_",CR_taxtable$o)
  CR_taxtable$f<-gsub(".*\\;f_","f_",CR_taxtable$f)
  rownames(CR_taxtable)<-CR_taxtable$X
  CR_taxtable<-subset(CR_taxtable, select = c("k","p","c","o","f"))
  
  return(CR_taxtable)
  
}



EV_taxtable<-phyloseq_table("EV")
MF_taxtable<-phyloseq_table("MF")
PB_taxtable<-phyloseq_table("PB")
TD_taxtable<-phyloseq_table("TD")

taxtable_all<-rbind(EV_taxtable, MF_taxtable, PB_taxtable, TD_taxtable)



taxtable_all<-distinct(taxtable_all)


EVOTU<-read.csv("EV_OTU_family.csv")
MFOTU<-read.csv("MF_OTU_family.csv")
PBOTU<-read.csv("PB_OTU_family.csv")
TDOTU<-read.csv("TD_OTU_family.csv")



OTU_all<-merge(EVOTU, MFOTU, by.x = "X", by.y="X", all = TRUE)
OTU_all<-merge(OTU_all, PBOTU, by.x = "X", by.y="X", all = TRUE)
OTU_all<-merge(OTU_all, TDOTU, by.x = "X", by.y="X", all = TRUE)
OTU_all[is.na(OTU_all)]<-0

rownames(OTU_all)<-OTU_all$X
OTU_all<-subset(OTU_all, select=c("EV1","EV2","EV3","MF1","MF2","MF3","PB1","PB2","PB3","TD1","TD2","TD3"))
write.csv(taxtable_all,"taxtable_all_Nov24.csv")
write.csv(OTU_all,"OTUtable_all_Nov24.csv")





library("phyloseq")
otu<-read.csv("OTUtable_all_Nov24.csv")
tax<-read.csv("taxtable_all_Nov24.csv")

rownames(otu)<-otu$X
rownames(tax)<-tax$X
otu$X<-NULL
tax$X<-NULL

otu_matrix<-as.matrix(otu)
tax_matrix<-as.matrix(tax)


OTU <- otu_table(otu_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(tax_matrix)

physeq = phyloseq(OTU, TAX)
# this is taxa abundance plot; but already get using sqmtools
#plot_bar(physeq, fill = "c")
#plot_bar(physeq, fill = "p")

plot_richness(physeq)
ggsave("alpha_diversity_phyloseq.pdf")

# alpha diversity



# phyloseq
# follow https://joey711.github.io/phyloseq/plot_ordination-examples.html
library("phyloseq")
library("ggplot2")
library("plyr")

sample<-read.table("metadata.txt", header=TRUE)
sample$sample<-row.names(sample)
sample<-sample_data(sample)


physeq2 = phyloseq(OTU, TAX, sample)
plot_heatmap(physeq2, sample.order = c("EV1","EV2","EV3","PB1","PB2","PB3","MF1","MF2","MF3","TD1","TD2","TD3"))



wh0 = genefilter_sample(physeq2, filterfun_sample(function(x) x > 5), A=0.5*nsamples(physeq2))
physeq3 = prune_taxa(wh0, physeq2)
physeq3 = transform_sample_counts(physeq3, function(x) 1E6 * x/sum(x))

phylum.sum = tapply(taxa_sums(physeq3), tax_table(physeq3)[, "p"], sum, na.rm=TRUE)
top9phyla = names(sort(phylum.sum, TRUE))[1:9]
physeq3 = prune_taxa((tax_table(physeq3)[, "p"] %in% top9phyla), physeq3)

physeq.ord <- ordinate(physeq3, "DCA", "bray")
p1 = plot_ordination(physeq3, physeq.ord, type="taxa", color="p", title="taxa")
print(p1)
p1 + facet_wrap(~p, 3)
plot_heatmap(physeq3, sample.order = c("EV1","EV2","EV3","PB1","PB2","PB3","MF1","MF2","MF3","TD1","TD2","TD3"))





mangrove = get_variable(physeq2, "diet") %in% c("omnivore","herbivore")
sample_data(physeq3)$mangrove <- factor(mangrove)

p2 = plot_ordination(physeq3, physeq.ord, type="sample", color="species", shape="mangrove") 
p2 + geom_polygon(aes(fill=species)) + geom_point(size=5) + ggtitle("samples")


