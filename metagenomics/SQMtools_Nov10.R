library('SQMtools')
library('ggplot2')

crab<-loadSQM("EV",engine="data.table")

#TD is run on SQMtools 0.6.3 (on Mac computer)
#EV/MF/PB are run on server's R; SQMtools v1.6.0


# see what is the read composition (how many host reads, how many bacterial reads)
crab$taxa$superkingdom$abund
crab$taxa$superkingdom$percent

subtax_bac<-subsetTax(crab, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)



phylum_table<-subtax_bac$taxa$phylum$abund
class_table<-subtax_bac$taxa$class$abund
order_table<-subtax_bac$taxa$order$abund

#write.csv(phylum_table, "subbac_phylum_table.csv")
#write.csv(class_table, "subbac_class_table.csv")
#write.csv(order_table, "subbac_order_table.csv")

summary_bac<-summary(subtax_bac)


plotTaxonomy(subtax_bac, rank='phylum', count='percent', rescale=T)
#ggsave("plotTax_phylum.pdf")
plotTaxonomy(subtax_bac, rank='class', count='percent', rescale=T)
#ggsave("plotTax_class.pdf")
plotTaxonomy(subtax_bac, rank='order', count='percent', rescale=T)
#ggsave("plotTax_order.pdf")



#---------------------------------------------------------------------

EV<-loadSQM("EV",engine="data.table")


# see what is the read composition (how many host reads, how many bacterial reads)
EV_superkingdom_abund<-EV$taxa$superkingdom$abund
#write.csv(EV_superkingdom_abund,"EV_taxa_superkingdom_abund.csv")
EV_superkingdom_percent<-EV$taxa$superkingdom$percent
#write.csv(EV_superkingdom_percent,"EV_taxa_superkingdom_percent.csv")



EV_subtax_bac<-subsetTax(EV, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)



EV_phylum_table<-EV_subtax_bac$taxa$phylum$abund
EV_class_table<-EV_subtax_bac$taxa$class$abund
EV_order_table<-EV_subtax_bac$taxa$order$abund

#write.csv(EV_phylum_table, "EV_subbac_phylum_table.csv")
#write.csv(EV_class_table, "EV_subbac_class_table.csv")
#write.csv(EV_order_table, "EV_subbac_order_table.csv")


plotTaxonomy(EV_subtax_bac, rank='phylum', count='percent', rescale=T)
#ggsave("EV_plotTax_phylum.pdf")
plotTaxonomy(EV_subtax_bac, rank='class', count='percent', rescale=T)
#ggsave("EV_plotTax_class.pdf")
plotTaxonomy(EV_subtax_bac, rank='order', count='percent', rescale=T)
#ggsave("EV_plotTax_order.pdf")



plotFunctions(EV_subtax_bac, fun_level = 'KEGG', count = 'tpm', N = 10)
#ggsave("EV_plotFun_KEGG_tpm.pdf")
plotFunctions(EV_subtax_bac, fun_level = 'cazy', count = 'tpm', N = 10)
#ggsave("EV_plotFun_cazy_tpm.pdf")





EV_cazy <- subsetFun(EV_subtax_bac, fun = cazy_list, rescale_copy_number = F) 
plotTaxonomy(EV_cazy, rank='phylum', count='percent', N = 10, rescale = T)












