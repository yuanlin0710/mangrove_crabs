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





#________________________________________________________________________________________________




# the bottom is the new code














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











PB<-loadSQM("PB",engine="data.table")


# see what is the read composition (how many host reads, how many bacterial reads)
PB_superkingdom_abund<-PB$taxa$superkingdom$abund
#write.csv(PB_superkingdom_abund,"PB_taxa_superkingdom_abund.csv")
PB_superkingdom_percent<-PB$taxa$superkingdom$percent
#write.csv(PB_superkingdom_percent,"PB_taxa_superkingdom_percent.csv")



PB_subtax_bac<-subsetTax(PB, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)



PB_phylum_table<-PB_subtax_bac$taxa$phylum$abund
PB_class_table<-PB_subtax_bac$taxa$class$abund
PB_order_table<-PB_subtax_bac$taxa$order$abund

#write.csv(PB_phylum_table, "PB_subbac_phylum_table.csv")
#write.csv(PB_class_table, "PB_subbac_class_table.csv")
#write.csv(PB_order_table, "PB_subbac_order_table.csv")


plotTaxonomy(PB_subtax_bac, rank='phylum', count='percent', rescale=T)
#ggsave("PB_plotTax_phylum.pdf")
plotTaxonomy(PB_subtax_bac, rank='class', count='percent', rescale=T)
#ggsave("PB_plotTax_class.pdf")
plotTaxonomy(PB_subtax_bac, rank='order', count='percent', rescale=T)
#ggsave("PB_plotTax_order.pdf")


plotFunctions(PB_subtax_bac, fun_level = 'KEGG', count = 'tpm', N = 10)
#ggsave("PB_plotFun_KEGG_tpm.pdf")
plotFunctions(PB_subtax_bac, fun_level = 'cazy', count = 'tpm', N = 10)
#ggsave("PB_plotFun_cazy_tpm.pdf")







MF<-loadSQM("MF",engine="data.table")


# see what is the read composition (how many host reads, how many bacterial reads)
MF_superkingdom_abund<-MF$taxa$superkingdom$abund
#write.csv(MF_superkingdom_abund,"MF_taxa_superkingdom_abund.csv")
MF_superkingdom_percent<-MF$taxa$superkingdom$percent
#write.csv(MF_superkingdom_percent,"MF_taxa_superkingdom_percent.csv")



MF_subtax_bac<-subsetTax(MF, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)



MF_phylum_table<-MF_subtax_bac$taxa$phylum$abund
MF_class_table<-MF_subtax_bac$taxa$class$abund
MF_order_table<-MF_subtax_bac$taxa$order$abund

#write.csv(MF_phylum_table, "MF_subbac_phylum_table.csv")
#write.csv(MF_class_table, "MF_subbac_class_table.csv")
#write.csv(MF_order_table, "MF_subbac_order_table.csv")


plotTaxonomy(MF_subtax_bac, rank='phylum', count='percent', rescale=T)
#ggsave("MF_plotTax_phylum.pdf")
plotTaxonomy(MF_subtax_bac, rank='class', count='percent', rescale=T)
#ggsave("MF_plotTax_class.pdf")
plotTaxonomy(MF_subtax_bac, rank='order', count='percent', rescale=T)
#ggsave("MF_plotTax_order.pdf")




plotFunctions(MF_subtax_bac, fun_level = 'KEGG', count = 'tpm', N = 10)
#ggsave("MF_plotFun_KEGG_tpm.pdf")
plotFunctions(MF_subtax_bac, fun_level = 'cazy', count = 'tpm', N = 10)
#ggsave("MF_plotFun_cazy_tpm.pdf")






TD<-loadSQM("TD1_3")


# see what is the read composition (how many host reads, how many bacterial reads)
TD_superkingdom_abund<-TD$taxa$superkingdom$abund
#write.csv(TD_superkingdom_abund,"TD_taxa_superkingdom_abund.csv")
TD_superkingdom_percent<-TD$taxa$superkingdom$percent
#write.csv(TD_superkingdom_percent,"TD_taxa_superkingdom_percent.csv")



TD_subtax_bac<-subsetTax(TD, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)



TD_phylum_table<-TD_subtax_bac$taxa$phylum$abund
TD_class_table<-TD_subtax_bac$taxa$class$abund
TD_order_table<-TD_subtax_bac$taxa$order$abund

#write.csv(TD_phylum_table, "TD_subbac_phylum_table.csv")
#write.csv(TD_class_table, "TD_subbac_class_table.csv")
#write.csv(TD_order_table, "TD_subbac_order_table.csv")


plotTaxonomy(TD_subtax_bac, rank='phylum', count='percent', rescale=T)
#ggsave("TD_plotTax_phylum.pdf")
plotTaxonomy(TD_subtax_bac, rank='class', count='percent', rescale=T)
#ggsave("TD_plotTax_class.pdf")
plotTaxonomy(TD_subtax_bac, rank='order', count='percent', rescale=T)
#ggsave("TD_plotTax_order.pdf")


plotFunctions(TD_subtax_bac, fun_level = 'KEGG', count = 'tpm', N = 10)
#ggsave("TD_plotFun_KEGG_tpm.pdf")
plotFunctions(TD_subtax_bac, fun_level = 'cazy', count = 'tpm', N = 10)
#ggsave("TD_plotFun_cazy_tpm.pdf")





