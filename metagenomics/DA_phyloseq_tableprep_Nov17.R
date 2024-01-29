library('SQMtools')
library('ggplot2')

#TD is run on SQMtools 0.6.3 (on Mac computer)
#EV/MF/PB are run on server's R; SQMtools v1.6.0




# cazy gene list for chordtable
tablecazy<-read.table("funlist.c.subfam.txt")
binding<-subset(tablecazy, Substrate == 'Binding')
cellulose<-subset(tablecazy, Substrate == 'Cellulose')
hemicellulose<-subset(tablecazy, Substrate == 'Hemicellulose')
lignin<-subset(tablecazy, Substrate == 'Lignin')
pectin<-subset(tablecazy, Substrate == 'Pectin')



get_genelist<-function(gene_table){
  str<-''
  for (thing in gene_table) {
    str<-paste(str, thing, sep="\\*?$|")
    # so that can keep cazy annotation with *
    # asterisk (*) means trusted function, which means in SqueezeMeta, the functional assignment was done by both best hit and best average scores, therefore is more reliable
  }
  str<-substring(str, 6)
  str<-paste(str, "\\*?$", sep='')
  return(str)
}

cazy_str<-get_genelist(tablecazy$Family)
binding_str<-get_genelist(binding$Family)
cellulose_str<-get_genelist(cellulose$Family)
hemicellulose_str<-get_genelist(hemicellulose$Family)
lignin_str<-get_genelist(lignin$Family)
pectin_str<-get_genelist(pectin$Family)











TD<-loadSQM("TD",engine="data.table")
TD_subtax_bac<-subsetTax(TD, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)
# TD is done on another computer


TD_b<-subsetFun(TD_subtax_bac,binding_str, columns = "cazy", rescale_tpm = TRUE)
TD_c<-subsetFun(TD_subtax_bac,cellulose_str, columns = "cazy", rescale_tpm = TRUE)
TD_h<-subsetFun(TD_subtax_bac,hemicellulose_str, columns = "cazy", rescale_tpm = TRUE)
TD_l<-subsetFun(TD_subtax_bac,lignin_str, columns = "cazy", rescale_tpm = TRUE)
TD_p<-subsetFun(TD_subtax_bac,pectin_str, columns = "cazy", rescale_tpm = TRUE)
TD_cazy<-subsetFun(TD_subtax_bac,cazy_str, columns = "cazy", rescale_tpm = TRUE)

table<-data.frame(TD_cazy$functions$cazy$tpm)
write.csv(table,"TD_cazyGOIsubset_cazy_tpm.csv")



table<-data.frame(TD_b$taxa$phylum$abund)
write.csv(table,"TD_bindingsubset_phylum_abund.csv")
table<-data.frame(TD_b$taxa$phylum$percent)
write.csv(table,"TD_bindingsubset_phylum_percent.csv")

table<-data.frame(TD_c$taxa$phylum$abund)
write.csv(table,"TD_cellulosesubset_phylum_abund.csv")
table<-data.frame(TD_c$taxa$phylum$percent)
write.csv(table,"TD_cellulosesubset_phylum_percent.csv")

table<-data.frame(TD_h$taxa$phylum$abund)
write.csv(table,"TD_hemicellulosesubset_phylum_abund.csv")
table<-data.frame(TD_h$taxa$phylum$percent)
write.csv(table,"TD_hemicellulosesubset_phylum_percent.csv")

table<-data.frame(TD_l$taxa$phylum$abund)
write.csv(table,"TD_ligninsubset_phylum_abund.csv")
table<-data.frame(TD_l$taxa$phylum$percent)
write.csv(table,"TD_ligninsubset_phylum_percent.csv")

table<-data.frame(TD_p$taxa$phylum$abund)
write.csv(table,"TD_pectinsubset_phylum_abund.csv")
table<-data.frame(TD_p$taxa$phylum$percent)
write.csv(table,"TD_pectinsubset_phylum_percent.csv")




# prepare for phyloseq
# get the OTU and tax table at family level
TD_family_abund<-data.frame(TD_subtax_bac$taxa$family$abund)
TD_tax_long<-TD_subtax_bac$misc$tax_names_long
TD_tax_family<-data.frame(TD_tax_long$family)



write.csv(TD_tax_family, "TD_tax_family.csv",row.names = TRUE)
write.csv(TD_family_abund, "TD_OTU_family.csv")




# prepare for DE of gene
table<-TD_subtax_bac$functions$cazy$tpm
write.csv(table,"TD_cazyall_DESeq2.csv")









EV<-loadSQM("EV",engine="data.table")
EV_subtax_bac<-subsetTax(EV, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)


EV_b<-subsetFun(EV_subtax_bac,binding_str, columns = "cazy", rescale_tpm = TRUE)
EV_c<-subsetFun(EV_subtax_bac,cellulose_str, columns = "cazy", rescale_tpm = TRUE)
EV_h<-subsetFun(EV_subtax_bac,hemicellulose_str, columns = "cazy", rescale_tpm = TRUE)
EV_l<-subsetFun(EV_subtax_bac,lignin_str, columns = "cazy", rescale_tpm = TRUE)
EV_p<-subsetFun(EV_subtax_bac,pectin_str, columns = "cazy", rescale_tpm = TRUE)
EV_cazy<-subsetFun(EV_subtax_bac,cazy_str, columns = "cazy", rescale_tpm = TRUE)

table<-data.frame(EV_cazy$functions$cazy$tpm)
write.csv(table,"EV_cazyGOIsubset_cazy_tpm.csv")



table<-data.frame(EV_b$taxa$phylum$abund)
write.csv(table,"EV_bindingsubset_phylum_abund.csv")
table<-data.frame(EV_b$taxa$phylum$percent)
write.csv(table,"EV_bindingsubset_phylum_percent.csv")

table<-data.frame(EV_c$taxa$phylum$abund)
write.csv(table,"EV_cellulosesubset_phylum_abund.csv")
table<-data.frame(EV_c$taxa$phylum$percent)
write.csv(table,"EV_cellulosesubset_phylum_percent.csv")

table<-data.frame(EV_h$taxa$phylum$abund)
write.csv(table,"EV_hemicellulosesubset_phylum_abund.csv")
table<-data.frame(EV_h$taxa$phylum$percent)
write.csv(table,"EV_hemicellulosesubset_phylum_percent.csv")

table<-data.frame(EV_l$taxa$phylum$abund)
write.csv(table,"EV_ligninsubset_phylum_abund.csv")
table<-data.frame(EV_l$taxa$phylum$percent)
write.csv(table,"EV_ligninsubset_phylum_percent.csv")

table<-data.frame(EV_p$taxa$phylum$abund)
write.csv(table,"EV_pectinsubset_phylum_abund.csv")
table<-data.frame(EV_p$taxa$phylum$percent)
write.csv(table,"EV_pectinsubset_phylum_percent.csv")

# prepare for DE
table<-data.frame(EV_subtax_bac$functions$cazy$abund)
write.csv(table,"EV_gene_cazy_all_abund.csv")



# prepare for phyloseq
# get the OTU and tax table at family lEVel
EV_family_abund<-data.frame(EV_subtax_bac$taxa$family$abund)
EV_tax_long<-EV_subtax_bac$misc$tax_names_long
EV_tax_family<-data.frame(EV_tax_long$family)



write.csv(EV_tax_family, "EV_tax_family.csv",row.names = TRUE)
write.csv(EV_family_abund, "EV_OTU_family.csv")

























MF<-loadSQM("MF",engine="data.table")
MF_subtax_bac<-subsetTax(MF, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)


MF_b<-subsetFun(MF_subtax_bac,binding_str, columns = "cazy", rescale_tpm = TRUE)
MF_c<-subsetFun(MF_subtax_bac,cellulose_str, columns = "cazy", rescale_tpm = TRUE)
MF_h<-subsetFun(MF_subtax_bac,hemicellulose_str, columns = "cazy", rescale_tpm = TRUE)
MF_l<-subsetFun(MF_subtax_bac,lignin_str, columns = "cazy", rescale_tpm = TRUE)
MF_p<-subsetFun(MF_subtax_bac,pectin_str, columns = "cazy", rescale_tpm = TRUE)
MF_cazy<-subsetFun(MF_subtax_bac,cazy_str, columns = "cazy", rescale_tpm = TRUE)

table<-data.frame(MF_cazy$functions$cazy$tpm)
write.csv(table,"MF_cazyGOIsubset_cazy_tpm.csv")



table<-data.frame(MF_b$taxa$phylum$abund)
write.csv(table,"MF_bindingsubset_phylum_abund.csv")
table<-data.frame(MF_b$taxa$phylum$percent)
write.csv(table,"MF_bindingsubset_phylum_percent.csv")

table<-data.frame(MF_c$taxa$phylum$abund)
write.csv(table,"MF_cellulosesubset_phylum_abund.csv")
table<-data.frame(MF_c$taxa$phylum$percent)
write.csv(table,"MF_cellulosesubset_phylum_percent.csv")

table<-data.frame(MF_h$taxa$phylum$abund)
write.csv(table,"MF_hemicellulosesubset_phylum_abund.csv")
table<-data.frame(MF_h$taxa$phylum$percent)
write.csv(table,"MF_hemicellulosesubset_phylum_percent.csv")

table<-data.frame(MF_l$taxa$phylum$abund)
write.csv(table,"MF_ligninsubset_phylum_abund.csv")
table<-data.frame(MF_l$taxa$phylum$percent)
write.csv(table,"MF_ligninsubset_phylum_percent.csv")

table<-data.frame(MF_p$taxa$phylum$abund)
write.csv(table,"MF_pectinsubset_phylum_abund.csv")
table<-data.frame(MF_p$taxa$phylum$percent)
write.csv(table,"MF_pectinsubset_phylum_percent.csv")




# prepare for phyloseq
# get the OTU and tax table at family lMFel
MF_family_abund<-data.frame(MF_subtax_bac$taxa$family$abund)
MF_tax_long<-MF_subtax_bac$misc$tax_names_long
MF_tax_family<-data.frame(MF_tax_long$family)



write.csv(MF_tax_family, "MF_tax_family.csv",row.names = TRUE)
write.csv(MF_family_abund, "MF_OTU_family.csv")






# prepare for DE of gene
table<-data.frame(MF_subtax_bac$functions$cazy$abund)
write.csv(table,"MF_gene_cazy_all_abund.csv")
























PB<-loadSQM("PB",engine="data.table")
PB_subtax_bac<-subsetTax(PB, rank = "superkingdom",tax="Bacteria",rescale_tpm = TRUE)


PB_b<-subsetFun(PB_subtax_bac,binding_str, columns = "cazy", rescale_tpm = TRUE)
PB_c<-subsetFun(PB_subtax_bac,cellulose_str, columns = "cazy", rescale_tpm = TRUE)
PB_h<-subsetFun(PB_subtax_bac,hemicellulose_str, columns = "cazy", rescale_tpm = TRUE)
PB_l<-subsetFun(PB_subtax_bac,lignin_str, columns = "cazy", rescale_tpm = TRUE)
PB_p<-subsetFun(PB_subtax_bac,pectin_str, columns = "cazy", rescale_tpm = TRUE)
PB_cazy<-subsetFun(PB_subtax_bac,cazy_str, columns = "cazy", rescale_tpm = TRUE)

table<-data.frame(PB_cazy$functions$cazy$tpm)
write.csv(table,"PB_cazyGOIsubset_cazy_tpm.csv")



table<-data.frame(PB_b$taxa$phylum$abund)
write.csv(table,"PB_bindingsubset_phylum_abund.csv")
table<-data.frame(PB_b$taxa$phylum$percent)
write.csv(table,"PB_bindingsubset_phylum_percent.csv")

table<-data.frame(PB_c$taxa$phylum$abund)
write.csv(table,"PB_cellulosesubset_phylum_abund.csv")
table<-data.frame(PB_c$taxa$phylum$percent)
write.csv(table,"PB_cellulosesubset_phylum_percent.csv")

table<-data.frame(PB_h$taxa$phylum$abund)
write.csv(table,"PB_hemicellulosesubset_phylum_abund.csv")
table<-data.frame(PB_h$taxa$phylum$percent)
write.csv(table,"PB_hemicellulosesubset_phylum_percent.csv")

table<-data.frame(PB_l$taxa$phylum$abund)
write.csv(table,"PB_ligninsubset_phylum_abund.csv")
table<-data.frame(PB_l$taxa$phylum$percent)
write.csv(table,"PB_ligninsubset_phylum_percent.csv")

table<-data.frame(PB_p$taxa$phylum$abund)
write.csv(table,"PB_pectinsubset_phylum_abund.csv")
table<-data.frame(PB_p$taxa$phylum$percent)
write.csv(table,"PB_pectinsubset_phylum_percent.csv")




# prepare for phyloseq
# get the OTU and tax table at family lPBel
PB_family_abund<-data.frame(PB_subtax_bac$taxa$family$abund)
PB_tax_long<-PB_subtax_bac$misc$tax_names_long
PB_tax_family<-data.frame(PB_tax_long$family)



write.csv(PB_tax_family, "PB_tax_family.csv",row.names = TRUE)
write.csv(PB_family_abund, "PB_OTU_family.csv")






# prepare for DE of gene
table<-data.frame(PB_subtax_bac$functions$cazy$abund)
write.csv(table,"PB_gene_cazy_all_abund.csv")
























