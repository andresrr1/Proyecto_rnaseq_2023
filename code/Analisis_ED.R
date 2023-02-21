library("edgeR")
library("ggplot2")
library("limma")
library("pheatmap")

rse_gene_ERP110066 <- expand_sra_attributes(rse_gene_ERP110066)

colData(rse_gene_ERP110066)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_ERP110066)))
]

## Se revisa que los atributos de interes tengan el formato adecuado
rse_gene_ERP110066$sra_attribute.stimulus<-factor(tolower((rse_gene_ERP110066$sra_attribute.stimulus)))
rse_gene_ERP110066$sra_attribute.growth_condition<-factor(tolower((rse_gene_ERP110066$sra_attribute.growth_condition)))