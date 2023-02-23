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

hist(rse_gene_ERP110066$assigned_gene_prop, xlab = "Assigned gene proportion", main = 'Assigned gene proportion')


summary(as.data.frame(colData(rse_gene_ERP110066)[
  ,
  grepl("^sra_attribute.[growt_condition|stimulus]", colnames(colData(rse_gene_ERP110066)))
]))


rse_gene_ERP110066$assigned_gene_prop <- rse_gene_ERP110066$recount_qc.gene_fc_count_all.assigned/rse_gene_ERP110066$recount_qc.gene_fc_count_all.total

summary(rse_gene_ERP110066$assigned_gene_prop)

promedios_genes <- rowMeans(assay(rse_gene_ERP110066, "counts"))
summary(promedios_genes)


rse_gene_ERP110066_completo <- rse_gene_ERP110066
rse_gene_ERP110066 <- rse_gene_ERP110066[promedios_genes > 0.1, ]


round(nrow(rse_gene_ERP110066) / nrow(rse_gene_ERP110066_completo) * 100, 2)

DGE <- DGEList(
  counts = assay(rse_gene_ERP110066, "counts"),
  genes = rowData(rse_gene_ERP110066)
)

DGE <- calcNormFactors(DGE)

