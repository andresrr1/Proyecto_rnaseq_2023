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

ggplot(as.data.frame(colData(rse_gene_ERP110066)), aes(y = assigned_gene_prop, x = sra_attribute.growth_condition)) +
  geom_boxplot() +
  theme_bw(base_size = 20) +
  ylab("Assigned Gene Prop") +
  xlab("Growth condition")

mod <- model.matrix(~sra_attribute.growth_condition + assigned_gene_prop,
                    data = colData(rse_gene_ERP110066))

colnames(mod)


mod

vGenes <- voom(DGE, mod, plot = T)

eb_results <- eBayes(lmFit(vGenes))

de_results <- topTable(
  eb_results,
  coef = 2,
  number = nrow(rse_gene_ERP110066),
  sort.by = "p"
)


head(de_results)
dim(de_results)

table(de_results$adj.P.Val<0.05)


plotMA(eb_results, coef = 2)
volcanoplot(eb_results, coef = 2, highlight = 3, names = de_results$gene_name)


exprs_heatmap <- vGenes$E[rank(de_results$adj.P.Val) <= 50,]

rownames(exprs_heatmap) <- de_results[row.names(exprs_heatmap),"gene_name"]


df <- as.data.frame(colData(rse_gene_ERP110066)[, c("sra_attribute.growth_condition", "sra_attribute.stimulus")])
colData(rse_gene_ERP110066)[, "sra_attribute.growth_condition","sra_attribute.stimulus" ]


colnames(df) <- c("growth_condition", "stimulus")

View(df)

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = df
)

plotMDS(vGenes$E, labels = df$growth_condition, col = c("blue","red"))

