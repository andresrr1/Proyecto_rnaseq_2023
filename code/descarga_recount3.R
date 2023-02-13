#Descargar los datos del proyecto
#Generar el objeto 'RangeSummarizedObject'
library(recount3)
library(iSEE)

#Verificar el URL actual
getOption(
  "recount3_url",
  "http://duffel.rail.bio/recount3"
)
#Cambiar el URL 
options(recount3_url = "https://recount-opendata.s3.amazonaws.com/recount3/release")
#Buscar los proyectos realizados con celulas humanas
human_projects <- available_projects(organism = 'human')
#Obtener el proyecto de interes
proj_info <- subset(
  human_projects,
  project == "ERP110066" & project_type == "data_sources"
)
#Crear el objeto RSE
rse_gene_ERP110066 <- create_rse(proj_info)
#Convertir las cuentas por n a cuentas por lectura

assay(rse_gene_ERP110066, "counts") <- compute_read_counts(rse_gene_ERP11006)
#Generar un data frame mas facil de acceder
rse_gene_ERP110066 <- expand_sra_attributes(rse_gene_ERP110066)
colData(rse_gene_ERP110066)[
  ,
  grepl("^sra_attribute", colnames(colData(rse_gene_ERP110066)))
]
#Se guarda el objeto RSE en un archivo
save(rse_gene_ERP110066, file="rse_gene_ERP110066.RData")
#Se explora la informacion de filas y columnas
colData(rse_gene_ERP110066)
rowData(rse_gene_ERP110066)
#Se exploran los datos de forma interactiva
iSee::iSee(rse_gene_ERP110066)

