# MATRIZ-DE-EXPRESION

#Cargar las librerías
library("TCGAbiolinks")
library("limma")
library("edgeR")
library("glmnet")
library("factoextra")
library("FactoMineR")
library("caret")
library("SummarizedExperiment")
library("gplots")
library("survival")
library("survminer")
library("RColorBrewer")
library("gProfileR")
library("genefilter")

#En primer intancia se deben conocer todos los proyectos a los que 
#se pueden tener acceso desde el TCGA mediante la función "getGDCprojects"
#la cual la contendremos en el objeto "GDCprojects".

GDCprojects = getGDCprojects()

#Se obtendrán los primeros 10 datos del objeto aplicandole la función head
#y también la función dim() para obtener la dimensión de los datos.

head(GDCprojects[c("project_id", "name")])

dim(GDCprojects)

#Con la función getProjectSummary() se obtendrá toda la información que 
#comprende el proyecto "TCGA-GBM".

TCGAbiolinks:::getProjectSummary("TCGA-GBM")

#Con la funcion GDCquery() se especifican los datos del proyecto con los
#que vamos a estar trabajando los cuales los vamos a poner en el objeto
#"query_TCGA"

query_TCGA = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", 
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts")

#Ahora necesitamos visualizar los datos, para ello, al objeto "query_TCGA" 
#se le aplicará la función getResults() para obtenerlos en forma de tabla 
#y se contendrán en el objeto "lihc_res"

lihc_res = getResults(query_TCGA)

#Con la función head() se visualizarán solo los datos de los primeros 6
#pacientes y con la función colnames() se visualizarán todas las columnas
#presentes en la tabla.

head(lihc_res)

colnames(lihc_res) 

#Ahora se visualizarán solo los primeros tipos de muestras que contiene nuestro
#objeto "lihc_res"

head(lihc_res$sample_type)

#En este paso se pretende obtener la información de todas las variables
#contenidas en la columna sample_type del objeto lihc_res, sin embargo, 
#al ejecutar la función solo me debuelve la cantidad de datos que contiene
#la columna y no todos los distintos tipos de variables

summary(lihc_res$sample_type)#hay que revisar este comando

#El siguiente paso consta en obtener un objeto el cual contenga solamente 
#la información correspondiente a las muestras y nuestros controles, ignorando
#cualquier otro dato que no se encuentre en los paámetros, esto se realizará 
#gracias a la función GDCquery()

query_TCGA = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "HTSeq - Counts",
  sample.type = c("Primary Tumor",
                  "Solid Tissue Normal",
                  "Recurrent Tumor"))

#Ahora, el siguiente paso es descargar los datos que hemos minado hasta este 
#punto con el fin de poder almacenarlos localmente, lo cual se logrará mediante
#la función GDCdownload()

GDCdownload(query = query_TCGA)

#En este paso se descargarán los datos reales de RNASeq al equipo mediante la 
#siguiente función.

tcga_data = GDCprepare(query_TCGA)

#Visualizar el tamaño del objeto con la función dim()

dim(tcga_data)
