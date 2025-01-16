#1. PREPROCESAMIENTO ########################################
#objetivo: Preparar datos fenotípicos y de metilación para el análisis
#Input: Tesis
#Output: qcReport.pdf , preprocesadosAdversidad.RData (datos_fenotípicos + Matriz_Met_Flt + mSetSqFlt + annEPIC)
#Author: karenchurba
#Fecha creación 2024/08/11
#Ultima modificación 2024/14/11
############################################################

#---1.1. Paquetes---------------------------------------------------------------

#if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install ("minfi")
#BiocManager::install ("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#BiocManager::install ("IlluminaHumanMethylationEPICmanifest")
#install.packages("RCurl")
#install.packages("GenomicRanges")
#install.packages("GenomeInfoDb")
#BiocManager::install("DMRcate", update = TRUE)
#BiocManager::install("Gviz")
#install.packages("futile.logger")               


#---1.2. Librerías--------------------------------------------------------------
library(minfi)
library(knitr)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
library(Gviz) 
library(DMRcate)
library(bumphunter)
library(GenomicRanges)
library(EpiDISH)
library(ggplot2)

#---2. Configuración inicial----------------------------------------------------
# Directorio de trabajo
setwd("C:/Users/Karen/Documents/Tesis")
dataDirectory <- "C:/Users/Karen/Documents/Tesis"
#list.files(dataDirectory, recursive = TRUE)

# Cargar anotaciones de EPIC array
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#annEPIC contiene la información de sitios de metilación de genoma humano

#---3.Cargar Datos--------------------------------------------------------------
## 3.1 Cargar datos fenotípicos
datos_fenotipicos <- read.csv("C:/Users/Karen/Documents/Tesis/Datos fenotipicos completa.csv")


##3.2 Leer el sample sheet y datos IDAT
# Cargar información del Sample Sheet
targets <- read.metharray.sheet(dataDirectory, pattern="Scorza_Project_001_Sample_Sheet.csv")

#Leer los datos IDAT y guardarlos en rgSet
rgSet <- read.metharray.exp(dataDirectory,recursive=TRUE)


## 3.3 Renombrar columnas de rgSet
# 1. Extraer los nombres de las muestras desde targets
sample_names <- targets$Sample_Name

# 2: crear los nombres de las columnas uniendo "Slide" y "Array"
column_names <- paste(targets$Slide, targets$Array, sep="_")

# 3: Renombrar las columnas de rgSet (si es que las columnas de rgSet coinciden con los nombres generados)

if (all(column_names %in% colnames(rgSet))) {
  # Ordenar los nombres de las muestras según el orden de las columnas en rgSet
  sample_names_ordered <- sample_names[match(colnames(rgSet), column_names)]
  # Renombrar las columnas de detP con los nombres de las muestras
  colnames(rgSet) <- sample_names_ordered
} else {
  stop("Los nombres generados no coinciden con rgSet. Revisa 'targets'.")
}

print(colnames(rgSet))

# Eliminar la muestra "NA10858_2"
rgSet <- rgSet[, sampleNames(rgSet) != "NA10858_2"]
sampleNames(rgSet)


#---4.0 Control de Calidad------------------------------------------------------
##4.1 P-valores de detección para cada sitio cg de cada muestra
detP <- detectionP(rgSet) #en detP cada columna es una muestra y cada fila un sitio cg
head(detP)

##4.2 Hacer promedio de p-valores de detección para cada muestra y ordenarlos decrecientemente
mean_detP <- colMeans(detP)
mean_detP_sorted <- sort(mean_detP, decreasing = TRUE)

##4.3 Visualizar en un gráfico los p-valores de detección promedio para cada muestra
barplot(mean_detP_sorted, las=2, cex.names=0.8, ylim=c(0,0.0005) , ylab="Mean detection p-values")

#abline(h=0.0005,col="red")# 

##4.4 Hacer el qcreport
qcReport(rgSet, sampNames=targets$Sample_Name, 
         pdf="qcReport.pdf") #chequear los warnings

#---5.0 Normalización-----------------------------------------------------------
mSetSq <- preprocessQuantile(rgSet)
mSetRaw <- preprocessRaw(rgSet)
#row.names(mSetSq)


# Visualizar los datos antes y después de la normalización
#par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Name,main="Raw", legend=FALSE) 
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Name, main="Normalized", legend=FALSE) 

rm(mSetRaw);gc()  

#text.col=brewer.pal(8,"Dark2")
#---6.Eliminar muestras duplicadas----------------------------------------------
muestras_duplicadas <- c("29_2", "33_2", "17_2")
mSetSq <- mSetSq[, !colnames(mSetSq) %in% muestras_duplicadas]
dim(mSetSq)
datos_fenotipicos <- subset(datos_fenotipicos, !Muestra.Metiloma %in% muestras_duplicadas)
dim(datos_fenotipicos)

#Eliminar duplicados de detP
detP <- detP[, !colnames(detP) %in% muestras_duplicadas]
dim(detP)

#---7. Exploración--------------------------------------------------------------
#7.1 MDS -------
Matriz_met<-getM(mSetSq) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado mSetSq)#
head(Matriz_met)

#definir paleta de colores
pal <- brewer.pal(8,"Dark2")

par(mfrow = c(1, 1))
plotMDS(Matriz_met, top=1000, gene.selection="common")  

# Para examinar otras dimensiones y buscar otras fuentes de variación#
#plotMDS(Matriz_metF, top=1000, gene.selection="common", dim=c(1,3))

##7.1.2 MDS con color por grupo adversidad
# Definir una paleta de colores para los grupos de ACE.score
colores_grupo <- c("A" = "red", "B" = "blue", "C" = "green")

# Crear un vector de colores asociando cada muestra a su grupo adversidad
muestra_colores_grupo <- colores_grupo[as.character(datos_fenotipicos$Grupo.según.tipo.de.adversidad)]
names(muestra_colores_grupo) <- datos_fenotipicos$Muestra.Metiloma

# Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores_grupo'
colores_muestras_usadas_grupo<- muestra_colores_grupo[colnames(Matriz_met)]

# Generar el gráfico MDS con las muestras coloreadas según su grupo ACE.score y añadir una leyenda para identificar los grupos
plotMDS(Matriz_met, top=1000, gene.selection="common", col=colores_muestras_usadas_grupo)
legend("topright", legend=names(colores_grupo), fill=colores_grupo, title="grupo.adversidad")

#7.2 PCA -------
#Análisis de PCA
pca_result <- prcomp(t(Matriz_met), scale. = TRUE)
# Convertir los resultados de PCA en un data frame
pca_df <- as.data.frame(pca_result$x)
#Agregar información de las muestras
pca_df$Sample_Name <- factor(colnames(Matriz_met))

#Crear el gráfico
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample_Name)) +
  geom_text(size = 3) +
  labs(title = "PCA of Methylation Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) 


#---8. Filtrado de datos--------------------------------------------------------
##8.1 Filtrado mantener sólo las sondas que también están presentes en mSetSq
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
#para verificar algunos valores de detP
head(detP) 
tail(detP)
set.seed(6) 
detP[sample(nrow(detP),6),]


##8.2 Filtrar CpGs con p-valores de detección > 0.01 en todas las muestras
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)
mSetSqFlt <- mSetSq[keep,]


##8.3 Eliminar CpGs asociados a SNPs
mSetSqFlt <- dropLociWithSnps(mSetSqFlt) 

##8.4 filtrado de sondas en cromosomas sexuales
keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

##8.5 Filtrar CpGs no variables
# Obtener los valores beta de la matriz de datos normalizados y filtrados
betaVals <- getBeta(mSetSqFlt)

# Calcular el rango interpercentil para cada CpG
interpercentile_range <- apply(betaVals, 1, function(x) {
  quantile(x, 0.90) - quantile(x, 0.10)})


#Visualizar cantidad de sitios por quantil en un histograma
hist(interpercentile_range , breaks = 100 ) 
abline(v=0.05, col="blue")

# Definir el umbral 
umbral <- 0.05
# Identificar los CpGs no variables (aquellos que varían entre las muestras menos menos que el umbral)
non_variable_cpgs <- rownames(betaVals)[which(interpercentile_range < umbral)]
length(non_variable_cpgs) 

# Filtrar para excluir los CpGs no variables
betaVals <- betaVals[!rownames(betaVals) %in% non_variable_cpgs, ]
dim(betaVals)
mSetSqFlt <- mSetSqFlt[!rownames(mSetSqFlt) %in% non_variable_cpgs,]
dim(mSetSqFlt)

#Eliminar objetos que ya no usaremos más
rm(interpercentile_range);gc() 
rm(non_variable_cpgs);gc() 



##8.6 filtrado Funcional (promotores)
table(annEPIC$Regulatory_Feature_Group) #Para ver los niveles que toma Regulatory_f_group

# selección de las sondas que están en promotores
sondas_promotores <- annEPIC$Name[annEPIC$Regulatory_Feature_Group %in% c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")]

# Filtrar Matriz_met_Flt y mSetSqFlt para mantener sólo sondas de promotores 
mSetSqFlt <- mSetSqFlt[rownames(mSetSqFlt) %in% sondas_promotores, ]

#eliminar objeto que ya no vamos a usar
rm(sondas_promotores);gc() 

#7.8 Verificar dimensiones y generar la matriz de valores M de metilación para las muestras filtradas
Matriz_met_Flt<-getM(mSetSqFlt)

dim(mSetSqFlt)
dim(Matriz_met_Flt)

#9. MDS ------------------------------------------------------------------------
Matriz_met_Flt<-getM(mSetSqFlt) 
head(Matriz_met_Flt)

colores_muestras_usadas_grupo<- muestra_colores_grupo[colnames(Matriz_met_Flt)]
colores_muestras_usadas_grupo

# Generar el gráfico MDS con las muestras coloreadas según su grupo ACE.score y añadir una leyenda para identificar los grupos
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common", col=colores_muestras_usadas_grupo)
legend("topright", legend=names(colores_grupo), fill=colores_grupo, title="grupo.adversidad")



#10. guardar resultados en RData------------------------------------------------
save(datos_fenotipicos, Matriz_met_Flt, mSetSqFlt, annEPIC, file = "preprocesadosAdversidad.RData")
