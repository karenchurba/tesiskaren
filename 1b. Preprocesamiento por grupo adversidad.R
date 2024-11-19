########################################
#objetivo: 
#Input: Tesis
#Output: preprocesados_2.RData (datos_fenotípicos + Matriz_Met_Flt + mSetSqFlt)
#Author: karenchurba
#Fecha creación 2024/08/11
#Ultima modificación 2024/14/11
#########################################

#---1.Paquetes-------

#install.packages("BiocManager", quietly =TRUE)
#BiocManager::install ("minfi")
#BiocManager::install ("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
#BiocManager::install ("IlluminaHumanMethylationEPICmanifest")
#install.packages("RCurl")
#install.packages("GenomicRanges")
#install.packages("GenomeInfoDb")#
#BiocManager::install("DMRcate", update = TRUE)
#if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager")

#BiocManager::install("Gviz")


library(minfi)
library(knitr)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
#library(missMethyl) #Este no lo pude descargar#
#library(minfiData) #Este no lo pude descargar#
library(Gviz) 
library(DMRcate)
#library(stringr)#
library(bumphunter)
library(GenomicRanges)


if  (!("devtools" %in% installed.packages()[, "Package"])) {
  install.packages("devtools")
}
devtools::install_github("markgene/maxprobes")


library(maxprobes) #No me funcionó

library(EpiDISH)

#-2. Import Data-------
setwd("C:/Users/Karen/Documents/Tesis")

dataDirectory <- "C:/Users/Karen/Documents/Tesis"
#list.files(dataDirectory, recursive = TRUE)

#Eliminar objeto del environment
#rm(mSetRaw);gc()

##Obtener la información de sitios de metilación de genoma humano y guardarla en annEPIC
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
#2.1 Carga de los datos fenotípicos----
datos_fenotipicos_completa <- read.csv("C:/Users/Karen/Documents/Tesis/Datos fenotípicos completa.csv")


#2.1 Leer el sample sheet y guardar la info en Targets-----

targets <- read.metharray.sheet(dataDirectory, pattern="Scorza_Project_001_Sample_Sheet.csv")

#sum(targets$Sample_Plate=="Tube") para contar la cantidad de filas que dicen "tube"#

#2.2 Lectura de datos de los IDAT y guardado en rgSet-------

rgSet <- read.metharray.exp(dataDirectory,recursive=TRUE)

#2.2.1 Nombres de las muestras a los datos de rgSet---- 
# 1. Extraigo los nombres de las muestras desde targets
sample_names <- targets$Sample_Name

# 2: creo los nombres de las columnas uniendo "Slide" y "Array" con un guion bajo
column_names <- paste(targets$Slide, targets$Array, sep="_")

# 3: Renombro las columnas de rgSet (si es que las columnas de rgSet coinciden con los nombres generados)

if (all(column_names %in% colnames(rgSet))) {
  # Ordenar los nombres de las muestras según el orden de las columnas en rgSet
  sample_names_ordered <- sample_names[match(colnames(rgSet), column_names)]
  # Renombrar las columnas de detP con los nombres de las muestras
  colnames(rgSet) <- sample_names_ordered
} else {
  stop("Los nombres de las columnas generados no coinciden con los nombres de las columnas en rgSet")
}

print(colnames(rgSet))

# Eliminar el "NA10858_2" del rgSet
rgSet <- rgSet[, sampleNames(rgSet) != "NA10858_2"]
sampleNames(rgSet)

#2.3 P-valores de detección para cada sitio cg de cada muestra -------
#en detP cada columna es una muestra y cada fila un sitio cg
detP <- detectionP(rgSet)
head(detP)
#ncol(detP)
#nrow(detP)
#dim(detP)#



#---4.0 NORMALIZACION-------############ 

mSetSq <- preprocessQuantile(rgSet)
mSetRaw <- preprocessRaw(rgSet)
#row.names(mSetSq)


# Visualizar cómo se ven los datos antes y después de la normalización
#par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Name,main="Raw", legend=FALSE) 
text.col=brewer.pal(8,"Dark2")
densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Name,
            main="Normalized", legend=FALSE) 
rm(mSetRaw);gc()     

#---5. EXPLORACIÓN-------
#5.1 MDS-------
Matriz_met<-getM(mSetSq) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado mSetSq)#
#Eliminar muestras duplicadas
Matriz_met <- Matriz_met[, !colnames(Matriz_met) %in% c("29_2", "33_2", "17_2")]
head(Matriz_met)

#definir paleta de colores
pal <- brewer.pal(8,"Dark2")

par(mfrow = c(1, 1))
plotMDS(Matriz_met, top=1000, gene.selection="common")  

# Para examinar otras dimensiones y buscar otras fuentes de variación#
#plotMDS(Matriz_metF, top=1000, gene.selection="common", dim=c(1,3))
#5.2 MDS: Agregar color por grupo adversidad-------
# Definir una paleta de colores para los grupos de ACE.score
colores_grupo <- c("A" = "red", "B" = "blue", "C" = "green")

# Crear un vector de colores asociando cada muestra a su grupo adversidad
muestra_colores_grupo <- colores_grupo[as.character(datos_fenotipicos_completa$Grupo.según.tipo.de.adversidad)]
names(muestra_colores_grupo) <- datos_fenotipicos_completa$Muestra.Metiloma
muestra_colores_grupo

# Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores_grupo'
colores_muestras_usadas_grupo<- muestra_colores_grupo[colnames(Matriz_met)]
colores_muestras_usadas_grupo
# Generar el gráfico MDS con las muestras coloreadas según su grupo ACE.score y añadir una leyenda para identificar los grupos
plotMDS(Matriz_met, top=1000, gene.selection="common", col=colores_muestras_usadas_grupo)
legend("topright", legend=names(colores_grupo), fill=colores_grupo, title="grupo.adversidad")

#---6. FILTRADO-------
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] #Mantiene en detP sólo la información de las sondas que están también en mSetSq (que es el archivo que estaba normalizado)

featureNames(mSetSq)
rownames(detP)
head(detP)
nrow(detP)


keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]


mSetSqFlt <- dropLociWithSnps(mSetSqFlt) #elimina las sondas que pueden tener SNPs comunes que afectan el CpG

###Prueba:filtrado de dobles reactivos usando la base de datos que saqué de: https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "1031_CpG_sites_removed_from_MethylationEPIC_15073387_v1-0_bpm.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID) #ver que extrae featurenames
table(keep) #no me cierra porque no obtuve ningún FALSE (osea no va a filtrar nada)#
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt
head(xReactiveProbes)
head(featureNames(mSetSqFlt))

rm(xReactiveProbes);gc() 
#6.1 filtrado de sondas en cromosomas sexuales----
# if your data includes males and females, remove probes on the sex chromosomes
keep <- !(featureNames(mSetSqFlt) %in% annEPIC$Name[annEPIC$chr %in% 
                                                      c("chrX","chrY")])
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]

#6.2 Filtrado por CpGs no-variables----
# Obtener los valores beta de la matriz de datos normalizados y filtrados
betaVals <- getBeta(mSetSqFlt)
dim(betaVals)
# Calcular el rango interpercentil (percentil 90 - percentil 10) para cada CpG
interpercentile_range <- apply(betaVals, 1, function(x) {
  quantile(x, 0.90) - quantile(x, 0.10)
})


#Visualizar cantidad de sitios por quantil en histograma
hist(interpercentile_range , breaks = 100 ) 
abline(v=0.05, col="blue")

# Definir el umbral del 5%
threshold <- 0.05

# Identificar los CpGs no variables
non_variable_cpgs <- rownames(betaVals)[which(interpercentile_range < threshold)]
length(non_variable_cpgs) # Ver cuántos CpGs son no variables

# Crear una matriz de datos filtrada que excluye los CpGs no variables
filtered_betaVals <- betaVals[!rownames(betaVals) %in% non_variable_cpgs, ]


# Verificar las dimensiones de la matriz filtrada
dim(filtered_betaVals)

# Filtrar el objeto mSetSqFlt para excluir los CpGs no variables
mSetSqFlt <- mSetSqFlt[!rownames(mSetSqFlt) %in% non_variable_cpgs,]

# Verificar las dimensiones del objeto filtrado
dim(mSetSqFlt)
Matriz_met_Flt<-getM(mSetSqFlt)
dim(Matriz_met_Flt)

rm(interpercentile_range);gc() 
rm(non_variable_cpgs , betaVals);gc() 



#6.3 filtrado Funcional: nos quedamos con sondas de promotores---- 
table(annEPIC$Regulatory_Feature_Group) #Para ver los niveles que toma Regulatory_f_group
# selección de las sondas que están en promotores
sondas_promotores <- annEPIC$Name[annEPIC$Regulatory_Feature_Group %in% c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")]

# Filtrar la matriz de metilación (Matriz_met_Flt) y mSetSqFlt para quedarte solo con las sondas en promotores
Matriz_met_Flt <- Matriz_met_Flt[rownames(Matriz_met_Flt) %in% sondas_promotores, ]
mSetSqFlt <- mSetSqFlt[rownames(mSetSqFlt) %in% sondas_promotores, ]
# Verificar las dimensiones después del filtrado
dim(Matriz_met_Flt)
dim(mSetSqFlt)
rm(sondas_promotores);gc() 

#6.4 Eliminar duplicados----
mSetSqFlt <- mSetSqFlt[, !colnames(mSetSqFlt) %in% c("29_2", "33_2", "17_2")]
dim(mSetSqFlt)
#Eliminar duplicados de datos fenotipicos
datos_fenotipicos_completa <- subset(datos_fenotipicos_completa, !Muestra.Metiloma %in% c("29_2", "33_2", "17_2"))


#7 MDS: Agregar color por grupo adversidad-------
Matriz_met_Flt<-getM(mSetSqFlt) 
#Eliminar muestras duplicadas
Matriz_met_Flt <- Matriz_met_Flt[, !colnames(Matriz_met_Flt) %in% c("29_2", "33_2", "17_2")]
head(Matriz_met_Flt)
# Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores_grupo'
colores_muestras_usadas_grupo<- muestra_colores_grupo[colnames(Matriz_met_Flt)]
colores_muestras_usadas_grupo
# Generar el gráfico MDS con las muestras coloreadas según su grupo ACE.score y añadir una leyenda para identificar los grupos
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common", col=colores_muestras_usadas_grupo)
legend("topright", legend=names(colores_grupo), fill=colores_grupo, title="grupo.adversidad")



#8. guardar archivos en RData----
save(datos_fenotipicos_completa, Matriz_met_Flt, mSetSqFlt, file = "preprocesadosAdversidad.RData")
