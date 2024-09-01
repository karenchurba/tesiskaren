########################################
#objetivo:
#Input: Tesis
#Output: 
#Author: karenchurba
#Fecha creación 2024/03/01
#Ultima modificación 2024/07/17
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

library(minfi)
library(knitr)
library(limma)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) 
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
#library(missMethyl) #Este no lo pude descargar#
#library(minfiData) #Este no lo pude descargar#
#library(Gviz) #Este no lo pude descargar# este quizás lo usemos
library(DMRcate)
#library(stringr)#

if (!("dif (!("dif (!("devtools" %in% installed.packages()[, "Package"])) {
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


#2.1 Leer el sample sheet y guardar la info en Targets-----

targets <- read.metharray.sheet(dataDirectory, pattern="Scorza_Project_001_Sample_Sheet.csv")

#sum(targets$Sample_Plate=="Tube") para contar la cantidad de filas que dicen "tube"#

#2.2 Lectura de datos de los IDAT y guardado en rgSet-------

rgSet <- read.metharray.exp(dataDirectory,recursive=TRUE)

#2.3 P-valores de detección para cada sitio cg de cada muestra -------
#en detP cada columna es una muestra y cada fila un sitio cg
detP <- detectionP(rgSet)
#head(detP)
#ncol(detP)
#nrow(detP)
#dim(detP)#

#2.3.1 Nombres de las muestras a las columnas de detP---- 
# Paso 1: Extraer los nombres de las muestras desde targets
sample_names <- targets$Sample_Name

# Paso 2: Crear los nombres de las columnas combinando "Slide" y "Array" con un guion bajo
column_names <- paste(targets$Slide, targets$Array, sep="_")

# Paso 3: Renombrar las columnas de detP
# Verificar que las columnas de detP coincidan con los nombres generados
if (all(column_names %in% colnames(detP))) {
  # Ordenar los nombres de las muestras según el orden de las columnas en detP
  sample_names_ordered <- sample_names[match(colnames(detP), column_names)]
  # Renombrar las columnas de detP con los nombres de las muestras
  colnames(detP) <- sample_names_ordered
} else {
  stop("Los nombres de las columnas generados no coinciden con los nombres de las columnas en detP")
}

print(colnames(detP))

# Eliminar el NA
index_to_remove <- which(sample_names == "NA10858_2")
detP <- detP[ , -index_to_remove]


#3.0 QUALITY CONTROL-------

# Hacer promedio de p-valores de detección para cada muestra y ordenarlos decrecientemente
mean_detP <- colMeans(detP)
mean_detP_sorted <- sort(mean_detP, decreasing = TRUE)

# Visualizar en un gráfico los p-valores de detección promedios para cada muestra
barplot(mean_detP_sorted, las=2, 
        cex.names=0.8, ylim=c(0,0.0005) , ylab="Mean detection p-values")

#abline(h=0.0005,col="red")# 

#Hacer el qcreport
qcReport(rgSet, sampNames=targets$Sample_Name, 
         pdf="qcReport3.pdf") #chequear los warnings


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

#--5.1.1 Nombres de las muestras a mSetSq----
original_column_names <- colnames(mSetSq)
# Verificar que las columnas de mSetSq coincidan con los nombres generados
if (all(column_names %in% colnames(mSetSq))) {
  # Ordenar los nombres de las muestras según el orden de las columnas en mSetSq
  sample_names_ordered <- sample_names[match(colnames(mSetSq), column_names)]
  # Renombrar las columnas de mSetSq con los nombres de las muestras
  colnames(mSetSq) <- sample_names_ordered
} else {
  stop("Los nombres de las columnas generados no coinciden con los nombres de las columnas en mSetSq")
}

# Comprobar los nombres de las columnas de mSetSq
print(colnames(mSetSq))

#tabla de comparación de los nombres antes y después
#comparison_table <- data.frame(Archivo = original_column_names, Muestra = colnames(mSetSq))

# Verifica los nombres de las muestras en el objeto mSetSq
#sample_names_mset <- sampleNames(mSetSq)
#print(sample_names_mset)

# ELIMINAR EL NA
#Encuentra el índice de la muestra que deseas eliminar
index_to_remove <- which(sample_names_mset == "NA10858_2")

# Elimina la muestra del objeto mSetSq
mSetSq <- mSetSq[ , -index_to_remove]


#---5. EXPLORACIÓN-------
#5.1 MDS-------
Matriz_met<-getM(mSetSq) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado mSetSq)#
head(Matriz_met)

#pal <- brewer.pal(8,"Dark2") define la paleta de colores#

plotMDS(Matriz_met, top=1000, gene.selection="common")  

# Para examinar otras dimensiones y buscar otras fuentes de variación#
#plotMDS(Matriz_metF, top=1000, gene.selection="common", dim=c(1,3))

#5.2 PCA-------

library(ggplot2)

#1:Extraer los valores de beta para PCA - NO HACE FALTA. podemos usar la matriz que creamos antes con los M values
#beta_values <- getBeta(mSetSq)
#2:Análisis de PCA
pca_result <- prcomp(t(Matriz_metF), scale. = TRUE)

# Convertir los resultados de PCA en un data frame
pca_df <- as.data.frame(pca_result$x)
#Agregar información de las muestras
pca_df$Sample_Name <- factor(colnames(Matriz_metF))


#Crear el gráfico
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample_Name)) +
  geom_text(size = 3) +
  labs(title = "PCA of Methylation Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) 


#Para ver los influence scores de PC1---

#Extraer los scores del PC1
#pc1_scores <- pca_result$x[, "PC1"]
#print(pc1_scores)# Mostrar los scores del PC1

###Para visualizar en un gráfico los scores del PC1

#Crear un data frame
#pc1_df <- data.frame(Sample = rownames(pca_result$x), PC1 = pc1_scores, Group = targets$Sample_Name)

# Ordenar por los scores de PC2
#pc1_df <- pc1_df[order(pc1_df$PC1), ]

# Gráfico de barras de los scores de PC2
#ggplot(pc1_df, aes(x = Group, y = PC1, fill = Group)) +
 # geom_bar(stat = "identity") + ylim(-250,250) +  labs(title = "Scores of PC1",
  #     x = "Group",
   #    y = "PC1 Score")



####Para ver los influence scores de PC2

#Extraer los scores del PC2
# pc2_scores <- pca_result$x[, "PC2"]
#print(pc2_scores)# Mostrar los scores del PC2

### Visualizar en un gráfico los scores del PC2

#Crear un data frame
#pc2_df <- data.frame(Sample = rownames(pca_result$x), PC2 = pc2_scores, Group = targets$Sample_Name)

# Ordenar por los scores de PC2
#pc2_df <- pc2_df[order(pc2_df$PC2), ]

# Gráfico de barras de los scores de PC2
#ggplot(pc2_df, aes(x = Group, y = PC2, fill = Group)) +
 # geom_bar(stat = "identity") +
  #labs(title = "Scores of PC2",
   #    x = "Group",
    #   y = "PC2 Score")





#---6. FILTRADO-------
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] #Mantiene en detP sólo la información de las sondas que están también en mSetSq (que es el archivo que estaba normalizado)

head(detP)
ncol(detP)
detP<- detP[,colnames(detP)!="NA10858_2"] #Eliminar el NA de detP#

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

#6.2 Prueba: filtrado por CpGs no-variables----
# Obtener los valores beta de la matriz de datos normalizados y filtrados
betaVals <- getBeta(mSetSqFlt)

# Calcular el rango interpercentil (percentil 90 - percentil 10) para cada CpG
interpercentile_range <- apply(betaVals, 1, function(x) {
  quantile(x, 0.90) - quantile(x, 0.10)
})


dim()
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
rm(non_variable_cpgs , filtered_betaVals);gc() 



#6.2.1 opcional: importar cpgs no variables de estudio previo---- 
#OJO: este estudio fue hecho cpara 450k
#library(RCurl)
#x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Placenta_CpGs.csv")
#y <- read.csv(text = x)
#nonvariable_placenta<-y$CpG

# Convertir ambas listas a conjuntos
#set_non_variable_cpgs <- as.character(non_variable_cpgs)
#set_nonvariable_placenta <- as.character(nonvariable_placenta)

# Encontrar los CpGs no variables en común
#common_non_variable_cpgs <- intersect(set_non_variable_cpgs, set_nonvariable_placenta) 

# Ver cuántos CpGs no variables hay en común
#length(common_non_variable_cpgs)
#head(common_non_variable_cpgs)

# Filtrar la matriz de datos excluyendo los CpGs no variables en común
#Matriz_met_Flt <- Matriz_met_Flt[!rownames(Matriz_met_Flt) %in% common_non_variable_cpgs, ]

#---9. Identificación de fracciones de tipos celulares-------
Matriz_met_Flt<-getM(mSetSqFlt)
data(centEpiFibIC.m)

# Estimación de las fracciones de tipos celulares
resultado <- epidish(beta.m = Matriz_met_Flt, ref.m = centEpiFibIC.m, method = "RPC")

# Fracciones estimadas de tipos celulares
print(resultado$estF)

# Dimensiones de la matriz de referencia usada
print(dim(resultado$ref))

# Dimensiones de la matriz de datos usada para la estimación
print(dim(resultado$dataREF))

data(centBloodSub.m)
frac.m <- hepidish(beta.m = Matriz_met_Flt, ref1.m = centEpiFibIC.m, ref2.m = centBloodSub.m, h.CT.idx = 3, method = 'RPC')
summary(frac.m)
apply(frac.m , 2 , sd)

#6.3 filtrado Funcional: nos quedamos con sondas de promotores---- 
#table(annEPIC$Regulatory_Feature_Group) #Para ver los niveles que toma Regulatory_f_group
# selección de las sondas que están en promotores
sondas_promotores <- annEPIC$Name[annEPIC$Regulatory_Feature_Group %in% c("Promoter_Associated", "Promoter_Associated_Cell_type_specific")]

# Filtrar la matriz de metilación (Matriz_met_Flt) para quedarte solo con las sondas en promotores
Matriz_met_Flt <- Matriz_met_Flt[rownames(Matriz_met_Flt) %in% sondas_promotores, ]

# Verificar las dimensiones después del filtrado
dim(Matriz_met_Flt)

rm(sondas_promotores);gc() 

#---7. EXPLORACION POST-FILTRADO-------
Matriz_met_Flt<-getM(mSetSqFlt) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado Y FILTRADO mSetSqFlt)#

plotMDS(Matriz_met_Flt, top=1000, gene.selection="common")  

rm(Matriz_met);gc()


#7.1 Carga de los datos fenotípicos----
datos_fenotipicos <- read.csv("C:/Users/Karen/Documents/Tesis/datos_fenotipicos.csv")

# Identificar las filas que tienen duplicados
filas_duplicadas <- datos_fenotipicos[!is.na(datos_fenotipicos$ID.Illumina.de.duplicados...Barcode) &
                                        !is.na(datos_fenotipicos$ID.Illumina.de.duplicados...Sentrix), ]

# Crear las filas duplicadas con las columnas actualizadas
filas_duplicadas_nuevas <- filas_duplicadas
filas_duplicadas_nuevas$ID.Illumina...Barcode <- filas_duplicadas$ID.Illumina.de.duplicados...Barcode
filas_duplicadas_nuevas$ID.Illumina...Sentrix <- filas_duplicadas$ID.Illumina.de.duplicados...Sentrix


# Agregar las filas duplicadas a la tabla original
datos_fenotipicos <- rbind(datos_fenotipicos, filas_duplicadas_nuevas)

###Agregado de la columna Sample_Names a datos_fenotipicos
# Renombrar columnas para facilitar el merge
names(datos_fenotipicos)[names(datos_fenotipicos) == "ID.Illumina...Barcode"] <- "Array"
names(datos_fenotipicos)[names(datos_fenotipicos) == "ID.Illumina...Sentrix"] <- "Slide"

# Hacer el merge para agregar el Sample_Name
merged_data <- merge(datos_fenotipicos, targets[, c("Array", "Slide", "Sample_Name")], by = c("Array", "Slide"))
datos_fenotipicos <- merged_data


rm(filas_duplicadas,filas_duplicadas_nuevas,merged_data);gc() 

#7.2 MDS: Agregar color por grupo MTR-------


#definir paleta de colores#
colores <- c("1" = "red", "0" = "blue")                          

#Crear un vector de colores asociando cada muestra a su grupo MTR#
muestra_colores <- colores[as.character(datos_fenotipicos$MTR)]
names(muestra_colores) <- datos_fenotipicos$Sample_Name
#Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores'#
colores_muestras_usadas <- muestra_colores[colnames(Matriz_met_Flt)]
#Genera el gráfico MDS con las muestras coloreadas según su grupo y añade una leyenda para identificar los grupos#
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common", col=colores_muestras_usadas)
legend("topright", legend=c("1", "0"), fill=c("red", "blue"), title="MTR")


#####versión 2 del gráfico: con la leyenda fuera del gráfico#
# Configurar el layout para tener un espacio superior para la leyenda
layout(matrix(c(1, 2), nrow = 2), heights = c(1, 4))
# Espacio para la leyenda
par(mar = c(0, 4, 0, 4)) # Margenes (bottom, left, top, right)
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "") # Gráfico vacío
legend("center", legend = c("Grupo 1", "Grupo 0"), fill = c("red", "blue"), title = "Grupos")
# Gráfico MDS
par(mar = c(5, 4, 4, 2) + 0.1) # Resetear margenes para el gráfico MDS
plotMDS(Matriz_met_Flt, top = 1000, gene.selection = "common", col = colores_muestras_usadas)

rm(colores_muestras_usadas , muestra_colores);gc()
#7.3 MDS: Agregar color por grupo ARRAY-------
colores_array <- c("R01C01" = "red", "R02C01" = "blue", "R03C01" = "green", 
                   "R04C01" = "purple", "R05C01" = "orange", "R06C01" = "brown", 
                   "R07C01" = "pink", "R08C01" = "cyan")

# Crear un vector de colores asociando cada muestra a su grupo Array #
muestra_colores_array <- colores_array[as.character(datos_fenotipicos$Array)]
names(muestra_colores_array) <- datos_fenotipicos$Sample_Name

# Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores_array' #
colores_muestras_usadas_array <- muestra_colores_array[colnames(Matriz_met_Flt)]

plotMDS(Matriz_met_Flt, top=1000, gene.selection="common", col=colores_muestras_usadas_array)

legend("topright", inset=c(-0.3, 0), legend=names(colores_array), fill=colores_array, title="Array", xpd=TRUE)

rm(colores_array , muestra_colores_array , colores_muestras_usadas_array);gc()
#7.3 MDS: Agregar color por grupo slide----
colores_slide <- c("206960650083" = "red", 
                   "206960650086" = "blue", 
                   "206960650112" = "green", 
                   "206960650116" = "purple", 
                   "206960650123" = "orange")

# Crear un vector de colores asociando cada muestra a su grupo Slide
muestra_colores_slide <- colores_slide[as.character(datos_fenotipicos$Slide)]
names(muestra_colores_slide) <- datos_fenotipicos$Sample_Name

# Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores_slide'
colores_muestras_usadas_slide <- muestra_colores_slide[colnames(Matriz_met_Flt)]

# Generar el gráfico MDS con las muestras coloreadas según su grupo de Slide y añadir una leyenda para identificar los grupos
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common", col=colores_muestras_usadas_slide)
legend("topright", legend=names(colores_slide), fill=colores_slide, title="Slide")
#7.3 MDS: Agregar color por ACE score----

# Definir una paleta de colores para los grupos de ACE.score
colores_ace <- c("0" = "red", "1" = "blue", "2" = "green", "3" = "purple", "4" = "orange", "5" = "pink", "6" = "brown")

# Crear un vector de colores asociando cada muestra a su grupo ACE.score
muestra_colores_ace <- colores_ace[as.character(datos_fenotipicos$ACE.score)]
names(muestra_colores_ace) <- datos_fenotipicos$Sample_Name

# Asegurarse de que los nombres de las muestras en 'Matriz_met_Flt' coincidan con 'muestra_colores_ace'
colores_muestras_usadas_ace <- muestra_colores_ace[colnames(Matriz_met_Flt)]

# Generar el gráfico MDS con las muestras coloreadas según su grupo ACE.score y añadir una leyenda para identificar los grupos
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common", col=colores_muestras_usadas_ace)
legend("topright", legend=names(colores_ace), fill=colores_ace, title="ACE.score")

rm(colores_ace,muestra_colores_ace,colores_muestras_usadas_ace);gc() 

#7.4 PCA post filtrado-------
pca_result <- prcomp(t(Matriz_met_Flt), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
#Agregar información de las muestras
pca_df$Sample_Name <- factor(colnames(Matriz_met_Flt))
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample_Name)) +
  geom_text(size = 3) +
  labs(title = "PCA of Methylation Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) 

rm(pca_result,pca_df);gc() 



#---10. ANÁLISIS DE METILACIÓN DIFERENCIAL POR SONDA-------


# Crear el factor de interés (MTR)
MTR <- factor(datos_fenotipicos$MTR)
levels(MTR) <- c("Nivel0", "Nivel1")
# Crear una matriz de diseño
design <- model.matrix(~0 + MTR)
colnames(design) <- c("Nivel0", "Nivel1")
print(design)
#print(dim(design))
dim(Matriz_met_Flt)
dim(design)
# Ajustar el modelo lineal
fit <- lmFit(Matriz_met_Flt, design)


print(fit)
# Crear una matriz de contrastes para la comparación de interés (por ejemplo, MTR 1 vs MTR 0)
contMatrix <- makeContrasts(Nivel1 - Nivel0, levels=design)
contMatrix

# Ajustar los contrastes
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)


# Resumen de los resultados
summary(decideTests(fit2))


# Obtener la tabla de resultados para el contraste especificado
annEPICSub <- annEPIC[match(rownames(Matriz_met_Flt), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
head(DMPs)


# Guardar los resultados en un archivo CSV
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# Graficar los 4 sitios CpG más significativamente diferenciados
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(getBeta(mSetSqFlt), cpg=cpg, pheno=datos_fenotipicos$MTR, ylab="Beta values")
})
View(datos_fenotipicos[,c("Sample_Name","MTR")])








 