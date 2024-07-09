########################################
#objetivo:
#Input: Tesis
#Output: 
#Author: karenchurba
#Fecha creación 2024/03/01
#Ultima modificación 2024/07/9
#########################################

#---1.Paquetes-------

install.packages("BiocManager", quietly =TRUE)

BiocManager::install ("minfi")

install.packages("minfi")
install.packages("RCurl")
install.packages("GenomicRanges")
#install.packages("GenomeInfoDb")#

library(minfi)
force = TRUE

install.packages("IlluminaHumanMethylationEPICmanifest")
library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) 
library(IlluminaHumanMethylationEPICmanifest)
library(RColorBrewer)
#library(missMethyl) #Este no lo pude descargar#
#library(minfiData) #Este no lo pude descargar#
#library(Gviz) #Este no lo pude descargar#
#library(DMRcate)Este no lo pude descargar# Este lo vamos a necesitar
#library(stringr)#


#---2. Import Data-------


dataDirectory <- "C:/Users/Karen/Documents/Tesis"
#list.files(dataDirectory, recursive = TRUE)

#Eliminar objeto del environment
#rm(mSetRaw);gc()

##Obtener la información de sitios de metilación de genoma humano y guardarla en ann450k
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19) #chequear si necesitamos cambiar este paquete por uno que se corresponda con las sondas que usaron#

#---2.1 Leer el sample sheet y guardar la info en Targets-----

targets <- read.metharray.sheet(dataDirectory, pattern="Scorza_Project_001_Sample_Sheet.csv")

#sum(targets$Sample_Plate=="Tube") para contar la cantidad de filas que dicen "tube"#

#---2.2 Leer los datos de los IDAT y guardarlos en rgSet-------

rgSet <- read.metharray.exp(dataDirectory,recursive=TRUE)


#---3.0 QUALITY CONTROL-------

#Calcular los p-valores de detección para cada sitio cg de cada muestra - en detP cada columna es una muestra y cada fila un sitio cg
detP <- detectionP(rgSet)
#head(detP)
#ncol(detP)
#nrow(detP)
#dim(detP)#

#Poner los nombres de las muestras a las columnas de detP 
targets$Slide_Array <- apply(targets[,c("Slide","Array")],MARGIN = 1,FUN = paste,collapse="_")
colnames(detP)<- targets$Sample_Name
#head(detP)#

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


#---5. EXPLORACIÓN-------
#---5.1 MDS-------
sample_names <- targets$Sample_Name

Matriz_met<-getM(mSetSq) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado mSetSq)#
head(Matriz_met)

colnames(Matriz_met) <- sample_names #¿Cómo sé que se asignó correctamente a cada muestra el nomrbe que le corresponde?#
#head(Matriz_met) 


Matriz_metF<- Matriz_met[,colnames(Matriz_met)!="NA10858_2"] #Eliminar el NA#

#pal <- brewer.pal(8,"Dark2") define la paleta de colores#

plotMDS(Matriz_metF, top=1000, gene.selection="common")  

# Para examinar otras dimensiones y buscar otras fuentes de variación#
plotMDS(Matriz_metF, top=1000, gene.selection="common", dim=c(1,3))

#---5.2 PCA-------

install.packages("ggplot2")
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


####Para ver los influence scores de PC1

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
#hacer nuevamente MDS y PCA después
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
head(detP)

keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]
#mSetSqFlt

mSetSqFlt <- dropLociWithSnps(mSetSqFlt) #elimina las sondas que pueden tener SNPs comunes que afectan el CpG

###Prueba:filtrado de dobles reactivos usando la base de datos que saqué de: https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "1031_CpG_sites_removed_from_MethylationEPIC_15073387_v1-0_bpm.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes$TargetID) #ver que extrae featurenames
table(keep) #no me cierra porque no obtuve ningún FALSE (osea no va a filtrar nada)#
mSetSqFlt <- mSetSqFlt[keep,] 
mSetSqFlt

#---7. EXPLORACION POST-FILTRADO-------
Matriz_met_Flt<-getM(mSetSqFlt) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado Y FILTRADO mSetSqFlt)#
colnames(Matriz_met_Flt) <- sample_names #¿Cómo sé que se asignó correctamente a cada muestra el nomrbe que le corresponde?#
Matriz_met_Flt<- Matriz_met_Flt[,colnames(Matriz_met_Flt)!="NA10858_2"] #Eliminar el NA#
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common")  

rm(Matriz_met);gc()

#Agregar color por grupo al MDS
####Carga de los datos fenotípicos#
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


rm(filas_duplicadas);gc() 
rm(filas_duplicadas_nuevas);gc() 
rm(merged_data);gc() 

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

#Hacer de nuevo el PCA
pca_result <- prcomp(t(Matriz_met_Flt), scale. = TRUE)
pca_df <- as.data.frame(pca_result$x)
#Agregar información de las muestras
pca_df$Sample_Name <- factor(colnames(Matriz_met_Flt))
ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample_Name)) +
  geom_text(size = 3) +
  labs(title = "PCA of Methylation Data",
       x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 2), "%)"),
       y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 2), "%)")) 


#---8. ANÁLISIS DE METILACIÓN DIFERENCIAL POR SONDA-------


# Crear el factor de interés (MTR)
MTR <- factor(datos_fenotipicos$MTR)
levels(MTR) <- c("Nivel0", "Nivel1")
# Crear una matriz de diseño
design <- model.matrix(~0 + MTR)
colnames(design) <- c("Nivel0", "Nivel1")
print(design)
#print(dim(design))

# Ajustar el modelo lineal
fit <- lmFit(Matriz_met_Flt, design)
colnames(Matriz_met_Flt)
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
ann450kSub <- ann450k[match(rownames(Matriz_met_Flt), ann450k$Name), c(1:4, 12:19, 24:ncol(ann450k))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=ann450kSub)
head(DMPs)

# Guardar los resultados en un archivo CSV
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# Graficar los 4 sitios CpG más significativamente diferenciados
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(getBeta(mSetSqFlt), cpg=cpg, pheno=datos_fenotipicos$MTR, ylab="Beta values")
})
View(datos_fenotipicos[,c("Sample_Name","MTR")])


AAAAA

