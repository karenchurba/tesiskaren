########################################
#objetivo: obtener DMPs, DMRs, y sus genes asociados
#Input: Preprocesamiento
#Output: 
#Author: karenchurba
#Fecha creación 2024/03/01
#Ultima modificación 2024/07/17
#########################################

#---1.Paquetes -------

install.packages("BiocManager")

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

#9. cargar archivos preprocesados----
setwd("C:/Users/Karen/Documents/Tesis")

load("preprocesados.RData")
pal <- brewer.pal(8,"Dark2")
##Obtener la información de sitios de metilación de genoma humano y guardarla en annEPIC
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

#---10. ANÁLISIS DE METILACIÓN DIFERENCIAL POR SONDA-------

#Eliminar duplicados de datos fenotipicos
datos_fenotipicos_completa <- subset(datos_fenotipicos_completa, !Muestra.Metiloma %in% c("29_2", "33_2", "17_2"))

# Crear el factor de interés (Adversidad)
Adversidad <- factor(datos_fenotipicos_completa$Grupo.según.tipo.de.adversidad)
#levels(Adversidad) <- c("A", "B","C")

# Crear una matriz de diseño
design <- model.matrix(~0 + Adversidad, data=datos_fenotipicos_completa)

colnames(design) <- c("A", "B", "C")
print(design)
dim(Matriz_met_Flt)
dim(design)
# Ajustar el modelo lineal
fit <- lmFit(Matriz_met_Flt, design)


print(fit)
# Crear una matriz de contrastes para la comparación de interés (por ejemplo, MTR 1 vs MTR 0)
contMatrix <- makeContrasts(
  A - B,  # Contraste entre grupo adversidad 3 y 1
  A - C,  # Contraste entre grupo adversidad 6 y 0
  B - C,  # Contraste entre grupo adversidad 4 y 0
  levels = design
)

contMatrix

# Ajustar los contrastes
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

# Obtener la tabla de resultados para el contraste especificado
annEPICSub <- annEPIC[match(rownames(Matriz_met_Flt), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
head(DMPs)
# Graficar los 4 sitios CpG más significativamente diferenciados
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(getBeta(mSetSqFlt), cpg=cpg, pheno=datos_fenotipicos_completa$Grupo.según.tipo.de.adversidad, ylab="Beta values")
})
head(Matriz_met_Flt)
#---11.1 Identificacion de DMRs----
#---11.1AC Adversidad A vs C----
myAnnotation <- cpg.annotate(
  object = Matriz_met_Flt, 
  datatype = "array", 
  what = "M", 
  analysis.type = "differential", 
  design = design, 
  contrasts = TRUE, 
  cont.matrix = contMatrix, 
  coef = "A - C",  
  arraytype = "EPIC")
#fdr=0.7)

str(myAnnotation)
# Identificar DMRs 
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.1)
#DMRs$results
str(DMRs)
length(DMRs@min_smoothed_fdr)
sort(DMRs@min_smoothed_fdr)
# Extraer los resultados
results.ranges <- extractRanges(DMRs)
results.ranges

#guardar los resultados en un excel
write.table(results.ranges, file="DMRs_GrupoAvsC.csv", sep=";", row.names=FALSE)
#---11.1BC Adversidad B vs C----
myAnnotationBC <- cpg.annotate(
  object = Matriz_met_Flt, 
  datatype = "array", 
  what = "M", 
  analysis.type = "differential", 
  design = design, 
  contrasts = TRUE, 
  cont.matrix = contMatrix, 
  coef = "B - C",  
  arraytype = "EPIC")
#fdr=0.7)

str(myAnnotationBC)
# Identificar DMRs 
DMRsBC <- dmrcate(myAnnotationBC, lambda = 1000, C = 2, pcutoff = 0.05)
#DMRs$results
str(DMRsBC)
length(DMRsBC@min_smoothed_fdr)
sort(DMRsBC@min_smoothed_fdr)
# Extraer los resultados
results.rangesBC <- extractRanges(DMRsBC)
results.rangesBC

#guardar los resultados en un excel
write.table(results.rangesBC, file="DMRs_GrupoBvsC.csv", sep=";", row.names=FALSE)


#---11.1BA Adversidad B vs A----
myAnnotationBA <- cpg.annotate(
  object = Matriz_met_Flt, 
  datatype = "array", 
  what = "M", 
  analysis.type = "differential", 
  design = design, 
  contrasts = TRUE, 
  cont.matrix = contMatrix, 
  coef = "A - B",  
  arraytype = "EPIC")
#fdr=0.7)

str(myAnnotationBA)
# Identificar DMRs 
DMRsBA <- dmrcate(myAnnotationBA, lambda = 1000, C = 2, pcutoff = 0.05)
#DMRs$results
str(DMRsBA)
length(DMRsBA@min_smoothed_fdr)
sort(DMRsBA@min_smoothed_fdr)
# Extraer los resultados
results.rangesBA <- extractRanges(DMRsBA)
results.rangesBA

#guardar los resultados en un excel
write.table(results.rangesBA, file="DMRs_GrupoAvsB.csv", sep=";", row.names=FALSE)
