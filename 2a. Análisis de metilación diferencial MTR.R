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

load("datos_preprocesadosMTR.RData")
pal <- brewer.pal(8,"Dark2")
##Obtener la información de sitios de metilación de genoma humano y guardarla en annEPIC
annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

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
str(DMPs)


# Guardar los resultados en un archivo CSV
write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# Graficar los 4 sitios CpG más significativamente diferenciados
par(mfrow=c(2,2))
sapply(rownames(DMPs)[1:4], function(cpg){
  plotCpg(getBeta(mSetSqFlt), cpg=cpg, pheno=datos_fenotipicos$MTR, ylab="Beta values")
})

#View(datos_fenotipicos[,c("Sample_Name","MTR")])

#---10. POR ACE SCORE: ANÁLISIS DE METILACIÓN DIFERENCIAL POR SONDA-------
# Crear el factor para ACE score
ACE_score <- factor(datos_fenotipicos$ACE.score, labels = paste0("ACE_", levels(factor(datos_fenotipicos$ACE.score))))
design <- model.matrix(~0 + ACE_score)
colnames(design) <- levels(ACE_score)
#head(design)

# Ajustar el modelo lineal para cada sonda usando los valores en 'Matriz_met_Flt'
fit <- lmFit(Matriz_met_Flt, design)
fit <- eBayes(fit)

# Crear la matriz de contraste 
contMatrix <- makeContrasts(
  #ACE_3 - ACE_1,  # Contraste entre ACE score 3 y 1
  #ACE_6 - ACE_0,  # Contraste entre ACE score 6 y 0
  #ACE_3 - ACE_0,  # Contraste entre ACE score 4 y 0
  ACE_5 - ACE_0,  # Contraste entre ACE score 5 y 0
  levels = design
)
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))

# Obtener la tabla de resultados para el contraste especificado
annEPICSub <- annEPIC[match(rownames(Matriz_met_Flt), annEPIC$Name), c(1:4, 12:19, 24:ncol(annEPIC))]
DMPs <- topTable(fit2, num=Inf, coef=1, genelist=annEPICSub)
str(DMPs)



# Identificar los DMPs a un nivel de significancia dado (por ejemplo, p < 0.01)
#DMPs <- topTable(fit1, coef = 1, adjust = "fdr", number = Inf)
head(DMPs)

# Guardar los resultados en un archivo CSV
write.table(DMPs, file="DMPs2.csv", sep=",", row.names=FALSE)

# Filtrar DMPs con p-valor ajustado menor a 0.1
significant_DMPs <- DMPs[DMPs$adj.P.Val < 0.1, ]
head(significant_DMPs)
# Seleccionar los 10 sitios con mayor diferencia en metilación (ordenados por logFC o adj.P.Val)
top_DMPs <- head(significant_DMPs[order(significant_DMPs$logFC, decreasing = TRUE), ], 10)

library(ggplot2)

# Extraer datos de metilación para los sitios seleccionados
top_DMP_data <- Matriz_met_Flt[rownames(top_DMPs), ]

# Transformar a formato largo para ggplot
plot_data <- reshape2::melt(top_DMP_data)
colnames(plot_data) <- c("CpG_site", "Sample", "Methylation_Value")

# Agregar los valores de ACE score
plot_data$ACE_score <- datos_fenotipicos$ACE_score[match(plot_data$Sample, rownames(datos_fenotipicos))]

# Graficar cada CpG
ggplot(plot_data, aes(x = factor(ACE_score), y = Methylation_Value)) +
  geom_boxplot(aes(fill = factor(ACE_score))) +
  facet_wrap(~ CpG_site, scales = "free_y") +
  labs(x = "ACE score", y = "Methylation Value", title = "DMPs by ACE score") +
  theme_minimal()


#---11.1 Identificacion de DMRs----
myAnnotation <- cpg.annotate(
  object = Matriz_met_Flt, 
  datatype = "array", 
  what = "M", 
  analysis.type = "differential", 
  design = design, 
  contrasts = TRUE, 
  cont.matrix = contMatrix, 
  coef = "Nivel1 - Nivel0",  
  arraytype = "EPIC")
#fdr=0.3)

str(myAnnotation)
# Identificar DMRs 
DMRs <- dmrcate(myAnnotation, lambda = 1000, C = 2, pcutoff = 0.05)
#DMRs$results
str(DMRs)
length(DMRs@min_smoothed_fdr)
sort(DMRs@min_smoothed_fdr)
# Extraer los resultados
results.ranges <- extractRanges(DMRs)


#guardar los resultados en un excel
write.table(results.ranges, file="DMRs2_results.csv", sep=";", row.names=FALSE)

#11.2 graficar DMRs----

library(rtracklayer)
ucscTables <- ucscTables(genome="hg19", track="NCBI RefSeq")
print(ucscTables)
session <- browserSession("UCSC")
genome(session) <- "hg19"
query <- ucscTableQuery(session, track="NCBI RefSeq")
availableTables <- tableNames(query)
print(availableTables)
query$table <- "refGene"
trackData <- getTable(query)





# Colores para las muestras según el grupo MTR
groups <- c("red", "blue") # rojo para "Nivel1", azul para "Nivel0"
names(groups) <- levels(factor(datos_fenotipicos$MTR))
cols <- groups[as.character(factor(datos_fenotipicos$MTR))]

betaVals <- getBeta(mSetSqFlt)


DMR.plot(ranges = results.ranges, dmr = 17, CpGs = Matriz_met_Flt, phen.col = cols, 
         what = "M", arraytype = "EPIC", genome = "hg19", group.means=TRUE)

DMR.plot(ranges = results.ranges, dmr = 3, CpGs = betaVals, phen.col = cols, 
         what = "Beta", arraytype = "EPIC", genome = "hg19", group.means=TRUE)

#?DMR.plot

results.ranges$no.cpgs
#11.3 Gráficos más facheros----
# indicar el genoma
gen <- "hg19"
# el número de DMR que graficaremos 
dmrIndex <- 7 

# extraer el número de cromosoma y la localización de DMR results 
chrom <- as.character(seqnames(results.ranges[dmrIndex]))
start <- as.numeric(start(results.ranges[dmrIndex]))
end <- as.numeric(end(results.ranges[dmrIndex]))

# establecer límites del gráfico 
minbase <- start - (7*(end-start))
maxbase <- end + (7*(end-start))

# cargar CpG islands
islandHMM <- read.csv(paste0(dataDirectory,
                             "/model-based-cpg-islands-hg19.txt"),
                     sep="\t", stringsAsFactors=FALSE, header=TRUE)


islandData <- GRanges(seqnames=Rle(islandHMM[,1]),
                     ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                    strand=Rle(strand(rep("*",nrow(islandHMM)))))


# cargar DNAseI hypersensitive sites
dnase <- read.csv(paste0(dataDirectory,"/sitioshipersensibles.csv"), 
                 sep="",stringsAsFactors=FALSE,header=TRUE)

dnaseData <- GRanges(seqnames=dnase[,2],
                    ranges=IRanges(start=dnase[,3], end=dnase[,4]),
                   strand=Rle(rep("*",nrow(dnase))),
                  data=dnase[,6])

# crear tracks de contexto genómico
iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
#rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
 #                   from=minbase, to=maxbase, trackType="GeneRegionTrack", 
  #                  rstarts="exonStarts", rends="exonEnds", gene="name", 
   #                 symbol="name2", transcript="name", strand="strand", 
    #                fill="darkblue",stacking="squish", name="RefSeq", 
     #               showId=TRUE, geneSymbol=TRUE)
rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                    from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                    rstarts="exonStarts", rends="exonEnds", gene="name", 
                    symbol="name2", transcript="name", strand="strand", 
                    fill="darkblue", stacking="squish", name="RefSeq", 
                    showId=TRUE, geneSymbol=TRUE, 
                   transcriptAnnotation="symbol")



#Ordenar los datos por cromosoma y posición de las bases
annEPICOrd <- annEPICSub[order(annEPICSub$chr,annEPICSub$pos),]

#alinear los bvals, para que las posiciones genómicas de las sondas y sus valores de metilación estén alineados en el gráfico
bValsOrd <- betaVals[match(annEPICOrd$Name,rownames(betaVals)),]

#Crear un GRanges con los CpGs con su posición genómica y valores beta asociados
cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                   ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                   strand=Rle(rep("*",nrow(annEPICOrd))),
                   betas=bValsOrd)

# extraer datos de CpGs en el DMR
cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])

# Crear el track de metilación, agrupando por MTR
methTrack <- DataTrack(range=cpgData, groups=MTR, genome=gen,
                       chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                       type=c("a","p"), name="DNA Meth.\n(beta value)",
                       background.panel="white", legend=TRUE, cex.title=0.8,
                       cex.axis=0.8, cex.legend=0.8)

# Track de CpG island 
                              
# Filtrar las islas CpG solo para el DMR elegido
islandDataSub <- subset(islandData, seqnames(islandData) == chrom & 
                          start(islandData) >= minbase & 
                          end(islandData) <= maxbase)

#crear el track 
islandTrack <- AnnotationTrack(range=islandDataSub, genome=gen, name="CpG Is.", 
                               chromosome=chrom, fill="darkgreen")


# DNaseI hypersensitive site data track
dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                        type="gradient", chromosome=chrom)
#Track de la posición del DMR
dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                            chromosome=chrom,fill="darkred")
#Generar gráfico
tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack, rTrack)
sizes <- c(0.5,2,5,0.5,0.5,3,5)  # ajustar los tamaños relativos de los tracks 
plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
           add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes)


#plotTracks(tracks, from=1230000, to=1244200, showTitle=TRUE, add53=TRUE, 
#           add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes)
#plotTracks(iTrack, from=minbase, to=maxbase)
#plotTracks(gTrack, from=minbase, to=maxbase)
#plotTracks(methTrack, from=minbase, to=maxbase)
#plotTracks(rTrack, from=minbase, to=maxbase)
#plotTracks(islandTrack, from=minbase, to=maxbase)

#crear función DMRplot----
DMRplot <- function(dmrIndex, zoomout1 , zoomout2) {
  # indicar el genoma
  gen <- "hg19"

  # extraer el número de cromosoma y la localización de DMR results 
  chrom <- as.character(seqnames(results.ranges[dmrIndex]))
  start <- as.numeric(start(results.ranges[dmrIndex]))
  end <- as.numeric(end(results.ranges[dmrIndex]))
  
  # establecer límites del gráfico 
  minbase <- start - (zoomout1 *(end-start))
  maxbase <- end + (zoomout2 *(end-start))
  
  # cargar CpG islands
  islandHMM <- read.csv(paste0(dataDirectory,
                               "/model-based-cpg-islands-hg19.txt"),
                        sep="\t", stringsAsFactors=FALSE, header=TRUE)
  
  
  islandData <- GRanges(seqnames=Rle(islandHMM[,1]),
                        ranges=IRanges(start=islandHMM[,2], end=islandHMM[,3]),
                        strand=Rle(strand(rep("*",nrow(islandHMM)))))
  
  
  # cargar DNAseI hypersensitive sites
  dnase <- read.csv(paste0(dataDirectory,"/sitioshipersensibles.csv"), 
                    sep="",stringsAsFactors=FALSE,header=TRUE)
  
  dnaseData <- GRanges(seqnames=dnase[,2],
                       ranges=IRanges(start=dnase[,3], end=dnase[,4]),
                       strand=Rle(rep("*",nrow(dnase))),
                       data=dnase[,6])
  
  # crear tracks de contexto genómico
  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name="")
  gTrack <- GenomeAxisTrack(col="black", cex=1, name="", fontcolor="black")
  
  rTrack <- UcscTrack(genome=gen, chromosome=chrom, track="NCBI RefSeq", 
                      from=minbase, to=maxbase, trackType="GeneRegionTrack", 
                      rstarts="exonStarts", rends="exonEnds", gene="name", 
                      symbol="name2", transcript="name", strand="strand", 
                      fill="darkblue", stacking="squish", name="RefSeq", 
                      showId=TRUE, geneSymbol=TRUE, 
                      transcriptAnnotation="symbol")
  
  #Ordenar los datos por cromosoma y posición de las bases
  annEPICOrd <- annEPICSub[order(annEPICSub$chr,annEPICSub$pos),]
  
  #alinear los bvals, para que las posiciones genómicas de las sondas y sus valores de metilación estén alineados en el gráfico
  bValsOrd <- betaVals[match(annEPICOrd$Name,rownames(betaVals)),]
  
  #Crear un GRanges con los CpGs con su posición genómica y valores beta asociados
  cpgData <- GRanges(seqnames=Rle(annEPICOrd$chr),
                     ranges=IRanges(start=annEPICOrd$pos, end=annEPICOrd$pos),
                     strand=Rle(rep("*",nrow(annEPICOrd))),
                     betas=bValsOrd)
  
  # extraer datos de CpGs en el DMR
  cpgData <- subsetByOverlaps(cpgData, results.ranges[dmrIndex])
  
  # Crear el track de metilación, agrupando por MTR
  methTrack <- DataTrack(range=cpgData, groups=MTR, genome=gen,
                         chromosome=chrom, ylim=c(-0.05,1.05), col=pal,
                         type=c("a","p"), name="DNA Meth.\n(beta value)",
                         background.panel="white", legend=TRUE, cex.title=0.8,
                         cex.axis=0.8, cex.legend=0.8)
  
  # Track de CpG island 
  islandDataSub <- subset(islandData, seqnames(islandData) == chrom & 
                            start(islandData) >= minbase & 
                            end(islandData) <= maxbase)
  islandTrack <- AnnotationTrack(range=islandDataSub, genome=gen, name="CpG Is.", 
                                 chromosome=chrom, fill="darkgreen")
  
  
  # DNaseI hypersensitive site data track
  dnaseTrack <- DataTrack(range=dnaseData, genome=gen, name="DNAseI", 
                          type="gradient", chromosome=chrom)
  #Track de la posición del DMR
  dmrTrack <- AnnotationTrack(start=start, end=end, genome=gen, name="DMR", 
                              chromosome=chrom,fill="darkred")
  #Generar gráfico
  tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack, rTrack)
  sizes <- c(0.5,2,5,0.5,0.5,3,5)  # ajustar los tamaños relativos de los tracks 
  plotTracks(tracks, from=minbase, to=maxbase, showTitle=TRUE, add53=TRUE, 
             add35=TRUE, grid=TRUE, lty.grid=3, sizes=sizes)
}
#11.4 Gráficos más facheros usando mi función----
DMRplot(dmrIndex = 26, zoomout1 = 2.5, zoomout2 = 0.5)

