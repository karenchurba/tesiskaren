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
library("VennDiagram") 

#9. cargar archivos preprocesados----
setwd("C:/Users/Karen/Documents/Tesis")

load("preprocesadosAdversidad.RData")
pal <- brewer.pal(8,"Dark2")

##Obtener la información de sitios de metilación de genoma humano y guardarla en annEPIC
#annEPIC <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

betaVals <- getBeta(mSetSqFlt)

#---10. ANÁLISIS DE METILACIÓN DIFERENCIAL POR SONDA-------

# Crear el factor de interés (Adversidad)
Adversidad <- factor(datos_fenotipicos$Grupo.según.tipo.de.adversidad)
#levels(Adversidad) <- c("A", "B","C")

# Crear una matriz de diseño
design <- model.matrix(~0 + Adversidad, data=datos_fenotipicos)
rownames(design) <- datos_fenotipicos$Muestra.Metiloma
colnames(design) <- c("A", "B", "C")
print(design)

#dim(Matriz_met_Flt)
#dim(design)

# Ajustar el modelo lineal
fit <- lmFit(Matriz_met_Flt, design)

print(fit)

# Crear una matriz de contrastes para la comparación de interés (por ejemplo, MTR 1 vs MTR 0)
contMatrix <- makeContrasts(
  A - B,  # Contraste entre grupo adversidad A y B
  A - C,  # Contraste entre grupo adversidad A y C
  B - C,  # Contraste entre grupo adversidad B y C
  levels = design
)


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
  plotCpg(getBeta(mSetSqFlt), cpg=cpg, pheno=datos_fenotipicos$Grupo.según.tipo.de.adversidad, ylab="Beta values")
})

# Guardar los resultados en un archivo CSV
write.table(DMPs, file="DMPs_adversidad.csv", sep=",", row.names=FALSE)


#---11.1 Identificacion de DMRs----
#Generar fx para asignar IDs a los DMRs siguiendo esta estructura: comparación_cromosoma_start_end
generate_ID <- function(gr, comparison) {
  paste0(comparison, "_", seqnames(gr), "_", start(gr), "_", end(gr))
}

#---11.1AC Adversidad A vs C----
myAnnotationAC <- cpg.annotate(
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

str(myAnnotationAC)
# Identificar DMRs 
DMRsAC <- dmrcate(myAnnotationAC, lambda = 1000, C = 2, pcutoff = 0.05)

length(DMRsAC@min_smoothed_fdr)
sort(DMRsAC@min_smoothed_fdr)

# Extraer los resultados
results.rangesAC <- extractRanges(DMRsAC)
length(results.rangesAC)

#Agrego ID a los DMRs 
mcols(results.rangesAC)$ID <- generate_ID(results.rangesAC, "AC")
head(results.rangesAC)

#guardar los resultados para los grupos A vs C en un archivo CSV. Cada fila representa un DMR#
write.table(results.rangesAC, file="DMRs_GrupoAvsC.csv", sep=";", row.names=FALSE) 


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

#str(myAnnotationBC)
# Identificar DMRs 
DMRsBC <- dmrcate(myAnnotationBC, lambda = 1000, C = 2, pcutoff = 0.05)
#DMRs$results
str(DMRsBC)
length(DMRsBC@min_smoothed_fdr)

# Extraer los resultados
results.rangesBC <- extractRanges(DMRsBC)
length(results.rangesBC)
#Agregar ID
mcols(results.rangesBC)$ID <- generate_ID(results.rangesBC, "BC")
head(results.rangesBC)

#guardar los resultados para los grupos B vs C en un archivo CSV. Cada fila representa un DMR#
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


# Identificar DMRs 
DMRsBA <- dmrcate(myAnnotationBA, lambda = 1000, C = 2, pcutoff = 0.05)

#str(DMRsBA)
length(DMRsBA@min_smoothed_fdr)
sort(DMRsBA@min_smoothed_fdr)

# Extraer los resultados
results.rangesBA <- extractRanges(DMRsBA)

mcols(results.rangesBA)$ID <- generate_ID(results.rangesBA, "BA")

#guardar los resultados para los grupos A vs B en un archivo CSV. Cada fila representa un DMR#
write.table(results.rangesBA, file="DMRs_GrupoAB.csv", sep=";", row.names=FALSE)

#11.2 graficar DMRs----
# indicar el genoma
gen <- "hg19"
# el número de DMR que graficaremos 
dmrIndex <- 1 

# extraer el número de cromosoma y la localización de DMR results 
chrom <- as.character(seqnames(results.rangesBC[dmrIndex]))
start <- as.numeric(start(results.rangesBC[dmrIndex]))
end <- as.numeric(end(results.rangesBC[dmrIndex]))

# establecer límites del gráfico 
minbase <- start - (10*(end-start))
maxbase <- end + (10*(end-start))

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
cpgData <- subsetByOverlaps(cpgData, results.rangesBC[dmrIndex])

# Crear el track de metilación, agrupando por MTR
methTrack <- DataTrack(range=cpgData, groups=Adversidad, genome=gen,
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

#11.2 crear función DMRplot----
DMRplot2 <- function(TablaDMRs, dmrIndex, zoomout1, zoomout2) {
  # Indicar el genoma
  gen <- "hg19"
  
  # Extraer el número de cromosoma y la localización del DMR
  chrom <- as.character(seqnames(TablaDMRs[dmrIndex]))
  start <- start(TablaDMRs[dmrIndex])
  end <- end(TablaDMRs[dmrIndex])
  
  # Establecer límites del gráfico
  minbase <- start - (zoomout1 * (end - start))
  maxbase <- end + (zoomout2 * (end - start))
  
  # Cargar CpG islands
  islandHMM <- read.csv(paste0(dataDirectory, "/model-based-cpg-islands-hg19.txt"),
                        sep = "\t", stringsAsFactors = FALSE, header = TRUE)
  
  islandData <- GRanges(seqnames = Rle(islandHMM[, 1]),
                        ranges = IRanges(start = islandHMM[, 2], end = islandHMM[, 3]),
                        strand = Rle(strand(rep("*", nrow(islandHMM)))))
  
  # Cargar DNAseI hypersensitive sites
  dnase <- read.csv(paste0(dataDirectory, "/sitioshipersensibles.csv"), 
                    sep = "", stringsAsFactors = FALSE, header = TRUE)
  
  dnaseData <- GRanges(seqnames = dnase[, 2],
                       ranges = IRanges(start = dnase[, 3], end = dnase[, 4]),
                       strand = Rle(rep("*", nrow(dnase))),
                       data = dnase[, 6])
  
  # Crear tracks de contexto genómico
  iTrack <- IdeogramTrack(genome = gen, chromosome = chrom, name = "")
  gTrack <- GenomeAxisTrack(col = "black", cex = 1, name = "", fontcolor = "black")
  
  rTrack <- UcscTrack(genome = gen, chromosome = chrom, track = "NCBI RefSeq", 
                      from = minbase, to = maxbase, trackType = "GeneRegionTrack", 
                      rstarts = "exonStarts", rends = "exonEnds", gene = "name", 
                      symbol = "name2", transcript = "name", strand = "strand", 
                      fill = "darkblue", stacking = "squish", name = "RefSeq", 
                      showId = TRUE, geneSymbol = TRUE, 
                      transcriptAnnotation = "symbol")
  
  # Ordenar los datos por cromosoma y posición de las bases
  annEPICOrd <- annEPICSub[order(annEPICSub$chr, annEPICSub$pos),]
  
  # Alinear los bvals
  bValsOrd <- betaVals[match(annEPICOrd$Name, rownames(betaVals)),]
  
  # Crear un GRanges con los CpGs con su posición genómica y valores beta asociados
  cpgData <- GRanges(seqnames = Rle(annEPICOrd$chr),
                     ranges = IRanges(start = annEPICOrd$pos, end = annEPICOrd$pos),
                     strand = Rle(rep("*", nrow(annEPICOrd))),
                     betas = bValsOrd)
  
  # Extraer datos de CpGs en el DMR
  cpgData <- subsetByOverlaps(cpgData, TablaDMRs[dmrIndex])
  
  # Crear el track de metilación, agrupando por MTR
  methTrack <- DataTrack(range = cpgData, groups = Adversidad, genome = gen,
                         chromosome = chrom, ylim = c(-0.05, 1.05), col = pal,
                         type = c("a", "p"), name = "DNA Meth.\n(beta value)",
                         background.panel = "white", legend = TRUE, cex.title = 0.8,
                         cex.axis = 0.8, cex.legend = 0.8)
  
  # Track de CpG island
  islandDataSub <- subset(islandData, seqnames(islandData) == chrom & 
                            start(islandData) >= minbase & 
                            end(islandData) <= maxbase)
  islandTrack <- AnnotationTrack(range = islandDataSub, genome = gen, name = "CpG Is.", 
                                 chromosome = chrom, fill = "darkgreen")
  
  # DNaseI hypersensitive site data track
  dnaseTrack <- DataTrack(range = dnaseData, genome = gen, name = "DNAseI", 
                          type = "gradient", chromosome = chrom)
  # Track de la posición del DMR
  dmrTrack <- AnnotationTrack(start = start, end = end, genome = gen, name = "DMR", 
                              chromosome = chrom, fill = "darkred")
  # Generar gráfico
  tracks <- list(iTrack, gTrack, methTrack, dmrTrack, islandTrack, dnaseTrack, rTrack)
  sizes <- c(0.5, 2, 5, 0.5, 0.5, 3, 5)  # ajustar los tamaños relativos de los tracks 
  plotTracks(tracks, from = minbase, to = maxbase, showTitle = TRUE, add53 = TRUE, 
             add35 = TRUE, grid = TRUE, lty.grid = 3, sizes = sizes)
}
#11.2 Graficar DMRs----
DMRplot2(TablaDMRs=results.rangesBA , dmrIndex = 67, zoomout1 = 1 , zoomout2 = 0.5)


#12. Comparación entre DMRs
#12.1 Diagramas de Venn----

genesAC <- unique(unlist(strsplit(as.character(results.rangesAC$overlapping.genes), ";")))
genesBC <- unique(unlist(strsplit(as.character(results.rangesBC$overlapping.genes), ";")))
genesBA <- unique(unlist(strsplit(as.character(results.rangesBA$overlapping.genes), ";")))
genesAC <- na.omit(genesAC)
genesBC <- na.omit(genesBC)
genesBA <- na.omit(genesBA)

grid.newpage() 
venn.plot <- draw.triple.venn(
  area1 = length(genesAC),
  area2 = length(genesBC),
  area3 = length(genesBA),
  n12 = length(intersect(genesAC, genesBC)),
  n23 = length(intersect(genesBC, genesBA)),
  n13 = length(intersect(genesAC, genesBA)),
  n123 = length(Reduce(intersect, list(genesAC, genesBC, genesBA))),
  category = c("AC", "BC", "BA"),
  col = c("red", "blue", "green"),
  fill = c("lightpink", "lightblue", "lightgreen"),
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20, -90)
)

#versión 2: separando los genes que corresponden a un mismo DMR----
genesAC <- unique(unlist(strsplit(as.character(results.rangesAC$overlapping.genes), split = ",\\s*")))
genesBC <- unique(unlist(strsplit(as.character(results.rangesBC$overlapping.genes), split = ",\\s*")))
genesBA <- unique(unlist(strsplit(as.character(results.rangesBA$overlapping.genes), split = ",\\s*")))
genesAC <- na.omit(genesAC)
genesBC <- na.omit(genesBC)
genesBA <- na.omit(genesBA)
grid.newpage()
venn.plot <- draw.triple.venn(
  area1 = length(genesAC),
  area2 = length(genesBC),
  area3 = length(genesBA),
  n12 = length(intersect(genesAC, genesBC)),
  n23 = length(intersect(genesBC, genesBA)),
  n13 = length(intersect(genesAC, genesBA)),
  n123 = length(Reduce(intersect, list(genesAC, genesBC, genesBA))),
  category = c("AC", "BC", "BA"),
  col = c("red", "blue", "green"),
  fill = c("lightpink", "lightblue", "lightgreen"),
  cex = 1.5,
  cat.cex = 1.2,
  cat.pos = c(-20, 20, -90)
)


#versión 3: con los porcentajes----

library(ggvenn)

# Crear un data frame con la lista completa de genes y en qué comparaciones aparecen 
all_genes <- unique(c(genesAC, genesBC, genesBA)) 
venn_data <- data.frame(
  Gene = all_genes,
  AC = all_genes %in% genesAC,
  BC = all_genes %in% genesBC,
  BA = all_genes %in% genesBA
)
#head(venn_data)

# Convertir a formato tibble para ggvenn
library(tibble)
venn_data_tibble <- as_tibble(venn_data)

# Crear el gráfico de Venn
ggvenn(
  data = venn_data_tibble,
  columns = c("AC", "BC", "BA"),
  show_percentage = TRUE,  
  fill_color = c("lightpink", "lightblue", "lightgreen"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 2.5
)
?ggvenn

###Obtener la tabla de los genes en cada comparación----

gene_comparison_table <- data.frame(
  Gene = all_genes,
  AC = all_genes %in% genesAC,
  BC = all_genes %in% genesBC,
  BA = all_genes %in% genesBA
)

# Crear una columna que indique en qué comparaciones está cada gen
gene_comparison_table$Comparaciones <- apply(gene_comparison_table[, 2:4], 1, function(x) {
  comparaciones <- names(x)[x]
  paste(comparaciones, collapse = "_")
})

head(gene_comparison_table)
#guardar los resultados en un excel
write.table(gene_comparison_table, file="comparacion genes por grupo.csv", sep=";", row.names=FALSE)


# Encontrar superposiciones de posiciones DMRs entre las tres tablas----
overlaps_AC_BA <- findOverlaps(results.rangesAC, results.rangesBA)
overlaps_AC_BC <- findOverlaps(results.rangesAC, results.rangesBC)
overlaps_BA_BC <- findOverlaps(results.rangesBA, results.rangesBC)
#Hits: superposiciones encontradas entre los DMRs 
#queryHits: Índices de los DMRs en el primer objeto que tienen superposiciones.
#subjectHits: Índices de los DMRs en el segundo objeto que se superponen con los del primer objeto.

# DMRs compartidos
dmrs_AC_BA <- results.rangesAC[queryHits(overlaps_AC_BA)]
dmrs_BA_AC <- results.rangesBA[subjectHits(overlaps_AC_BA)]
length(dmrs_AC_BA@seqnames)
length(dmrs_BA_AC@seqnames)

dmrs_AC_BC <- results.rangesAC[queryHits(overlaps_AC_BC)]
dmrs_BC_AC <- results.rangesBC[subjectHits(overlaps_AC_BC)]
length(dmrs_AC_BC@seqnames)
length(dmrs_BC_AC@seqnames)

dmrs_BA_BC <- results.rangesBA[queryHits(overlaps_BA_BC)]
dmrs_BC_BA <- results.rangesBC[subjectHits(overlaps_BA_BC)]
length(dmrs_BA_BC@seqnames)
length(dmrs_BC_BA@seqnames)

# Unifico en un solo GRanges con todos los DMRs compartidos
dmrs_comparacion_red <- reduce(c(dmrs_AC_BA, dmrs_BA_AC, dmrs_AC_BC, dmrs_BC_AC, dmrs_BA_BC, dmrs_BC_BA))

# resultados
dmrs_comparacion
write.table(dmrs_comparacion_red, file="DMRs_comparacion_red.csv", sep=";", row.names=FALSE)


#13.metilación-apego----
#datos_fenotipicos$Grupo.según.tipo.de.adversidad <- as.factor(datos_fenotipicos$Grupo.según.tipo.de.adversidad)
#datos_fenotipicos$ATTACHMENT <- as.factor(datos_fenotipicos$ATTACHMENT)

design <- model.matrix(~ ATTACHMENT + Adversidad + ATTACHMENT:Adversidad, data = datos_fenotipicos)
rownames(design) <- datos_fenotipicos$Muestra.Metiloma
Matriz_met_Flt <- Matriz_met_Flt[, rownames(design)]

fit <- lmFit(Matriz_met_Flt, design)
fit <- eBayes(fit)
topTable(fit, coef = "AdversidadC", adjust = "fdr")
summary(fit)
fit
# Crear una matriz de contrastes para la comparación de interés (por ejemplo, MTR 1 vs MTR 0)
contMatrix <- makeContrasts(
  ATTACHMENT0 - AdversidadC,  # Contraste entre grupo adversidad 3 y 1
  ATTACHMENT1 - C,  # Contraste entre grupo adversidad 6 y 0
  B - C,  # Contraste entre grupo adversidad 4 y 0
  levels = design
)


# Ajustar los contrastes
fit2 <- contrasts.fit(fit, contMatrix)
fit2 <- eBayes(fit2)

summary(decideTests(fit2))


