#1.1 IMPORT DATA y lectura de IDAT----

dataDirectory <- "C:/Users/Karen/Documents/Tesis"
targets <- read.metharray.sheet(dataDirectory, pattern="Scorza_Project_001_Sample_Sheet.csv")
rgSet <- read.metharray.exp(dataDirectory,recursive=TRUE)
detP <- detectionP(rgSet) 

#1.2 Poner a detP los nombres de las muestras----
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

# Comprobar los nombres de las columnas de detP
print(colnames(detP))


#2. NORMALICACIÓN----
mSetSq <- preprocessQuantile(rgSet)

#2.1 Renombrar las columnas de mSetSq para que sean los nombres de las muestras----

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
comparison_table <- data.frame(Archivo = original_column_names, Muestra = colnames(mSetSq))


#3. Filtrado----
detP <- detP[match(featureNames(mSetSq),rownames(detP)),] 
keep <- rowSums(detP < 0.01) == ncol(mSetSq) 
table(keep)

mSetSqFlt <- mSetSq[keep,]

mSetSqFlt <- dropLociWithSnps(mSetSqFlt) #elimina las sondas que pueden tener SNPs comunes que afectan el CpG

#3.1 MDS----

Matriz_met_Flt<-getM(mSetSqFlt) #devuelve una matriz que tiene los datos de metilación M para todas las sondas y muestras (usa el conjunto de datos normalizado Y FILTRADO mSetSqFlt)#
Matriz_met_Flt<- Matriz_met_Flt[,colnames(Matriz_met_Flt)!="NA10858_2"] #Eliminar el NA#
plotMDS(Matriz_met_Flt, top=1000, gene.selection="common")  
#3.2 MDS: Agregar color por grupo ARRAY-------
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
