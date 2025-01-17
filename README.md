## README: Análisis de Metilación Diferencial
README: Análisis de Metilación Diferencial - Proyecto Metiloma Placentas

Este repositorio contiene dos scripts principales desarrollados para el análisis de metilación diferencial en datos de metilomas humanos obtenidos mediante arrays EPIC. Los scripts deben ejecutarse en orden y realizan el preprocesamiento de datos y el análisis diferencial de metilación, respectivamente.

Contexto del análisis

Estos scripts fueron desarrollados para analizar muestras humanas de placentas con el objetivo de identificar patrones de metilación diferencial asociados con experiencias de adversidad en la infancia. Utilizan datos obtenidos mediante la tecnología Illumina EPIC array, enfocándose en CpGs vinculados a genes y secuencias regulatorias. Las muestras se clasifican en grupos de adversidad que se definen de acuerdo al tipo de adversidades reportadas:
* Grupo A: Adversidades directas como abuso emocional, físico o sexual, y negligencia emocional o física.
* Grupo B: Adversidades indirectas como divorcio, violencia contra la madre, consumo de sustancias, enfermedad mental de los padres o encarcelamiento.
* Grupo C: Sin adversidades reportadas.


**Script 1. Preprocesamiento**

Archivo: 01_preprocesamiento.R
Objetivo: Preparar los datos fenotípicos y de metilación para el análisis diferencial. Esto incluye normalización de datos, eliminación de muestras duplicadas, filtrado de CpGs y exploración inicial de los datos.

Inputs:
* Archivos IDAT provenientes de arrays EPIC.
* Archivo CSV con datos fenotípicos que incluyen información del grupo de adversidad de cada muestra.

Outputs:
Archivo RData: preprocesado.RData, que contiene:
* datos_fenotipicos: Datos fenotípicos procesados.
* Matriz_met_Flt: Matriz de valores de metilación filtrada y normalizada.
* mSetSqFlt: Objeto con datos normalizados y filtrados.

Ejecución:
1. Configurar el directorio de trabajo.
2. Asegurarse de que los archivos IDAT y el archivo CSV de datos fenotípicos estén disponibles en las rutas especificadas en el script.
3. Correr el script para generar los datos preprocesados.




**Script 2. Análisis de Metilación Diferencial**

Archivo: 02_analisis_de_metilacion_diferencial_adversidad.R
Objetivo:
Identificar posiciones CpG y regiones diferencialmente metiladas (DMPs y DMRs) entre diferentes grupos de adversidad, y asociar estas regiones con genes relevantes.

Inputs:
* Archivo RData generado en el preprocesamiento: preprocesado.RData.

Outputs:
* Archivo CSV con DMPs: DMPs_adversidad.csv.
* Archivos CSV con DMRs por comparación de grupos (e.g., DMRs_GrupoAvsC.csv).
* Gráficos de resultados, incluyendo:
* Gráficos de CpGs más significativos.
* Visualizaciones de DMRs en contexto genómico.

Ejecución:
1. Asegurarse de que el archivo preprocesado.RData esté disponible en el directorio de trabajo.
2. Configurar el directorio de trabajo en el script.
3. Correr el script para realizar el análisis diferencial y generar los resultados.




Si tenés preguntas o comentarios, por favor contactá a: karenchurba@gmail.com
