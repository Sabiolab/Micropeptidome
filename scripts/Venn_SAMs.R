rm(list = ls())

library(ggplot2)
library(dplyr)
library(readr)
library(VennDiagram)
library(grid)

# Leer los archivos CSV
datos_visceral <- read_csv("Visceral.locus_summary.all_loci.csv")
datos_liver <- read_csv("Liver.locus_summary.all_loci.csv")

# Filtrar por probabilidad
prob_threshold <- 0.9

visceral_filtrado <- datos_visceral %>%
  filter(max_prob >= prob_threshold)

liver_filtrado <- datos_liver %>%
  filter(max_prob >= prob_threshold)

# Obtener los conjuntos de SAMs (locus)
sams_visceral <- unique(visceral_filtrado$locus)
sams_liver <- unique(liver_filtrado$locus)

# Calcular intersección y diferencias
sams_comunes <- intersect(sams_visceral, sams_liver)
sams_solo_visceral <- setdiff(sams_visceral, sams_liver)
sams_solo_liver <- setdiff(sams_liver, sams_visceral)


# Colores personalizados (con transparencia)
colores_personalizados <- c("#1571ae", "#e8781b")

# Diagrama de Venn
venn_plot <- venn.diagram(
  x = list(
    Visceral = sams_visceral,
    Liver = sams_liver
  ),
  category.names = c("Visceral fat", "Liver"),
  filename = NULL,
  output = TRUE,
  
  # Colores
  fill = colores_personalizados,
  alpha = 0.5,
  
  # Bordes
  lwd = 2,
  col = colores_personalizados,
  
  # Tamaños de texto
  cex = 1.2,
  fontface = "bold",
  fontfamily = "sans",
  
  # Categorías
  cat.cex = 1.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27),
  cat.dist = c(0.055, 0.055),
  cat.fontfamily = "sans",
  cat.col = colores_personalizados,
  
  # Hacer círculos proporcionales al tamaño
  scaled = TRUE,
  
  # Espaciado
  ext.text = TRUE,
  ext.line.lwd = 2,
  ext.dist = -0.15,
  ext.length = 0.9,
  ext.pos = c(0, 180),
  
  # Formatear números con separador de miles
  print.mode = "raw"
)

# Formatear los números en el diagrama con separador de miles
for (i in 1:length(venn_plot)) {
  if (grepl("text", venn_plot[[i]]$name)) {
    if (!is.null(venn_plot[[i]]$label)) {
      # Intentar convertir a número y formatear
      num_value <- suppressWarnings(as.numeric(venn_plot[[i]]$label))
      if (!is.na(num_value)) {
        venn_plot[[i]]$label <- format(num_value, big.mark = ",", scientific = FALSE)
      }
    }
  }
}

# Mostrar el diagrama
grid.newpage()
pushViewport(viewport(x = 0.5, y = 0.45, width = 0.92, height = 0.9))
grid.draw(venn_plot)
popViewport()
grid.text("Microproteins overlap", 
          x = 0.5, y = 0.92, 
          gp = gpar(fontsize = 16, fontface = "bold"))
