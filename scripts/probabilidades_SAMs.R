rm(list = ls())

library(ggplot2)
library(dplyr)
library(readr)
library(ggridges) # para el Ridge plot

# Leer los archivos CSV
datos_visceral <- read_csv("Visceral.locus_summary.all_loci.csv")
datos_liver <- read_csv("Liver.locus_summary.all_loci.csv")

# Añadir columna de tejido
datos_visceral$tissue <- "Visceral fat"
datos_liver$tissue <- "Liver"

# Combinar datos
datos <- rbind(datos_visceral, datos_liver)

# Factor para el orden
datos$tissue <- factor(datos$tissue, levels = c("Visceral fat", "Liver"))

# Colores personalizados
colores_personalizados <- c(
  "Visceral fat" = "#1571ae",
  "Liver" = "#e8781b"
)

# Calcular estadísticas para cada tejido
stats <- datos %>%
  group_by(tissue) %>%
  summarise(
    media = mean(max_prob, na.rm = TRUE),
    mediana = median(max_prob, na.rm = TRUE),
    n_total = n(),
    n_above_09 = sum(max_prob >= 0.9, na.rm = TRUE),
    perc_above_09 = (n_above_09 / n_total) * 100
  )

# Mostrar estadísticas
print("Estadísticas de probabilidad por tejido:")
print(stats)


p_ridge <- ggplot(datos, aes(x = max_prob, y = tissue, fill = tissue)) +
  geom_density_ridges(alpha = 0.7, scale = 0.9, linewidth = 1.2) +
  geom_vline(xintercept = 0.9, linetype = "dashed", color = "black", 
             linewidth = 0.8) +
  scale_fill_manual(values = colores_personalizados) +
  labs(
    title = "Distribution of probabilities",
    x = "Maximum probability",
    y = ""
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1),
    expand = c(0.01, 0)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold", size = 11),
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_blank(),
    panel.border = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.y = element_text(face = "bold")
  ) +
  annotate("text", x = 0.9, y = 2.4, 
           label = "Threshold = 0.9", 
           hjust = 1.1, vjust = 0, size = 3.5, fontface = "italic")

print(p_ridge)


# Crear tabla resumen con diferentes umbrales para ayudarte a decidir
umbrales <- seq(0.5, 0.9, by = 0.05)
tabla_umbrales <- data.frame()

for (umbral in umbrales) {
  temp <- datos %>%
    group_by(tissue) %>%
    summarise(
      umbral = umbral,
      n_sams = sum(max_prob >= umbral, na.rm = TRUE),
      porcentaje = (n_sams / n()) * 100
    )
  tabla_umbrales <- rbind(tabla_umbrales, temp)
}

print("\nNúmero y porcentaje de SAMs según diferentes umbrales:")
print(tabla_umbrales %>% arrange(tissue, umbral))

# Guardar gráficos
ggsave("SAM_probability_density.png", plot = p, width = 10, height = 6, dpi = 300)
ggsave("SAM_probability_histogram.png", plot = p_hist, width = 10, height = 6, dpi = 300)

# Guardar tabla de umbrales
write.csv(tabla_umbrales, "SAM_probability_thresholds.csv", row.names = FALSE)