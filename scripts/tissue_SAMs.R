rm(list = ls())

library(ggplot2)
library(dplyr)

# Datos de distribución de frecuencias para Visceral
datos_visceral <- data.frame(
  n_patients = c(1:91),
  n_sams = c(41735, 10800, 5743, 3883, 2737, 2092, 1760, 1449, 1198, 999,
             864, 817, 647, 638, 585, 535, 465, 436, 441, 404,
             335, 355, 283, 306, 273, 234, 267, 234, 228, 208,
             209, 201, 204, 186, 192, 171, 183, 151, 162, 166,
             147, 136, 136, 130, 119, 121, 126, 103, 99, 102,
             93, 109, 99, 90, 89, 89, 90, 94, 85, 88,
             100, 73, 72, 88, 88, 82, 90, 85, 77, 91,
             99, 72, 90, 87, 86, 101, 100, 111, 88, 85,
             103, 106, 125, 138, 141, 164, 145, 216, 261, 448, 4459)
)

# Datos de distribución de frecuencias para Liver
datos_liver <- data.frame(
  n_patients = c(1:91),
  n_sams = c(40830, 10273, 5468, 3544, 2593, 2016, 1585, 1242, 1120, 974,
             835, 754, 686, 574, 541, 483, 451, 420, 430, 366,
             309, 293, 288, 268, 276, 238, 232, 237, 216, 208,
             205, 161, 134, 151, 205, 178, 151, 134, 134, 122,
             123, 129, 116, 138, 113, 106, 112, 82, 98, 110,
             101, 100, 89, 97, 91, 90, 92, 75, 65, 99,
             68, 77, 94, 75, 84, 74, 73, 81, 73, 76,
             87, 90, 103, 69, 66, 85, 72, 87, 103, 107,
             112, 108, 117, 110, 125, 128, 174, 200, 323, 528, 3681)
)

# Añadir columna de tejido
datos_visceral$tissue <- "Visceral"
datos_liver$tissue <- "Liver"

# Combinar datos
datos <- rbind(datos_visceral, datos_liver)

# Factor para el orden
datos$tissue <- factor(datos$tissue, levels = c("Visceral", "Liver"))

# Colores personalizados (consistente con tu estilo)
colores_personalizados <- c(
  "Visceral" = "#1571ae",
  "Liver" = "#e8781b"
)

# Crear gráfico de líneas
p <- ggplot(datos, aes(x = n_patients, y = n_sams, color = tissue)) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  geom_point(data = datos %>% filter(n_patients == 91), 
             size = 4, alpha = 0.85) +
  scale_color_manual(values = colores_personalizados) +
  labs(
    title = "SAM distribution across patient frequencies",
    x = "Number of patients",
    y = "Number of SAMs",
    color = "Tissue"
  ) +
  scale_x_continuous(
    breaks = c(1, seq(10, 90, by = 10), 91),
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 45000),
    breaks = seq(0, 45000, by = 5000),
    expand = c(0, 0),
    labels = scales::comma
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold", size = 11),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold", size = 11),
    axis.line = element_line(color = "black"),
    panel.border = element_blank(),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = 91, y = 5500, label = "Core SAMs", 
           hjust = -0.1, vjust = 0, size = 3.5, fontface = "italic")

# Mostrar gráfico
print(p)

# Opcional: Versión con escala logarítmica para mejor visualización
p_log <- ggplot(datos, aes(x = n_patients, y = n_sams, color = tissue)) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  geom_point(data = datos %>% filter(n_patients == 91), 
             size = 4, alpha = 0.85) +
  scale_color_manual(values = colores_personalizados) +
  labs(
    title = "Microprotein distribution",
    x = "Number of patients",
    y = "No. of microproteins (log)",
    color = "Tissue"
  ) +
  scale_x_continuous(
    breaks = c(1, seq(10, 100, by = 10)),
    expand = c(0.01, 0)
  ) +
  scale_y_log10(
    breaks = c(10, 50, 100, 500, 1000, 5000, 10000, 50000),
    labels = scales::comma
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.85, 0.80),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
    axis.title.y = element_text(margin = margin(r = 10), face = "bold", size = 11),
    axis.title.x = element_text(margin = margin(t = 10), face = "bold", size = 11),
    axis.line = element_line(color = "black"),
    panel.border = element_blank(),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = 91, y = 4000, label = "Core SAMs", 
           hjust = -0.1, vjust = 0, size = 3.5, fontface = "italic")

# Mostrar gráfico con escala logarítmica
print(p_log)

# Guardar gráficos
ggsave("SAM_distribution_linear.png", plot = p, width = 10, height = 6, dpi = 300)
ggsave("SAM_distribution_log.png", plot = p_log, width = 10, height = 6, dpi = 300)