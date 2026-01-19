rm(list = ls())

library(ggplot2)
library(dplyr)

# Crear datos
datos <- data.frame(
  categoria = c("RAW", "SAX", "SEERa", "SEERb"),
  valor = c(123, 265, 386, 193)
)

# Asegurar el orden de las categorías
datos$categoria <- factor(datos$categoria, 
                          levels = c("RAW", "SAX", "SEERa", "SEERb"))

# Colores personalizados
colores_personalizados <- c(
  "RAW" = "#35943a",
  "SAX" = "#1571ae",
  "SEERa" = "#e8781b",
  "SEERb" = "#d62728"
)

# Crear gráfico de barras
p <- ggplot(datos, aes(x = categoria, y = valor, fill = categoria)) +
  geom_bar(stat = "identity", alpha = 0.70, width = 0.7) +
  scale_fill_manual(values = colores_personalizados) +
  labs(
    title = "Identified microproteins",
    x = "",
    y = "Microprotein peptides (no.)"
  ) +
  scale_y_continuous(
    limits = c(0, 400),
    breaks = seq(0, 400, by = 50),
    expand = c(0, 0)
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.line = element_line(color = "black"),
    axis.ticks.x = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 11, face = "bold")
  )

# Mostrar gráfico
print(p)