global_palette <- list(teal = "#008080",
                       coral = "#F08080")

gl_condition_colors <- c(
  "Naive_EZH2i" = "#76c6c7", # "#5F9EA0",
  "Naive_Untreated" = "#278b8b",
  "Primed_EZH2i" = "#f48c8c", #  "#f47770",
  "Primed_Untreated" = "#f44b34"
)

theme_default <- function(base_size = 12) {
  theme_classic(base_size = base_size)
}


my_svg <- function(file, width, height) {
  library(svglite)
  svglite(file = file, width = width, height = height, bg = "white")
}
