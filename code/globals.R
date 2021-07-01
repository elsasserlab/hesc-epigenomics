global_palette <- list(teal = "#008080",
                       coral = "#F08080")

gl_condition_colors <- c(
  "Naive_EZH2i" = "#76c6c7", # "#5F9EA0",
  "Naive_Untreated" = "#278b8b",
  "Primed_EZH2i" = "#f48c8c", #  "#f47770",
  "Primed_Untreated" = "#f44b34"
)

gl_mark_colors <- list(
  "H2Aub" = "#8140c1",
  "H3K27m3" = "#3e5aa8",
  "H3K4m3" = "#b64c28"
)

theme_default <- function(base_size = 12) {
  theme_classic(base_size = base_size)
}


my_svg <- function(file, width, height) {
  library(svglite)
  svglite(file = file, width = width, height = height, bg = "white")
}

bwdir <- file.path("./data/bw/Kumar_2020")
bwfiles_rep <-
  list(
    k4_naive = list.files(bwdir, pattern = "H3K4m3_H9_Ni_rep[1-3].hg38.scaled.bw", full.names = T),
    k4_naive_ezh2i = list.files(bwdir, pattern = "H3K4m3_H9_Ni-EZH2i_rep[1-3].hg38.scaled.bw", full.names = T),
    k4_primed = list.files(bwdir, pattern = "H3K4m3_H9_Pr_rep[1-3].hg38.scaled.bw", full.names = T),
    k4_primed_ezh2i = list.files(bwdir, pattern = "H3K4m3_H9_Pr-EZH2i_rep[1-3].hg38.scaled.bw", full.names = T),
    k27_naive = list.files(bwdir, pattern = "H3K27m3_H9_Ni_rep[1-3].hg38.scaled.bw", full.names = T),
    k27_primed = list.files(bwdir, pattern = "H3K27m3_H9_Pr_rep[1-3].hg38.scaled.bw", full.names = T),
    ub_naive = list.files(bwdir, pattern = "H2Aub_H9_Ni_rep[1-3].hg38.scaled.bw", full.names = T),
    ub_naive_ezh2i = list.files(bwdir, pattern = "H2Aub_H9_Ni-EZH2i_rep[1-3].hg38.scaled.bw", full.names = T),
    ub_primed = list.files(bwdir, pattern = "H2Aub_H9_Pr_rep[1-3].hg38.scaled.bw", full.names = T),
    ub_primed_ezh2i = list.files(bwdir, pattern = "H2Aub_H9_Pr-EZH2i_rep[1-3].hg38.scaled.bw", full.names = T),
    in_naive = list.files(bwdir, pattern = "IN_H9_Ni.*rep[1-3].hg38.*.bw", full.names = T),
    in_naive_ezh2i = list.files(bwdir, pattern = "IN_H9_Ni-EZH2i.*rep[1-3].hg38.*.bw", full.names = T),
    in_primed = list.files(bwdir, pattern = "IN_H9_Pr_rep[1-3].hg38.*.bw", full.names = T),
    in_primed_ezh2i = list.files(bwdir, pattern = "IN_H9_Pr-EZH2i.*rep[1-3].hg38.*.bw", full.names = T)
  )

bwfiles <-
  list(
    k4 = list.files(bwdir, pattern = "H3K4m3.*pooled.hg38.scaled.*", full.names = T),
    k27 = list.files(bwdir, pattern = "H3K27m3.*pooled.hg38.scaled.*", full.names = T),
    ub = list.files(bwdir, pattern = "H2Aub.*pooled.hg38.scaled.*", full.names = T),
    input = list.files(bwdir, pattern = "IN.*pooled.hg38.*", full.names = T)
  )

bwfiles_unscaled <-
  list(
    k4 = list.files(bwdir, pattern = "H3K4m3.*pooled.hg38.unscaled.*", full.names = T),
    k27 = list.files(bwdir, pattern = "H3K27m3.*pooled.hg38.unscaled.*", full.names = T),
    ub = list.files(bwdir, pattern = "H2Aub.*pooled.hg38.unscaled.*", full.names = T),
    input = list.files(bwdir, pattern = "IN.*pooled.hg38.*", full.names = T)
  )

