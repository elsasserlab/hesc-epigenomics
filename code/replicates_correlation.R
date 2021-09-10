library(wigglescout)

replicates_scatter <- function(flist) {
  # sapply(flist, function(x) sapply(flist, function(y) plot_bw_bins_scatter(x,y, genome = "hg38", bin_size = 10000)))
  p1 <- plot_bw_bins_scatter(flist[[1]], flist[[2]], genome = "hg38", bin_size = 10000, verbose = F)
  p2 <- plot_bw_bins_scatter(flist[[1]], flist[[3]], genome = "hg38", bin_size = 10000, verbose = F)
  p3 <- plot_bw_bins_scatter(flist[[2]], flist[[3]], genome = "hg38", bin_size = 10000, verbose = F)

  list(p1,p2,p3)
}

# K4-Ni
scatter_list <- replicates_scatter(bwfiles_rep$k4_naive)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k4_naive_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k4_naive_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k4_naive_03.png", width = 7, height = 7)

scatter_list <- replicates_scatter(bwfiles_rep$k4_naive_ezh2i)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k4_naive_ezh2i_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k4_naive_ezh2i_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k4_naive_ezh2i_03.png", width = 7, height = 7)

scatter_list <- replicates_scatter(bwfiles_rep$k4_primed)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k4_primed_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k4_primed_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k4_primed_03.png", width = 7, height = 7)


scatter_list <- replicates_scatter(bwfiles_rep$k4_primed_ezh2i)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k4_primed_ezh2i_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k4_primed_ezh2i_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k4_primed_ezh2i_03.png", width = 7, height = 7)



# k27-Ni
scatter_list <- replicates_scatter(bwfiles_rep$k27_naive)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k27_naive_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k27_naive_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k27_naive_03.png", width = 7, height = 7)

scatter_list <- replicates_scatter(bwfiles_rep$k27_naive_ezh2i)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k27_naive_ezh2i_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k27_naive_ezh2i_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k27_naive_ezh2i_03.png", width = 7, height = 7)

scatter_list <- replicates_scatter(bwfiles_rep$k27_primed)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k27_primed_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k27_primed_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k27_primed_03.png", width = 7, height = 7)


scatter_list <- replicates_scatter(bwfiles_rep$k27_primed_ezh2i)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_k27_primed_ezh2i_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_k27_primed_ezh2i_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_k27_primed_ezh2i_03.png", width = 7, height = 7)




# k27-Ni
scatter_list <- replicates_scatter(bwfiles_rep$ub_naive)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_h2aub_naive_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_h2aub_naive_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_h2aub_naive_03.png", width = 7, height = 7)

scatter_list <- replicates_scatter(bwfiles_rep$ub_naive_ezh2i)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_h2aub_naive_ezh2i_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_h2aub_naive_ezh2i_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_h2aub_naive_ezh2i_03.png", width = 7, height = 7)

scatter_list <- replicates_scatter(bwfiles_rep$ub_primed)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_h2aub_primed_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_h2aub_primed_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_h2aub_primed_03.png", width = 7, height = 7)


scatter_list <- replicates_scatter(bwfiles_rep$ub_primed_ezh2i)
ggsave(plot = scatter_list[[1]], filename = "./output/sup_01_correlate_h2aub_primed_ezh2i_01.png", width = 7, height = 7)
ggsave(plot = scatter_list[[2]], filename = "./output/sup_01_correlate_h2aub_primed_ezh2i_02.png", width = 7, height = 7)
ggsave(plot = scatter_list[[3]], filename = "./output/sup_01_correlate_h2aub_primed_ezh2i_03.png", width = 7, height = 7)

