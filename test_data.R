library(openxlsx)
library(ggplot2)


setwd("//research/antom/Projects/PeakCluster/peak_clustR")
source('functions.R')
source("visualizations.R")

X <- read.xlsx("HILIC_NEG_TEST.xlsx", sheet = 1)
P <- read.xlsx("HILIC_NEG_TEST.xlsx", sheet = 2)


conn <- find_connections(data = X, peaks = P, corr_thresh = 0.85, rt_window = 1/60, name_col = "Name", mz_col = "mz", rt_col = "rt")

clusters <- find_clusters(conn)

pulled <- pull_peaks(clusters, data= X, peaks = P, name_col = "Name")
cdata <- pulled$cdata
cpeaks <- pulled$cpeaks

multisheet_xlsx(dfs = pulled, filename = "HILIC_NEG_output2.xlsx")


cluster <- clusters[[2]]
length(cluster$peaks)

P$MPA <- sapply(X[P$Name], median, na.rm = TRUE)
plot_cluster(peaks = P, cluster, name_col = "Name", mz_col = "mz", rt_col = "rt")


net <- igraph_to_networkD3(cluster$graph)
saveNetwork(forceNetwork(Links = net$links, Nodes = net$nodes, NodeID = "name", Group = "name", zoom = TRUE),
            "net.html")

source("visualizations.R")
plot_intensities(P, cluster, name_col = "Name", mz_col = "mz", rt_col = "rt")

