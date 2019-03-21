library(openxlsx)
library(ggplot2)
library(ggrepel)


setwd("//research/antom/Projects/PeakCluster/peak_clustR")
source('functions.R')
source("visualizations.R")

X <- read.xlsx("HILIC_NEG_TEST.xlsx", sheet = 1)
P <- read.xlsx("HILIC_NEG_TEST.xlsx", sheet = 2)


conn <- find_connections(data = X, features = P, corr_thresh = 0.85, rt_window = 1/60, name_col = "Name", mz_col = "mz", rt_col = "rt")

clusters <- find_clusters(conn)

pulled <- pull_features(clusters, data= X, features = P, name_col = "Name")
cdata <- pulled$cdata
cpeaks <- pulled$cfeatures

multisheet_xlsx(dfs = pulled, filename = "HILIC_NEG_output2.xlsx")


cluster <- clusters[[2]]
length(cluster$features)

P$MPA <- sapply(X[P$Name], median, na.rm = TRUE)
plot_cluster(features = P, cluster, name_col = "Name", mz_col = "mz", rt_col = "rt")


net <- igraph_to_networkD3(cluster$graph)
saveNetwork(forceNetwork(Links = net$links, Nodes = net$nodes, NodeID = "name", Group = "name", zoom = TRUE),
            "net.html")

source("visualizations.R")
plot_features(P, cluster, name_col = "Name", mz_col = "mz", rt_col = "rt", rt_window = 1/60)

plot_heatmaps(data = X, features = P, name_col = "Name", mz_col = "mz", rt_col = "rt")

visualize_clusters(data = X, features = P, clusters = clusters, rt_window = 1/60, name_col = "Name", mz_col = "mz", rt_col = "rt",
                   file_path = "figures/")
