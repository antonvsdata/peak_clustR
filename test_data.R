library(igraph)
library(openxlsx)


setwd("//research/antom/Projects/PeakCluster/peak_clustR")
source('functions.R')

n <- 100

set.seed(38)
x1 <- runif(n, 1000, 2000)
x2 <- x1 + runif(n, -200, 200)
x3 <- x1 + runif(n, -200, 200)


x4 <- x3 + runif(n, -300, 300)

x5 <- x4 + runif(n, -200, 200)
x6 <- x5

x7 <- runif(n, 3000, 3500)
x8 <- x7 + runif(n, 100, 200)

x9 <- runif(n, 10000, 12000)

X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9)
#cor(X)

colnames(X) <- paste0("COMPOUND_", 1:ncol(X))

P <- data.frame(Name = colnames(X), mzmin = 100, mz = 100, mzmax = 100,
                rtmin = 1, rt = 1, rtmax = 1, stringsAsFactors = FALSE)
X$Group <- rep(c("A", "B"), each = 50)

conn <- find_connections(data = X, peaks = P, corr_thresh = 0.85, rt_window = 2, name_col = "Name", mz_col = "mz", rt_col = "rt")


clusters <- find_clusters(conn)

pulled <- pull_peaks(clusters, data= X, peaks = P, name_col = "Name")
cdata <- pulled$cdata
cpeaks <- pulled$cpeaks

multisheet_xlsx(dfs = pulled, filename = "output.xlsx")




