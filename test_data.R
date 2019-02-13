library(igraph)
library(openxlsx)


setwd("//research/antom/Projects/PeakCluster")
source('R_PeakCluster.R')

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

X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8)
cor(X)

colnames(X) <- paste0("COMPOUND_", 1:ncol(X))

P <- data.frame(Name = colnames(X), mzmin = 100, mz = 100, mzmax = 100,
                rtmin = 1, rt = 1, rtmax = 1)
P

conn <- find_connections(data = X, peaks = P, corr_thresh = 0.85, rt_window = 2, name_col = "Name", mz_col = "mz", rt_col = "rt")


conn








