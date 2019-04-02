library(testthat)
library(igraph)
library(doParallel)

# Load functions
source("//research/antom/Projects/PeakCluster/peak_clustR/functions.R")

# Create toy data
n <- 100

set.seed(38)
# First cluster
x1 <- runif(n, 1000, 2000)
x2 <- x1 + runif(n, -200, 200)
x3 <- x1 + runif(n, -200, 200)
# Second cluster connected to the first
x4 <- x3 + runif(n, -300, 300)
x5 <- x4 + runif(n, -200, 200)
x6 <- x5
# Third cluster
x7 <- runif(n, 3000, 3500)
x8 <- x7 + runif(n, 100, 200)
# A loner signal
x9 <- runif(n, 10000, 12000)

X <- data.frame(x1, x2, x3, x4, x5, x6, x7, x8, x9)
colnames(X) <- paste0("F_", 1:ncol(X))

# Peak info
P <- data.frame(Name = colnames(X), mz = runif(ncol(X), 100, 400), rt = c(runif(ncol(X), 1, 1.02)),
                stringsAsFactors = FALSE)
# Add groups for samples
X$Group <- rep(c("A", "B"), each = 50)

first <- data.frame(x = c("F_1", "F_1", "F_2"), y = c("F_2", "F_3", "F_3"),
                    cor = c(cor(x1,x2), cor(x1,x3), cor(x2,x3)),
                    rt_diff = c(P$rt[2]- P$rt[1], P$rt[3]- P$rt[1], P$rt[3]- P$rt[2]),
                    mz_diff = c(P$mz[2]- P$mz[1], P$mz[3]- P$mz[1], P$mz[3]- P$mz[2]))

second <- data.frame(x = c("F_4", "F_4", "F_5"), y = c("F_5", "F_6", "F_6"),
                     cor = c(cor(x4,x5), cor(x4,x6), cor(x5,x6)),
                     rt_diff = c(P$rt[5]- P$rt[4], P$rt[6]- P$rt[4], P$rt[6]- P$rt[5]),
                     mz_diff = c(P$mz[5]- P$mz[4], P$mz[6]- P$mz[4], P$mz[6]- P$mz[5]))

third <- data.frame(x = c("F_7"), y = c("F_8"),
                    cor = cor(x7,x8),
                    rt_diff = P$rt[8]- P$rt[7],
                    mz_diff = P$mz[8]- P$mz[7])

test_that("Connections found correctly", {
  # First cluster
  
  expect_equal(find_connections(X[, 1:3], P[1:3,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               first)
  
  # Second cluster
  
  expect_equal(find_connections(X[, 4:6], P[4:6,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               second)
  
  # Third cluster
  expect_equal(find_connections(X[, 7:8], P[7:8,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               third)
  
  # Loner should produce an error
  expect_error(find_connections(X[9], P[9,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               "Need at least 2")
  
  # All together
  link <- data.frame(x = c("F_3"), y = c("F_4"),
                     cor = cor(x3,x4),
                     rt_diff = P$rt[4]- P$rt[3],
                     mz_diff = P$mz[4]- P$mz[3])
  expect_equal(find_connections(X, P, corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               rbind(first, link, second, third))
  
})



test_that("Clusters contain the correct peaks", {
  conn <- find_connections(data = X, features = P, corr_thresh = 0.85, rt_window = 0.03, name_col = "Name", mz_col = "mz", rt_col = "rt")
  
  clusters <- find_clusters(conn)
  
  expect_equal(clusters[[1]]$features, c("F_4", "F_5", "F_6"))
  expect_equal(clusters[[2]]$features, c("F_7", "F_8"))
  expect_equal(clusters[[3]]$features, c("F_1", "F_2", "F_3"))
})

graph_equal <- function(g, h){
  expect_equal(names(V(g)), names(V(h)))
  expect_equal(E(g)$weight, E(h)$weight)
}

test_that("Cluster subgraphs are formed correctly", {
  
  conn <- find_connections(data = X, features = P, corr_thresh = 0.85, rt_window = 0.03, name_col = "Name", mz_col = "mz", rt_col = "rt")
  
  # Expected graphs
  g <- graph_from_edgelist(as.matrix(conn[1:2]), directed = FALSE)
  g <- set_edge_attr(graph = g, name = "weight", value = conn[, "cor"])
  g1 <- delete.vertices(g, paste0("F_", c(1:3, 7:8)))
  g2 <- delete.vertices(g, paste0("F_", 1:6))
  g3 <- delete.vertices(g, paste0("F_", 4:8))
  
  # The clusters returned
  clusters <- find_clusters(conn)
  cg1 <- clusters[[1]]$graph
  cg2 <- clusters[[2]]$graph
  cg3 <- clusters[[3]]$graph
  
  # Test equality
  graph_equal(g1, cg1)
  graph_equal(g2, cg2)
  graph_equal(g3, cg3)
})


test_that("feature pulling works", {
  
  # Expected
  cdata <- data.frame(Group = X$Group,
                      Cluster_F_3 = X$F_3,
                      Cluster_F_5 = X$F_5,
                      Cluster_F_8 = X$F_8,
                      F_9 = X$F_9,
                      stringsAsFactors = FALSE)
  cfeatures <- data.frame(Cluster_ID = c("Cluster_F_3", "Cluster_F_5", "Cluster_F_8", "F_9"),
                       n_features = c(3, 3, 2, 1),
                       Features = c("F_1;F_2;F_3", "F_4;F_5;F_6", "F_7;F_8", "F_9"),
                       P[P$Name %in% c("F_3", "F_5", "F_8", "F_9"), ],
                       MPA = sapply(X[c("F_3", "F_5", "F_8", "F_9")], median),
                       stringsAsFactors = FALSE)
  rownames(cfeatures) <- 1:nrow(cfeatures)
  # Returned
  conn <- find_connections(data = X, features = P, corr_thresh = 0.85, rt_window = 0.03, name_col = "Name", mz_col = "mz", rt_col = "rt")
  clusters <- find_clusters(conn)
  pulled <- pull_features(clusters, data= X, features = P, name_col = "Name")
  
  expect_identical(pulled$cdata, cdata)
  expect_identical(pulled$cfeatures, cfeatures)
})
