library(testthat)

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

test_that("Connections found correctly", {
  # First cluster
  first <- data.frame(x = c("F_1", "F_1", "F_2"), y = c("F_2", "F_3", "F_3"),
                      cor = c(cor(x1,x2), cor(x1,x3), cor(x2,x3)),
                      rt_diff = c(P$rt[2]- P$rt[1], P$rt[3]- P$rt[1], P$rt[3]- P$rt[2]),
                      mz_diff = c(P$mz[2]- P$mz[1], P$mz[3]- P$mz[1], P$mz[3]- P$mz[2]))
  expect_equal(find_connections(X[, 1:3], P[1:3,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               first)
  
  # Second cluster
  second <- data.frame(x = c("F_4", "F_4", "F_5"), y = c("F_5", "F_6", "F_6"),
                      cor = c(cor(x4,x5), cor(x4,x6), cor(x5,x6)),
                      rt_diff = c(P$rt[5]- P$rt[4], P$rt[6]- P$rt[4], P$rt[6]- P$rt[5]),
                      mz_diff = c(P$mz[5]- P$mz[4], P$mz[6]- P$mz[4], P$mz[6]- P$mz[5]))
  expect_equal(find_connections(X[, 4:6], P[4:6,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               second)
  
  third <- data.frame(x = c("F_7"), y = c("F_8"),
                       cor = cor(x7,x8),
                       rt_diff = P$rt[8]- P$rt[7],
                       mz_diff = P$mz[8]- P$mz[7])
  expect_equal(find_connections(X[, 7:8], P[7:8,], corr_thresh = 0.85,
                                name_col = "Name", mz_col = "mz", rt_col = "rt"),
               third)
  
  # Loner
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

