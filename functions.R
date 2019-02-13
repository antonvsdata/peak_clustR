








# Find out which peaks are correlated within a specified retention time window

# Parameters:
# data: data frame with the measurements and possibly some sample information
#- size samples x (peaks + sample info columns)
# peaks: data frame of peak information: at least name, retention time and mass-to-charge ratio
#- size peaks x peak info columns
# corr_thresh: the threshhold of correlation to use in linking signals
# rt_window: the retention time window to use in linking signals
# name_col: a string, name of the column in peaks that contains signal names
# mz_col: a string, name of the column in peaks that contains mass-to-charge ratios
# rt_col: a string, name of the column in peaks that contains retention times
#
# Returns:
# a data frame of pairs of signals that are linked together
#   - x & y: indexes and names of the signals
#   - cor: correlation coefficient
#   - mz_diff & rt_diff: mass and retention time difference
find_connections <- function(data, peaks, corr_thresh = 0.9, rt_window = 1/60,
                         name_col, mz_col, rt_col) {
  
  D <- data[peaks[, name_col]]
  C <- cor(D)
  n <- nrow(peaks)
  connections <- data.frame()
  for(i in 1:(n-1)){
    if (i %% 100 == 0){
      print(i)
    }
    for (j in (i+1):n){
      rt_diff <- abs(peaks[i, rt_col] - peaks[j, rt_col])
      mz_diff <- abs(peaks[i, mz_col] - peaks[j, mz_col])
      if (rt_diff < rt_window & C[i,j] > corr_thresh){
        connections = rbind(connections, data.frame(x = peaks[i, name_col], y = peaks[j, name_col],
                                                    cor = C[i,j], rt_diff = rt_diff, mz_diff = mz_diff))
      }
    }
    
  }
  connections
}


# Extract the densely connected clusters
#
# Parameters:
# connections:  data frame of pairs of signals that are linked together,
#               output of find_connections
# d_thresh:     numeric, the minimum degree required for each signal in a cluster
#               expressed as a percentage of the maximum degree in the cluster
#
# Returns:
# a list of clusters, each a list of:
#   - features: character vector of the names of the features included in the cluster
#   - graph: an igraph object of the cluster
find_clusters <- function(connections, d_thresh = 0.8){
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("The igraph package is required for this function")
  }
  
  # Construct graph from the given edges
  g <- igraph::graph_from_edgelist(as.matrix(conn[1:2]), directed = FALSE)
  
  # Initialize list of clusters
  clusters <- list()
  k <- 1
  
  # Repeatedly extract densely connected clusters from the graph
  while(length(igraph::V(g))){
    # Connected components of the remaining graph
    comp <- igraph::decompose(g)
    
    # Only keep the densely connected part of each component (subgraph)
    for(subg in comp){
      
      n_nodes <- length(V(subg))
      d <- igraph::degree(subg)
      # The limit of the degree a node needs to be kept
      d_lim <- round(d_thresh * (n_nodes-1))
      
      # 
      if(n_nodes >= 3) {
        # Remove the node with the smallest degree until all nodes in the cluster have
        # a degree above the limit
        while(any(d < d_lim)){
          idx <- which(d == min(d))[1]
          subg <- igraph::delete.vertices(subg, v = V(subg)[idx])
          d <- igraph::degree(subg)
          n_nodes <- n_nodes - 1
          d_lim <- round(d_thresh * (n_nodes-1))
        }
      }
      
      # Record the final cluster and remove the nodes from the main graph
      clusters[[k]] <- list(features = names(V(subg)),
                            graph = subg)
      k <- k + 1
      g <- igraph::delete.vertices(g, v = names(igraph::V(subg)))
    }
    
  }
  clusters
}


pull_peaks <- function(clusters, data, peaks,
                       name_col, mz_col, rt_col){
  
}


# Write .xlsx file with multiple sheets
#
# parameters:
# dfs:        a list of dataframes, one per sheet
# filename:   character, the name of the file
# sheetnames: charcater vector, names of the sheets
multisheet_xlsx <- function(dfs, filename, sheetnames){
  
  if (length(dfs) != length(sheetnames)) {
    stop("The number of dataframes and the number of sheet names does not match!")
  }
  
  wb <- openxlsx::createWorkbook()
  
  for (i in 1:length(sheetnames)) {
    openxlsx::addWorksheet(wb, sheetName = sheetnames[i])
    openxlsx::writeData(wb, sheet = i, x = dfs[[i]])
  }
  
  openxlsx::saveWorkbook(wb, filename, overwrite = TRUE)
}

save_clusters <- function(data, peaks, clusters){
  "hi"
}


# A helper function for plotting, scales the values in X
# between new min and max
rescale <- function(x, new_min, new_max){
  (new_max - new_min) * (x - min(x)) / (max(x) - min(x)) + new_min
}


# UNFINISHED!!
plot_cluster <- function(peaks, cluster, name_col, mz_col, rt_col){
  
  # Ensure a correct order of the rows
  g <- cluster$graph
  vertices <- data.frame(Name = igraph::V(g)$name, stringsAsFactors = FALSE)
  colnames(vertices) <- name_col
  P <- dplyr::left_join(vertices, peaks)
  
  # Scaling of MPA to correct size
  # Square root to scale area, not radius
  size <- sqrt(rescale(P$MPA, new_min = 15^2, new_max = 40^2))
  
  V(g)$label <- P$mz
  V(g)$size <- size
  V(g)$color <- c(1,2,3)
  
  g$palette <- brewer.pal(n = max(3, length(unique(degree(g)))), name = "Blues")
  pdf("lol.pdf", width = 6, height = 6)
  plot(g, vertex.label.dist = 0.105*V(g)$size, vertex.label.degree = -pi/2)
  dev.off()
}

