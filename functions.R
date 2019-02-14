
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
  
  D <- data[peaks[,name_col]]
  if (ncol(D) < 2) {
    stop("Need at least 2 signals to do any clustering!")
  }
  C <- cor(D)
  n <- nrow(peaks)
  connections <- data.frame()
  for(i in 1:(n-1)){
    if (i %% 100 == 0){
      print(i)
    }
    for (j in (i+1):n){
      rt_diff <- peaks[j, rt_col] - peaks[i, rt_col]
      if (abs(rt_diff) < rt_window & C[i,j] > corr_thresh){
        mz_diff <- peaks[j, mz_col] - peaks[i, mz_col]
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
#   - peaks: character vector of the names of the peaks included in the cluster
#   - graph: an igraph object of the cluster
find_clusters <- function(connections, peaks, name_col, d_thresh = 0.8){
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("The igraph package is required for this function")
  }
  
  # Construct graph from the given edges
  g <- igraph::graph_from_edgelist(as.matrix(connections[1:2]), directed = FALSE)
  
  plot(g)
  # Initialize list of clusters
  clusters <- list()
  k <- 1
  
  # Repeatedly extract densely connected clusters from the graph
  while(length(igraph::V(g))){
    # Connected components of the remaining graph
    comp <- igraph::decompose(g)
    
    # Only keep the densely connected part of each component (subgraph)
    for(subg in comp){
      
      n_nodes <- length(igraph::V(subg))
      d <- igraph::degree(subg)
      # The limit of the degree a node needs to be kept
      d_lim <- round(d_thresh * (n_nodes-1))
      
      # 
      if(n_nodes >= 3) {
        # Remove the node with the smallest degree until all nodes in the cluster have
        # a degree above the limit
        while(any(d < d_lim)){
          idx <- which(d == min(d))[1]
          subg <- igraph::delete.vertices(subg, v = igraph::V(subg)[idx])
          d <- igraph::degree(subg)
          n_nodes <- n_nodes - 1
          d_lim <- round(d_thresh * (n_nodes-1))
        }
      }
      
      # Record the final cluster and remove the nodes from the main graph
      clusters[[k]] <- list(peaks = names(igraph::V(subg)),
                            graph = subg)
      k <- k + 1
      g <- igraph::delete.vertices(g, v = names(igraph::V(subg)))
    }
    
  }
  clusters
}

# Extract information of the peaks of all clusters
# The LC-MS data of the peak with largest median peak area is retained,
# all the peaks in every cluster are recorded
#
# Parameters:
# - clusters: list of cluster information, as returned by find_clusters
# - data: data frame of the original LC-MS data
# - peaks: data frame holding the peak information
# - name_col: a string, name of the column in peaks that contains signal names
#
# Returns:
# List of two items:
#   - cdata: a new data frame with the combined LC-MS data
#   - cpeaks: data frame, peak information per cluster
pull_peaks <- function(clusters, data, peaks,
                       name_col){
  
  # Median peak area
  peaks$MPA <- sapply(data[peaks[, name_col]], median, na.rm = TRUE)
  
  cpeaks <- data.frame()
  sample_cols <- setdiff(colnames(data), peaks[, name_col])
  cdata <- data[sample_cols]
  handled_peaks <- c()
  
  # Retain the strongest signal (MPA) from each cluster 
  for (cluster in clusters) {
    peaks_tmp <- peaks[peaks[, name_col] %in% cluster$peaks, ]
    
    # Find the peak with maximal MPA
    max_mpa_idx <- which(peaks_tmp$MPA == max(peaks_tmp$MPA, na.rm = TRUE))[1]
    cluster_row <- peaks_tmp[max_mpa_idx, ]
    # Record all the peaks in the cluster
    cluster_row$Peaks <- paste(sort(cluster$peaks), collapse = ";")
    # Create cluster ID
    cluster_row$Cluster_ID <- paste0("Cluster_", cluster_row[, name_col])
    cpeaks <- rbind(cpeaks, cluster_row)
    
    # Take the LC-MS data of the largest peak
    cdata_col <- data[peaks_tmp[max_mpa_idx, name_col]]
    colnames(cdata_col) <- cluster_row$Cluster_ID
    cdata <- cbind(cdata, cdata_col)
    
    handled_peaks <- c(handled_peaks, cluster$peaks)
  }
  
  # All the peaks that were not in the clusters are retained unchanged
  missed_peaks <- peaks[!peaks[, name_col] %in% handled_peaks, ]
  missed_peaks$Cluster_ID <- missed_peaks[, name_col]
  missed_peaks$Peaks <- missed_peaks[, name_col]
  
  cpeaks <- rbind(missed_peaks, cpeaks)
  cdata <- cbind(cdata, data[missed_peaks[, name_col]])
  
  # Reorder columns
  cpeaks <- dplyr::select(cpeaks, "Cluster_ID", "Peaks", name_col, dplyr::everything())
  
  list(cdata = cdata, cpeaks = cpeaks)
}


# Write .xlsx file with multiple sheets
#
# parameters:
# dfs:        a list of dataframes, one per sheet
# filename:   character, the name of the file
# sheetnames: charcater vector, names of the sheets
multisheet_xlsx <- function(dfs, filename, sheetnames = c("Data", "Peaks")){
  
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
