
# Find out which features are correlated within a specified retention time window

# Parameters:
# data: data frame with the measurements and possibly some sample information
#- size samples x (features + sample info columns)
# features: data frame of feature information: at least name, retention time and mass-to-charge ratio
#- size features x feature info columns
# corr_thresh: the threshhold of correlation to use in linking signals
# rt_window: the retention time window to use in linking signals
# name_col: a string, name of the column in features that contains signal names
# mz_col: a string, name of the column in features that contains mass-to-charge ratios
# rt_col: a string, name of the column in features that contains retention times
#
# Returns:
# a data frame of pairs of signals that are linked together
#   - x & y: indexes and names of the signals
#   - cor: correlation coefficient
#   - mz_diff & rt_diff: mass and retention time difference
find_connections <- function(data, features, corr_thresh = 0.9, rt_window = 1/60,
                         name_col, mz_col, rt_col) {
  
  D <- data[features[,name_col]]
  if (ncol(D) < 2) {
    stop("Need at least 2 signals to do any clustering!")
  }
  n <- nrow(features)
  connections <- foreach::foreach(i = seq_len(n-1), .combine = rbind) %dopar% {
    if (i %% 100 == 0){
      print(i)
    }
    connections_tmp <- data.frame()
    for (j in (i+1):n){
      rt_diff <- features[j, rt_col] - features[i, rt_col]
      cor_coef <- cor(D[, i], D[, j])
      if (abs(rt_diff) < rt_window & cor_coef > corr_thresh){
        mz_diff <- features[j, mz_col] - features[i, mz_col]
        connections_tmp <- rbind(connections_tmp, data.frame(x = features[i, name_col], y = features[j, name_col],
                                                    cor = cor_coef, rt_diff = rt_diff, mz_diff = mz_diff))
      }
    }
    connections_tmp
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
find_clusters <- function(connections, features, name_col, d_thresh = 0.8){
  if(!requireNamespace("igraph", quietly = TRUE)){
    stop("The igraph package is required for this function")
  }
  
  # Construct graph from the given edges
  g <- igraph::graph_from_edgelist(as.matrix(connections[1:2]), directed = FALSE)
  g <- igraph::set.edge.attribute(graph = g, name = "weight", value = connections$cor)
  
  # Initialize list of clusters
  clusters <- list()
  k <- 1
  
  # Repeatedly extract densely connected clusters from the graph
  while(length(igraph::V(g))){
    # Connected components of the remaining graph
    comp <- igraph::decompose(g)
    n_comp <- length(comp)
    cat(paste(n_comp, "components_found\n\n"))
    
    # Only keep the densely connected part of each component (subgraph)
    for(i in seq_len(n_comp)) {
      
      if (i %% 100 == 0) {
        cat(paste("Component", i, "/", n_comp,"\n"))
      }
      subg <- comp[[i]]
      
      n_nodes <- length(igraph::V(subg))
      d <- igraph::degree(subg)
      # The limit of the degree a node needs to be kept
      d_lim <- round(d_thresh * (n_nodes-1))
      
      # 
      if(n_nodes >= 3) {
        # Remove the node with the smallest degree until all nodes in the cluster have
        # a degree above the limit
        while(any(d < d_lim)){
          idx <- which(d == min(d))
          if (length(idx) > 1) {
            edgesums <- sapply(igraph::V(subg)$name[idx], function(x) sum(igraph::E(subg)[from(igraph::V(subg)[x])]$weight))
            idx <- idx[which(edgesums == min(edgesums))[1]]
          }
          subg <- igraph::delete.vertices(subg, v = igraph::V(subg)[idx])
          d <- igraph::degree(subg)
          n_nodes <- n_nodes - 1
          d_lim <- round(d_thresh * (n_nodes-1))
        }
      }
      
      # Record the final cluster and remove the nodes from the main graph
      clusters[[k]] <- list(features = names(igraph::V(subg)),
                            graph = subg)
      k <- k + 1
      g <- igraph::delete.vertices(g, v = names(igraph::V(subg)))
    }
    
  }
  clusters
}

# Extract information of the features of all clusters
# The LC-MS data of the feature with largest median peak area is retained,
# all the features in every cluster are recorded
#
# Parameters:
# - clusters: list of cluster information, as returned by find_clusters
# - data: data frame of the original LC-MS data
# - features: data frame holding the feature information
# - name_col: a string, name of the column in features that contains signal names
#
# Returns:
# List of two items:
#   - cdata: a new data frame with the combined LC-MS data
#   - cfeatures: data frame, feature information per cluster
pull_features <- function(clusters, data, features,
                       name_col){
  
  # Median peak area
  features$MPA <- sapply(data[features[, name_col]], median, na.rm = TRUE)
  
  cfeatures <- data.frame()
  sample_cols <- setdiff(colnames(data), features[, name_col])
  cdata <- data[sample_cols]
  handled_features <- c()
  
  n_clusters <- length(clusters)
  # Retain the strongest signal (MPA) from each cluster 
  for (i in seq_along(clusters)) {
    if (i %% 100 == 0) {
      cat(paste("Cluster", i, "/", n_clusters, "\n"))
    }
    
    cluster <- clusters[[i]]
    if (length(cluster$features) > 1) {
      features_tmp <- features[features[, name_col] %in% cluster$features, ]
      
      # Find the feature with maximal MPA
      max_mpa_idx <- which(features_tmp$MPA == max(features_tmp$MPA, na.rm = TRUE))[1]
      cluster_row <- features_tmp[max_mpa_idx, ]
      # Record all the features in the cluster
      cluster_row$Features <- paste(sort(cluster$features), collapse = ";")
      cluster_row$n_features <- length(cluster$features)
      # Create cluster ID
      cluster_row$Cluster_ID <- paste0("Cluster_", cluster_row[, name_col])
      cfeatures <- rbind(cfeatures, cluster_row)
      
      # Take the LC-MS data of the largest feature
      cdata_col <- data[features_tmp[max_mpa_idx, name_col]]
      colnames(cdata_col) <- cluster_row$Cluster_ID
      cdata <- cbind(cdata, cdata_col)
      
      handled_features <- c(handled_features, cluster$features)
    }
  }
  
  # Reorganise
  cfeatures <- dplyr::arrange(cfeatures, Cluster_ID)
  cdata <- cdata[c(sample_cols, cfeatures$Cluster_ID)]
  
  # All the features that were not in the clusters are retained unchanged
  missed_features <- features[!features[, name_col] %in% handled_features, ]
  missed_features$Features <- missed_features[, name_col]
  missed_features$n_features <- 1
  missed_features$Cluster_ID <- missed_features[, name_col]
  
  
  cfeatures <- rbind(cfeatures, missed_features)
  cdata <- cbind(cdata, data[missed_features[, name_col]])
  
  # Reorder columns
  cfeatures <- dplyr::select(cfeatures, "Cluster_ID", "n_features", "Features", name_col, dplyr::everything())
  rownames(cfeatures) <- 1:nrow(cfeatures)
  
  list(cdata = cdata, cfeatures = cfeatures)
}


# Write .xlsx file with multiple sheets
#
# parameters:
# dfs:        a list of dataframes, one per sheet
# filename:   character, the name of the file
# sheetnames: charcater vector, names of the sheets
multisheet_xlsx <- function(dfs, filename, sheetnames = c("Data", "Features")){
  
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
