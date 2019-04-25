
# A helper function for plotting, scales the values in X
# between new min and max
rescale <- function(x, new_min, new_max) {
  (new_max - new_min) * (x - min(x)) / (max(x) - min(x)) + new_min
}


# UNFINISHED!!
plot_graph <- function(features, cluster, name_col, mz_col, rt_col) {
  
  # Ensure a correct order of the rows
  g <- cluster$graph
  vertices <- data.frame(Name = igraph::V(g)$name, stringsAsFactors = FALSE)
  colnames(vertices) <- name_col
  features_tmp <- dplyr::left_join(vertices, features, by = name_col)
  
  # Scaling of MPA to correct size
  # Square root to scale area, not radius
  size <- sqrt(rescale(features_tmp$MPA, new_min = 15^2, new_max = 40^2))
  
  if (length(igraph::E(g)) <= 20) {
    edge_labels <- as.character(round(igraph::E(g)$weight, digits = 2))
  } else {
    edge_labels <- NULL
  }
  
  
  g$palette <- RColorBrewer::brewer.pal(n = max(3, length(unique(igraph::degree(g)))), name = "Blues")
  plot(g, vertex.label = as.character(features_tmp[, mz_col]), vertex.size = size, vertex.label.dist = 0.1*size, vertex.label.degree = -pi/2,
       vertex.color = igraph::degree(g), edge.label = edge_labels)
}

plot_features <- function(features, cluster, name_col, mz_col, rt_col, rt_window) {
  
  features_tmp <- features[features[, name_col] %in% cluster$features, ]
  
  p1 <- ggplot(features_tmp, aes_string(mz_col, "MPA")) +
    geom_point(size = 3, color = "steelblue4") +
    geom_segment(aes_string(x = mz_col, yend = "MPA", xend = mz_col), y = 0, color = "steelblue4") +
    geom_label_repel(aes_string(label = mz_col), color = "steelblue4") +
    theme_minimal() +
    xlim(0.9*min(features_tmp[, mz_col], na.rm = TRUE),
         1.15*max(features_tmp[, mz_col], na.rm = FALSE)) +
    expand_limits(y=0) +
    labs(x = "Mass-to-charge ratio", y = "Median Peak Area")
  
  features_tmp$rtmin <- features_tmp[, rt_col] - rt_window
  features_tmp$rtmax <- features_tmp[, rt_col] + rt_window
  
  p2 <- ggplot(features_tmp, aes_string(rt_col, mz_col)) +
    geom_point(size = 3, color = "steelblue4") +
    geom_errorbarh(aes(xmin = rtmin, xmax = rtmax), color = "steelblue4") +
    theme_minimal() +
    labs(x = "Retention time", y = "Mass-to-charge ratio", title = "Retention time & tolerance")
  
  plot(p1)
  plot(p2)
}

plot_heatmaps <- function(data, features, cluster, name_col, mz_col, rt_col) {
  
  D <- data[, cluster$features]
  
  features_tmp <- features[features[, name_col] %in% cluster$features, ]
  
  n <- length(cluster$features)
  mz_rt <- data.frame()
  
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      mz_rt <- rbind(mz_rt, data.frame(x = features_tmp[i, name_col],
                                       y = features_tmp[j, name_col],
                                       mz_diff = features_tmp[i, mz_col] - features_tmp[j, mz_col],
                                       rt_diff = features_tmp[i, rt_col] - features_tmp[j, rt_col],
                                       stringsAsFactors = FALSE))
    }
  }
  
  mz_ord <- features_tmp[, name_col][order(features_tmp[, mz_col])]
  mz_rt$x <- factor(mz_rt$x, levels = mz_ord)
  mz_rt$y <- factor(mz_rt$y, levels = rev(mz_ord))
  
  p1 <- ggplot(mz_rt, aes(x = x, y = y, fill = mz_diff)) +
    geom_tile(color = "grey80") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
    scale_fill_gradient2()
  
  if (nrow(mz_rt) <= 10) {
    p <- p + geom_text(aes(label = round(mz_diff, digits = 2)))
  }
  
  # rt_ord <- features_tmp[, name_col][order(features_tmp[, rt_col])]
  # mz_rt$x <- factor(mz_rt$x, levels = rt_ord)
  # mz_rt$y <- factor(mz_rt$y, levels = rev(rt_ord))
  # 
  # ggplot(mz_rt, aes(x = x, y = y, fill = rt_diff)) +
  #   geom_tile(color = "grey80") +
  #   geom_text(aes(label = round(rt_diff, digits = 5))) +
  #   theme_minimal() +
  #   theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  #   scale_fill_gradient2()
  
  plot(p1)
}

visualize_clusters <- function(data, features, clusters, min_size, rt_window, name_col, mz_col, rt_col, file_path) {
  
  for (i in seq_along(clusters)) {
    if (i %% 100 == 0) {
      print(paste(i, "/", length(clusters)))
    }
    cluster <- clusters[[i]]
    if (length(cluster$features) >= min_size) {
      features_tmp <- features[features[, name_col] %in% cluster$features, ]
      cluster_id <- features_tmp[, name_col][which(features_tmp$MPA == max(features_tmp$MPA, na.rm = TRUE))[1]]
      
      pdf(paste0(file_path, "Cluster_", i, "_", cluster_id, ".pdf"), width = 10, height = 10)
      plot_heatmaps(data, features, cluster, name_col, mz_col, rt_col)
      plot_features(features, cluster, name_col, mz_col, rt_col, rt_window)
      plot_graph(features, cluster, name_col, mz_col, rt_col)
      dev.off()
    }
  }
  
  
}