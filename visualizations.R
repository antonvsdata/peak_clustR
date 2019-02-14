
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

