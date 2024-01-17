### Authors: Andrea Tizzani (EUI) and Marlene Thomas (EUI)

#' This function computes hierarchical clustering based on dissimilarity matrix coming from
#' pair-wise distances of nodes within a network based on road connectivity and blending
#' of geographical units
#'
#' @param roads Linestrings of roads
#' @param points Points to include in the network (generally centroids of polygons)
#' @param tol This defines the tolarance level (in meters) to blend points
#' inside the network. Default is 5km.
#' @param idvar ID identifier for points (e.g. ISTAT code for italian municipalities)
#' @param weight_met Weight method for computing shortest paths: 1) distance; 2) random draws
#' @param hierarchical_method Hierarchical_method: Default is 'ward.D2'
#' @param num_clusters Number of clusters
#' This algorithm works for not directed networks. Applications to directed network would require
#' a careful editing of the algorithm below
#'
RoadHierCluster <- function(roads, points, tol = 5000, idvar, weight_met = c('1', '2'),
                          hierarchical_method = 'ward.D2', num_clusters){
  id <- idvar

  ## Make sure about the CRS
  points     <- points %>%
    st_transform(crs(roads))

  ## Construct the network and blend points
  net <- as_sfnetwork(roads, directed=F)
  net_blended <-  st_network_blend(net,points, tolerance= tol)
  ## Plotting
  ## Extract joined points
  pts_net <- net_blended %>%
    activate("nodes") %>%
    st_as_sf()
  ## Extract joined lines
  lines_net <- net_blended %>%
    activate("edges") %>%
    st_as_sf()

  ## Plot the network: coral points are blended, blue points are original
  ggplot()+
    geom_sf(data=lines_net,alpha=1, linewidth = 1.5, color = "gray82")+
    geom_sf(data=pts_net%>% filter(!is.na(id)),color="coral", size=2)+
    geom_sf(data=points,color="lightblue", size = 0.8)
    theme_light()

  ## Compute weights: default are 1) distances 2) random draws from Ber.

  net_blended <- net_blended %>%
    activate("edges") %>% # activate edges to create weights
    dplyr::mutate(weight1 = edge_length(), # length of road piece
                  weight2=runif(nrow(lines_blended))) ## plot generated weights

  ## Plot distance weights
  net_blended %>%
    st_as_sf() %>%
    select(from, contains("weight")) %>%
    ggplot()+
    geom_sf(aes(color=units::drop_units(weight1)))+
    geom_sf(data=pts_net%>% filter(!is.na(id)),color="black", size=0.5)+
    scale_color_viridis_c()+
    theme_light()+
    theme(legend.position = "right")+
    labs(color="Weight")

  ## Compute shortest distances
  orig <- with(pts_net,which(!is.na(id)))
  dest <- with(pts_net,which(!is.na(id)))

  if (weight_met == '1') {
    cost_matrix = st_network_cost(net_blended, from = orig, to= dest, weights = "weight1")
  }
  if (weight_met == '2') {
    cost_matrix = st_network_cost(net_blended, from = orig, to= dest, weights = "weight2")
  }
  ## Replace Inf with a large number
  cost_matrix[is.infinite(cost_matrix)] <- max(cost_matrix[is.finite(cost_matrix)], na.rm = TRUE)*10

  ## Compute distance matrix
  d <- as.dist(cost_matrix)

  ## Perform hierarchical clustering
  hc <- hclust(d, method = hierarchical_method)

  # Cut the tree into clusters
  clusters <- cutree(hc, k = num_clusters)

  # Assign cluster to succesfully clustered points
  pts_net <- pts_net %>% select(id)
  pts_net$cluster <- NA  # Initialize cluster column
  pts_net$cluster[orig] <- clusters  # Assign cluster IDs to nodes with idvar

  cluster_frame <- pts_net %>% st_drop_geometry() %>%
    as.data.frame()  %>% filter(!is.na(id))

  points_cluster <- points %>%
    left_join(cluster_frame, by = "id") %>% filter(!is.na(cluster))

  return(points_cluster)
}
