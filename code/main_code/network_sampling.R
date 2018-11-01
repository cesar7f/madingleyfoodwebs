# This function returns a sampled version in number of links and edge property
sampling.poisson <- function(n_links_sampled, G, lambda_rate = "weight"){

    mat <- as_adjacency_matrix(G, attr=lambda_rate) 
    S   <- nrow(mat)
    L   <- sum(mat>0)
  
    if(L > n_links_sampled ) {
      
        #print(paste("S",S,"L", n_links_sampled))
        
        # required sampling effort (observation time) in order for the expected number of links to be equal to the required
        sampling_tao <- uniroot(function(x) n_links_sampled - sum(1-exp(-x*mat)), c(-2,20000000000), tol=0.0000000000001)$root
        
        # probability of each link giving the appropiate tao
        poisson_prob_link <- 1-exp(-sampling_tao * mat)
        
        # The sampled adjacency matrix that will be used for the final sample graph
        sample_matrix <- matrix(runif(S*S), ncol=S) < poisson_prob_link
        
        return(graph.adjacency(as.matrix(1*sample_matrix)))
    } else {
      #print(paste("SS",S,"LL", n_links_sampled))
      G
    }
}

# This funtion just return a sampling in wich the top links (based in edge_property)
# are sampled. Therefore it is a deterministic result.
sampling.top <- function(n_links_sampled, madingley_aggregated, edge_property = "weight"){
    
    dt_nodes <- madingley_aggregated@nodes
    dt_edges <- madingley_aggregated@edges
    n_nodes <- nrow(madingley_aggregated@nodes)
    n_edges <- nrow(madingley_aggregated@edges)
    
    dt_edges_sampled = dt_edges[order(-get(edge_property))][c(1:n_links_sampled)]
    dt_edges_sampled$time = c(1:nrow(dt_edges_sampled))
    
    graph <- graph.data.frame(dt_edges_sampled,directed = T,vertices = dt_nodes)  
    
    madingley_object <- new("Madingley", nodes = dt_nodes, edges = dt_edges_sampled, 
                            prop = madingley_aggregated@prop, G = graph)
    
    return(madingley_object)
}

# This functions uses sampling using roulette wheel selection based in the values of the passed edge_property
sampling.roulette.wheel  <- function(n_links_sampled, madingley_aggregated, edge_property = "weight"){
    
    dt_nodes <- madingley_aggregated@nodes
    dt_edges <- copy(madingley_aggregated@edges)
    n_nodes <- nrow(madingley_aggregated@nodes)
    n_edges <- nrow(madingley_aggregated@edges)
    
    dt_edges$weight <- dt_edges[,get(edge_property)]
    weights = dt_edges$weight
    weights <- weights/sum(weights)
    
    #start with an empty list of links
    sampled_links_idx <- rep(0,n_links_sampled)
    dt_edges$time = 0
    
    t <- 1
    while(n_links_sampled > 0){
        idx <- sample(c(1:n_edges), 1,
                      prob=weights)
        
        sampled_links_idx[n_links_sampled] = idx
        
        weights[idx] = 0
        weights <- weights/sum(weights)
        
        n_links_sampled <- n_links_sampled - 1
        dt_edges$time[idx] = t
        t = t+1
    }
    
    dt_edges_sampled <- dt_edges[sampled_links_idx,]
    #dt_edges_sampled$time <- rev(c(1:nrow(dt_edges_sampled)))
    
    graph <- graph.data.frame(dt_edges_sampled,directed = T,vertices = dt_nodes)  
    
    madingley_object <- new("Madingley", nodes = dt_nodes, edges = dt_edges_sampled, 
                            prop = madingley_aggregated@prop, G = graph)
    
    return(madingley_object)
}


# This function maps the edge_property to a poisson process in which the edge_proerty is used as the
# lambda (rate) poisson parameter.
sampling.poisson.effort <- function(sampling_tao=1, G, lambda_rate = "weight"){
  
  mat <- as_adjacency_matrix(G, attr=lambda_rate) 
  S   <- nrow(mat)

  # probability of each link giving the sampling effort tao (time)
  poisson_prob_link <- 1-exp(-sampling_tao * mat)
  
  poisson_prob_link <- as.matrix(1-exp(-sampling_tao * mat))
  
  # The sampled adjacency matrix that will be used for the final sample graph
  sample_matrix <- matrix(runif(S*S), ncol=S) < poisson_prob_link
  
  return(graph.adjacency(as.matrix(1*sample_matrix)))
}

# This function maps the edge_property to a poisson process in which the edge_proerty is used as the
# lambda (rate) poisson parameter.
sampling.poisson.effort.prob_matrix <- function(sampling_tao=1, G, lambda_rate = "weight"){
  
  mat <- as_adjacency_matrix(G, attr=lambda_rate) 
  S   <- nrow(mat)
  
  # probability of each link giving the sampling effort tao (time)
  poisson_prob_link <- 1-exp(-sampling_tao * mat)
  
  return(poisson_prob_link)
}
