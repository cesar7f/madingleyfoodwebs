# This file includes the requiered functions to perform cohort aggregation and renormalization of properties (i.e
# juvenille mass, biomass.flow) for both edges and nodes.

# It performs the aggregation of the network, by first 
cohort.aggregation <- function(madingley_object,S,node_aggregation_function=node.kmeans.aggregation)
{
    dt_nodes <- node_aggregation_function(madingley_object, S)
    
    my_dts <- edge.aggregation(dt_nodes, madingley_object@edges)
    
    graph <- graph.data.frame(my_dts$edges,directed = T,vertices = my_dts$nodes)  
    
    madingley_object <- new("Madingley", nodes = my_dts$nodes, edges = my_dts$edges, 
                            prop = madingley_object@prop, G = graph)
    
    return(madingley_object)
}

# It performs the aggregation of nodes by using kmeans. Kmeans will assign a different cluster to each cohort,
# which later on in the function edges.aggregation will be used to tenormalize edge and node properties
node.kmeans.aggregation <- function(madingley_object,k)
{
    dt_nodes = copy(madingley_object@nodes)
    nrows = nrow(dt_nodes)
    data1 <- dt_nodes[node!=0,list(Mj = log10(Mj),
                                   Ma = log10(Ma), Fg=10000*as.integer(Fg) )]
    data1[,rand_cord:=runif(nrows)]
    data.p <- as.matrix(data1)
    
    Fgs_size = length(unique(data1[,Fg]))
    stopifnot(k>=Fgs_size) #you can not put inside the same cluster different functional groups
    
    #Will try kmeans until I found a cluster assigment that satisfy that no different functional groups are assigned the same cluster
    n_trials = 10
    repeat{
        n_trials = n_trials-1
        #print(k)
        km <- kmeans(data.p, centers=k, nstart = 100, iter.max = 300)
        success <- TRUE
        for (i in 1:k)
        {
            if(length(unique(data.p[which(km$cluster==i),3]))!=1){
                success <- FALSE
                break
            }
        }
        if(n_trials < 1 || success == TRUE) break
    }
    stopifnot(success == TRUE)
    
    dt_nodes[node!=0,cluster:=km$cluster]
    dt_nodes[node==0,cluster:=0]
    #dt_nodes[node!=0]$cluster=km$cluster
    #dt_nodes[node==0]$cluster=0

    return (dt_nodes)
}

# Another way of aggregating nodes. In this case, we aggregate the nodes using jaccard in the
# trophic link similarity
node.jaccard.aggregation <- function(madingley_object, k_groups = 60, mode = "all")
{
    dt_nodes <- copy(madingley_object@nodes)
    G <- madingley_object@G
    # we need jaccard distance
    jac_dist <- 1-similarity.jaccard(G, mode = mode, loops = FALSE) 
    # Avoid merge inside same cluster
    fun_dist = 100*(as.matrix(dist(dt_nodes[,Fg]),upper=T))
    
    dist_mat = as.dist(jac_dist + fun_dist)
    
    hc = hclust(dist_mat)
    clust = cutree(hc,k=k_groups)
    
    #dt_nodes[node!=0]$cluster=clust
    #dt_nodes[node==0]$cluster=0
    
    dt_nodes[node!=0,cluster:=clust]
    dt_nodes[node==0,cluster:=0]
    
    return (dt_nodes)
}

# This function perform the aggregation of the edges. It aggregates the properties of the edges, and renormalize them
# to different potential edge properties that may be used later for sampling.
edge.aggregation <- function(dt_nodes,dt_edges)
{
    a <- 3 #useless, just for breaking point
    # Match edges with cluster nodes
    #dt_nodes$cluster <- paste("c",dt_nodes$cluster,sep="")  
    dt_edges$prey <- dt_nodes$cluster[ match(dt_edges$prey, dt_nodes$node) ]
    dt_edges$pred   <- dt_nodes$cluster[ match(dt_edges$pred,   dt_nodes$node) ]
    #setnames(dt_edges,old=c("from","to","biomass.flow"),new=c("prey","pred","biomass"))
    
    # Aggregation of nodes and edges based in cluster assgiment
    dt_nodes <- dt_nodes[,list(Fg=mean(Fg),Mj=sum(N*Mj)/sum(N), Ma=sum(N*Ma)/sum(N), M = sum(N*M)/sum(N), N=sum(N), B = sum(N*M)),by=cluster]
    setkey(dt_nodes,cluster)
    setkey(dt_edges,prey)
    M_prey = dt_nodes[dt_edges,M]
    Nq_prey=dt_nodes[dt_edges,N]
    setkey(dt_edges,pred)
    M_pred = dt_nodes[dt_edges,M]
    Nq_pred=dt_nodes[dt_edges,N]
    
    dt_edges_temp <- data.table(prey=dt_edges$prey,pred=dt_edges$pred,biomass.flow=dt_edges$biomass.flow,
                                #Mp=dt_nodes[dt_edges$prey,M],
                                M_prey=M_prey,
                                Nq_pred=Nq_pred,
                                #Nq_prey=dt_nodes[dt_edges$prey,N])
                                Nq_prey=Nq_prey)
    dt_edges_temp <- dt_edges_temp[,list(biomass.flow=sum(biomass.flow), 
                                         individual_transfer_rate_prey = sum(biomass.flow)/mean(M_prey),
                                         individual_transfer_rate_pred = sum(biomass.flow)/mean(M_pred),
                                         individual_transfer_rate_prey_by_predator = sum(biomass.flow)/(mean(M_prey)*mean(Nq_pred)),
                                         individual_transfer_rate_pred_by_predator = sum(biomass.flow)/(mean(M_pred)*mean(Nq_pred))
    ), 
    by=list(prey,pred)]
    
    dt_edges <- dt_edges_temp[order(prey,pred)]
    dt_edges[,weight:=log10(1+biomass.flow)]# individual_transfer_rate_prey_by_predator]
    
    dt_nodes[,x:=log10(Mj)]
    dt_nodes[,y:=log10(Ma)]
    #dt_nodes[,color:=rgb(COHORT_COLORS[Fg==dt_nodes[,Fg], list(R,G,B)], maxColorValue = 255),on = c("Fg" = "Fg")]
    dt_nodes[COHORT_COLORS,color:=rgb(i.R,i.G,i.B, maxColorValue = 255), on = "Fg"] 
    #dt_nodes[,size:=log2(B)]
    
    setnames(dt_nodes,old=c("cluster"),new=c("node"))
    
    return(list(nodes=dt_nodes,edges=dt_edges))
}