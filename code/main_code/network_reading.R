# This is the main functions to deal with the original network data. The network data will be encapsulated in a class
# called Madingley. This class will is used during all the analysis, and it contains a data.table for both nodes
# and edges, a graph object, and a property list


COHORT_COLORS = data.table(Fg = c(seq(10,18),0), 
                           R=c(0,255,0,180,255,180,0,100,0,0), G = c(100,0,0,255,180,180,255,0,0,0), B = c(0,0,100,180,180,255,0,0,255,0))
setkey(COHORT_COLORS,Fg)


setOldClass("igraph")
setClass("Madingley", representation(nodes = "data.table", edges = "data.table", G = "igraph", prop = "list"))

read.madingley.data <- function(file_nodes, file_edges, include_source_node = FALSE, timeStep = 1200)
{
  #read edges file
  dt_edges <- fread(input=file_edges, header = T, sep = '\t')
  dt_edges <- dt_edges[dt_edges$time_step == timeStep]
  
  #remove sources (or not)
  if(include_source_node == FALSE)
  { 
    dt_edges <- dt_edges[!(dt_edges$Prey_ID %in% c("S1","S2"))]
  } else{
    dt_edges[dt_edges$Prey_ID %in% c("S1","S2"),Prey_ID:="0"]
  }
  
  #convert Prey_ID to integer
  dt_edges <- transform(dt_edges, Prey_ID = as.integer(Prey_ID))
  #no need for time if you do not plan to use it/also reorder for analysis in igraph
  dt_edges <- subset(dt_edges, select=c("Prey_ID","Pred_ID","Biomass_Assimilated"))
  setnames(dt_edges, c("prey","pred","biomass.flow"))
  dt_edges[]
  dt_edges[,log.biomass.flow:=log10(1+biomass.flow)]
  
  #read nodes file
  dt_nodes <- fread(input=file_nodes, header = T, sep = '\t')
  
  #only use one time
  dt_nodes <- dt_nodes[dt_nodes$TimeStep==timeStep]
  
  #remove sources (or not)
  if(include_source_node == FALSE)
  { 
    dt_nodes <- dt_nodes[!(dt_nodes$FunctionalGroup %in% c("S1","S2"))]
  } else{
    dt_nodes <- dt_nodes[!(dt_nodes$FunctionalGroup == "S2")]
    dt_nodes[dt_nodes$FunctionalGroup %in% c("S1","S2"),ID:=0]
    dt_nodes[dt_nodes$FunctionalGroup %in% c("S1","S2"),FunctionalGroup:="0"]
    dt_nodes[ID==0, AdultMass := min(dt_nodes[ID!=0,AdultMass],na.rm = T)/15]
    dt_nodes[ID==0, JuvenileMass := min(dt_nodes[ID!=0,JuvenileMass],na.rm = T)/15]
    dt_nodes[ID==0, AdultMass := min(dt_nodes[ID!=0,AdultMass],na.rm = T)/15]
    dt_nodes[ID==0, CohortAbundance := 1]
  }
  
  #convert to appropiete format
  dt_nodes <- transform(dt_nodes, FunctionalGroup = as.integer(FunctionalGroup))
  dt_nodes <- transform(dt_nodes, ID = as.character(ID))
  
  properties <- list()
  properties <- c(properties, title = "Madingley Model", habitat = "Virtual World", M.units = "kg", N.units = 'km^2')
  properties$lat <- dt_nodes[1,Latitude]
  properties$long <- dt_nodes[1,Longitude]
  
  #remove columns that are nor neccesary for the analysis/also move ID to first column for analysis in igraph
  dt_nodes <- subset(dt_nodes, select=c("ID","FunctionalGroup","JuvenileMass","AdultMass","IndividualBodyMass","CohortAbundance"))
  setnames(dt_nodes,c("node","Fg","Mj","Ma","M","N" ))
  
  dt_nodes[, B:=N*M] # Total biomass
  #dt_nodes[,color:=rgb(COHORT_COLORS[Fg %in% as.numeric(dt_nodes[,Fg]), list(R,G,B)], maxColorValue = 255)] # color
  dt_nodes[COHORT_COLORS,color:=rgb(i.R,i.G,i.B, maxColorValue = 255), on = "Fg"] # color
  
  dt_nodes[,x:=log10(Ma)]
  dt_nodes[,y:=log10(Mj)]
  
  graph <- graph.data.frame(dt_edges,directed = T,vertices = dt_nodes) 
  
  madingley_object <- new("Madingley", nodes = dt_nodes, edges = dt_edges,
                          prop = properties, G = graph)
  
  
  return (madingley_object)
}

get.graph <- function(dt_nodes, dt_edges)
{
  g <- graph.data.frame(dt_edges,directed = T,vertices = dt_nodes) 
  return(g)
}


madingley.to.cheddar <- function(madingley_object)
{
  prop <- madingley_object@prop
  #prop$title <- list(title = "Madingley Model", habitat = "Virtual World", M.units = "kg", N.units = 'km^2')
  df_nodes <- as.data.frame(madingley_object@nodes)
  df_edges <- as.data.frame(madingley_object@edges)
  
  df_nodes$node <- paste("n",df_nodes$node,sep="")
  setnames(df_edges,old=c("prey","pred"),new=c("consumer","resource"))
  df_edges$resource <- paste("n",df_edges$resource,sep="")
  df_edges$consumer <- paste("n",df_edges$consumer,sep="")
  ched <- Community(df_nodes, properties=prop, trophic.links=df_edges)
}



# ___________________________________________________________________________________________
# -------------------------------------------------------------------------------------------
# CENTROIDS INITIALIZATION
# ------------------------
# Centroids Initialization is according to functional groups proportion
# ratio <- rbind((k*table(Fgs)/length(Fgs)),rep(1,Fgs_size))
# assigned_clusters = round(apply(ratio,2,max))
# 
# # If the number of initial centroids is different than the k desired number of clusters....
# n_sum_clust <- sum(assigned_clusters)
# n_max_clust <- max(assigned_clusters)
# sel_fgs <- sample(which(assigned_clusters == n_max_clust), abs(k-n_sum_clust), replace = T)
# if(n_sum_clust < k)
# { 
#   for(i in sel_fgs) assigned_clusters[i] = assigned_clusters[i] + 1
# } else if(n_sum_clust > k){
#   for(i in sel_fgs) assigned_clusters[i] = assigned_clusters[i] - 1
# }
# 
# func_temp <- function(dt_edges,dt_nodes)
# {
#   dt_edges_temp <- data.table(prey=dt_edges$prey,pred=dt_edges$pred,biomass=dt_edges$biomass,Mp=dt_nodes1[dt_edges$pred,M],Nq=dt_nodes1[dt_edges$prey,N])
#   dt_edges_temp <- dt_edges_temp[,list(biomass=sum(biomass), lambda_indvidual_rate = sum(biomass)/(mean(Mp)*mean(Nq))), by=list(prey,pred)]
#   return(dt_edges_temp[order(prey,pred)])
# }
# 
# func_temp_ver2 <- function(dt_edges,dt_nodes1)
# {
#   dt_nodes1 <- dt_nodes1[order(cluster)]
#   dt_nodes2 <- data.table(sqldf("select cluster, sum(N) as N, sum(N*Mj)/sum(N) as Mj, sum(N*Ma)/sum(N) as Ma, sum(N*M)/sum(N) as M, sum(N*M) as B from dt_nodes group by cluster order by cluster"))
#   
#   dt_e<-as.data.table(sqldf("SELECT dt_edges.prey,dt_edges.pred, sum(dt_edges.biomass), sum(dn1.M), sum(dn2.N) 
#          FROM dt_edges 
#          INNER JOIN dt_nodes1 dn1 ON dn1.cluster = dt_edges.prey
#          INNER JOIN dt_nodes1 dn2 ON dn2.cluster = dt_edges.pred
#          GROUP BY dt_edges.prey,dt_edges.pred
#          ORDER by dt_edges.prey,dt_edges.pred
#   ")) 
#   return(as.data.table(dt_e))
# }
# 
# ss <- func_temp(dt_edges,dt_nodes1)
# rr <- func_temp_ver2(dt_edges,dt_nodes1)
# 
# 
# microbenchmark(ss <- func_temp(dt_edges,dt_nodes1),rr <- func_temp_ver2(dt_edges,dt_nodes1))