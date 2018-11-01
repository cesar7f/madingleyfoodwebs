#Adapted from package fwga
calculateGraphMetrics <-function (g) 
{
  A = as.matrix(as_adjacency_matrix(g))
  if ((sum(A == 0) + sum(A == 1)) != length(A)) {
    fw <- ifelse(A > 0, 1, 0)
  }  else {
    fw <- A
  }
  A.graph <- graph.adjacency(A, mode = "directed", diag = FALSE)
  A.graph2 <- graph.adjacency(A, mode = "undirected", diag = FALSE)
  S <- ncol(fw)
  L <- sum(fw)
  #diag(fw) <- 0
  connect <- S/(L * L)
  basal <- sum(apply(fw, 2, sum, na.rm = TRUE) == 0, na.rm = TRUE)
  int <- sum((apply(fw, 1, sum, na.rm = TRUE) != 0) & (apply(fw, 
                                                             2, sum, na.rm = TRUE) != 0))
  top <- sum(apply(fw, 1, sum, na.rm = TRUE) == 0, na.rm = TRUE)
  gen <- L/(top + int) 
  vul <- L/(basal + int)
  
  mean.pl <- average.path.length(A.graph2)
  mean.tl <- mean(TrophInd(fw)$TL, na.rm = TRUE)
  mean.oi <- mean(TrophInd(fw)$OI, na.rm = TRUE)
  
  #le.mod <- modularity(leading.eigenvector.community(A.graph2))
  le.mod <- modularity(walktrap.community(A.graph2))
  
  transglobal <-transitivity(A.graph,type = 'global')
  
  variables = c("basal","int","top","gen","vul",
                "mean.pl", "mean.tl", "mean.oi",
                "le.mod", "transglobal")
  
  return(list(S=S,L=L,connect = connect, basal = basal, int = int, 
              top = top, gen = gen, vul = vul,
              mean.pl = mean.pl, mean.tl = mean.tl, mean.oi = mean.oi, 
              le.mod = le.mod,
              transglobal = transglobal))
}

null_bernoulli_graph <- function(S, L){
  p = L/(S*S)
  return(graph.adjacency(replicate(S, rbinom(n=S,prob=p,size=1)), mode = "directed"))
}
