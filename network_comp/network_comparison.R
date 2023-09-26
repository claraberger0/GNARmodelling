load_libraries <- function(){
  if (!require("dplyr"))
    install.packages("dplyr"); library(dplyr) 
  if (!require("readr"))
    install.packages("readr"); library(readr) 
  if (!require("igraph"))
    install.packages("igraph"); library(igraph) 
  if (!require("ggplot2"))
    install.packages("ggplot2"); library(ggplot2)
  if (!require("sf"))
    install.packages("sf"); library(sf)
  if (!require("maps"))
    install.packages("maps"); library(maps)
  if (!require("ggmap"))
    install.packages("ggmap"); library(ggmap)
  if (!require("osmdata"))
    install.packages("osmdata"); library(osmdata)
  if (!require("sfnetworks"))
    install.packages("sfnetworks"); library(sfnetworks)
  if (!require("tidygraph"))
    install.packages("tidygraph"); library(tidygraph)
}


load_libraries()

###############################################################################################################
# ************** CHANGE VALUES HERE ******************
###############################################################################################################

# set working directory
wkdr = "/Users/claraberger/Library/CloudStorage/OneDrive-Nexus365/Documents/Dissertation/GNARmodelling/network_comp/"

edges.filename = "fake_edges2.csv"



###############################################################################################################
# Build both networks
###############################################################################################################

##### GIS Network
# read in the node names and locations 
nodes.path <- paste0(wkdr,"nodes.csv")
nodes<- read_csv(nodes.path) 

# read in adjacency matrix of substations
adj.path <- paste0(wkdr,"grid_adj.csv")
adj<- read_csv(adj.path) 
adj <- adj[, -1]    # remove the first variable
# convert to matrix
adj.matrix <- as.matrix(adj)  
# convert adjacency matrix into a network
grid <- graph_from_adjacency_matrix(adj.matrix, mode="max", weighted=NULL) 


##### Auroras Network
# read in the node names and locations from Aurora 
aurnodes.path <- paste0(wkdr,"aurnodes.csv")
aurnodes<- read_csv(aurnodes.path) 
# make the node names into strings 
aurnodes$NODE_ID <- as.character(aurnodes$NODE_ID)

# read in Aurora's connections 
auradj.path <- paste0(wkdr,edges.filename)
auradj <- read_csv(auradj.path) 
# convert all numeric columns to characters
auradj <- auradj %>% mutate_all(as.character)
# convert into a network
aurgrid <- graph_from_data_frame(auradj, directed=FALSE)


##### Matches between the two within 5km 
# read in Aurora's connections 
matched.path <- paste0(wkdr,"matched_tot.csv")
matched<- read_csv(matched.path) 
# convert all numeric columns to characters
matched <- matched %>% mutate_all(as.character)

# Start writing to an output file
sink(paste0(wkdr,'analysis-output.txt'), type="output")
cat("Summary of public data network")
summary(grid)
cat("Summary of Aurora's network")
summary(aurgrid)


###############################################################################################################
# Extract various network statistics
###############################################################################################################
# Compare various statistics for the two networks
compare_statistics <- function(graph1, graph2){

  # Get the value of each statistic for both graphs to compare
  #density
  density_graph1 <- graph.density(graph1, loops=TRUE) 
  density_graph2 <- graph.density(graph2, loops=TRUE) 
  # print output to log
  print(paste("Density of", deparse(substitute(graph1)),"is", round(density_graph1,5)))
  print(paste("Density of", deparse(substitute(graph2)),"is", round(density_graph2,5)))
  
  #average path length
  avg.path1 <- average.path.length(graph1)
  avg.path2 <- average.path.length(graph2)
  # print output to log
  print(paste("average path length of", deparse(substitute(graph1)),"is", round(avg.path1,5)))
  print(paste("average path length of", deparse(substitute(graph2)),"is", round(avg.path2,5)))
  
  #clustering coefficient
  gcc1 <- transitivity(graph1, type=c("undirected"), vids=NULL)
  gcc2 <- transitivity(graph2, type=c("undirected"), vids=NULL)
  # print output to log
  print(paste("global clustering coefficient of", deparse(substitute(graph1)),"is", round(gcc1,5)))
  print(paste("global clustering coefficient of", deparse(substitute(graph2)),"is", round(gcc2,5)))

}

#run the comparison of statistics
compare_statistics(grid, aurgrid)

# Compare the betweenness of nodes and which nodes have the highest betweenness
compare_betweenness <- function(graph1, nodes1, graph2, nodes2){
  nodes1 <- nodes1[,c("NODE_ID", "LATITUDE", "LONGITUDE")]
  nodes2 <- nodes2[,c("NODE_ID", "LATITUDE", "LONGITUDE")]
  # Find the vertex betweenness for the two graphs and compare their distributions
  #compute betweenness
  betw1 <- as.data.frame(betweenness(graph1, v=V(graph1), directed = FALSE))
  colnames(betw1) <- 'betweenness' #rename value column
  betw2 <- as.data.frame(betweenness(graph2, v=V(graph2), directed = FALSE))
  colnames(betw2) <- 'betweenness'
  #add column for node names
  betw1$NODE_ID <- rownames(betw1) 
  betw2$NODE_ID <- rownames(betw2) 
  #how is the betweenness distributed for each graph
  print(paste("For", deparse(substitute(graph1))))
  print(summary(betw1))
  print(paste("For", deparse(substitute(graph2))))
  print(summary(betw2))
  
  # Plot the top betweenness nodes on a map
  #get top5 nodes with highest betweenness
  top.betw1 <- betw1[c(order(betw1[,1] , decreasing = TRUE)),][1:15,]
  top.betw1 <- inner_join(top.betw1, nodes1, by='NODE_ID') # also bring in node locations
  top.betw2 <- betw2[c(order(betw2[,1] , decreasing = TRUE)),][1:15,]
  top.betw2 <- inner_join(top.betw2, nodes2, by='NODE_ID') # also bring in node locations
  
  top.betw <-rbind(transform(top.betw1, Source="GIS nodes"), transform(top.betw2, Source="Aurora nodes"))
  #plot a map of Texas and add a point on the map for highest betweenness nodes
  texas <- get_map( getbb('texas'), source="stamen")
  
  betw_fig <- ggmap(texas) +
    geom_point(top.betw, mapping=aes(x = LONGITUDE, y = LATITUDE, color = Source, shape = Source), size = 4) +
    scale_color_manual(values = c("GIS nodes" = "red", "Aurora nodes" = "blue")) +
    labs(title = "Nodes with the highest betweenness") 
  png(filename = paste0(wkdr,"betweenness_nodes.png"))
  print(betw_fig)
  dev.off()
}

#run the comparison of betweenness
compare_betweenness(grid, nodes, aurgrid, aurnodes)


###############################################################################################################
# Look at two well connected nodes we know are matched
###############################################################################################################

one_to_one_comp <- function(graph1, nodes1, station1, graph2, nodes2, station2){
  # Find the neighbors of the node of interest
  neighbors1 <- neighbors(graph1, station1)
  neighbors2 <- neighbors(graph2, station2)
  
  # Get the coordinates for this sample in each case
  coords1 <- nodes1[nodes1$NODE_ID %in% c(station1,names(unlist(neighbors1))), ] %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
  coords2 <- nodes2[nodes2$NODE_ID %in% c(station2,names(unlist(neighbors2))), ] %>%
    st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)
  
  # Plot the two subnetworks on top of eachother and compare the edges
  png(filename = paste0(wkdr,station1,"-",station2,"_one_to_one.png"))
  induced_subgraph(graph1, c(station1,names(unlist(neighbors1)))) %>%
    as_tbl_graph() %>%
    left_join(coords1, by = c('name'='NODE_ID')) %>% 
    as_sfnetwork(directed = FALSE, edges_as_lines = TRUE) %>%
    plot(col='black', cex = 4,  
    )
  
  induced_subgraph(graph2, c(station2,names(unlist(neighbors2)))) %>%
    as_tbl_graph() %>%
    left_join(coords2, by = c('name'='NODE_ID')) %>% 
    as_sfnetwork(directed = FALSE, edges_as_lines = TRUE) %>%
    plot(col='firebrick', cex = 1.5, add=TRUE, 
    ) 
  legend("topleft", legend=c("HIFLD nodes", "Aurora nodes"),
         col=c("black", "firebrick"), lty=1)
  dev.off()
}

# We know that STATION357 and node 40700 are matched
one_to_one_comp(grid, nodes, 'STATION673', aurgrid, aurnodes, '1126')


###############################################################################################################
# Look at only the matched nodes
###############################################################################################################

create_matched_networks <- function(){
  # Number of unique GIS stations when matched based on the matched_tot file
  print(paste(length(unique(matched$NODE_ID_left)), "unique nodes from Aurora were matched with", 
              length(unique(matched$NODE_ID_right)), "unique nodes from GIS"))
  
  #make sure the matched nodes actually have edges present
  matched_nodes <- matched$NODE_ID_right[matched$NODE_ID_left %in% unlist(auradj)]
  aurmatched_nodes <- matched$NODE_ID_left[matched$NODE_ID_left %in% unlist(auradj)]
  
  # Create subgraphs with only the nodes that have been matched between the two
  grid.matched <- induced_subgraph(grid, matched_nodes) %>% simplify()
  aurgrid.matchedfull <- induced_subgraph(aurgrid, aurmatched_nodes) %>% simplify()
  
  
  # Now we want to collapse the many to one nodes where there are multiple Aurora nodes at one point
  #find the matched GIS station name "STATION123" for each Aurora nods
  V(aurgrid.matchedfull)$match <- sapply(V(aurgrid.matchedfull)$name, function(x) matched$NODE_ID_right[matched$NODE_ID_left== x])
  #collapse the nodes that are all associated with the same station
  aurgrid.matched <- contract(aurgrid.matchedfull, factor(V(aurgrid.matchedfull)$match),
                              vertex.attr.comb=list(match="first",name="ignore"))
  #give the network the station names to be able to compare the two grids 1 to 1
  V(aurgrid.matched)$name <- V(aurgrid.matched)$match
  
  return(list(grid.matched,aurgrid.matched))
}


#get the overlapping networks
matched_networks <- create_matched_networks()
grid.matched <- matched_networks[[1]]
aurgrid.matched <- matched_networks[[2]]

cat("Summary of matched public data network")
summary(grid.matched)
cat("Summary of matched Aurora's network")
summary(aurgrid.matched)

#run the comparison of statistics of the matched grids
compare_statistics(grid.matched, aurgrid.matched)

# How many edges are different between the two
print(paste("there are", length(E(aurgrid.matched %m% grid.matched)),"edges present in the Aurora grid that are not in the GIS grid"))
print(paste("there are", length(E(grid.matched %m% aurgrid.matched)),"edges present in the GIS grid that are not in the Aurora grid"))

# NOTE: now the matched grids have the same node names 
one_to_one_comp(grid.matched, nodes, 'STATION357', aurgrid.matched, nodes, 'STATION357')



cat(".............Finished")

sink()



