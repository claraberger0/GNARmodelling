library("igraph")
library("GNAR")
library("dplyr")
library("readr")
library('geosphere')
library('ggplot2')
library("sf")
library("maps")
library("ggmap")
library("osmdata")
library("sf")
library("sfnetworks")
library("tidygraph")


# set working directory
wkdr = "/Users/claraberger/Library/CloudStorage/OneDrive-Nexus365/Documents/Dissertation/tables/"
fig.path = "/Users/claraberger/Library/CloudStorage/OneDrive-Nexus365/Documents/Dissertation/figures/"

###############################################################################################################
# Create network from grid adjacency matrix
###############################################################################################################

# read in adjacency matrix of substations
adj.path <- paste0(wkdr,"grid_adj.csv")
adj<- read_csv(adj.path) 
adj <- adj[, -1]    ## remove the first variable
# convert to matrix
adjm <- as.matrix(adj)  


# make symmetric
# element-wise maximum between 'adj.matrix' and its transpose
make_symmetric <- function(matrix){
  sym_mat <- pmax(matrix, t(matrix))
  return(sym_mat)
}
adj.matrix <- make_symmetric(adjm)

adj.matrix[adj.matrix[,'STATION2'] == 1 , 'STATION2']

# read in the node names and locations 
nodes.path <- paste0(wkdr,"nodes.csv")
nodes<- read_csv(nodes.path) 

# read in the node names and locations from Aurora 
aurnodes.path <- paste0(wkdr,"aurnodes.csv")
aurnodes<- read_csv(aurnodes.path) 


# convert adjacency matrix into a network
grid <- graph.adjacency(adj.matrix, mode="undirected", weighted=NULL) 
summary(grid)
con.comps <-components(grid)

# Visualise network on map of Texas
# make the color and size of each node correspond to the node degree
V(grid)$color <- scales::dscale(degree(grid) %>% cut(c(0,3,6,9,12,17)), sequential_pal)
sizes <- c(0.5,1,1.5,2,2.5,3)
coords <- as.matrix(nodes[c('LONGITUDE','LATITUDE')]) # coordinates of each node
plot(
    grid,
    vertex.size = degree(grid)/2 + 0.5,  # set vertex size by the degree
    #vertex.color = colors,
    edge.arrow.size = .25,
    vertex.label = NA,
    layout = coords
  )
legend("topleft", title="Node Degree", legend = levels(degree(grid) %>% cut(c(0,3,6,9,12,17))), 
       pch = 21, pt.bg=sequential_pal(5), col='black', cex = 0.5, pt.cex = sizes)



# Visualise the different connected components of the network
plot(
  grid,
  vertex.size = 3,  # set vertex size by the degree
  vertex.color = con.comps$membership,
  edge.arrow.size = .25,
  vertex.label = NA,
  layout = coords
)



###############################################################################################################
# Extract various network statistics
###############################################################################################################

# get the degree sequence of the nodes
degs <- degree(grid)
head(sort(degs, decreasing = TRUE))
head(sort(degs, decreasing = FALSE))
# STATION357  STATION287 STATION1829 STATION3328  STATION428  STATION476 
#      21          18          17          17          16          16 


#graph density
density <- graph.density(grid, loops=TRUE) # self loops are possible



#average path length
avg.path <- average.path.length(grid, directed=TRUE, unconnected=TRUE)

#clustering coefficient
gcc <- transitivity(grid, type=c("undirected"), vids=NULL)

#diameter
diam <- diameter(grid, directed = TRUE, unconnected = TRUE, weights = NULL)


#vertex betweenness
betw <- as.data.frame(betweenness(grid, v=V(grid), directed = TRUE))
colnames(betw) <- 'betweenness'
betw$NODE_ID <- rownames(betw) # add column for the node names
summary(betw)
# get top5 nodes with highest betweenness
top.betw <- head(betw[c(order(betw[,1] , decreasing = TRUE)),])
top.betw <- inner_join(top.betw, nodes, by='NODE_ID') # also bring in node locations

#### take a look at the high betweenness vertex
# Plot a map of Texas and add a point on the map for highest betweenness nodes
texas <- get_map( getbb('texas'), source="stamen")

betw_fig <- ggmap(texas) +
  geom_point(top.betw, mapping=aes(x = longitude, y = latitude, color = betweenness), size = 4) +
  scale_colour_gradient(low = "orange", high = "red") +
  ggtitle("6 Nodes with the highest betweenness")
png(paste0(fig.path,"betweenness_nodes.png"))
print(betw_fig)
dev.off()






