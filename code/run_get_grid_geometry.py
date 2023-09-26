import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiLineString
from os import listdir
from os.path import isfile, join
import contextily as cx

########################################################################################
#
#  Constructing network topology of the ERCOT grid using available GIS transmission
#  line data
#  Clara Berger July 2023
#
########################################################################################


from get_grid_geometry import *


#### Set key directories and paths
wkdr = "/Users/claraberger/Library/CloudStorage/OneDrive-Nexus365/Documents/Dissertation/"

# Use shapefile of US transmission lines here:
# https://www.arcgis.com/home/item.html?id=d4090758322c4d32a4cd002ffaa0aa12&view=list&sortOrder=desc&sortField=defaultFSOrder
grid_path = wkdr + "grid_shape/us_grid.shp"

#shapefile of US states from US Census TIGER API https://www.census.gov/geographies/mapping-files/time-series/
# geo/cartographic-boundary.html
states_path = wkdr + "cb_2022_us_state_5m/cb_2022_us_state_5m.shp"
# set desired epsg
epsg_proj = 'epsg:3857' # CRS from shapefile metadata
epsg = 'epsg:4326'   # desired final CRS

########################################################################################
#  Get grid of Texas and extract nodes to create an adjacency matrix
########################################################################################

#### Get grid of Texas
tx_grid = get_texas_grid(grid_path=grid_path, states_path=states_path, epsg_grid=epsg_proj, epsg=epsg)

tx_grid['geometry'].to_file(wkdr+'writeup/figures/tx_grid.shp') 
# list of columns -- subset the ones we are interested in for ease of reading tables
#   'FID', 'OBJECTID', 'ID', 'TYPE', 'STATUS', 'NAICS_CODE', 'NAICS_DESC',
#   'SOURCE', 'SOURCEDATE', 'VAL_METHOD', 'VAL_DATE', 'OWNER', 'VOLTAGE',
#   'VOLT_CLASS', 'INFERRED', 'SUB_1', 'SUB_2', 'SHAPE__Len', 'Shape__L_1',
#   'GlobalID', 'geometry', 'index_right', 'STATEFP', 'STATENS', 'AFFGEOID',
#   'GEOID', 'STUSPS', 'NAME', 'LSAD', 'ALAND', 'AWATER',
tx_grid = tx_grid[['OWNER', 'VOLTAGE', 'VOLT_CLASS', 'SUB_1', 'SUB_2','SHAPE__Len', 'Shape__L_1', 'geometry',]]
# explode MultiLineStrings to LineStrings 
tx_grid = tx_grid.explode(index_parts=True)


#### Get Nodes
# get list of unique nodes and their locations
# also get a dataframe of which nodes are connected to which
nodes, nodes_lines = get_lat_long(df=tx_grid, epsg=epsg, start_node='SUB_1', end_node='SUB_2')

# make an adjacency matrix from the node connections
adj_matrix = make_adj_matrix(nodes_lines, 'start_node', 'end_node', nodes, 'NODE_ID')


#### Save tables as csv
adj_path = wkdr+"tables/grid_adj.csv" 
adj_matrix.to_csv(adj_path)
# stations, voltages, start and end locations
stations_path = wkdr+"tables/tx_stations.csv" 
nodes_lines.to_csv(stations_path)
# list of nodes with locations
nodes_path = wkdr+"tables/nodes.csv" 
nodes.to_csv(nodes_path)


########################################################################################
#  Merge nodes from GIS with those from Aurora
########################################################################################

aurora_dir = wkdr+"Aurora_data/" # set directory path

# Node locations from Aurora
aurora_nodes_path = aurora_dir +"node_line_data/node_and_location.csv"
aurnode_ids = pd.read_csv(aurora_nodes_path, header=[1])
# put into a GeoDataFrame
aurnodes = gpd.GeoDataFrame(aurnode_ids, geometry=gpd.points_from_xy(aurnode_ids.LONGITUDE, aurnode_ids.LATITUDE, crs=epsg))
aurnodes.to_csv(wkdr+"/tables/aurnodes.csv") # save GDF
# merge all years data from the csv files into one table for the congestion data, the demand data, 
# and the generation data
variables = ("cong","demand","gen") # which values we're interested in
hist_dict = read_aurora_files(dir=aurora_dir, variables=variables)
# and set each table to its own variable
hist_cong = hist_dict['cong']
hist_demand = hist_dict['demand']
hist_gen = hist_dict['gen']

# which nodes have congestion, demand, and generation data
congnodes_ids = list(hist_cong.columns[1:])      
demnodes_ids = list(hist_demand.columns[2:])
gennodes_ids = list(hist_gen.columns[1:])

# indicate whether each node is has data on demand, price, or generation or none
aurnodes['congestion'] = aurnodes['NODE_ID'].apply(lambda x: 1 if str(x) in congnodes_ids else 0)
aurnodes['demand'] = aurnodes['NODE_ID'].apply(lambda x: 1 if str(x) in demnodes_ids else 0)
aurnodes['gen'] = aurnodes['NODE_ID'].apply(lambda x: 1 if str(x) in gennodes_ids else 0)

#### Merge the nodes with congestion costs with our GIS nodes
congnodes = aurnodes[aurnodes.congestion == 1]
congnodes = congnodes.to_crs(epsg_proj)
congnodes.crs

# perform a nearest spatial join between the GIS nodes and Auroras nodes
# then we can do a nearest spatial join with a max distance of 5 kilometers to find another node
cong_merged = find_node_matches(df_l=congnodes, df_r=nodes, epsg_proj=epsg_proj, threshold=5000, lsuffix='aur', rsuffix='gis')
print("number of congestion nodes:", len(congnodes), "number of matches from Aurora:", 
    len(cong_merged.NODE_ID_aur.unique()), "number of matches from GIS nodes:", len(cong_merged.NODE_ID_gis.unique()))
# 193 / 527 can be matched within 5 meters of our nodes
# 281 / 527 can be matched within 50 meters
# 454 / 527 can be matched within 5000 meters

# do the same spatial join for all the nodes, not just congestion nodes
merged = find_node_matches(df_l=aurnodes, df_r=nodes, epsg_proj=epsg_proj, threshold=5000, lsuffix='aur', rsuffix='gis')
len(merged.NODE_ID_aur.unique())

# save tables in home directory
merged.to_csv(wkdr+"matched_tot.csv")
not_matched = lmpnodes[~lmpnodes.NODE_ID.isin(lmp_merged.NODE_ID)]
not_matched.to_csv(wkdr+"notmatched.csv")

