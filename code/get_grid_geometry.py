import geopandas as gpd
import pandas as pd
import matplotlib.pyplot as plt
from shapely.geometry import Point, LineString, MultiLineString
from os import listdir
from os.path import isfile, join

########################################################################################
#
#  Constructing network topology of the ERCOT grid using available GIS transmission
#  line data
#  Clara Berger July 2023
#
########################################################################################


#### Get grid of Texas using transmission line data and the shape of texas
def get_texas_grid(grid_path, states_path, epsg_grid, epsg):
    # Read in US grid shapefile
    grid_full = gpd.read_file(grid_path)
    # assign the CRS
    grid_full = grid_full.set_crs(epsg_grid)
    grid_full = grid_full.to_crs(epsg)

    # Cut out Texas 
    # read in shapefile of US states 
    us_states = gpd.read_file(states_path)
    # align the CRS of the states shapefile with the grid shapefile 
    us_states = us_states.to_crs(epsg)

    # take only Texas
    texas = us_states[us_states['NAME']=='Texas']

    # perform spatial join between the grid and the shape Texas 
    # to approximate the ERCOT grid
    tx_grid = grid_full.sjoin(texas, how="inner", predicate='intersects').reset_index(drop=True)

    return tx_grid


#### Find the location of each node
def get_lat_long(df, epsg, start_node, end_node):
    # find the start and end point of each transmission line
    for index, row in df.iterrows():
        coords = [(coords) for coords in list(row['geometry'].coords)]  #.boundary.geoms)]
        first_coord, last_coord = [ coords[i] for i in (0, -1) ]
        df.at[index,'start'] = Point(first_coord)
        df.at[index,'end'] = Point(last_coord)
    df_bounds = df.set_crs(epsg)

    # to find a unique list of nodes get all of the starting stations
    start_df = df[[start_node, 'start']]
    start_df.columns = ['SUB', 'location']
    # and the ending stations
    end_df = df[[end_node, 'end']]
    end_df.columns = ['SUB', 'location']
    
    # and set concatenate these lists so that we can find unique points
    subs = pd.concat([start_df, end_df]).set_geometry('location')
    subs = subs.set_crs(epsg)

    ########
    # for now just take all of the unique points - can investigate names as well

    # create dataframe with unique name for each node
    node_names = ['STATION{}'.format(i) for i in range(1,len(subs['location'].unique())+1)]
    node_df = gpd.GeoDataFrame({'NODE_ID':node_names, 'geometry':subs['location'].unique()})

    # finally add lat and long columns to the nodes table
    node_df['LATITUDE'] = node_df.geometry.y
    node_df['LONGITUDE'] = node_df.geometry.x

    # merge start and end points with nodes to get start_node end_node
    # these are unique by point so we will use these instead of the SUB names
    nodes_lines = pd.merge(df_bounds, node_df, left_on='start', right_on='geometry').rename(columns={'NODE_ID':'start_node'}) 
    nodes_lines = pd.merge(nodes_lines, node_df, left_on='end', right_on='geometry').rename(columns={'NODE_ID':'end_node'}) 
    nodes_lines = nodes_lines.drop(['geometry_y','geometry'], axis='columns')

    return node_df, nodes_lines


#### Make an adjacency matrix from this table showing which nodes are connected to which
def make_adj_matrix(df, index, columns, nodes, node_col):
    # use crosstab to get which nodes are connected
    adj = pd.crosstab(df[index], df[columns])

    # since the network is undirected, make sure opposite edges exist

    # reindex such that all nodes are along the rows and columns
    idx = nodes[node_col]  #indices are list of nodes
    adj = adj.reindex(index = idx, columns=idx, fill_value=0)

    # since the network is undirected, make sure opposite edges exist

    return adj


#### Read in other relevant dataframes from Aurors and combine all of the years togetehr
def read_aurora_files(dir, variables):
    # dictionary for concatenated dfs
    dict_of_dfs = dict()
    for var in variables:
        # get path to the price data, the demand data, and the generation data
        paths = [join(dir,p) for p in listdir(dir) if p.startswith(var)]
        # combine multiple years
        for path in paths:
            # get all files in subdirectory
            files = [f for f in listdir(path) if f.endswith(".csv")] 
            # combine multiple years
            df_csv_concat = pd.concat([pd.read_csv(join(path, fl), header=[1]) for fl in files ], ignore_index=True)
            dict_of_dfs[var] = df_csv_concat

    return dict_of_dfs


#### Do a spatial join between the nodes from GIS and Aurora's nodes
def find_node_matches(df_l, df_r, epsg_proj, threshold, lsuffix, rsuffix):
    # first we have to project the nodes to a projected CRS if they are not already in one
    df_l_proj = df_l.to_crs(epsg_proj)
    df_r_proj = df_r.to_crs(epsg_proj)

    # then we can do a nearest spatial join with a max distance threshold to find another node
    merged = gpd.sjoin_nearest(df_l_proj, df_r_proj, max_distance=threshold, lsuffix=lsuffix, rsuffix=rsuffix) 
    return merged


