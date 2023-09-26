library("dplyr")
library("readr")
library("xts")
library("GNAR")
library("ggplot2")
library("igraph")
library("tseries")
library("Metrics")
library("forecast")
library("MASS")
library("scales")
library("lubridate")

# set working directory
wkdr = "/Users/claraberger/Library/CloudStorage/OneDrive-Nexus365/Documents/Dissertation/"
fig.path = "/Users/claraberger/Library/CloudStorage/OneDrive-Nexus365/Documents/Dissertation/figures/"

source(paste0(wkdr,"code/GNAR.R"))


######################################################################################################
##### Create Network
######################################################################################################
# Create grid network using adjacency matrix
adj.path <- paste0(wkdr,"tables/grid_adj.csv")
grid <- create_network(wkdr, adj.path)

# Read in the node names and locations 
nodes.path <- paste0(wkdr,"tables/nodes.csv")
nodes<- read_csv(nodes.path) 

# read in Aurora's connections 
matched.path <- paste0(wkdr,"tables/matched_tot.csv")
matched<- read_csv(matched.path) %>% mutate_all(as.character) #convert all numeric columns to characters


######################################################################################################
##### Manipulate Historical Congestion Data
######################################################################################################

# Get historical congestion data
cong.path <- paste0(wkdr,"tables/hist_cong.csv")
hist.congraw<- read_csv(cong.path) 
times <- hist.congraw$timeutcminus6

y2016 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2016) # 213 nodes
y2017 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2017) # 221 nodes
y2018 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2018) # 230 nodes
y2019 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2019) # 233 nodes
y2020 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2020) # 251 nodes
y2021 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2021) # 281 nodes
y2022 <- which(format(hist.congraw$timeutcminus6, format = "%Y") == 2022) # 303 nodes

year.start <- 2017
year.end <- 2017

hist.cong<- manipulate_congestion_ts(raw.input = hist.congraw,
                                     matches.aur = matched$NODE_ID_aur,
                                     matches.gis = matched$NODE_ID_gis,
                                     log.bool = FALSE,
                                     clip.bool = TRUE,
                                     year.idx = c(y2017),
                                     rmpartialNA.bool = TRUE
)


######################################################################################################
##### Get One Particular Year and a Connected Component
######################################################################################################

# get largest connected component after connecting 2nd or 3rd stage neighbors
con.grid <- get_connected_component(grid.input = grid,
                                    connect.stage = 1,
                                    hist.cong = hist.cong)

summary(con.grid)

# Make time series with only those nodes in the connected connected component
cong.ts <- ts(hist.cong[,colnames(hist.cong) %in% V(con.grid)$name])

# Test taking differences
take_differences <- function(ts.input, lags=1, seasonal.option, seasonal.lags=24){
  if(seasonal.option == 'season'){
    # Take the first difference from the seasonal lag
    diff.ts <- ts(sapply(ts.input, function(x) diff(as.numeric(x), lag = seasonal.lags, differences = 1))) 
    return(diff.ts)
  }
  if(seasonal.option == 'seasonfirst'){
    # Take the first difference from the seasonal lag
    diffseason.ts <- ts(sapply(ts.input, function(x) diff(as.numeric(x), lag = seasonal.lags, differences = 1))) 
    # then take the first different for the first lag
    diff.ts <- ts(sapply(diffseason.ts, function(x) diff(as.numeric(x), lag =  lags, differences = 1))) 
    return(diff.ts)
  }
  if(seasonal.option == 'firstseason'){
    # then take the first different for the first lag
    diffone.ts <- ts(sapply(ts.input, function(x) diff(as.numeric(x), lag = lags, differences = 1))) 
    # then take the first difference for the seasonal lag
    diff.ts <- ts(sapply(diffone.ts, function(x) diff(as.numeric(x), lag = seasonal.lags, differences = 1))) 
    return(diff.ts)
  }
  if(seasonal.option == 'first'){
    # Take first difference  
    diff.ts <- ts(sapply(ts.input, function(x) diff(as.numeric(x), lag = lags, differences = 1))) 
    return(diff.ts)
  }
}

# Try different combinations of differences 
diff.s <- take_differences(ts(cong.ts), seasonal.option = 'season', seasonal.lags=24)
diff.sf <- take_differences(ts(cong.ts), seasonal.option = 'seasonfirst', seasonal.lags=24)
diff.fs <- take_differences(ts(cong.ts), seasonal.option = 'firstseason', seasonal.lags=24)
diff.f <- take_differences(ts(cong.ts), seasonal.option = 'first')


######################################################################################################
##### Prepare the Data for Modeling
######################################################################################################

# Convert graph from igraph to GNAR
gridnet <- igraphtoGNAR(con.grid)
plot(gridnet)

# Reorder the time series to match up nodes and time series columns
name.order <- V(con.grid)$name
diff.s <- diff.s[, name.order]
diff.sf <- diff.sf[, name.order] 
diff.fs <- diff.fs[, name.order] 
diff.f <- diff.f[, name.order]
cong.ts <- cong.ts[, name.order]


colnames(diff.s)
V(con.grid)

######################################################################################################
##### Fit GNAR Model
######################################################################################################
# Set up the combinations of parameters we would like to try
alphas <- c(rep(1,2), rep(2,6) , rep(24,3))
betas <- list(c(0), c(1), 
              c(0,0), c(1,0), c(1,1), c(2,0), c(2,1), c(2,2), 
              c(rep(0,24)), c(1,rep(0,23)), c(2,rep(0,23)))

# Specify the date range for predicting
# Here we use data through 30th November to predict 1st December
start <- as.POSIXct(paste0(year.start,"-01-01 00:00:00"), tz = "UTC")
end <- as.POSIXct(paste0(year.end,"-12-31 23:00:00"), tz = "UTC")
predstart <- as.POSIXct(paste0(year.end,"-12-01 00:00:00"), tz = "UTC")


# Fit the GNAR mdoel for each combination of parameters and calculate the BIC of the model fit, and the 
# SSE, RMSE, and MASE of the 1 hour forecast
model.params <- find_best_params(network = gridnet, 
                                 diff.ts = diff.s, 
                                 alphas = alphas, 
                                 betas = betas,
                                 times = times,
                                 start.date = start,
                                 end.date = end,
                                 predstart.date = predstart,
                                 n.pred = 1
)

# Get the best model parameters based on the lowest forecast RMSE
bestalpha <- as.numeric(model.params$alpha[which.min(model.params$RMSE)])
bestbeta <- eval(parse(text=model.params$beta[which.min(model.params$RMSE)]))
bestglobalalpha <- as.logical(model.params$globalalpha[which.min(model.params$RMSE)])


# Use best results from the parameters and retrieve the fitted data
gnar.model <- GNARfit(vts=diff.s, net = gridnet, alphaOrder = bestalpha, betaOrder = bestbeta, globalalpha = bestglobalalpha)
gnar.model
fit <- fitted(gnar.model)
gnar.coeffs <- cbind(as.data.frame(gnar.model$mod$coefficients), confint(gnar.model))

# Undifference data 
rebuilt <- undo_differences(cong.ts, fit, seasonal.option = 'season')


######################################################################################################
##### Find the Best Model for each Year and Connection Stage
######################################################################################################
# which years to test and the ids associated with that year
test.years <- list(c(start = 2016, end = 2016, ids = list(y2016)), 
                   c(start = 2017, end = 2017, ids = list(y2017)),
                   c(start = 2018, end = 2018, ids = list(y2018)),
                   c(start = 2019, end = 2019, ids = list(y2019)),
                   c(start = 2020, end = 2020, ids = list(y2020)),
                   c(start = 2021, end = 2021, ids = list(y2021)),
                   c(start = 2022, end = 2022, ids = list(y2022)),
                   c(start = 2016, end = 2022, ids = list(c(y2016,y2017,y2018,y2019,y2020,y2021,y2022)))
)
# which stage connections to test
test.conn <- c(1,2)

# create dataframe to store the results of the best model each year 
test_all_years_conn <- function(){
  # Set up dataframe to check multiple years
  year.params <- data.frame(year = character(0),
                            conn_stage = numeric(0),
                            nnodes = numeric(0),
                            RMSE = numeric(0),
                            MASE = numeric(0),
                            best_model = character(0)
  )
  j <- 1
  for(year in test.years){
    for(conn in test.conn){
      year.params.row <- best_model_by_year(year.start = year$start, year.end = year$end, year.idx = year$ids, connect.stage = conn, n.pred = 1)
      year.params[j, ] <- year.params.row
      j <- j +1
    }
  }
  return(year.params)
}
year.params <- test_all_years_conn()


######################################################################################################
##### Assess Model Performance
######################################################################################################

# Look at model fit residuals
# convert the residuals to matrix form
model.resids <- as.data.frame(residToMat(gnar.model, nnodes=length(V(con.grid)))$resid)
se.residuals <- as.data.frame(sapply(model.resids, function(x) sqrt(var(x))))

rmses.diff <- calculate_rmse(diff.s[(nrow(diff.s)-nrow(fit)+1):nrow(diff.s),], fit)
rmses <- calculate_rmse(cong.ts[(nrow(cong.ts)-nrow(rebuilt)+1):nrow(cong.ts),], rebuilt)
min.rmse <- which.min(rmses$RMSE)
max.rmse <- which.max(rmses$RMSE)

# Make plots of the model fit, the PDC plots, and the residuals
make_performance_plots(fit = fit,
                       year = paste0(as.character(year.start), "-",
                                     as.character(year.end)),
                       model.resids = model.resids,
                       name.order = name.order,
                       diff.ts = diff.s,
                       cong.ts = cong.ts,
                       rebuilt = rebuilt
)

######################################################################################################
##### Forecast
######################################################################################################

# Specify the date range
start <- as.POSIXct(paste0(year.start,"-01-01 00:00:00"), tz = "UTC")
end <- as.POSIXct(paste0(year.end,"-12-31 23:00:00"), tz = "UTC")
predstart <- as.POSIXct(paste0(year.end,"-12-01 00:00:00"), tz = "UTC")

# Set number of hours ahead we want to predict
n.pred <- 1
forecast.sse <- forecast_congestion(diff.ts = diff.s,
                                    times = times,
                                    start.date = start,
                                    end.date = end,
                                    predstart.date = predstart,
                                    n.pred = n.pred,
                                    year = paste0(substr(as.character(year.start),3,4), "-",
                                                  substr(as.character(year.end),3,4)),
                                    alpha = bestalpha,
                                    beta = bestbeta,
                                    globalalpha = bestglobalalpha,
                                    gridnet = gridnet,
                                    cong.ts = cong.ts,
                                    fit = fit,
                                    rebuilt = rebuilt,
                                    name.order = name.order
)


