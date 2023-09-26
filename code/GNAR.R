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

######################################################################################################
##### Create Network
######################################################################################################
# Create grid network using adjacency matrix
create_network <- function(wkdr, adj.path){
  # Read in adjacency matrix of substations
  adj<- read_csv(adj.path) 
  adj <- adj[, -1]    # remove the first variable
  # convert to matrix
  adj.matrix <- as.matrix(adj)  
  # make symmetric by element-wise maximum between 'adj.matrix' and its transpose
  adj.matrix <- pmax(adj.matrix, t(adj.matrix))
  # convert adjacency matrix into a network
  grid <- graph_from_adjacency_matrix(adj.matrix, mode="max", weighted=NULL) 
  return(grid)
}

######################################################################################################
##### Manipulate Historical Congestion Data
######################################################################################################

# Average the costs across duplicated nodes
average_cong_costs <- function(hist.input, matches.aur, matches.gis){
  # Only keep the columns of Aurora nodes that had a match in the GIS data
  hist.input <- hist.input %>% 
    mutate_at(vars(-timeutcminus6), as.numeric) %>% 
    select_if(names(.) %in% matches.aur)
  
  # Use recode() to rename the column names with the names of the corresponding GIS nodes
  #this will cause multiple columns with the same name if multiple Aurora points were at the same location
  colnames(hist.input) <- recode(colnames(hist.input), !!!setNames(matches.gis, matches.aur))
  
  # For all columns with duplicate names, take the average of the congestion costs 
  hist.avg <- as.data.frame(
    sapply(unique(names(hist.input)), # for each unique column name
           function(col) rowMeans(hist.input[names(hist.input) == col], na.rm = TRUE) # calculate row means
    ))
  return(hist.avg)
}
# Take the log of the costs adding a constant to shrink the spread of the data
log_costs <- function(hist.input, mins){
  # Take the log of all the data adding constant based on the values of each node
  hist.log <- as.data.frame(mapply(function(x,y) log(x + 1.01*abs(y)), hist.input, mins))
}
# Undo the log transform by exponentiating the costs and subtracting the constant
exp_costs <- function(hist.log, mins){
  # Take the log of all the data adding constant based on the values of each node
  hist.output <- as.data.frame(mapply(function(x,y){exp(x) - 1.01*abs(y)}, hist.log, mins))
}
# Since we don't care about the magnitude of the big spikes, cut off values outside of 3 stdevs from the mean
clip_spikes <- function(vec.input){
  # Get the 97.7 and 2.3 percentiles to be the upper and lower bounds of data we want
  upper <- mean(vec.input, na.rm = TRUE) + 2*sd(vec.input, na.rm = TRUE)
  lower <- mean(vec.input, na.rm = TRUE) - 2*sd(vec.input, na.rm = TRUE)
  # If a value is outside 3 stdevs, set it to the upper and lower limit we set 
  clipped <- ifelse(vec.input>upper, upper, vec.input)
  clipped <- ifelse(clipped<lower, lower, clipped)
  return(clipped)
}
# We want each column to be a timeseries object
create_ts <- function(hist.input, times){
  # Add back in a time variable
  hist.input$time <- times
  # Create time series object from each column of input 
  cong.xts <- xts(hist.input, order.by = time)
  cong.ts <- ts(hist.input)#, frequency = 24)
  #cong.ts <- cong.ts[, -ncol(cong.ts)]
}

manipulate_congestion_ts <- function(raw.input, matches.aur, matches.gis, log.bool, clip.bool, year.idx, rmpartialNA.bool){
  # average congestion costs over nodes at the same location
  hist.cong <- average_cong_costs(raw.input, matches.aur=matches.aur, matches.gis=matches.gis)

  # do we want to log the data? 
  if(log.bool == TRUE){
    # to take the long of the data we want to make the data positive by adding the minimum value 
    # at each node 
    mins <- lapply(hist.cong, min, na.rm=TRUE)
    hist.conglog <- log_costs(hist.cong, mins)
    # reassign to the historical congestion costs
    hist.cong <- hist.conglog
  }
 
  if(clip.bool == TRUE){
    hist.congclip <- as.data.frame(sapply(hist.cong, clip_spikes))
    hist.cong <- hist.congclip
  }

  # create a multivariate time series from the the costs
  cong.ts <- ts(hist.cong)
  
  if(rmpartialNA.bool == TRUE){
    # take just the year(s) you are interested in 
    # and years where all data is available - i.e. remove columns with any NAs
    cong.year <- cong.ts[year.idx, colSums(!is.na(hist.cong[year.idx,])) == nrow(hist.cong[year.idx,])] 
    return(cong.year)
  }
  else{
    # else just take just the year(s) you are interested in but only remove columns where every value is NA 
    cong.year <- cong.ts[year.idx, colSums(!is.na(hist.cong[year.idx,])) != 0] 
    return(cong.year)
  }
}

######################################################################################################
##### Get One Particular Year and a Connected Component
######################################################################################################

# get largest connected component after connecting 2nd or 3rd stage neighbors
get_connected_component <- function(grid.input, connect.stage, hist.cong){
  # connected nodes with their X-stage neighbors
  grid <- connect(grid.input, connect.stage)
  
  # Create subgraph with only those nodes we're using with cost data
  sub.grid <- induced_subgraph(grid, colnames(hist.cong))
  # Find the components and get those nodes that are a part of the largest connected component
  con.comps <- components(sub.grid)
  connected.nodes <- V(sub.grid)[con.comps$membership == which.max(con.comps$csize)]
  # Induce another subgraph for just the largest connected component 
  con.grid <- induced_subgraph(sub.grid, connected.nodes) %>% simplify()
  
  return(con.grid)
}

# Take differences
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


######################################################################################################
##### Fit GNAR Model
######################################################################################################

# Find the best model parameters based on BIC
find_best_params <- function(network, diff.ts, alphas, betas, times, start.date, end.date, predstart.date, n.pred){
  
  params.df <- data.frame(alpha = numeric(0),
                       beta = character(0),
                       globalalpha = logical(0),
                       BIC = numeric(0),
                       SSE = numeric(0),
                       RMSE = numeric(0),
                       MASE = numeric(0))
  j <- 1 
  for(i in 1:length(alphas)){
    print(paste("lags=",alphas[i],"neighbors=",as.character(betas[i])))
    # Fit the GNAR model with input parameters 
    gnar.model <- GNARfit(vts=diff.ts, net = network, alphaOrder = alphas[i], betaOrder = betas[[i]], globalalpha = TRUE)
    fit <- fitted(gnar.model)
    
    bic <- BIC(gnar.model)
    
    predend.date <- predstart.date + hours(n.pred-1)
    
    # List of relevant times based on all the data we have and then just the prediction region
    all.times <- times[times >= start.date & times <= end.date]
    pred.times <- times[times >= start.date & times <= predend.date] # times used to make prediction + predicted indices
    
    # Index within the differenced dataframe of the end of values used for prediction
    pred.idx <- length(diff.ts[,1]) - 1 - as.numeric(difftime(end.date, predstart.date, units="hours")) # back up the appropriate number of days 
    
    # Use the data up until the predicted region to forecast ahead the set number of values
    pred <- predict(GNARfit(vts=diff.ts[1:pred.idx,], net = network, alphaOrder = alphas[i], betaOrder = betas[[i]], globalalpha = TRUE), 
                    n.ahead = n.pred)
    
    sse.d <- round(sum((diff.ts[(pred.idx+1):(pred.idx+n.pred), ] - pred)^2),4)
    rmse.d <- round(rmse(diff.ts[(pred.idx+1):(pred.idx+n.pred), ], pred),4)
    mase.d <- round(computeMASE(forecast = pred, train = diff.ts[1:(pred.idx), ], test = diff.ts[(pred.idx+1):(pred.idx+n.pred), ], period = 24), 4)
    
    
    params.row <- c(alpha=alphas[i], beta=as.character(betas[i]), globalalpha='TRUE' , BIC=round(bic,4), SSE = sse.d, RMSE = rmse.d, MASE = mase.d)
    params.df[j,] <- params.row
    j = j+1
    
  }
  return(params.df)
}

# Undifference data 
undo_differences <- function(ts.input, fit, lags=1, seasonal.option, seasonal.lags=24){
  if(seasonal.option == 'season'){
    # Get the initial values to build off of
    xi24 <- ts.input[1:24,]
    # Undifference values starting with the initial 
    rebuilt <- diffinv(fit, xi=xi24, lag=seasonal.lags,differences=1)
    return(rebuilt)
  }
  if(seasonal.option == 'seasonfirst'){
    # Get the initial values to build off of
    xi24 <- ts.input[1:24,]
    xi1 <- t(diff.s[1,])
    # Undifference values starting with the initial 
    rebuilt <- diffinv(diffinv(fit, xi=xi1, lag=lags, differences=1),xi=xi24, lag=seasonal.lags, differences=1)
    return(rebuilt)
  }
  if(seasonal.option == 'firstseason'){
    # Get the initial values to build off of
    xi24 <- diff.f[1:24,]
    xi1 <- t(ts.input[1,]) 
    # Undifference values starting with the initial 
    rebuilt <- diffinv(diffinv(fit, xi=xi24, lag=seasonal.lags, differences=1),xi=xi1, lag=lags, differences=1)
    return(rebuilt)
  }
  if(seasonal.option == 'first'){
    # Get the initial values to build off of
    xi1 <- t(ts.input[1,]) 
    # Undifference values starting with the initial 
    rebuilt <- diffinv(fit, xi=xi1, lag=lags, differences=1)
    return(rebuilt)
  }
}


######################################################################################################
##### Find the Best Model for each Year and Connection Stage
######################################################################################################

# Fit the GNAR model for each year and stage connection 
# Uses seasonally differenced data
best_model_by_year <- function(year.start, year.end, year.idx, connect.stage, n.pred){
  ##### Manipulate Historical Congestion Data 
  hist.cong<- manipulate_congestion_ts(raw.input = hist.congraw,
                                       matches.aur = matched$NODE_ID_aur,
                                       matches.gis = matched$NODE_ID_gis,
                                       log.bool = FALSE,
                                       clip.bool = TRUE,
                                       year.idx = year.idx,
                                       rmpartialNA.bool = TRUE
  )
  ##### Get One Particular Year and a Connected Component
  con.grid <- get_connected_component(grid.input = grid,
                                      connect.stage = connect.stage,
                                      hist.cong = hist.cong)
  # Make time series with only those nodes in the connected connected component
  cong.ts <- ts(hist.cong[,colnames(hist.cong) %in% V(con.grid)$name])
  
  # Take seasonal differences 
  diff.s <- take_differences(ts(cong.ts), seasonal.option = 'season', seasonal.lags=24)
  
  # Convert graph from igraph to GNAR
  gridnet <- igraphtoGNAR(con.grid)
  
  # Reorder the time series to match up nodes and time series columns
  name.order <- V(con.grid)$name
  print(name.order)
  diff.s <- diff.s[, name.order]
  cong.ts <- cong.ts[, name.order]
  
  # Set up the combinations of parameters we would like to try
  alphas <- c(rep(1,2), rep(2,6) , rep(24,3))
  betas <- list(c(0), c(1), 
                c(0,0), c(1,0), c(1,1), c(2,0), c(2,1), c(2,2), 
                c(rep(0,24)), c(1,rep(0,23)), c(2,rep(0,23)))
  
  # Specify the date range
  start <- as.POSIXct(paste0(year.start,"-01-01 00:00:00"), tz = "UTC")
  end <- as.POSIXct(paste0(year.end,"-12-31 23:00:00"), tz = "UTC")
  predstart <- as.POSIXct(paste0(year.end,"-12-01 00:00:00"), tz = "UTC")
  #predend <- as.POSIXct(paste0(year.end,"-12-02 23:00:00"), tz = "UTC")
  
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
  # Extract the best parameters based on SSE 
  bestalpha <- as.numeric(model.params$alpha[which.min(model.params$RMSE)])
  bestbeta <- eval(parse(text=model.params$beta[which.min(model.params$RMSE)]))
  bestglobalalpha <- as.logical(model.params$globalalpha[which.min(model.params$RMSE)])
  bestrmse <- as.numeric(model.params$RMSE[which.min(model.params$RMSE)])
  bestmase <- as.numeric(model.params$MASE[which.min(model.params$RMSE)])
  
  # Output the year and the best GNAR model based on SSE
  year.params.row <- c(year = paste(year.start,"to",year.end),
                       conn_stage = connect.stage,
                       nnodes = vcount(con.grid),
                       RMSE = bestrmse, 
                       MASE = bestmase,
                       best_model = paste("GNAR(",bestalpha,",[",toString(bestbeta),"])"))
  
  return(year.params.row)
}


######################################################################################################
##### Assess Model Performance
######################################################################################################

# Function to calculate Root Mean Squared Error (RMSE)
calculate_rmse <- function(actual, predicted) {
  # Make into DFs to deal with NAs
  actualv <- as.data.frame(actual) %>% mutate_all(as.numeric)
  predictedv <- as.data.frame(predicted) %>% mutate_all(as.numeric)
  
  # Calculate RMSE for each corresponding column
  rmse_values <- lapply(1:ncol(actualv), function(i) {
    round(rmse(na.omit(actualv[,i]), na.omit(predictedv[,i])),4)}
  )
  
  # Create a dataframe to store the RMSE values
  rmse_df <- data.frame(
    Node = names(actualv),
    RMSE = unlist(rmse_values)
  )
  mean_row <- data.frame(Node = "average", RMSE = mean(rmse_df$RMSE))
  rmse_df <- rbind(rmse_df, mean_row)
  return(rmse_df)
}

# Get a scaled error performance measyre
computeMASE <- function(forecast,train,test,period){
  
  # forecast - forecasted values
  # train - data used for forecasting .. used to find scaling factor
  # test - actual data used for finding MASE.. same length as forecast
  # period - in case of seasonal data.. if not, use 1
  
  forecast <- as.vector(forecast)
  train <- as.vector(train)
  test <- as.vector(test)
  
  n <- length(train)
  scalingFactor <- sum(abs(train[(period+1):n] - train[1:(n-period)])) / (n-period)
  
  et <- abs(test-forecast)
  qt <- et/scalingFactor
  meanMASE <- mean(qt)
  return(meanMASE)
}

# Make performance plots
# Must have 3 subdirectiories in the directory given by "fig.path" named model_fit, PDC, and resids 
make_performance_plots <- function(year, fit, model.resids, name.order, diff.ts, cong.ts, rebuilt ){
  for(i in 1:ncol(fit)){
    ########### Visualise the model fit, isolating one node to see how the fit looks through time
    # On the differenced scale
    png(filename = paste0(fig.path,"model_fit/",name.order[i],"_",year,"fit_diff.png"))
    par(cex =1.5)
    plot(diff.ts[(nrow(diff.ts)-nrow(fit)):nrow(diff.ts),i], type = "l", ylab = "Cost Differences", xlab = "Time")
    lines( fit[,i], col = "red")
    legend("topright", legend = c("Actual", "Fitted"), col = c("black", "red"), lwd = 2)
    title(paste(name.order[i],"for",year))
    dev.off()
    # On the original undifferenced scale
    png(filename = paste0(fig.path,"model_fit/",name.order[i],"_",year,"fit_undiff.png"))
    par(cex =1.5)
    plot(cong.ts[(nrow(cong.ts)-nrow(rebuilt)):nrow(cong.ts),i], type = "l", ylab = "Costs", xlab = "Time")
    lines(rebuilt[, i], col = 'red')
    legend("topright", legend = c("Actual", "Fitted"), col = c("black", "red"), lwd = 2)
    title(paste(name.order[i],"for",year))
    dev.off()
    
    ########### Create PDC plots
    # On the differenced scale
    png(filename = paste0(fig.path,"PDC/",name.order[i],"_",year,"_PDC_diff.png"))
    par(cex =1.5)
    plot(sort(diff.ts[,i], decreasing = TRUE), type='l',  ylab = "Cost Differences", xlab = "Rank")
    lines(sort(fit[,i], decreasing = TRUE), type='l',  col="red")
    legend("topright", legend = c("Actual", "Fitted"), col = c("black", "red"), lwd = 2)
    title(paste(name.order[i],"for",year))
    dev.off()
    # On the original undifferenced scale
    png(filename = paste0(fig.path,"PDC/",name.order[i],"_",year,"_PDC_undiff.png"))
    par(cex =1.5)
    plot(sort(cong.ts[,i], decreasing = TRUE), type='l',  ylab = "Cost", xlab = "Rank")
    lines(sort(rebuilt[,i], decreasing = TRUE), type='l',  col="red")
    legend("topright", legend = c("Actual", "Fitted"), col = c("black", "red"), lwd = 2)
    title(paste(name.order[i],"for",year))
    dev.off()
    
    # Plot the residuals over time and the histogram of the residuals
    png(filename = paste0(fig.path,"resids/",name.order[i],"_",year,"_diff_resids.png"))
    par(mfrow = c(2, 1), par(cex =1.5))
    plot(ts(model.resids[,i]), ylab="")
    title(paste(name.order[i],"GNAR model residuals"))
    hist(model.resids[,i], main = "", breaks=50, xlab=paste(name.order[i],"model residuals"))
    par(mfrow = c(1, 1), par(cex =1.5))
    dev.off()
    
    # Plot the residuals using the check residuals function
    png(filename = paste0(fig.path,"resids/",name.order[i],"_",year,"_diff_checkresids.png"))
    par(cex =1.5)
    p <- ggtsdisplay(model.resids[,i], plot.type = "histogram", xlab="Time", main = paste(name.order[i],"Residuals"), theme=theme(text = element_text(size = 16)))
    p <- p + labs(y = "Frequency", x = "Residuals", cex=1.6)   
    print(p, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 2))
    dev.off()
    
    # Plot QQPlot of residuals
    png(filename = paste0(fig.path,"resids/",name.order[i],"_",year,"_qqplot.png"))
    par(cex =1.5)
    qqnorm(model.resids[,i], main = paste("Normal Q-Q Plot for",name.order[i]))
    qqline(model.resids[,i])
    dev.off()
    
  }
  
}


######################################################################################################
##### Forecast
######################################################################################################

# Forecast future data points
# Must have subdirectiories in the directory given by "fig.path" named forecasts
forecast_congestion <- function(diff.ts, times, start.date, end.date, predstart.date, n.pred, year, alpha, beta, globalalpha, gridnet, cong.ts, fit, rebuilt, name.order){
  predend.date <- predstart.date + hours(n.pred-1)
  
  # List of relevant times based on all the data we have and then just the prediction region
  all.times <- times[times >= start.date & times <= end.date]
  pred.times <- times[times >= start.date & times <= predend.date] # times used to make prediction + predicted indices
  
  # Index within the differenced dataframe of the end of values used for prediction
  pred.idx <- length(diff.ts[,1]) - 1 - as.numeric(difftime(end.date, predstart.date, units="hours")) # back up the appropriate number of days 
  
  # Use the data up until the predicted region to forecast ahead the set number of values
  pred <- predict(GNARfit(vts=diff.ts[1:pred.idx,], net = gridnet, alphaOrder = alpha, betaOrder = beta, globalalpha = globalalpha), 
                  n.ahead = n.pred)
  old.plus.pred <- rbind(diff.ts[1:(pred.idx), ], pred) # combine predicted values with the differenced data that came before
  rebuilt.old.plus.pred <- undo_differences(cong.ts, old.plus.pred, seasonal.option = 'season') # undifference the data
  # for plotting we want to show the model and the forecast together
  fit.plus.pred <- rbind(fit[1:(pred.idx-24), ], pred) # combine the fit with the predicted differenced values
  rebuilt.fit.plus.pred <- rbind(rebuilt[1:(pred.idx +1), ], rebuilt.old.plus.pred[(nrow(rebuilt.old.plus.pred) - (n.pred-1)):nrow(rebuilt.old.plus.pred), ])
  
  if(n.pred==1){
    print(paste("n.pred = ", n.pred))
    ############## Create plots of the forecasted values
    # Differenced Time Scale
    diff.df <- diff.ts %>%
      as.data.frame()
    fit.df <- fit.plus.pred %>%
      as.data.frame()  %>%
      setNames(name.order) %>%
      mutate(pred = ifelse(row_number() > (n() - n.pred), 1, 0)) # add a variable for whether it is a predicted value
    # Create plots
    for(i in 1:ncol(diff.ts)){
      preddiffplot <- ggplot(data = diff.df, aes(x = c(1:nrow(diff.df)), y = eval(parse(text=colnames(diff.df)[i])))) +
        geom_line() +
        geom_point(data=pred, aes(x = (pred.idx+1), y = pred[,i], color = "red"), size = 2)+
        geom_line(data=fit.df, aes(x = c(24:(nrow(fit.df)+23)), y = eval(parse(text=colnames(fit.df)[i])), linetype = factor(pred)), color='red' ) +
        xlim((pred.idx-100),(pred.idx+n.pred)) +
        scale_linetype_manual(name="", labels= c("Model Fit",""), values=c("solid","blank")) +
        scale_color_manual(name="", labels= c("Forecasted"), values = c("red")) +
        labs(x = "Index", y = paste(colnames(diff.df)[i], " Differenced Congestion Costs"), title=paste("Cost Diff Forecasts for",date(predend.date),"using", date(start.date),"to",date(predstart.date-1)))+
        theme(text = element_text(size = 24))
      ggsave(filename = paste0(fig.path,"/forecasts/",colnames(diff.df)[i],"_diff_",n.pred,"n_",year,".png"))
    }
    
    # Original Time Scale
    # Add times and dates for plotting to the rebuilt dataset
    range.start <- predstart.date - days(3) # range of dates to display on plot
    range.end <- predend.date
    rebuilt.df <- rebuilt.fit.plus.pred %>%
      as.data.frame() %>%
      setNames(name.order) %>%
      mutate(time = pred.times[24:(length(pred.times))]) %>% # add datetime column
      mutate(pred = ifelse(row_number() > (n() - n.pred), 1, 0)) %>% # add a variable for whether it is a predicted value
      filter(time >= range.start & time <= range.end) # subset based on the range we want to plot
    cong.df <- cong.ts %>%
      as.data.frame() %>%
      mutate(time = all.times) %>%
      filter(time >= range.start & time <= range.end)
    
    for(i in 1:ncol(cong.df)){
      predplot <- ggplot(data = cong.df, aes(x = time, y = eval(parse(text=colnames(cong.df)[i])) )) +
        geom_line() +
        geom_line(data=rebuilt.df, aes(x = time, y = eval(parse(text=colnames(rebuilt.df)[i])), linetype = factor(pred)), color='red' ) +
        geom_point(data=rebuilt.df, aes(x = predstart.date, y = rebuilt.df[nrow(rebuilt.df),i], color = "red"), size = 2) +
        scale_linetype_manual(name="", labels= c("Model Fit",""), values=c("solid","blank")) +
        scale_color_manual(name="", labels= c("Forecasted"), values = c("red")) +
        labs(x = "Day", y = paste(colnames(cong.df)[i],"Congestion Costs"), title=paste("Forecasts for",date(predend.date),"using", date(start.date),"to",date(predstart.date-1))) +
        theme(text = element_text(size = 24))
      ggsave(filename = paste0(fig.path,"/forecasts/",colnames(cong.df)[i],"_",n.pred,"n_",year,".png"))
    }
    
    # Find sum of squared differences between predicted day and actual
    sse.df <- data.frame(Station = character(0),
                         SSE_diff = numeric(0),
                         RMSE_diff = numeric(0),
                         SSE = numeric(0),
                         RMSE = numeric(0))
    for(i in 1:ncol(diff.ts)){
      sse.d <- round(sum((diff.ts[(pred.idx+1):(pred.idx+n.pred),i] - pred[,i])^2),6)
      rmse.d <- round(rmse(diff.ts[(pred.idx+1):(pred.idx+n.pred),i], pred[,i]),4)
      sse <- round(sum((cong.df[((nrow(cong.df)-n.pred):nrow(cong.df)),i] - rebuilt.df[((nrow(rebuilt.df)-n.pred):nrow(rebuilt.df)),i])^2),6)
      rmse <- round(rmse(cong.ts[((nrow(cong.df)-n.pred):nrow(cong.df)),i], rebuilt.df[((nrow(rebuilt.df)-n.pred):nrow(rebuilt.df)),i]),4)
      sse.row <- c(Station=colnames(diff.ts)[i], SSE_diff=sse.d, RMSE_diff=rmse.d, SSE=sse, RMSE=rmse)
      sse.df[i,] <- sse.row
    }
    return(sse.df)
  }
  else{
    print(paste("n.pred = ", n.pred))
    ############## Create plots of the forecasted values
    # Differenced Time Scale
    diff.df <- diff.ts %>%
      as.data.frame()
    fit.df <- fit.plus.pred %>%
      as.data.frame()  %>%
      setNames(name.order) %>%
      mutate(pred = ifelse(row_number() > (n() - n.pred), 1, 0)) # add a variable for whether it is a predicted value
    # Create plot
    for(i in 1:ncol(diff.ts)){
      preddiffplot <- ggplot(data = diff.df, aes(x = c(1:nrow(diff.df)), y = eval(parse(text=colnames(diff.df)[i])))) +
        geom_line() +
        geom_line(data=fit.df, aes(x = c(24:(nrow(fit.df)+23)), y = eval(parse(text=colnames(fit.df)[i])), linetype = factor(pred)), color='red' ) +
        xlim((pred.idx-100),pred.idx+n.pred) +
        scale_linetype_manual(name="", labels= c("Model Fit","Forecasted"), values=c("solid","dashed")) +
        labs(x = "Index", y = paste(colnames(diff.df)[i], " Differenced Congestion Costs"), title=paste("Cost Diff Forecasts for",date(predend.date),"using", date(start.date),"to",date(predstart.date-1))) +
        theme(text = element_text(size = 24))
      ggsave(filename = paste0(fig.path,"/forecasts/",colnames(diff.df)[i],"_diff_",n.pred,"n_",year,".png"))
    }
    
    # Original Time Scale
    # Add times and dates for plotting to the rebuilt dataset
    range.start <- predstart.date - days(3) # range of dates to display on plot
    range.end <- predend.date
    rebuilt.df <- rebuilt.fit.plus.pred %>%
      as.data.frame() %>%
      setNames(name.order) %>%
      mutate(time = pred.times[24:(length(pred.times))]) %>% # add datetime column
      mutate(pred = ifelse(row_number() > (n() - n.pred), 1, 0)) %>% # add a variable for whether it is a predicted value
      filter(time >= range.start & time <= range.end) # subset based on the range we want to plot
    cong.df <- cong.ts %>%
      as.data.frame() %>%
      mutate(time = all.times) %>%
      filter(time >= range.start & time <= range.end)
    
    for(i in 1:ncol(cong.df)){
      predplot <- ggplot(data = cong.df, aes(x = time, y = eval(parse(text=colnames(cong.df)[i])) )) +
        geom_line() +
        geom_line(data=rebuilt.df, aes(x = time, y = eval(parse(text=colnames(rebuilt.df)[i])), linetype = factor(pred)), color='red' ) +
        scale_linetype_manual(name="", labels= c("Model Fit","Forecasted"), values=c("solid","dashed")) +
        labs(x = "Day", y = paste(colnames(cong.df)[i],"Congestion Costs"), title=paste("Forecasts for",date(predend.date),"using", date(start.date),"to",date(predstart.date-1))) +
        theme(text = element_text(size = 24))
      ggsave(filename = paste0(fig.path,"/forecasts/",colnames(cong.df)[i],"_",n.pred,"n_",year,".png"))
    }
    
    # Find sum of squared differences between predicted day and actual
    sse.df <- data.frame(Station = character(0),
                         SSE_diff = numeric(0),
                         RMSE_diff = numeric(0),
                         SSE = numeric(0),
                         RMSE = numeric(0))
    for(i in 1:ncol(diff.ts)){
      sse.d <- round(sum((diff.ts[(pred.idx+1):(pred.idx+n.pred),i] - pred[,i])^2),6)
      rmse.d <- round(rmse(diff.ts[(pred.idx+1):(pred.idx+n.pred),i], pred[,i]),4)
      sse <- round(sum((cong.df[((nrow(cong.df)-n.pred):nrow(cong.df)),i] - rebuilt.df[((nrow(rebuilt.df)-n.pred):nrow(rebuilt.df)),i])^2),6)
      rmse <- round(rmse(cong.ts[((nrow(cong.df)-n.pred):nrow(cong.df)),i], rebuilt.df[((nrow(rebuilt.df)-n.pred):nrow(rebuilt.df)),i]),4)
      sse.row <- c(Station=colnames(diff.ts)[i], SSE_diff=sse.d, RMSE_diff=rmse.d, SSE=sse, RMSE=rmse)
      sse.df[i,] <- sse.row
    }
    return(sse.df)
  }
  
  
}

