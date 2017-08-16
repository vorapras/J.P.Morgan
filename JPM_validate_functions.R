#####################################################################################
## These are the functions needed in order to run script "JPM_validate.R"
#####################################################################################
## Writen by: Tom Parry, Thanachai Voraprasertsak, Wenqian Cheng and Maria Shevchenko

## These functions are the same as in the file JPM_functions except that they
## omit any lines with print strings or produce plots
## The purpose of these functions is to be able to run and test models on 
## multiple days data

## format_data formats excel data into evenly spaced time series with one price per minute based on median price
# excel_data is the data to be formated
# name is the variable to format, i.e. "forward_bid_price"
format_data = function(excel_data, name){
  
  # Select first two column of data (timestamp and price)
  excel_data = excel_data[,1:2]
  
  # Format Timestamp and numbers
  raw_timestamp = rev(strptime(excel_data[,1], format="%F %T"))
  price = rev(as.numeric(gsub(",", "", excel_data[,2])))
  
  # increase timestamps by 5 hours to change EST timezone to GMT
  raw_timestamp = raw_timestamp+(60*60*5)+1
  
  # Take price for each minute as median value of previous minute
  trunc_timestamp = trunc(raw_timestamp, units="mins")
  df_temp = data.frame(trunc_timestamp, price)
  df = df_temp %>% group_by(trunc_timestamp) %>% summarise(price = median(price))
  colnames(df) = c("timestamp", name)
  
  return(df)
}

## count_data counts the number of data points per day
# price_data is the testing data
# name is the variable that you are working with
count_data = function(price_data, name){
  
  # Create dataframe with only date of timestamp
  trun_timestamp = strptime(as.matrix(price_data[,1]), format="%F")
  df_temp = data.frame(trun_timestamp, price_data[,2])
  
  # Count the number of data points per day
  num = df_temp %>% group_by(trun_timestamp) %>% summarize(count=n())
  colnames(num) = c("timestamp", paste("count",name,sep=""))
  
  return(num)
}

## choose_data filters out days where the number of data points is below the threshold value
# df_spot_bid is a dataframe containing the spot_bid values, similar for other arguments
choose_data = function(df_spot_bid, df_spot_ask, df_future_bid, df_future_ask, threshold){
  
  num_spot_bid = count_data(df_spot_bid,"_spot_bid")
  num_spot_ask = count_data(df_spot_ask,"_spot_ask")
  num_future_bid = count_data(df_future_bid,"_future_bid")
  num_future_ask = count_data(df_future_ask,"_future_ask")
  
  num_time_stamp_spot = merge(rbind(num_spot_bid)[,1], rbind(num_spot_ask)[,1])
  num_spot = cbind(num_time_stamp_spot, num_spot_bid[,2], num_spot_ask[,2])
  num_time_stamp_future = merge(rbind(num_future_bid)[,1], rbind(num_future_ask)[,1])
  num_future = cbind(num_time_stamp_future, num_future_bid[,2], num_future_ask[,2])
  num_data = merge(num_future, num_spot, by="timestamp")
  
  # Create array of days to be used based on days where number of data points is above threshold value
  condition = (num_data$count_future_bid > threshold)&(num_data$count_future_ask > threshold)&(num_data$count_spot_bid > threshold)&(num_data$count_spot_ask > threshold)
  condition_data = num_data[condition,]
  days = condition_data[,1] # Days to be used
  num_days = length(days)
  
  # Remove days where the number of data points is less than threshold value for data array
  df_spot_bid$temp.timestamp=strptime(as.matrix(df_spot_bid[,1]), format="%F")
  df_spot_ask$temp.timestamp=strptime(as.matrix(df_spot_ask[,1]), format="%F")
  df_future_bid$temp.timestamp=strptime(as.matrix(df_future_bid[,1]), format="%F")
  df_future_ask$temp.timestamp=strptime(as.matrix(df_future_ask[,1]), format="%F")
  
  final_spot_bid = vector("list",num_days)
  final_spot_ask = vector("list",num_days)
  final_future_bid = vector("list",num_days)
  final_future_ask = vector("list",num_days)
  final_data = vector("list",num_days)
  
  # Only include days that we want
  for (i in 1:num_days){
    final_spot_bid[[i]] = df_spot_bid[df_spot_bid$temp.timestamp == days[i],][,1:2]
    final_spot_ask[[i]] = df_spot_ask[df_spot_ask$temp.timestamp == days[i],][,1:2]
    final_future_bid[[i]] = df_future_bid[df_future_bid$temp.timestamp == days[i],][,1:2]
    final_future_ask[[i]] = df_future_ask[df_future_ask$temp.timestamp == days[i],][,1:2]
    future_merge = merge(final_future_bid[[i]], final_future_ask[[i]], by="timestamp")
    spot_merge = merge(final_spot_bid[[i]], final_spot_ask[[i]], by="timestamp")
    final_data[[i]] = merge(future_merge, spot_merge, by="timestamp")
  }
  
  return_objects = list("final_data"=final_data, "num_days"=num_days)
  return(return_objects)
}

## lag_1 lags a matrix (or vector) of time series data, where time runs down the rows.
# An NA is inserted in the front of each time series so the length remains the same.
lag_1 = function(data){
  # Test if data is a vector
  if (is.null(dim(data))){
    new_data = data
    new_data[1] = NA
    for (i in 2:length(data)){
      new_data[i] = data[i-1]
    }
    return(new_data)
  } else {
    new_data = data
    new_data[1,] = NA
    for (i in 2:dim(data)[1]){
      new_data[i,] = data[i-1,]
    }
    return(new_data)
  }
}

## Normalise a matrix of column vectors, i.e. each column of output is a unit vector
norm_vec = function(x){
  out = matrix(NA,nrow=dim(x)[1],ncol=dim(x)[2])
  for (i in 1:dim(x)[2]){
    out[,i] = x[,i]/sqrt(sum(x[,i]^2))
  }
  return(out)
}

## Build a function to predicted Values from the AR model
AR_predict = function(coe,test_real,q,X_test){
  #coe <- append(coe,rep(0,nparam-length(coe)))   
  fitdata <- matrix(data =NA,nrow = length(test_real),ncol=1)
  fitdata[1:q,1] <- NA
  for(j in (q+1):length(test_real)){
    fitdata[j,1] = sum((coe*rev(test_real[(j-q):(j-1)]))) 
  }
  lagX1 = lag_1(X_test)
  lagX1 = lagX1[2:length(lagX1)]
  fitdata = fitdata+lagX1
  return (fitdata)
}

# Coefficient Matrix
ARIMA_coef = function(X_train,p){
  coef = matrix(0,nrow = p,ncol=2*p+1)
  # Create the matrix to store AIC value
  aic = matrix(NA, nrow = 1 , ncol = p)
  for(j in 1:4){
    index = c(1:j,(p+1):(p+j),9)
    fitARIMA = arima(X_train, order = c(j,0,j), method = "CSS-ML")
    Info = AIC(fitARIMA)
    aic[1,j] = Info
    for(i in 1:length(fitARIMA$coef)){
      coef[j,index[i]] = fitARIMA$coef[i]
    }
  }
  return_object = list("coef"=coef, "aic"=aic)
  return (return_object)
}

## Predicted Values from the ARIMA model
ARIMA_predict = function(coef,test_real,X_test,p){
  # Calculate the residual term
  resid = test_real-lag_1(test_real)
  
  # Data Matrix Transformation
  Xlag = as.matrix(test_real)
  Residlag = as.matrix(resid)
  for (i in 1:(p-1)){
    Xlag1 = lag_1(Xlag[,i])
    Xlag = cbind(Xlag,Xlag1)
    Residlag1 = lag_1(Residlag[,i])
    Residlag = cbind(Residlag,Residlag1)
  }
  Constant = rep(1,length(test_real))
  new_X = cbind(Xlag,Residlag,Constant)
  new_X = t(new_X)
  
  # One-Step Ahead Forecasting
  Diff_X = coef%*%new_X
  Diff_X = t(Diff_X)
  lagX1 = lag_1(X_test)
  lagX1 = lagX1[2:length(lagX1)]
  fitdata = Diff_X+lagX1
  return (fitdata)
}

## Create predicted price series based on a VAR model
# model is the VAR model to be used
# DX_test is the differenced time series to predict on
# X_test is the orignal time series to predict on
# p is the lag order of the model
build_VAR = function(model, DX_test, X_test, p){
  predict_d = matrix(0,nrow=4,ncol=dim(DX_test)[1])
  lag_data = DX_test
  for (i in 1:p){
    coef = Acoef(model)[[i]]
    lag_data = lag_1(lag_data)
    predict_d = predict_d+(coef%*%t(lag_data))
  }
  predict = rbind(NA,t(predict_d))+lag_1(X_test)
  return(predict)
}

## ecm_model returns coefficients and AIC value for ecm models built up to a particular lag
# Y is the variable you want to model
# X is the variable used to model variable Y
# numlag is the number of lags to include in the model
ecm_model = function(Y,X,numlag){
  # Run Regression
  reg.ols = lm(Y~X)
  reg_coef = reg.ols$coefficients
  resid.ols = residuals(reg.ols) 
  
  # Build ECM Model 
  n = length(resid.ols)
  name = c("dY4","dX4","dY3","dX3","dY2","dX2","dY1","dX1","hatZ","dY") 
  # Take difference for our data
  dY  = diff(Y) # Y
  dX  = diff(X) # X
  Info = matrix(NA, ncol=numlag, nrow=1)
  
  # Preallocating space
  ecm_coef = matrix(0,nrow=4,ncol=10)
  index = seq(6,0,by=-2)
  
  # Create new lag matrix for ECM Model 
  for(i in 2:(numlag+1)){
    hatZ  = resid.ols[(i+1):n]  # leave out number of window
    mat_X = cbind(embed(cbind(dY,dX),i))
    new_X = cbind(mat_X[,1:(ncol(mat_X)-2)],hatZ)
    colnames(new_X) = name[(length(name)-(2*i-1)):(length(name)-1)]
    new_Y = mat_X[,(2*i-1)]
    reg_ecm = lm(new_Y~new_X)
    recoef = reg_ecm$coefficients
    # Store coefficient every time in the right position
    for(j in 1:length(recoef)){
      ecm_coef[(i-1),index[i-1]+j]=recoef[j]
    }
    # Model Selection Criterion
    aic      = AIC(reg_ecm)
    rownames(Info) = c("aic")
    colnames(Info) = c("lag1","lag2","lag3","lag4")
    Info[,(i-1)] = rbind(aic)
  }
  return_object = list("coef"=ecm_coef, "aic"=Info)
  return (return_object)
}

## ecm_predict produces a predicted price series for a specified ecm model
# ecm_coef is the coefficients of the model
# test_Y is the actual price series of the predicted variable
# test_X is the actual price series of the variable used to predict
# numlag is the number of lags to include in the model
ecm_predict = function(ecm_coef,test_Y,test_X,numlag){
  ndate   = as.numeric(nrow(test_Y))
  # Run Regression
  reg.ols = lm(test_Y~test_X)
  hatZ    = residuals(reg.ols)  # Calculate Zhat
  hatZ    = hatZ[(numlag+2):length(hatZ)]
  # Take difference for our data
  dX  = diff(test_X) # X
  dY  = diff(test_Y) # Y
  # Take lag1 function
  lagY1 = lag_1(test_Y)
  lagY1 = lagY1[(numlag+2):length(lagY1)]
  # Create new lag matrix for ECM Model 
  mat_X = cbind(embed(cbind(dY,dX),numlag+1))
  new_X = cbind(mat_X[,1:(ncol(mat_X)-2)],hatZ)
  # Matrix Transformation
  transposeX = t(new_X)
  constantX  = rep(1,ncol(transposeX))
  X_new = rbind(constantX,transposeX)
  rownames(X_new) = c("C","dY4","dX4","dY3","dX3","dY2","dX2","dY1","dX1","hatZ")
  # One-Step Ahead Forecasting 
  diffY  = ecm_coef%*%X_new
  diffY  = t(diffY)
  fitdata = diffY+lagY1
}

#### This function calculates the standard errors for
#### the fitted ECM model by cajorls
####
#### Written by Dr Bernhard Pfaff, received on 22 April 2014 by email
####
setGeneric("abStdErr", function(x, r) standardGeneric("abStdErr"))
setMethod("abStdErr", signature = c(x = "ca.jo", r = "numeric"), 
          function(x, r){
            r <- abs(as.integer(r))
            xrestr <- cajorls(z = x, r = r)
            alpha <- coef(xrestr$rlm)[1:r, ]
            beta <- xrestr$beta
            res <- resid(xrestr$rlm)
            N <- nrow(res)
            Sigma <- crossprod(res) / N
            ## Standard Errors of Loading Matrix
            alphaSE <- matrix(sqrt(kronecker(
              diag(solve(crossprod(cbind(x@ZK %*% beta, x@Z1)))[1:r, 1:r]), 
              diag(Sigma))), nrow = x@P, ncol = r)
            colnames(alphaSE) <- rownames(alpha)
            rownames(alphaSE) <- colnames(x@Z0)
            ## Standard Errors of Cointegration relation(s)
            semat <- matrix(sqrt(diag(kronecker(
              solve(crossprod(x@RK[, -c(1:r)])), 
              solve(alpha %*% solve(Sigma) %*% t(alpha))))), ncol = r, byrow = TRUE)
            zemat <- matrix(0, nrow = r, ncol = r)
            betaSE <- rbind(zemat, semat)
            dimnames(betaSE) <- dimnames(beta)
            return(list(alphaSE = alphaSE, betaSE = betaSE))
})

## Create predicted price series based on VECM model
# y is a variable containing the cointegrated time series created from the model
# name is the variable you want to model, i.e. "forward_bid_price"
# cmodels is number of cointegration models to be used
# K-1 is lag order of the vecm
build_vecm = function(vecm, y, X, name, cmodels, K){
  
  # Find index from column name
  index = which(colnames(X)==name)
  
  # Get coefficients from model
  vecm_coef = coef(vecm$rlm)[,index]
  
  # Get lags for ecm
  price_diff = rbind(NA,diff(X))
  price_dl = matrix(NA, nrow=length(X[,1]), ncol=1)
  for (i in 1:(K-1)){
    price_diff = lag_1(price_diff)
    price_dl = cbind(price_dl,price_diff)
  }
  price_dl = price_dl[,-1] # remove first column of NA's
  
  # Vector of ones for constant term
  vec_ones = rep(1,length(X[,1]))
  
  variables = cbind(lag_1(y[,1:cmodels]), vec_ones, price_dl)
  prediction_d = t(vecm_coef%*%t(variables))
  
  prediction = prediction_d + lag_1(X[,index])
  return(prediction)
}

## mm_strategy calculates profits from market making strategy
# predicted_bid and predicted_ask are time series of predicted prices
# actual_bid and actual_ask are time series of actual prices
# delta is a measure of how competitive the market is
mm_strategy = function(predicted_bid, predicted_ask, actual_bid, actual_ask, delta){
  # delta is the amount you increase/decrease your quote by to make it more/less competitive
  # in theory the smaller the delta, the more profit you make
  
  bid_profit = array(NA,length(predicted_bid)-1)
  ask_profit = array(NA,length(predicted_bid)-1)
  
  for (i in 1:(length(predicted_bid)-1)){
    
    if (actual_bid[i]<predicted_bid[i+1]){ # predict market up, make bid competitive
      bid_quote = actual_bid[i]+delta
      bid_profit[i] = actual_bid[i+1]-bid_quote
    } else { # predict market down, make bid uncompetitive
      bid_quote = actual_bid-delta
      bid_profit[i] = 0
    }
    
    if (actual_ask[i]<predicted_ask[i+1]){ # predict market up, make ask uncompetitive
      ask_quote = actual_ask[i]+delta
      ask_profit[i] = 0
    } else { # predict market down, make bid uncompetitive
      ask_quote = actual_ask[i]-delta
      ask_profit[i] = ask_quote-actual_ask[i+1]
    }
  }
  
  profit_ar = bid_profit+ask_profit
  return(profit_ar)
}

## test_forecast calculates various statistics of a predicted price series
# fit_data is the predicted price series
# actual_data is the actual price series
# data_for_spread is the corresponding predicted price series used to calculate spread, i.e. bid if measuring ask data
# data_for_profit is the true time series of data_for_spread
test_forecast = function(fit_data,actual_data,data_for_spread,data_for_profit,timestamp){
  
  # Remove NA's from beggining of prediction
  fit_data = na.omit(fit_data)
  new_len = length(fit_data)
  actual_data = tail(actual_data,new_len)
  data_for_spread = tail(data_for_spread,new_len)
  data_for_profit = tail(data_for_profit,new_len)
  timestamp = tail(timestamp,new_len)
  
  # Calculate MAPE and MSE
  diff_1    = abs((fit_data-actual_data)/actual_data)
  diff_2    = fit_data-actual_data
  diffsqr = diff_2^2
  MSE     = mean(diffsqr)
  MAPE    = 100*mean(diff_1)
  
  # Calculate Up/Down Signals
  count_up = 0
  count_down = 0
  for(i in 2:length(actual_data)){
    # Count Correct Predicted Upward
    if(fit_data[i]-actual_data[i-1] > 0 & actual_data[i]-actual_data[i-1] > 0){
      count_up = count_up+1
    }
    # Count Correct Predicted Downward
    if(fit_data[i]-actual_data[i-1] < 0 & actual_data[i]-actual_data[i-1] < 0){
      count_down = count_down+1
    }
  }
  Prob = (count_down+count_up)/(length(actual_data)-1)
  
  # Calculate standardeviation (volitility) of data
  sd = sd(fit_data)
  
  # Identify bid from ask data
  if (fit_data[1]<data_for_spread[1]){
    fit_bid = fit_data
    fit_ask = data_for_spread
  } else {
    fit_bid = data_for_spread
    fit_ask = fit_data
  }
  if (actual_data[1]<data_for_profit[1]){
    actual_bid = actual_data
    actual_ask = data_for_profit
  } else {
    actual_bid = data_for_profit
    actual_ask = actual_data
  }
  
  # Calculate average, min and max spread of quotes
  spread = fit_ask - fit_bid
  avg_spread = mean(spread)
  min_spread = min(spread)
  max_spread = max(spread)
  
  # Calculate profit from market making strategy
  # set delta, which is a measure of how competitive the market is
  delta = 0.001
  # position is the size of the quotes posted in USD
  position = 1000
  profit_ar = mm_strategy(fit_bid, fit_ask, actual_bid, actual_ask, delta)
  # Calculate profit in USD based on position size
  profit = (sum(profit_ar)*position)/tail(fit_bid,1)

  
  return_object = list("MSE"=MSE, "MAPE"=MAPE, "Prob"=Prob, "sd"=sd, "average spread"=avg_spread, "min spread"=min_spread, "max spread"=max_spread, "profit"=profit)
  return (return_object)
}

## add_to_table adds the results from testing a model to the results table
# results is the table of results
# model_result is the model results to add
# type is the type of model
# name is the variable being tested
add_to_table = function(results, model_result, type, lag, name, count){
  
  results$Model.type[count] = type
  results$lags.used[count] = lag
  results$Modelled.on[count] = name
  results$MSE[count] = model_result$MSE
  results$MAPE[count] = model_result$MAPE
  results$sd[count] = model_result$sd
  results$avg.spread[count] = model_result$`average spread`
  results$min.spread[count] = model_result$`min spread`
  results$max.spread[count] = model_result$`max spread`
  results$profit[count] = model_result$profit
  results$strat.accuracy[count] = model_result$Prob
  
  return(results)
}