#####################################################################################
## JPM_tests builds models and tests your final chosen model on 20% of the data
#####################################################################################
## Writen by: Tom Parry, Thanachai Voraprasertsak, Wenqian Cheng and Maria Shevchenko

## Note that all models are built in the script but only the model for final testing
## is included in the graphs. This is so if you change your model choice then you do
## not have to rewrite the whole script

# Clear current workspace
rm(list = ls()) # Remove all variables
graphics.off() # Remove all plots

# Set directory to where your data is saved
setwd("/Users/ZerWinner/Desktop/JPMorgan Challenge") 

# Import library, you will need to install these packages first if you haven't already
library(readxl) # Read excel tables into R
library(tseries) # Univariant time series analysis package
library(vars) # Multivariant time series analysis package
library(urca) # For unit root and cointegration tests
library(dplyr) # For dataframe manipulation
library(smooth) # For moving average data smoothing
library(Matrix) # For complex matrix construction
library(data.table) # For dataframe manipulation
library(tsDyn) # Vector Error Correction Model package 

# JPM_test_functions.R contains all functions needed for this script
source("/Users/ZerWinner/Desktop/JPMorgan Challenge/JPM_test_functions.R")

######///////////////////////////////////////////////////////////######
# Part 1 Cleaning Data
######///////////////////////////////////////////////////////////######

# Import and format data
spot_bid_data = as.matrix(read_excel("USDKRW_1311-0902_Bid.xlsx"))
spot_ask_data = as.matrix(read_excel("USDKRW_1311-0902_Ask.xlsx"))
future_bid_data = as.matrix(read_excel("KUX6_KUZ6_KUF7_KUG7_1311-0902_Bid.xlsx"))
future_ask_data = as.matrix(read_excel("KUX6_KUZ6_KUF7_KUG7_1311-0902_Ask.xlsx"))

# Format data to one price point per minute based on median value for previous minute
df_spot_bid = format_data(spot_bid_data, "spot_bid_price")
df_spot_ask = format_data(spot_ask_data, "spot_ask_price")
df_future_bid = format_data(future_bid_data, "future_bid_price")
df_future_ask = format_data(future_ask_data, "future_ask_price")

# Only include days where number of data points exceeds threshold value, should be 390 points per day
threshold = 380
# This function takes about 20 secs to run
final_data = choose_data(df_spot_bid, df_spot_ask, df_future_bid, df_future_ask, threshold)
num_days = final_data$"num_days"
final_data = final_data$"final_data"
  
# Creat list to store results from testing
all_results = vector("list",num_days-1)

# option decides which method you use to test the models
# option = 1, models are trained on day i and test on day i+1
# option = 2, models are trained on first 70% of day i and test on last 30%
option = 2

validation_end = floor(num_days*0.8)
# Loop through the last 20% of the data to test our final model choices
for (z in validation_end:(num_days-1)){
  
  if (option == 1){
  
    # Train model on day i and test model on day i+1
    # Data to build models with
    today_data = final_data[[z]]
    # Data to test models with
    tomorrow_data = final_data[[z+1]]
    # Create matrix for time series analysis
    X_train = as.matrix(today_data[-1])
    timestamp_train = as.matrix(today_data[1])
    X_test = as.matrix(tomorrow_data[-1])
    timestamp_test = as.matrix(tomorrow_data[1])
    
  } else if (option == 2){
  
    # Train model on first 70% of day i, test on last 30% of day i
    # Split X up into train and test data
    X = as.matrix(final_data[[z]][-1])
    timestamp = as.matrix(final_data[[z]][1])
    train_prop = 70
    test_prop = 100-train_prop
    train_num = floor(dim(X)[1]*(train_prop/100))
    test_num = dim(X)[1]-train_num
    X_train = head(X,train_num)
    timestamp_train = head(timestamp,train_num)
    X_test = tail(X,test_num)
    timestamp_test = tail(timestamp,test_num)
  }
  
  ######///////////////////////////////////////////////////////////######
  # PART 2 AR Model, ARIMA Model
  ######///////////////////////////////////////////////////////////######
  
  # Take the difference of the data
  DX_train = diff(X_train)
  DX_test = diff(X_test)
  
  ### AR Model ###
  ################
  
  AR_models = vector("list",4)
  empty_model = matrix(NA,nrow=length(timestamp_test), ncol=4)
  
  # Create the matrix to store AIC value
  Info = matrix(NA, nrow = 4 , ncol = 1)
  aic_ar = matrix(NA, nrow = 16 , ncol = 1)
  seq = c(1,5,9,13)
  # test lags 1 to 4
  for (q in 1:4){
    AR_models[[q]] = empty_model
    for (i in 1:4){
      #Find out the optimal order for AR model
      #fit the model and produce parameters
      fitAR <- ar(DX_train[,i],aic = FALSE,order.max = q,method = "mle",na.action = na.fail) 
      coe <- array(data = fitAR$ar,dim=q)
      Info[i] = fitAR$aic
      #fit the model
      AR_M <- AR_predict(coe,DX_test[,i],q,X_test[,i])
      AR_M = append(NA,AR_M)
      AR_models[[q]][,i] = AR_M
    }
    aic_ar[seq[q]:(seq[q]+3)]=Info[c(3,4,1,2)]
  }
  
  ### ARIMA Model ###
  ###################
  
  # future_bid_coefficients
  model = ARIMA_coef(DX_train[,1],4)
  future_bid_ARIMA_coef = model$coef
  future_bid_aic        = model$aic
  # future_ask_coefficients
  model = ARIMA_coef(DX_train[,2],4)
  future_ask_ARIMA_coef = model$coef
  future_ask_aic        = model$aic
  # spot_bid_coefficients
  model = ARIMA_coef(DX_train[,3],4)
  spot_bid_ARIMA_coef = model$coef
  spot_bid_aic        = model$aic
  # spot_ask_coefficients
  model = ARIMA_coef(DX_train[,4],4)
  spot_ask_ARIMA_coef = model$coef
  spot_ask_aic        = model$aic
  
  # Combine AIC results
  new_AIC = rbind(spot_bid_aic,spot_ask_aic,future_bid_aic,future_ask_aic)
  A = append(new_AIC[,1],new_AIC[,2])
  B = append(new_AIC[,3],new_AIC[,4])
  aic_arima = as.matrix(append(A,B))
  
  # future_bid_predicted
  future_bid_ARIMA = ARIMA_predict(future_bid_ARIMA_coef,DX_test[,1],X_test[,1],4)
  # future_ask_predicted
  future_ask_ARIMA = ARIMA_predict(future_ask_ARIMA_coef,DX_test[,2],X_test[,2],4)
  # spot_bid_predicted
  spot_bid_ARIMA = ARIMA_predict(spot_bid_ARIMA_coef,DX_test[,3],X_test[,3],4)
  # spot_ask_predicted
  spot_ask_ARIMA = ARIMA_predict(spot_ask_ARIMA_coef,DX_test[,4],X_test[,4],4)
  
  ######///////////////////////////////////////////////////////////######
  # PART 4 VAR model 
  ######///////////////////////////////////////////////////////////######
  
  # Create Matrix to store AIC Value
  var_aic = matrix(NA,ncol=1,nrow=4)
  aic_var = matrix(NA,ncol=1,nrow=16)

  # Build VAR using different lags from 1 to 4
  VAR_models = vector("list",4)
  mv = vector("list",4)
  for (p in 1:4){
    
    m = VAR(DX_train, p=p)
    # Restrict model to only include significant terms
    m = restrict(m)
    mv[[p]] = m
    var_aic[p] = AIC(m)
    aic_var[(4*p-3):(4*p)]=rep(var_aic[p],4)
    
    # Predict price series using VAR model
    VAR_models[[p]] = build_VAR(m, DX_test, X_test, p)
    colnames(VAR_models[[p]])<- c("future_bid_price", "future_ask_price", "spot_bid_price", "spot_ask_price")
    
  }

  ######///////////////////////////////////////////////////////////######
  # PART 4 ECM Model 
  ######///////////////////////////////////////////////////////////######
  
  # model based on future_bid_price leads spot_bid_price
  model = ecm_model(X_train[,3],X_train[,1],4)
  lmfwd_spot_bid = model$coef
  spot_bid_ecm   = model$aic
  
  # model based on future_ask_price leads spot_ask_price
  model = ecm_model(X_train[,4],X_train[,2],4)
  lmfwd_spot_ask = model$coef
  spot_ask_ecm   = model$aic
  
  # model based on spot_bid_price leads future_bid_price
  model = ecm_model(X_train[,1],X_train[,3],4)
  lmspot_fwd_bid = model$coef
  fwd_bid_ecm   = model$aic
  
  # model based on spot_ask_price leads future_ask_price
  model = ecm_model(X_train[,2],X_train[,4],4)
  lmspot_fwd_ask = model$coef
  fwd_ask_ecm   = model$aic
  
  # Combine AIC results
  new_AIC = rbind(spot_bid_ecm,spot_ask_ecm,fwd_bid_ecm,fwd_ask_ecm)
  A = append(new_AIC[,1],new_AIC[,2])
  B = append(new_AIC[,3],new_AIC[,4])
  aic_ecm = as.matrix(append(A,B))
  
  # Find Predicted Value
  spot_bid_ecm_predict=ecm_predict(lmfwd_spot_bid,X_test[,3],X_test[,1],4)
  spot_ask_ecm_predict=ecm_predict(lmfwd_spot_ask,X_test[,4],X_test[,2],4)
  future_bid_ecm_predict=ecm_predict(lmspot_fwd_bid,X_test[,1],X_test[,3],4)
  future_ask_ecm_predict=ecm_predict(lmspot_fwd_ask,X_test[,2],X_test[,4],4)
  
  ######///////////////////////////////////////////////////////////######
  # PART 5 VECM Model
  ######///////////////////////////////////////////////////////////######
  
  # Build VECM using different lags 1 to 4
  vecm_models = vector("list",4)
  # Create Matrix to store AIC Value
  vecm_aic = matrix(NA,ncol=1,nrow=4)
  aic_vecm = matrix(NA,ncol=1,nrow=16)
  
  for (K in 2:5){
    
    # Model VECM
    m = ca.jo(X_train, type="trace", ecdet="none", K=K, spec="transitory")
    y_test = as.matrix(X_test)%*%as.matrix(m@V)
    cmodels = 1
    vecm = cajorls(m, r=cmodels)
    
    # Calculate AIC Value
    vecm_est = VECM(X_train, lag = K-1, include = c("none"), estim = c("ML"))
    vecm_aic[K-1] = AIC(vecm_est)
    aic_vecm[(4*(K-1)-3):(4*(K-1))]=rep(vecm_aic[K-1],4)
    vecm_prediction_temp = data.frame("timestamp"=timestamp_test,"future_bid_prediction"=0,"future_ask_prediction"=0,"spot_bid_prediction"=0,"spot_ask_prediction"=0)
    
    # Calculate predictions
    for (i in 1:4){
      vecm_prediction_temp[i+1] = build_vecm(vecm, y_test, X_test, colnames(X_test)[i], cmodels, K)
    }
    vecm_models[[K-1]] = vecm_prediction_temp
  }
  
  ######///////////////////////////////////////////////////////////######
  # PART 6 Model Comparison
  ######///////////////////////////////////////////////////////////######
  
  # Combine All AIC results
  aic = rbind(as.matrix(rep(NA,4)),aic_ar,aic_arima,aic_var,aic_ecm,aic_vecm)
  
  # Create dummy models using last price as predicted price
  dummy_strat_future_bid = data.frame(timestamp_test, lag_1(X_test[,"future_bid_price"]))
  dummy_strat_future_ask = data.frame(timestamp_test, lag_1(X_test[,"future_ask_price"]))
  dummy_strat_spot_bid = data.frame(timestamp_test, lag_1(X_test[,"spot_bid_price"]))
  dummy_strat_spot_ask = data.frame(timestamp_test, lag_1(X_test[,"spot_ask_price"]))
  
  # Create data frame to store results from model comparrison
  # Specify number of models built to pre set dataframe size
  num_models = 84
  results = data.frame("Model type"=numeric(num_models), "lags used"=numeric(num_models), "Modelled on"=numeric(num_models),"AIC" = numeric(num_models),"MSE"=numeric(num_models), "MAPE"=numeric(num_models), "sd"=numeric(num_models), "avg spread"=numeric(num_models), "min spread"=numeric(num_models), "max spread"=numeric(num_models), "profit"=numeric(num_models), "strat accuracy"=numeric(num_models))
  count = 1
  
  #=================
  # DUMMY MODEL TEST
  #=================
  
  # Test spot_bid model
  spot_bid_result = test_forecast(dummy_strat_spot_bid[,2],X_test[,3],dummy_strat_spot_ask[,2],X_test[,4],timestamp_test)
  results = add_to_table(results, spot_bid_result, "Dummy", i, "spot_bid", count)
  count = count+1
  
  # Test spot_ask model
  spot_ask_result = test_forecast(dummy_strat_spot_ask[,2],X_test[,4],dummy_strat_spot_bid[,2],X_test[,3],timestamp_test)
  results = add_to_table(results, spot_ask_result, "Dummy", i, "spot_ask", count)
  count = count+1
  
  # Test future_bid model
  future_bid_result = test_forecast(dummy_strat_future_bid[,2],X_test[,1],dummy_strat_future_ask[,2],X_test[,2],timestamp_test)
  results = add_to_table(results, future_bid_result, "Dummy", i, "future_bid", count)
  count = count+1
  
  # Test future_ask_model
  future_ask_result = test_forecast(dummy_strat_future_ask[,2],X_test[,2],dummy_strat_future_bid[,2],X_test[,1],timestamp_test)
  results = add_to_table(results, future_ask_result, "Dummy", i, "future_ask", count)
  count = count+1
  
  #================
  # AR MODEL TEST
  #================
  
  for (i in 1:4){
    
    # Test spot_bid model
    spot_bid_result = test_forecast(AR_models[[i]][,3],X_test[,3],AR_models[[i]][,4],X_test[,4],timestamp_test)
    results = add_to_table(results, spot_bid_result, "AR", i, "spot_bid", count)
    count = count+1
    
    # Test spot_ask model
    spot_ask_result = test_forecast(AR_models[[i]][,4],X_test[,4],AR_models[[i]][,3],X_test[,3],timestamp_test)
    results = add_to_table(results, spot_ask_result, "AR", i, "spot_ask", count)
    count = count+1
    
    # Test future_bid model
    future_bid_result = test_forecast(AR_models[[i]][,1],X_test[,1],AR_models[[i]][,2],X_test[,2],timestamp_test)
    results = add_to_table(results, future_bid_result, "AR", i, "future_bid", count)
    count = count+1
    
    # Test future_ask model
    future_ask_result = test_forecast(AR_models[[i]][,2],X_test[,2],AR_models[[i]][,1],X_test[,1],timestamp_test)
    results = add_to_table(results, future_ask_result, "AR", i, "future_ask", count)
    count = count+1
  }
  
  #=================
  # ARIMA MODEL TEST
  #=================
  
  for (i in 1:4){
    
    # Test spot_bid model
    spot_bid_result = test_forecast(spot_bid_ARIMA[,i],X_test[,3],spot_ask_ARIMA[,i],X_test[,4],timestamp_test)
    results = add_to_table(results, spot_bid_result, "ARIMA", i, "spot_bid", count)
    count = count+1
    
    # Test spot_ask model
    spot_ask_result = test_forecast(spot_ask_ARIMA[,i],X_test[,4],spot_bid_ARIMA[,i],X_test[,3],timestamp_test)
    results = add_to_table(results, spot_ask_result, "ARIMA", i, "spot_ask", count)
    count = count+1
    
    # Test future_bid model
    future_bid_result = test_forecast(future_bid_ARIMA[,i],X_test[,1],future_ask_ARIMA[,i],X_test[,2],timestamp_test)
    results = add_to_table(results, future_bid_result, "ARIMA", i, "future_bid", count)
    count = count+1
    
    # Test future_ask model
    future_ask_result = test_forecast(future_ask_ARIMA[,i],X_test[,2],future_bid_ARIMA[,i],X_test[,1],timestamp_test)
    results = add_to_table(results, future_ask_result, "ARIMA", i, "future_ask", count)
    count = count+1
  }
  
  #================
  # VAR MODEL TEST
  #================
  
  # Test models of lag 1 and 4
  for (i in 1:4){
    
    # Test spot_bid model
    spot_bid_result = test_forecast(VAR_models[[i]][,3],X_test[,3],VAR_models[[i]][,4], X_test[,4],timestamp_test)
    results = add_to_table(results, spot_bid_result, "VAR", i, "spot_bid", count)
    count = count+1
    
    # Test spot_ask model
    spot_ask_result = test_forecast(VAR_models[[i]][,4],X_test[,4],VAR_models[[i]][,3], X_test[,3],timestamp_test)
    results = add_to_table(results, spot_ask_result, "VAR", i, "spot_ask", count)
    count = count+1
    
    # Test future_bid model
    future_bid_result = test_forecast(VAR_models[[i]][,1],X_test[,1],VAR_models[[i]][,2], X_test[,2],timestamp_test)
    results = add_to_table(results, future_bid_result, "VAR", i, "future_bid", count)
    count = count+1
    
    # Test future_ask model
    future_ask_result = test_forecast(VAR_models[[i]][,2],X_test[,2],VAR_models[[i]][,1], X_test[,1], timestamp_test)
    results = add_to_table(results, future_ask_result, "VAR", i, "future_ask", count)
    count = count+1
  }
  
  #================
  # ECM MODEL TEST
  #================
  
  # Test models of lag 1 to 4
  for (i in 1:4){
    
    # Test spot_bid model
    spot_bid_result = test_forecast(spot_bid_ecm_predict[,i],X_test[,3],spot_ask_ecm_predict[,i],X_test[,4],timestamp_test)
    results = add_to_table(results, spot_bid_result, "ECM", i, "spot_bid", count)
    count = count+1
    
    # Test spot_ask model
    spot_ask_result = test_forecast(spot_ask_ecm_predict[,i],X_test[,4],spot_bid_ecm_predict[,i],X_test[,3],timestamp_test)
    results = add_to_table(results, spot_ask_result, "ECM", i, "spot_ask", count)
    count = count+1
    
    # Test future_bid model
    future_bid_result = test_forecast(future_bid_ecm_predict[,i],X_test[,1],future_ask_ecm_predict[,i],X_test[,2],timestamp_test)
    results = add_to_table(results, future_bid_result, "ECM", i, "future_bid", count)
    count = count+1
    
    # Test future_ask model
    future_ask_result = test_forecast(future_ask_ecm_predict[,i],X_test[,2],future_bid_ecm_predict[,i],X_test[,1],timestamp_test)
    results = add_to_table(results, future_ask_result, "ECM", i, "future_ask", count)
    count = count+1
  }
  
  #================
  # VECM MODEL TEST
  #================
  
  # Test models of lag 1 to 4
  for (i in 1:4){
    
    # Test spot_bid model
    spot_bid_result = test_forecast(vecm_models[[i]][,4],X_test[,3],vecm_models[[i]][,5],X_test[,4],timestamp_test)
    results = add_to_table(results, spot_bid_result, "VECM", i, "spot_bid", count)
    count = count+1
    
    # Test spot_ask model
    spot_ask_result = test_forecast(vecm_models[[i]][,5],X_test[,4],vecm_models[[i]][,4],X_test[,3],timestamp_test)
    results = add_to_table(results, spot_ask_result, "VECM", i, "spot_ask", count)
    count = count+1
    
    # Test future_bid model
    future_bid_result = test_forecast(vecm_models[[i]][,2],X_test[,1],vecm_models[[i]][,3],X_test[,2],timestamp_test)
    results = add_to_table(results, future_bid_result, "VECM", i, "future_bid", count)
    count = count+1
    
    # Test future_ask model
    future_ask_result = test_forecast(vecm_models[[i]][,3],X_test[,2],vecm_models[[i]][,2],X_test[,1],timestamp_test)
    results = add_to_table(results, future_ask_result, "VECM", i, "future_ask", count)
    count = count+1
  }
  
  # Add AIC to table
  results$AIC = aic
  
  all_results[[z]] = results
}
 
# Average the results in each days table into one table 
avg_results = rbindlist(all_results)[,lapply(.SD,mean), c("Model.type","lags.used","Modelled.on")]
View(avg_results)

par(mfrow=c(1,2))
# Plot graphs of lag vs. MSE for each model type
for (i in c("spot_bid", "spot_ask", "future_bid", "future_ask")){
  
  Dummy_results = avg_results[((avg_results$Modelled.on==i)&(avg_results$Model.type=="Dummy")),]
  #AR_results = avg_results[((avg_results$Modelled.on==i)&(avg_results$Model.type=="AR")),]
  #ARIMA_results = avg_results[((avg_results$Modelled.on==i)&(avg_results$Model.type=="ARIMA")),]
  VAR_results = avg_results[((avg_results$Modelled.on==i)&(avg_results$Model.type=="VAR")),]
  #ECM_results = avg_results[((avg_results$Modelled.on==i)&(avg_results$Model.type=="ECM")),]
  #VECM_results = avg_results[((avg_results$Modelled.on==i)&(avg_results$Model.type=="VECM")),]
  #plot_limits = avg_results[(avg_results$Modelled.on==i),]$MSE
  #ymax = max(plot_limits)
  #ymin = min(plot_limits)
  lags=c(1,2,3,4)
  Dummy_plot = rep(Dummy_results$MSE,4)
  ymax = max(append(Dummy_plot,VAR_results$MSE))
  ymin = min(append(Dummy_plot,VAR_results$MSE))
  plot(lags, VAR_results$MSE, type="l", col="green", xlab="Lags", ylab="MSE", main=paste(i," model MSE",sep=""), xaxt='n', ylim=c(ymin,ymax))
  axis(1, at=1:4)
  #lines(lags, ARIMA_results$MSE, col="cyan")
  #lines(lags, VAR_results$MSE, col="green")
  #lines(lags, ECM_results$MSE, col="orange")
  #lines(lags, VECM_results$MSE, col="blue")
  lines(lags, Dummy_plot, col="black", lty=2)
  legend(x=2.5, y=ymax-0.001,c("VAR", "Dummy"), fill=c("green", "black"), cex=0.5, lty=c(1,2))
}

# Plot graphs of lag vs. profit for each model type
for (i in c("spot", "future")){
  
  name = paste(i,"_bid", sep="")
  Dummy_results = avg_results[((avg_results$Modelled.on==name)&(avg_results$Model.type=="Dummy")),]
  #AR_results = avg_results[((avg_results$Modelled.on==name)&(avg_results$Model.type=="AR")),]
  #ARIMA_results = avg_results[((avg_results$Modelled.on==name)&(avg_results$Model.type=="ARIMA")),]
  VAR_results = avg_results[((avg_results$Modelled.on==name)&(avg_results$Model.type=="VAR")),]
  #ECM_results = avg_results[((avg_results$Modelled.on==name)&(avg_results$Model.type=="ECM")),]
  #VECM_results = avg_results[((avg_results$Modelled.on==name)&(avg_results$Model.type=="VECM")),]
  #plot_limits = avg_results[(avg_results$Modelled.on==name),]$profit
  lags=c(1,2,3,4)
  Dummy_plot = rep(Dummy_results$profit,4)
  ymax = max(append(Dummy_plot,VAR_results$profit))
  ymin = min(append(Dummy_plot,VAR_results$profit))
  plot(lags, VAR_results$profit, type="l", col="green", xlab="Lags", ylab="Profit (USD) per trading period", main=paste(i," model profit",sep=""), xaxt='n', ylim=c(ymin,ymax))
  axis(1, at=1:4)
  #lines(lags, ARIMA_results$profit, col="cyan")
  #lines(lags, VAR_results$profit, col="green")
  #lines(lags, ECM_results$profit, col="orange")
  #lines(lags, VECM_results$profit, col="blue")
  lines(lags, Dummy_plot, col="black", lty=2)
  legend("topleft",c("VAR", "Dummy"), fill=c("green", "black"), cex=0.5, lty=c(1,2))
}

# For the final day, plot predicted series vs. actual series for visual comparison
timestamp_test = strptime(timestamp_test, format="%F %T")
par(mfrow=c(1,1))
plot(timestamp_test, X_test[,3],type='l', ylab="Price (Won)", xlab="Timestamp", col="red", main="Spot bid predicted vs. actual prices")
lines(timestamp_test, VAR_models[[3]][,3], col="green")
legend("topleft",c("Actual price","Predicted price"), fill=c("red","green"), cex=0.5)


# Export the results in excel table
# write.table(avg_results, "C:/Users/REG USER/Dropbox/JP Morgan Challenge/New Data/avg_result.csv",sep=",")

