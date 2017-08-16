#####################################################################################################
## JPM_main runs anaylsis for one days worth of data, it produces plots and prints some model outputs
#####################################################################################################
## Writen by: Tom Parry, Thanachai Voraprasertsak, Wenqian Cheng and Maria Shevchenko


# Clear current workspace
rm(list = ls()) # Remove all variables
graphics.off() # Remove all plots

# Set directory to where your data is saved
setwd("/Users/ZerWinner/Desktop/JPMorgan Challenge")

# Import library, you will need to install these packages first if you haven't already
library(readxl) # read in excel files
library(tseries) # Univariant time series analysis package
library(vars) # Multivariant time series analysis package
library(urca) # For unit root and cointegration tests
library(dplyr) # For dataframe manipulation
library(smooth) # For moving average data smoothing
library(Matrix) # For complex matrix construction
library(tsDyn) # Vector Error Correction Model package 

# JPM_main_functions.R contains all functions needed for this script
source("/Users/ZerWinner/Desktop/JPMorgan Challenge/JPM_main_functions.R")

######///////////////////////////////////////////////////////////######
# Part 1 Cleaning Data
######///////////////////////////////////////////////////////////######

# Import and format data
spot_bid_data = as.matrix(read_excel("USDKRW_0112_Bid.xlsx"))
spot_ask_data = as.matrix(read_excel("USDKRW_0112_Ask.xlsx"))
future_bid_data = as.matrix(read_excel("KUZ6_0112_Bid.xlsx"))
future_ask_data = as.matrix(read_excel("KUZ6_0112_Ask.xlsx"))

df_spot_bid = format_data(spot_bid_data, "spot_bid_price")
df_spot_ask = format_data(spot_ask_data, "spot_ask_price")
df_future_bid = format_data(future_bid_data, "future_bid_price")
df_future_ask = format_data(future_ask_data, "future_ask_price")

# Merge data to ensure timestamps start and end at same time
df_spot_price = merge(df_spot_bid, df_spot_ask, by="timestamp")
df_future_price = merge(df_future_bid, df_future_ask, by="timestamp")
df_price = merge(df_spot_price, df_future_price, by="timestamp")

# add 1 minute to each price point so price reperesents median price for previous minute
df_price[,1] = df_price[,1]+60

# Split prices up into vectors for analysis
timestamp = strptime(as.matrix(df_price)[,1], format="%F %T")
spot_bid_price = as.numeric(as.matrix(df_price)[,2])
spot_ask_price = as.numeric(as.matrix(df_price)[,3])
future_bid_price = as.numeric(as.matrix(df_price)[,4])
future_ask_price = as.numeric(as.matrix(df_price)[,5])

# Plot evenly spaced data
plot(timestamp, spot_bid_price, type='l', ylab="Price (Won)", xlab="Timestamp", col="blue", main="Plot of transformed (even) data")
lines(timestamp, spot_ask_price, col="red")
lines(timestamp, future_bid_price, col="green")
lines(timestamp, future_ask_price, col="yellow")
legend("topright",c("USDKRW Bid price","USDKW Ask price", "KUZ6 Bid price", "KUZ6 Ask price"), fill=c("blue","red", "green", "yellow"), cex=0.5)

# Create matrix for time series analysis
X = cbind(future_bid_price, future_ask_price, spot_bid_price, spot_ask_price)
log_X = log(X) # log prices

# Split X up into train and test data
train_prop = 70
test_prop = 100-train_prop
train_num = floor(dim(X)[1]*(train_prop/100))
test_num = dim(X)[1]-train_num
X_train = head(X,train_num)
DX_train = diff(X_train)
timestamp_train = head(timestamp,train_num)
X_test = tail(X,test_num)
DX_test = diff(X_test)
timestamp_test = tail(timestamp,test_num)

######////////////////////////////////////////////////////////////////////######
# PART 2 Test for stationary/non-stationary processes and finding leading series
######////////////////////////////////////////////////////////////////////######

# Difference the data
X_diff = diff(X)
DX_train = diff(X_train)
DX_test = diff(X_test)

# Plot evenly spaced data
par(mfrow=c(2,1))
plot(timestamp[-1], X_diff[,1], type='l', ylab="Differneced Price (Won)", xlab="Timestamp", col="blue", main="Plot of spot bid differenced data")
plot(timestamp[-1], X_diff[,2], type='l', ylab="Differneced Price (Won)", xlab="Timestamp", col="red", main="Plot of spot ask differenced data")
par(mfrow=c(2,1))
plot(timestamp[-1], X_diff[,3], type='l', ylab="Differneced Price (Won)", xlab="Timestamp", col="green", main="Plot of future bid differenced data")
plot(timestamp[-1], X_diff[,4], type='l', ylab="Differneced Price (Won)", xlab="Timestamp", col="cyan", main="Plot of future ask differenced data")

# Take the difference and calculate the auto correlation to see if it's stationary
acf(X_diff,NULL,type = c("correlation"),plot = TRUE,na.action = na.pass,demean = TRUE)

# Granger Casuality and find covariance and correlation
fit_train = VAR(DX_train, ic="AIC") 

# Plot impulse response functions
irfX=irf(fit_train) 
plot(irfX)
print("We can observe from the graphs that future is leading time series")

# create matrix of p-value of granger-causality test to support IRF
results_granger_causality <- matrix(NA, ncol=4, nrow=4)
rownames(results_granger_causality) <- c("future_bid_price", "future_ask_price", "spot_bid_price", "spot_ask_price")
colnames(results_granger_causality) <- c("future_bid_price", "future_ask_price", "spot_bid_price", "spot_ask_price")
# perform granger-causality test 
for (i in 1:4){
  for (j in 1:4){
    if(i==j){}
    else {
      MM=VAR(cbind(DX_train[,i],DX_train[,j]),ic="AIC")
      test = causality(MM)
      results_granger_causality[i,j]=test$Granger$p.value
    }}}
print("Formal test also support the idea that future is leading time series (p-values are below):")
results_granger_causality

######////////////////////////######
# PART 3 AR Model, ARIMA Model
######////////////////////////######

### AR Model ###
################

AR_models = vector("list",4)
empty_model = matrix(NA,nrow=length(timestamp_test), ncol=4)

# Create the matrix to store AIC value
Info = matrix(NA, nrow = 4 , ncol = 1)
aic_ar = matrix(NA, nrow = 16 , ncol = 1)
seq = c(1,5,9,13)
# Test lags 1 to 4
for (q in 1:4){
  AR_models[[q]] = empty_model
  for (i in 1:4){
    #Find out the optimal order for AR model
    #Fit the model using maximum likelihood method and return coefficients
    fitAR <- ar(DX_train[,i],aic = FALSE,order.max = q,method = "mle",na.action = na.fail) 
    coe <- array(data = fitAR$ar,dim=q)
    Info[i] = fitAR$aic
    #Apply the optimal coefficients into the AR prediction model
    AR_M <- AR_predict(coe,DX_test[,i],q,X_test[,i])
    AR_M = append(NA,AR_M)
    AR_models[[q]][,i] = AR_M
  }
  aic_ar[seq[q]:(seq[q]+3)]=Info[c(3,4,1,2)]
}


### ARIMA Model 1st model ###
###################
#ARIMA_predict = function(coe,test_real,p,q,resdidual,X_test){
#  fitdata <- matrix(data =NA,nrow = length(test_real),ncol=1)
#  fitdata[1:q,1] <- NA
#  for(j in (p+1):length(test_real)){
#    fitdata[j,1] = sum((coe[1:p]*rev(test_real[(j-p):(j-1)]))) + 
#      sum(coe[(p+1):(p+q)]*rev(append(res[(j-q+j-p):(j-1)],array((fitdata[p,1]-test_real(p)):(fitdata[j-1,1]-test_real(j-1)),dim = j-p))))  
#    
#  } 
#  lagX1 = lag_1(X_test)
#  lagX1 = lagX1[2:length(lagX1)]
#  fitdata = fitdata+lagX1
#  return (fitdata)
#}


#Find out the optimal orders for ARIMA model
#for (k in 1:1){
#  aic <- matrix(data = NA, nrow = 2, ncol = 2)
#  for (i in 1:4)  {
#    fitARIMA <- arima(DX_train[, k], order = c(i,0,i), method = "ML")
#    aic[i,i] = fitARIMA$aic
#    }
#  Ind <- arrayInd(which.max(aic), dim(aic))
#  p <- Ind[1]
#  q <- Ind[2]
  
  #fit the model and produce parameters
#  fitARIMA <- arima(DX_train[, k], order = c(p,0,q), method = "ML")
#  summary(fitARIMA)
#  M = fitARIMA$aic
#  coe = fitARIMA$coef
#  res = fitARIMA$residuals
  #plot residual to prove it's white noise
#  plot(res,type = "l") 
  #acf(res)
  
  # Use the ARIMA_Predict model
#  ARIMA_M <- ARIMA_predict(coe,DX_test[,k],p,q,res,X_test) 
#  Aplot(ARIMA_M,type = "l", col = 'red')
#  par (new = TRUE)
#  plot(X_test[,k],type = "l", col = 'black')
#} 



### ARIMA Model 2nd model ###
###################
# Future_bid_coefficients
# Write a function to return ARIMA model and store coeffiencts and AIC
model = ARIMA_coef(DX_train[,1],4)
future_bid_ARIMA_coef = model$coef
future_bid_aic        = model$aic
# Future_ask_coefficients
model = ARIMA_coef(DX_train[,2],4)
future_ask_ARIMA_coef = model$coef
future_ask_aic        = model$aic
# Spot_bid_coefficients
model = ARIMA_coef(DX_train[,3],4)
spot_bid_ARIMA_coef = model$coef
spot_bid_aic        = model$aic
# Spot_ask_coefficients
model = ARIMA_coef(DX_train[,4],4)
spot_ask_ARIMA_coef = model$coef
spot_ask_aic        = model$aic

# Combine AIC results
new_AIC = rbind(spot_bid_aic,spot_ask_aic,future_bid_aic,future_ask_aic)
A = append(new_AIC[,1],new_AIC[,2])
B = append(new_AIC[,3],new_AIC[,4])
aic_arima = as.matrix(append(A,B))

#Apply the prediction function for ARIMA Model
# Future_bid_predicted
future_bid_ARIMA = rbind(NA,ARIMA_predict(future_bid_ARIMA_coef,DX_test[,1],X_test[,1],4))
# Future_ask_predicted
future_ask_ARIMA = rbind(NA,ARIMA_predict(future_ask_ARIMA_coef,DX_test[,2],X_test[,2],4))
# Spot_bid_predicted
spot_bid_ARIMA = rbind(NA,ARIMA_predict(spot_bid_ARIMA_coef,DX_test[,3],X_test[,3],4))
# Spot_ask_predicted
spot_ask_ARIMA = rbind(NA,ARIMA_predict(spot_ask_ARIMA_coef,DX_test[,4],X_test[,4],4))


######///////////////////////////////////////////////////////////######
# PART 4 VAR Models
######///////////////////////////////////////////////////////////######

# Create Matrix to store AIC Value
var_aic = matrix(NA,ncol=1,nrow=4)
aic_var = matrix(NA,ncol=1,nrow=16)

# Build VAR using different lags from 1 to 4
VAR_models = vector("list",4)

for (p in 1:4){

  m = VAR(DX_train, p=p)
  # summary(m)
  # Restrict model to only include significant terms
  m = restrict(m)
  # serial.test(m, lags.pt=6, type="PT.adjusted")
  var_aic[p] = AIC(m)
  aic_var[(4*p-3):(4*p)]=rep(var_aic[p],4)
  
  # Predict price series using VAR model
  VAR_models[[p]] = build_VAR(m, DX_test, X_test, p)
  colnames(VAR_models[[p]])<- c("future_bid_price", "future_ask_price", "spot_bid_price", "spot_ask_price")
  
}

######///////////////////////////////////////////////////////////######
# PART 5 ECM Models 
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
# PART 6 VECM Models - Johansen's Likelihood Method
######///////////////////////////////////////////////////////////######

# Decide lag order to be used for initial analysis (lag used is 1 less than this numner, aka 4 means 3 lags used)
K = 5

# Fit data with model and apply the trace test

# Johansen test with trace test statistic
m1 = ca.jo(X, type="trace", ecdet="none", K=K, spec="transitory")

print("The critical values for this test are:")
m1@cval
print("The test statistics for this test are:")
m1@teststat
print("We cannot reject the null hypothesis for r=3")
print("This suggests there are 3 cointegrated relationships")

# Johansen test with eigen test statistic
m2 = ca.jo(X, type="eigen", ecdet="none", K=K, spec="transitory")

print("The critical values for this test are:")
m2@cval
print("The test statistics for this test are:")
m2@teststat
print("We cannot reject the null hypothesis for r=3")
print("This also suggests there are 3 cointegrated relationships")

# Calculate new timeseries using cointegrated weights
y_train = as.matrix(X_train)%*%as.matrix(m1@V)

# Plot of each new series along with their respective acf
par(mfrow=c(2,2))
plot(timestamp_train, y_train[,1], type="l", ylab="1st series", xlab="Timestamp", col="blue")
acf(y_train[,1], main="")
plot(timestamp_train, y_train[,2], type="l", ylab="2nd series", xlab="Timestamp", col="blue")
acf(y_train[,2], main="")
par(mfrow=c(2,2))
plot(timestamp_train, y_train[,3], type="l", ylab="3rd series", xlab="Timestamp", col="blue")
acf(y_train[,3], main="")
plot(timestamp_train, y_train[,4], type="l", ylab="4th series", xlab="Timestamp", col="blue")
acf(y_train[,4], main="")

print("Upon inspection of the graphs, only the first one appears to have no significant higher lags, suggesting only 1 cointegrated relationship")
print("Infact upon further analysis this is the case, including only 1 relationship yields the most accurate VECM")

# Decide how many cointegration relationships to include
cmodels = 1

# Build VECM including lags 1 to 4
vecm_models = vector("list",4)
# Create Matrix to store AIC Value
vecm_aic = matrix(NA,ncol=1,nrow=4)
aic_vecm = matrix(NA,ncol=1,nrow=16)

for (K in 2:5){
  
  # Model VECM
  m = ca.jo(X_train, type="trace", ecdet="none", K=K, spec="transitory")
  y_test = as.matrix(X_test)%*%as.matrix(m@V)
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

# Combine All AIC results
aic = rbind(as.matrix(rep(NA,4)),aic_ar,aic_arima,aic_var,aic_ecm,aic_vecm)

######///////////////////////////////////////////////////////////######
# PART 7 Model Comparison
######///////////////////////////////////////////////////////////######

# Create dummy models using last price as predicted price
dummy_strat_future_bid = data.frame(timestamp_test, lag_1(X_test[,"future_bid_price"]))
dummy_strat_future_ask = data.frame(timestamp_test, lag_1(X_test[,"future_ask_price"]))
dummy_strat_spot_bid = data.frame(timestamp_test, lag_1(X_test[,"spot_bid_price"]))
dummy_strat_spot_ask = data.frame(timestamp_test, lag_1(X_test[,"spot_ask_price"]))

# Create data frame to store results from model comparrison
# Specify number of models built to pre set dataframe size
num_models = 84
results = data.frame("Model type"=numeric(num_models), "lags used"=numeric(num_models), "Modelled on"=numeric(num_models), "MSE"=numeric(num_models), "MAPE"=numeric(num_models), "sd"=numeric(num_models), "avg spread"=numeric(num_models), "min spread"=numeric(num_models), "max spread"=numeric(num_models), "profit"=numeric(num_models), "strat accuracy"=numeric(num_models))
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

#================
# ARIMA MODEL TEST
#================

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

View(results)

par(mfrow=c(1,2))
# Plot graphs of lag vs. MSE for each model type
for (i in c("spot_bid", "spot_ask", "future_bid", "future_ask")){
  
  Dummy_results = results[((results$Modelled.on==i)&(results$Model.type=="Dummy")),]
  AR_results = results[((results$Modelled.on==i)&(results$Model.type=="AR")),]
  ARIMA_results = results[((results$Modelled.on==i)&(results$Model.type=="ARIMA")),]
  VAR_results = results[((results$Modelled.on==i)&(results$Model.type=="VAR")),]
  ECM_results = results[((results$Modelled.on==i)&(results$Model.type=="ECM")),]
  VECM_results = results[((results$Modelled.on==i)&(results$Model.type=="VECM")),]
  plot_limits = results[(results$Modelled.on==i),]$MSE
  ymax = max(plot_limits)
  ymin = min(plot_limits)
  lags=c(1,2,3,4)
  plot(lags, AR_results$MSE, type="l", col="red", xlab="Lags", ylab="MSE", main=paste(i,"model MSE"), xaxt='n', ylim=c(ymin,ymax))
  axis(1, at=1:4)
  lines(lags, ARIMA_results$MSE, col="cyan")
  lines(lags, VAR_results$MSE, col="green")
  lines(lags, ECM_results$MSE, col="orange")
  lines(lags, VECM_results$MSE, col="blue")
  Dummy_plot = rep(Dummy_results$MSE,4)
  lines(lags, Dummy_plot, col="black", lty=2)
  legend("topleft",c("AR","ARIMA" ,"VAR", "ECM","VECM", "Dummy"), fill=c("red","cyan","green","orange","blue","black"), cex=0.5, lty=c(1,1,1,1,1,2))
}

# Plot graphs of lag vs. profit for each model type
for (i in c("spot", "future")){
  
  name = paste(i,"_bid", sep="")
  Dummy_results = results[((results$Modelled.on==name)&(results$Model.type=="Dummy")),]
  AR_results = results[((results$Modelled.on==name)&(results$Model.type=="AR")),]
  ARIMA_results = results[((results$Modelled.on==name)&(results$Model.type=="ARIMA")),]
  VAR_results = results[((results$Modelled.on==name)&(results$Model.type=="VAR")),]
  ECM_results = results[((results$Modelled.on==name)&(results$Model.type=="ECM")),]
  VECM_results = results[((results$Modelled.on==name)&(results$Model.type=="VECM")),]
  plot_limits = results[(results$Modelled.on==name),]$profit
  ymax = max(plot_limits)
  ymin = min(plot_limits)
  lags=c(1,2,3,4)
  plot(lags, AR_results$profit, type="l", col="red", xlab="Lags", ylab="Profit (USD) per trading period", main=paste(i,"model profit"), xaxt='n', ylim=c(ymin,ymax))
  axis(1, at=1:4)
  lines(lags, ARIMA_results$profit, col="cyan")
  lines(lags, VAR_results$profit, col="green")
  lines(lags, ECM_results$profit, col="orange")
  lines(lags, VECM_results$profit, col="blue")
  Dummy_plot = rep(Dummy_results$profit,4)
  lines(lags, Dummy_plot, col="black", lty=2)
  legend("topleft",c("AR","ARIMA", "VAR", "ECM","VECM", "Dummy"), fill=c("red","cyan","green","orange","blue","black"), cex=0.5, lty=c(1,1,1,1,1,2))
}