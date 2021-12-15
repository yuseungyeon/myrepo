#setwd("C:\\Users\\km960\\Desktop\\?���?")
setwd('/Volumes/GoogleDrive/My Drive/Research-student/Kyeongmin/program')

rm(list=ls())
source('mLRD-library.R')
source("HAR_GARCH_outlier1.R")
#source("HAR-oxfordRV-outlier1.R")
library(forecast)
library(tsoutliers)
library(rugarch)
library(MASS)
result1 = result2 = result3 = result4 = data.frame()


for (i in 1:500) {
  
  mu=2.2;
  ######################### GENERATE DATA #########################
  
  y = farima_sim(300, 0.4)
  id = 150
  plot(y, type = 'l')
  
  ##### AO ####
  
  if (sign(y[id]) == -1){ y[id] = y[id] - rnorm(1, mu, .1)
  } else y[id] = y[id] + rnorm(1, mu, .1);
  
  
  #### TC ####
  
  a = rnorm(1, mu, .1)
  for (l in seq(id, length(y))) y[l] = y[l] + a * 0.7^(l+1-100)
  
  #### LS ####
  
  # b = rnorm(1, mu, .1)
  # for (k in id:length(y)) y[k] = y[k] + b
  
  y = y - mean(y)
  plot(y, type = 'l')
  
  ########################### TSOUTLIERS ###########################
  
  out1 = ARMAoutlier(y);
  out1
  
  out2 = HARoutlier(y);
  out2
  
  out.har = HAR(y);
  #    out2 = HARGARCH(y)
  resid = ts(as.numeric(ts(out.har$residuals)))
  fit2 =  arima(y, order=c(2,0,0))
  pars = out.har$pars
  types = c("AO", "LS", "TC"); delta = 0.7;
  sigma <- 1.483 * quantile(abs(resid - quantile(resid, probs = 0.5, 
                                                 na.rm = TRUE)), probs = 0.5, na.rm = TRUE)
  tmp <- outliers.tstatistics(pars = pars, resid = resid, types = types, 
                              sigma = sigma, delta = delta)
  ao = tmp[,1,2];
  ls = tmp[,2,2];
  tc = tmp[,3,2];
  
  
  plot.ts(tmp[,1,2])
  
  
#  library(fracdiff)
#  out4 =  fracdiff(y)
#   
#  tmp2 = outliers.tstatistics(pars, resid=out4$residuals, types=types, sigma=sigma, delta=delta)
#  plot.ts(tmp2[,1,2])
  
  # len = c(1:300)
  
  # TP1 = length(mo31$outliers[mo31$outliers$ind %in% id, ]$ind) # correctly detected as outliers
  # FP1 = length(mo31$outliers[!mo31$outliers$ind %in% id, ]$ind) # outliers?ƴѵ? detected as outliers
  # FN1 = length(id[!id %in% mo31$outliers$ind]) # outliers?ε? not detected as outliers
  # TN1 = length(len[-c(mo31$outliers$ind, id)])
  
  # specificity1 = TN1 / (TN1 + FP1)
  # sensitivity1 = TP1 / (TP1 + FN1) 
  
  # TP2 = length(mo32$outliers[mo32$outliers$ind %in% id, ]$ind) # correctly detected as outliers
  # FP2 = length(mo32$outliers[!mo32$outliers$ind %in% id, ]$ind) # outliers?ƴѵ? detected as outliers
  # FN2 = length(id[!id %in% mo32$outliers$ind]) # outliers?ε? not detected as outliers
  # TN2 = length(len[-c(mo32$outliers$ind, id)])
  
  # specificity2 = TN2 / (TN2 + FP2)
  # sensitivity2 = TP2 / (TP2 + FN2)
  
  #####################################################################
  
  # result[i, "specificity(ARIMA)"] = specificity1
  # result[i, "specificity(HAR)"] = specificity2
  # result[i, "sensitivity(ARIMA)"] = sensitivity1
  # result[i, "sensitivity(HAR)"] = sensitivity2
  # result[i, "TN(ARIMA)"] = TN1
  # result[i, "TP(ARIMA)"] = TP1
  # result[i, "FN(ARIMA)"] = FN1
  # result[i, "FP(ARIMA)"] = FP1 
  # result[i, "TN(HAR)"] = TN2
  # result[i, "TP(HAR)"] = TP2
  # result[i, "FN(HAR)"] = FN2
  # result[i, "FP(HAR)"] = FP2
  # result[i, "num(ARIMA)"] = TP1
  # result[i, "num(HAR)"] = TP2
  
  if(isFALSE(nrow(mo31$outliers) == 0)) {
    result1[i, 1:nrow(mo31$outliers)] = mo31$outliers$ind; colnames(result1) = rep('time', nrow(mo31$outliers))
    result2[i, 1:nrow(mo31$outliers)] = as.character(mo31$outliers$type); colnames(result2) = rep('type', nrow(mo31$outliers))
  }
  
  else {
    result1[i, ] = NA
    result2[i, ] = NA
  }
  
  if(isFALSE(nrow(ho32$outliers) == 0)) {
    result3[i, 1:nrow(ho32$outliers)] = ho32$outliers$ind; colnames(result3) = rep('time', nrow(ho32$outliers))
    result4[i, 1:nrow(ho32$outliers)] = as.character(ho32$outliers$type); colnames(result4) = rep('type', nrow(ho32$outliers))
  }
  
  else {
    result3[i, ] = NA
    result4[i, ] = NA
  }
  
  print(i)
  
}


warnings()

colnames(result1) = rep('time.AR', ncol(result1)); colnames(result2) = rep('type.AR', ncol(result2)); colnames(result3) = rep('time.HAR', ncol(result3)); colnames(result4) = rep('type.HAR', ncol(result4))


result11 = c(result1[, 1], result1[, 2])
result33 = c(result3[, 1], result3[, 2], result3[, 3], result3[, 4])


par(mfrow = c(1, 2))
boxplot(result11, main = "AR.index")
boxplot(result33, main = 'HAR.index')
par(mfrow = c(1, 1))

pd1 = 150 - log(length(y))
pd2 = 150 + log(length(y))

a = c()
for(i in 1:200) {
  a[i] = length(which(result1[i, ] >= pd1 & result1[i, ] <= pd2)) / sum(!is.na(result1[i, ]))
}

b = c()
for(i in 1:200) {
  b[i] = length(which(result3[i, ] >= pd1 & result3[i, ] <= pd2)) / sum(!is.na(result3[i, ]))
}

a[is.na(a)] = 0
mean(a)

b[is.na(b)] = 0
mean(b)

result = cbind(result1, result2, result3, result4)

nrow(result)

write.csv(result, "C:\\Users\\km960\\Desktop\\simulation_result\\SINGLE\\simul_210903_300_0.4_TC_FARIMA_SINGLE.csv", row.names = FALSE)
