install.packages("plotrix")
install.packages("kableExtra")
install.packages("dplyr")
install.packages("reshape2")
install.packages("car")
install.packages("lubridate")
install.packages("moments")
install.packages("xts")
install.packages("zoo")
install.packages("scales")
install.packages("dgof")
install.packages("tseries")
install.packages("rugarch")
install.packages("fpp2")
install.packages("RColorBrewer")
install.packages("QRM")
install.packages("quantmod")
install.packages("PerformanceAnalytics")

library(plotrix)
library(dplyr)
library(knitr)
library(reshape2)
library(car)
library(lubridate)
library(moments)
library(xts)
library(zoo)
library(scales)
library(dgof)
library(tseries)
library(rugarch)
library(tseries)
library(fpp2)
library(RColorBrewer)
library(QRM)
library(quantmod)
library(PerformanceAnalytics)


################## Amazon #########################3

rm(list=ls())

#Exporting the data and calculating returns

data <- read.csv("amzn.csv")
y <- diff(log(data$Adj.Close))
dates <- ymd(data$Date)[2:dim(data)[1]]
plot(dates, y, type = "l", xlab = "Year", ylab = "Returns")

#Performing Descriptive statistics

sd(y)
min(y)
max(y)
skewness(y)
kurtosis(y)
acf(y,1)
acf(y^2,1)
jarque.bera.test(y)
Box.test(y, lag = 20, type = c("Ljung-Box"))
Box.test(y^2, lag = 20, type = c("Ljung-Box"))
crisis_2008 <- y[year(dates) >= 2008 & year(dates) < 2010]
crisis_2020 <- y[year(dates) >= 2020 & year(dates) < 2022]
crisis_2022 <- y[year(dates) >= 2022]

qqPlot(y)
qqPlot(y,distribution="t", df=3, envelope=F)

df <- data.frame(dates, y)

#Calculating volatility for each financial crisis

vol <- sd(y)
vol_2008 <- sd(crisis_2008)
vol_2020 <- sd(crisis_2020)
vol_2022 <- sd(crisis_2022)

p <- 0.05

# Assume we have a portfolio value of 1000 USD
portfolio <- 1000

##### Backtesting VaR ##########

# Function that creates a GARCH forecast
DoGARCH <- function(y, spec, probability = 0.05, portfolio_value = 1, WE = 1000){
  old <- Sys.time()
  
  cat("Doing GARCH VaR forecast", "\n",
      "Estimation window:", WE, "\n",
      "Number of observations:", length(y), "\n",
      "VaR probability:", probability, "\n",
      "Portfolio value:", portfolio_value)
  
  # Number of observations
  n <- length(y)
  
  VaR <- rep(NA, n)
  
  for (i in 1:(n-WE)){
    
    window <- y[i:(i+WE-1)]
    
    res <- ugarchfit(spec = spec, data = window, solver = "hybrid")
    
    omega <- coef(res)['omega']
    alpha <- coef(res)['alpha1']
    beta <- coef(res)['beta1']
    
    sigma2 <- omega + alpha*tail(window,1)^2 + beta*tail(res@fit$var,1)
    
    VaR[i+WE] <- -sqrt(sigma2) * qnorm(probability) * portfolio_value
  }
  
  time <- difftime(Sys.time(), old, units = "secs")
  cat("\n", "Elapsed time:", round(time,4), "seconds")
  
  return(VaR)
}

#Defining the GARCH peculation (Not including the mean)
spec1 <- ugarchspec(
  variance.model = list(garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE)
)

#having window estimation for 500 and 1000 days
GARCH500 <- DoGARCH(y, spec = spec1, probability = 0.05, portfolio_value = 1000, WE = 500)
GARCH1000 <- DoGARCH(y, spec = spec1, probability = 0.05, portfolio_value = 1000, WE = 1000)

#speculation for GARCH
spec3 <- ugarchspec(
  variance.model = list(garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "std"
)

tGARCH <- ugarchfit(spec = spec3, data = y)

DotGARCH <- function(y, spec, probability = 0.05, portfolio_value = 1, WE = 1000){
  old <- Sys.time()
  
  cat("Doing tGARCH VaR forecast", "\n",
      "Estimation window:", WE, "\n",
      "Number of observations:", length(y), "\n",
      "VaR probability:", probability, "\n",
      "Portfolio value:", portfolio_value)
  
  n <- length(y)
  
  VaR <- rep(NA, n)
  
  for (i in 1:(n-WE)){
    
    # Subset the dataset to the estimation window
    window <- y[i:(i+WE-1)]
    
    res <- ugarchfit(spec = spec3, data = window, solver = "hybrid")
    
    omega <- coef(res)['omega']
    alpha <- coef(res)['alpha1']
    beta <- coef(res)['beta1']
    
    sigma2 <- omega + alpha*tail(window,1)^2 + beta*tail(res@fit$var,1)
    
    VaR[i+WE] <- -sqrt(sigma2) * qnorm(probability) * portfolio_value
  }
  
  time <- difftime(Sys.time(), old, units = "secs")
  cat("\n", "Elapsed time:", round(time,4), "seconds")
  
  return(VaR)
}

#tGARCH for 1000 days
tGARCH1000 <- DotGARCH(y, spec = spec3, probability = 0.05, portfolio_value = 1000, WE = 1000)

#Historical Simulation Backtesting
DoHS <- function(y, probability = 0.05, portfolio_value = 1000, WE = 500){

  old <- Sys.time()
  
  cat("Doing Historical simulation VaR forecast", "\n",
      "Estimation window:", WE, "\n",
      "Number of observations:", length(y), "\n",
      "VaR probability:", probability, "\n",
      "Portfolio value:", portfolio_value)
  
  n <- length(y)
  
  VaR <- rep(NA, n)
  
  for (i in 1:(n-WE)){
    
    ys <- sort(y[i:(i+WE-1)])
    
    quant <- ceiling(probability * length(ys))
    
    VaR[i+WE] <- -ys[quant]*portfolio_value
  }
  
  time <- Sys.time() - old
  cat("\n", "Elapsed time:", round(time,4), "seconds")
  
  return(VaR)
}

HS500 <- DoHS(y, probability = p, portfolio_value = 1000, WE = 500)
HS1000 <- DoHS(y, probability = p, portfolio_value = 1000, WE = 1000)

# EWMA Backtesting
#estimation window 300

p <- 0.05
portfolio_value <- 1000

portfolio_value = 1000
lambda <- 0.94
n <- length(y)
BurnTime <- 1000

EWMA_Variance <- rep(NA, length = n)

EWMA_Variance[1] <- var(y)

head(EWMA_Variance)

for (i in 2:n) {
  EWMA_Variance[i] <- lambda * EWMA_Variance[i-1] + (1-lambda) * y[i-1]^2
}

EWMA_Variance[1:BurnTime] <- NA

EWMA_cond_volatility <- sqrt(EWMA_Variance)


EWMA <- -qnorm(p) * EWMA_cond_volatility * portfolio_value

plot(dates, EWMA, type = "l", main = "EWMA VaR",
     las = 1, col = "red", xlab = "Date", ylab = "USD")


#VaR for each model
portfolio_value = 1000

VaR <- cbind(HS500, HS1000, EWMA, GARCH500, GARCH1000, tGARCH1000)

windows <- colSums(is.na(VaR))

start <- max(windows) + 1
end <- length(dates)

VaR <- as.data.frame(VaR)

#Violations for each risk models
Violations <- VaR
Violations[] <- NA

dim(Violations)

for(i in 1:dim(VaR)[2]){
  Violations[,i] <- (y*portfolio_value < -VaR[,i])
}
Violations[1:(start-1),] <- NA

w <- apply(Violations, 1, all)

sum(w, na.rm = TRUE)

# Counting Violations by model
colSums(Violations, na.rm = TRUE)

# Creating a Violation Ratio object

Violations <- Violations[!is.na(Violations[,1]),]

# Get the column sums
V <- colSums(Violations)

EV <- dim(Violations)[1]*p

# Violation Ratios
VR <- V/EV
round(VR,3)

#model assessment
model_assessment <- function(VR) {
  if (VR > 0.8 & VR < 1.2) {
    paste0(names(VR), "Model is good")
  } else if ((VR > 0.5 & VR <= 0.8) | (VR > 1.2 & VR <= 1.5)) {
    paste0(names(VR), "Model is acceptable")
  } else if ((VR > 0.3 & VR <= 0.5) | (VR > 1.5 & VR <= 2)) {
    paste0(names(VR), "Model is bad")
  } else {
    paste0(names(VR), "Model is useless")
  }
}

sapply(VR, model_assessment)

sort(round(abs(VR-1),3))


######### STRESS TESTING for each financial crisis ############

#for 2008

VaR_crisis_2008 <- VaR[year(dates) >= 2008 & year(dates) < 2010, ]

Violations_crisis_2008 <- VaR_crisis_2008
Violations_crisis_2008[] <- NA

for(i in 1:dim(VaR_crisis_2008)[2]){
  Violations_crisis_2008[,i] <- crisis_2008*portfolio_value < -VaR_crisis_2008[,i]
}

# Remove the rows with NA
Violations_crisis_2008 <- Violations_crisis_2008[!is.na(Violations_crisis_2008[,1]),]

# Get the column sums
V_crisis_2008 <- colSums(Violations_crisis_2008)

# Calculate expected violations
EV_crisis_2008 <- dim(Violations_crisis_2008)[1]*p

# Violation Ratios
VR_crisis_2008 <- V_crisis_2008/EV_crisis_2008

# Call object, rounding to 3 decimals
round(VR_crisis_2008,3)

model_assessment <- function(VR) {
  if (VR > 0.8 & VR < 1.2) {
    paste0(names(VR), "Model is good")
  } else if ((VR > 0.5 & VR <= 0.8) | (VR > 1.2 & VR <= 1.5)) {
    paste0(names(VR), "Model is acceptable")
  } else if ((VR > 0.3 & VR <= 0.5) | (VR > 1.5 & VR <= 2)) {
    paste0(names(VR), "Model is bad")
  } else {
    paste0(names(VR), "Model is useless")
  }
}

sapply(VR_crisis_2008, model_assessment)

sort(round(abs(VR_crisis_2008-1),3))

# 2020

VaR_crisis_2020 <- VaR[year(dates) >= 2020 & year(dates) < 2022, ]

Violations_crisis_2020 <- VaR_crisis_2020
Violations_crisis_2020[] <- NA

for(i in 1:dim(VaR_crisis_2020)[2]){
  Violations_crisis_2020[,i] <- crisis_2020*portfolio_value < -VaR_crisis_2020[,i]
}

Violations_crisis_2020 <- Violations_crisis_2020[!is.na(Violations_crisis_2020[,1]),]

V_crisis_2020 <- colSums(Violations_crisis_2020)

EV_crisis_2020 <- dim(Violations_crisis_2020)[1]*p

# Violation Ratios
VR_crisis_2020 <- V_crisis_2020/EV_crisis_2020

round(VR_crisis_2020,3)

model_assessment <- function(VR) {
  if (VR > 0.8 & VR < 1.2) {
    paste0(names(VR), "Model is good")
  } else if ((VR > 0.5 & VR <= 0.8) | (VR > 1.2 & VR <= 1.5)) {
    paste0(names(VR), "Model is acceptable")
  } else if ((VR > 0.3 & VR <= 0.5) | (VR > 1.5 & VR <= 2)) {
    paste0(names(VR), "Model is bad")
  } else {
    paste0(names(VR), "Model is useless")
  }
}

sapply(VR_crisis_2020, model_assessment)

sort(round(abs(VR_crisis_2020-1),3))


# 2022

VaR_crisis_2022 <- VaR[year(dates) >= 2022, ]

Violations_crisis_2022 <- VaR_crisis_2022
Violations_crisis_2022[] <- NA

for(i in 1:dim(VaR_crisis_2022)[2]){
  Violations_crisis_2022[,i] <- crisis_2022*portfolio_value < -VaR_crisis_2022[,i]
}


Violations_crisis_2022 <- Violations_crisis_2022[!is.na(Violations_crisis_2022[,1]),]

V_crisis_2022 <- colSums(Violations_crisis_2022)

EV_crisis_2022 <- dim(Violations_crisis_2022)[1]*p

# Violation Ratios
VR_crisis_2022 <- V_crisis_2022/EV_crisis_2022

round(VR_crisis_2022,3)

model_assessment <- function(VR) {
  if (VR > 0.8 & VR < 1.2) {
    paste0(names(VR), "Model is good")
  } else if ((VR > 0.5 & VR <= 0.8) | (VR > 1.2 & VR <= 1.5)) {
    paste0(names(VR), "Model is acceptable")
  } else if ((VR > 0.3 & VR <= 0.5) | (VR > 1.5 & VR <= 2)) {
    paste0(names(VR), "Model is bad")
  } else {
    paste0(names(VR), "Model is useless")
  }
}

sapply(VR_crisis_2022, model_assessment)

sort(round(abs(VR_crisis_2022-1),3))

#saving and loading the models 

save(GARCH1000, file = "GARCH1000AM.RData")
save(GARCH500, file = "GARCH500AM.RData")
save(tGARCH1000, file = "tGARCH1000AM.RData")
save(HS500, file = "HS500AM.RData")
save(HS1000, file = "HS1000AM.RData")
save(EWMA, file = "EWMAAM.RData")

load("GARCH1000AM.RData")
load("GARCH500AM.RData")
load("tGARCH1000AM.RData")
load("HS500AM.RData")
load("HS1000AM.RData")
load("EWMAAM.RData")

##### VaR  and ES###########

p <- 0.05

# Assume we have a portfolio value of 1000 USD
portfolio <- 1000

# Sort the values in y using sort()
crisis_2008_ys <- sort(crisis_2008)
crisis_2020_ys <- sort(crisis_2020)
crisis_2022_ys <- sort(crisis_2022)

quant_2008 <- ceiling(length(crisis_2008)*p)
quant_2020 <- ceiling(length(crisis_2020)*p)
quant_2022 <- ceiling(length(crisis_2022)*p)

# Calculating VaR for each crisis

VaR_2008 <- -crisis_2008_ys[quant_2008] * portfolio
VaR_2020 <- -crisis_2020_ys[quant_2020] * portfolio
VaR_2022 <- -crisis_2022_ys[quant_2022] * portfolio

#Calculating ES for each crisis

ES_2008 <- -mean(crisis_2008_ys[1:quant_2008]) * portfolio
ES_2020 <- -mean(crisis_2020_ys[1:quant_2020]) * portfolio
ES_2022 <- -mean(crisis_2022_ys[1:quant_2022]) * portfolio

#Plotting all the graphs for each financial crisis

#Focus on 2008-09
par(mar = c(5, 4, 4, 4) + 0.3)
matplot(dates[year(dates) >= 2008 & year(dates) < 2010], VaR_crisis_2008, type = "l", lty = 1, col = 1:6, xaxt = "n", xlab = "Date", ylab = "VaR USD")
par(new = TRUE)
plot(dates[year(dates) >= 2008 & year(dates) < 2010], y[year(dates) >= 2008 & year(dates) < 2010], type = "l", lty = 1, col = "grey", axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(y)))
axis.Date(1, at = seq(min(dates[year(dates) >= 2008 & year(dates) < 2010]), max(dates[year(dates) >= 2008 & year(dates) < 2010]), by = "years"))
mtext("Volatility", side = 4, line = 3) 
legend("topleft", legend = colnames(VaR), lty = 1, col = 1:6, cex = 0.65)

#Focus on 2020
par(mar = c(5, 4, 4, 4) + 0.3)
matplot(dates[year(dates) >= 2020 & year(dates) < 2022], VaR_crisis_2020, type = "l", lty = 1, col = 1:6, xaxt = "n", xlab = "Date", ylab = "VaR USD")
par(new = TRUE)
plot(dates[year(dates) >= 2020 & year(dates) < 2022], y[year(dates) >= 2020 & year(dates) < 2022], type = "l", lty = 1, col = "grey", axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(y)))
axis.Date(1, at = seq(min(dates[year(dates) >= 2020 & year(dates) < 2022]), max(dates[year(dates) >= 2020 & year(dates) < 2022]), by = "years"))
mtext("Volatility", side = 4, line = 3) 
legend("topright", legend = colnames(VaR), lty = 1, col = 1:6, cex = 0.65)

#Focus on 2020
par(mar = c(5, 4, 4, 4) + 0.3)
matplot(dates[year(dates) >= 2022], VaR_crisis_2022, type = "l", lty = 1, col = 1:6, xaxt = "n", xlab = "Date", ylab = "VaR USD")
par(new = TRUE)
plot(dates[year(dates) >= 2022], y[year(dates) >= 2022], type = "l", lty = 1, col = "grey", axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(y)))
axis.Date(1, at = seq(min(dates[year(dates) >= 2022]), max(dates[year(dates) >= 2022]), by = "years"))
mtext("Volatility", side = 4, line = 3) 
legend("topright", legend = colnames(VaR), lty = 1, col = 1:6, cex = 0.65)