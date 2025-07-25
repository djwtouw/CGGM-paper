rm(list=ls())


#### Function to estimate the HAR model of Corsi ####
estimate.HAR <- function(y)
{
  # Function : HAR model with daily, weekly and monthly components for RV/RD

  # INPUT
  # y: time series matrix of dimension Txp containing daily RV/RD in the columns

  # OUTPUT
  # resid.HAR : (T-22)xp containing the residuals in the columns of the
  #             univariate HAR regressions ran on each column of y

  # START CODE

  HARdata <- embed(y, 22+1)
  XHAR <- HARdata[, -1]

  YHAR <- HARdata[,1]
  X.D <- as.matrix(XHAR[,1])
  X.W <- as.matrix(apply(XHAR[,1:5]/5,1,sum))
  X.M <- as.matrix(apply(XHAR[,]/22,1,sum))
  X.HAR <- cbind(1, X.D,X.W,X.M)
  beta.HAR <- solve(t(X.HAR)%*%X.HAR)%*%t(X.HAR)%*%YHAR
  resid.HAR <- YHAR - X.HAR%*%beta.HAR

  return(resid.HAR)
}


data_select <- function(df, price_type)
{
  # Function : Extract the data for a given price type (High or Low) for each of
  #            the stocks. The column names contain the type, the first row of
  #            the input contains the tickers.

  # Select prices
  idx = which(sapply(colnames(df), function (cn) {
    return(grepl(price_type, cn, fixed = TRUE))
  }))
  data = df[, idx]

  # First row contains tickers, so set column names to the first row
  colnames(data) = data[1, ]

  # Drop first and second row (second is empty due to format)
  data = data[-c(1, 2), ]

  # Next issue, numbers are stored as char
  data = apply(data, c(1, 2), as.numeric)

  return(data)
}


#### Data ####
# Load raw price data
data = read.csv("applications/finance/data/sp100_raw.csv")

# Column price contains the dates
data$Price = NULL

# Get high and low prices
data_hi = data_select(data, "High")
data_lo = data_select(data, "Low")

# Compute daily range
data_dr = (log(data_hi) - log(data_lo))^2 / (4 * log(2))

# Use HAR model to obtain log realized daily range residuals
rd_clean = apply(data_dr, 2, estimate.HAR)

# Get the correlations with the first lag
test = apply(rd_clean, 2, function(X) {acf(X)$acf[2]})

# Plot some figures to inspect autocorrelation
round(test, 2)
par(mfrow = c(1, 2))
hist(test, main = "First-order autocorrelations")
plot(seq(0, 0.20, 0.005),
     sapply(seq(0, 0.20, 0.005), function (c) {
       sum(abs(test) < c) / length(test)
     }),
     type = "l",
     xlab = "c",
     ylab = "fraction within [-c, c]")

# Save data
save(rd_clean, file = "applications/finance/data/rd_clean.RData")
