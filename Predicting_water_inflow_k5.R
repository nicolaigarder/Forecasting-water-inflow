## -----------------------------------------------------------------------------------------------------------------
# Load the package
library(onlineforecast)
library(dplyr)


## -----------------------------------------------------------------------------------------------------------------
# Load necessary libraries
library(lubridate)

# Load the dataset
D <- read.csv("/Users/nicolaigarderhansen/Desktop/Bachelorprojekt/Forecasting-water-inflow/data_cleaned_and_interpolated.csv")

# Convert 'time' column to POSIXct format
D$t <- ymd(D$time)

# Assign the target variable (KarupBro) as y
D$y <- D$KarupBro

# Assign class for better structure
class(D) <- "data.list"

# Check structure
str(D)

make_tday <- function(t, kseq) {
  doy <- yday(t)  # Extract day of year (1-365 or 1-366)
  max_doy <- ifelse(leap_year(year(t)), 366, 365)  # Ensure correct max DOY per year
  
  # Compute future DOYs while correctly shifting
  out <- outer(doy, kseq, function(d, k) (d + k - 1) %% max_doy + 1)  
  
  as.data.frame(out)  # Convert to a dataframe
}



# Apply to Your Data
D$tday <- make_tday(D$t, 1:5)
colnames(D$tday) <- paste0("k", 1:ncol(D$tday))

horizons <- 1:5  # Define all required horizons

# SKAL FIKSES INDEN
# library(zoo)

# Fill missing values by carrying forward/backward
# D$y <- na.locf(D$y, na.rm = FALSE)  # Forward fill
# D$y <- na.locf(D$y, fromLast = TRUE, na.rm = FALSE)  # Backward fill

# Double-check if any NA values remain
# sum(is.na(D$y))  # Should return 0



## -----------------------------------------------------------------------------------------------------------------
# Set the score period
D$scoreperiod <- D$t >= as.POSIXct("2021-01-01", tz="UTC") & D$t <= as.POSIXct("2023-12-31", tz="UTC")


## ----output.lines=10----------------------------------------------------------------------------------------------
# Define a new model with low-pass filtering of the Ta input
model <- forecastmodel$new()

# Define output variable
model$output <- "y"

# Add raw inputs without low-pass filtering and a bias term
model$add_inputs(temp_mean_daily = "temp_mean_daily",
                 precip_total_daily = "precip_total_daily",
                 mu = "one()")

# Remove parameter bounds for `a1` since `lp(...)` is no longer used
model$add_prmbounds(lambda = c(0.85, 0.99, 0.9999))

# Add recursive least squares (RLS) regularization parameter
model$add_regprm("rls_prm(lambda=0.9)")

# Optimize parameters only on two horizons
kseqopt <- c(1,3,max(horizons))

# Convert input variables to forecast matrices
D$temp_mean_daily <- make_input(D$temp_mean_daily, horizons)
D$precip_total_daily <- make_input(D$precip_total_daily, horizons)

# Check if the transformation was successful
str(D$temp_mean_daily)  # Should be a dataframe with columns like k3, k18
str(D$precip_total_daily)  # Should also have k3, k18
rls_optim(model, D, kseqopt)

## -----------------------------------------------------------------------------------------------------------------
# Forecast for all horizons
model$kseq <- 1:max(horizons)
# Fit with RLS
fit1 <- rls_fit(model$prm, model, D)
# See the summary of the fit
summary(fit1)


## ----output.lines=10----------------------------------------------------------------------------------------------
# Add a diurnal curve using fourier series
model$add_inputs(mu_tday = "fs(tday/365, nharmonics=4)")
# Optimize the parameters
rls_optim(model, D, kseq=kseqopt)


## -----------------------------------------------------------------------------------------------------------------
# Fit with RLS
fit2 <- rls_fit(model$prm, model, D)
# Check the fit
summary(fit2)


## -----------------------------------------------------------------------------------------------------------------
# Keep the forecasts from each model by just inserting them in the data.list
D$Yhat1 <- fit1$Yhat
D$Yhat2 <- fit2$Yhat


## ----fig.height=figheight2----------------------------------------------------------------------------------------
# Plot to see the forecasts for the shortest and the longest horizon
plot_ts(subset(D,D$scoreperiod), c("^y|^Yhat1","^y|^Yhat2"), kseq = c(1,max(horizons)))


## ----fig.height=figheight2----------------------------------------------------------------------------------------
# Plot to see the forecasts for the shortest and the longest horizon
plot_ts(subset(D,which(D$scoreperiod)[1:30]), c("^y|^Yhat1","^y|^Yhat2"), kseq = c(1,max(horizons)))


## -----------------------------------------------------------------------------------------------------------------
# The simple persistence (forecast for same horizons as the model)
D$YhatP <- persistence(D$y, model$kseq)

# Ensure t is a properly formatted POSIXct vector inside the data.list
D$t <- as.POSIXct(D$t, origin = "1970-01-01", tz = "UTC")

plot_ts(subset(D, which(D$scoreperiod)[1:30]), c("^y$|YhatP$"), kseq = c(1, max(horizons)))


## -----------------------------------------------------------------------------------------------------------------
D$YhatP[1:4, 1:max(horizons)]
str(D$YhatP)


## -----------------------------------------------------------------------------------------------------------------
# Use the argument perlen to set the period length
D$YhatDP <- persistence(D$y, model$kseq, perlen=365)
# Plot a few horizons
plot_ts(D, c("^y$|YhatDP$"), kseq=c(1,max(horizons)))

# acf(D$y, lag.max=400)  # Look for a strong peak at lag 365


## -----------------------------------------------------------------------------------------------------------------
# Find the forecasts in D
nms <- grep("^Yhat", names(D), value=TRUE)
nms


## -----------------------------------------------------------------------------------------------------------------
# Find all complete cases for all forecasts and horizons
ok <- complete_cases(D[nms])


## -----------------------------------------------------------------------------------------------------------------
sum(ok)
length(ok)


## -----------------------------------------------------------------------------------------------------------------
D$Yhat1[1:max(horizons), 1:4]


## -----------------------------------------------------------------------------------------------------------------
D$y[59:72]
D$YhatP[59:72, 1]


## -----------------------------------------------------------------------------------------------------------------
ok <- ok & D$scoreperiod


## -----------------------------------------------------------------------------------------------------------------
sum(ok)
sum(D$scoreperiod)


## -----------------------------------------------------------------------------------------------------------------
# The score as a function of the horizon
R <- residuals(D$Yhat1, D$y)
score(R, ok & D$scoreperiod)


## -----------------------------------------------------------------------------------------------------------------
# Only complete cases are used per default
score(R, D$scoreperiod) == score(R, ok & D$scoreperiod)


## -----------------------------------------------------------------------------------------------------------------
# The score as a function of the horizon
score(R, usecomplete=FALSE) == score(R)


## -----------------------------------------------------------------------------------------------------------------
RMSE <- score(residuals(D[nms], D$y), D$scoreperiod)


## ----include=FALSE------------------------------------------------------------------------------------------------
# sapply(kseq, function(k){
#     rmse(y - lagdf(YhatDM[ ,pst("k",k)], k))
#     # hej det er vilfred jeg er peders sÃ¸n og jeg elsker min far go jeg god til matematik og jeg elsker ogsÃ¥ min mor 
# })


## ----fig.height=figheight2----------------------------------------------------------------------------------------
plot(0, type="n", xlim=range(model$kseq), ylim=range(RMSE), xlab="Horizon k", ylab="RMSE")
for(i in 1:ncol(RMSE)){
  points(model$kseq, RMSE[ ,i], type="b", col=i)
}
legend("topleft", nms, lty=1, col=1:length(nms))

