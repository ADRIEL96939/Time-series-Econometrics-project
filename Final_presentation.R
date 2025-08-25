# Load library
library(fpp2)
library(portes)
library(urca)
library(dynlm)
library(ggplot2)
library(tidyr)
library(vars)
library(tsDyn)

##### Preliminaries start
emissions <- read.csv("~/Desktop/USA data.csv")
head(emissions) 

df <- ts(emissions[,-1], frequency = 1, start = 1974)
head(df, 6)
tail(df, 6)

co2 <- df[,1]
gdp <- df[,2]
fd <- df[,3]

lco2 <- log(co2)
lgdp <- log(gdp)
lfd <- log(fd)

ldf <- cbind(lco2, lgdp, lfd)

head(ldf, 6)
tail(ldf, 6)

#Take first difference of whole time series dataframe 
dldf <- diff(ldf)

#Print first 6 observations of ddf
head(dldf)

#Change the column names
colnames(dldf) <- c("dlco2", "dlgdp", "dlfd")

#Print first 6 observations of ddf again to check labels
head(dldf)

dlco2 <- dldf[, "dlco2"]
dlgdp <- dldf[, "dlgdp"]
dlfd <- dldf[, "dlfd"]


##### Preliminaries end

# Our research question is: What are the effects of economic growth (GDP per capita), and financial development (domestic credit) on CO2 emissions in the USA, and how do these variables interact over the long term and short term?
# To observe the relationship, firstly, we will look at the trend and behaviour of the three variables.

##### Data description start
# Combine the time series into a data frame to plot all variables in one facet.
data <- data.frame(
  time = time(lco2),
  lco2 = as.numeric(lco2),
  lgdp = as.numeric(lgdp),
  lfd = as.numeric(lfd)
)

# Reshape the data into long format
data_long <- pivot_longer(data, cols = -time, names_to = "variable", values_to = "value")

# Plot using ggplot with facet_wrap
ggplot(data_long, aes(x = time, y = value)) +
  geom_line() +
  facet_wrap(~ variable, scales = "free_y", ncol = 1) +
  labs(x = "Time", y = "Value", title = "Time series plots of lco2, lgdp, and lfd") +
  theme_minimal()

##### Data description end

#### After looking at the general trends and behaviour of the variables, we want to check whether they are I(1) processess
### as I(1) processes are modelled differently from standard I(0) processess.

##### Unit root test (ADF) for lc02

#Choose lag truncation parameter
adf_selectP_fit <- ur.df(lco2, lags = 10, selectlags = "AIC", type = "drift")

#Obtain result (only to see number of lags required)
summary(adf_selectP_fit)

# 3 lags are included in the regression

#Re-estimate regression with lags set at 3
adf_fit <- ur.df(lco2, lags = 3, selectlags = "Fixed", type = "drift")

#Print estimation output
summary(adf_fit)

#Extract residuals
adf_residuals <- adf_fit@res

#Perform autocorrelation tests on regression residuals, with 5 parameters estimated
LjungBox(adf_residuals, lags = seq(1,15,1), fitdf = 5) 
LjungBox(adf_residuals, lags = 15, fitdf = 5) 

# Since p-value = 0.7782781 > 0.05, we conclude sufficient lags have been included

#the estimated regression equation given in pdf

adf_fit@teststat
# tau2 t stats = -1.368865
# phi1 t stats = 0.9873211

adf_fit@cval
# tau2 crit = -2.93
# phi1 crit = 4.86
# tstats larger than tau2 crit, no reject null of unit root but low power can affect this

############ Sequential Procedure (only for verification)

#Create variable dlco2 by first difference
drate <- diff(lco2)
#Estimate unrestricted model by OLS
lco2_UR <- dynlm(dlco2 ~ L(lco2,1) + L(dlco2,1:3))
#Obtain RSSU
rssu = deviance(lco2_UR)
rssu

#Estimate restricted model by OLS
lco2_R <- dynlm(dlco2 ~ L(dlco2,1:3) -1) #the '-1' is to not estimate an intercept
#Obtain RSSR
rssr = deviance(lco2_R)
rssr

#Calculate F1 stat
T = 38
K = 5
q = 2
F1 = ((rssr - rssu)/rssu) * (T-K)/q
print(F1)

# 0.9873 < (fcrit at around T = 50 for 5% level of significance which is 4.86)
# do not reject H0 : a = γ = 0 and conclude that the intercept is not significant in the regression equation
# follows non-standard distribution
# Prev tstats is correct --> fail to reject null of unit root, lco2 has presumably one unit root as we differenced once dlco2 becomes stationary.

#####################################

#### Unit root test (ADF) for dlc02

dlco2 <- diff(lco2, 1)
t.test(dlco2, mu = 0)
#p-value = 0.5366 > 0.05, fail to reject null of zero mean, so mean is zero.
# zero mean no deterministic trend

#Choose lag truncation parameter
adf_selectP_fit <- ur.df(dlco2, lags = 10, selectlags = "AIC", type = "none")

#Obtain result
summary(adf_selectP_fit) # this is the regression if sufficient lags included, dont forget the tstats

###### Ljung box test to see if sufficient lags are included

#Re-estimate regression with lags set at 2
adf_fit <- ur.df(dlco2, lags = 2, selectlags = "Fixed", type = "none")

#Print estimation output
summary(adf_fit)

#Extract residuals
adf_residuals <- adf_fit@res

#Perform autocorrelation tests on regression residuals
LjungBox(adf_residuals, lags = seq(1,15,1), fitdf = 3)
LjungBox(adf_residuals, lags = 15, fitdf = 3)

###Perform actual test for unit root
adf_fit@teststat
# tau1 t-stats = -2.9545
adf_fit@cval
# tau 1 crit = -1.95
# tstats < crit for tau1 --> reject null of unit root --> no unit root
# interpret with caution

#####################################

#### Unit root test (ADF) for lgdp

#Choose lag truncation parameter
adf_selectP_fit <- ur.df(lgdp, lags = 10, selectlags = "AIC", type = "trend")

#Obtain result
summary(adf_selectP_fit) # this is the regression if sufficient lags included, dont forget the tstats

###### Ljung box test to see if sufficient lags are included

#Re-estimate regression with lags set at 1
adf_fit <- ur.df(lgdp, lags = 1, selectlags = "Fixed", type = "trend")

#Print estimation output
summary(adf_fit)

#Extract residuals
adf_residuals <- adf_fit@res

#Perform autocorrelation tests on regression residuals
LjungBox(adf_residuals, lags = seq(1,15,1), fitdf = 4)
LjungBox(adf_residuals, lags = 15, fitdf = 4)
#sufficient lags

###Perform actual test for unit root
adf_fit@teststat
# tau3 t-stats = -2.192574
adf_fit@cval
# tau3 crit = -3.50
# t-stats > crit --> no reject null of unit root, at least one unit root, but differencing once unit root --> one unit root
# consider low power --> interpret with caution

#############################

#### Unit root test for dlgdp

dlgdp <- diff(lgdp, 1)
t.test(dlgdp, mu = 0)
# p-value = 7.701e-07 < 0.05, reject null that mean is 0
# mean is non-zero, and no deterministic trend

#Choose lag truncation parameter
adf_selectP_fit <- ur.df(dlgdp, lags = 10, selectlags = "AIC", type = "drift")

#Obtain result
summary(adf_selectP_fit) # this is the regression if sufficient lags included, dont forget the tstats

###### Ljung box test to see if sufficient lags are included

#Re-estimate regression with lags set at 1
adf_fit <- ur.df(dlgdp, lags = 1, selectlags = "Fixed", type = "drift")

#Print estimation output
summary(adf_fit)

#Extract residuals
adf_residuals <- adf_fit@res

#Perform autocorrelation tests on regression residuals
LjungBox(adf_residuals, lags = seq(1,15,1), fitdf = 3)
# 45h lags above sufficient; p-value > 0.05 using same hypothesis stated in question 3(a)
LjungBox(adf_residuals, lags = 15, fitdf = 3)

###Perform actual test for unit root
adf_fit@teststat
# tau1 tstats = -3.840307
adf_fit@cval
# tau1 crit = -2.93
# tstats < crit for tau1 --> reject null of unit root --> no unit root
# interpret with caution

################################

#### Unit root test for lfd

#Choose lag truncation parameter
adf_selectP_fit <- ur.df(lfd, lags = 10, selectlags = "AIC", type = "trend")

#Obtain result
summary(adf_selectP_fit) # this is the regression if sufficient lags included, dont forget the tstats

###### Ljung box test to see if sufficient lags are included

#Re-estimate regression with lags set at 1
adf_fit <- ur.df(lfd, lags = 8, selectlags = "Fixed", type = "trend")

#Print estimation output
summary(adf_fit)

#Extract residuals
adf_residuals <- adf_fit@res

#Perform autocorrelation tests on regression residuals
LjungBox(adf_residuals, lags = seq(1,20,1), fitdf = 11)
LjungBox(adf_residuals, lags = 15, fitdf = 11)
#sufficient lags

###Perform actual test for unit root
adf_fit@teststat
# tau3 tstats = -0.7440711
adf_fit@cval
# tau3 crit = -3.50
# t-stats > crit --> no reject null of unit root, at least one unit root, but differencing once unit root --> one unit root
# consider low power --> interpret with caution

#### Unit root test for dlfd

dlfd <- diff(lfd, 1)
t.test(dlfd, mu = 0)
# p-value = 0.02046 < 0.05, reject null of zero mean
# non-zero mean, no deterministic trend

#Choose lag truncation parameter
adf_selectP_fit <- ur.df(dlfd, lags = 10, selectlags = "AIC", type = "drift")

#Obtain result
summary(adf_selectP_fit) # this is the regression if sufficient lags included, dont forget the tstats

###### Ljung box test to see if sufficient lags are included

#Re-estimate regression with lags set at 4
adf_fit <- ur.df(dlfd, lags = 4, selectlags = "Fixed", type = "drift")

#Print estimation output
summary(adf_fit)

#Extract residuals
adf_residuals <- adf_fit@res

#Perform autocorrelation tests on regression residuals
LjungBox(adf_residuals, lags = seq(1,15,1), fitdf = 6)

# lags 7 above sufficient; p-value > 0.05 using same hypothesis stated in question 3(a)
LjungBox(adf_residuals, lags = 15, fitdf = 6)

###Perform actual test for unit root
adf_fit@teststat
# tau1 tstats = -3.350144
adf_fit@cval
# tau1 crit = -2.93
# tstats < crit for tau1 --> reject null of unit root --> no unit root
# interpret with caution

######################################### We notice that all of the original log transforemed variables are I(1) processes, with their first differenced form to be I(0) processes.

##### Multivariate model selection and justification

### Since we know all three variables are I(1) processes, we will not proceed with the VAR framework using the original log transformed variables.
## Instead we will perform a cointegration test to determine whether there is cointegration or not.
# If there is no cointegration, we will proceed with modelling the relationship b/w the 3 variables using the VAR model for the log transformed variables in their first difference.
# If there is cointegration, we will proceed with modelling the relationship b/w the 3 variables using the VECM model.
# However, in order to estimate the VECM model, we need to know which order it is, so we will derive it from the VAR model.

ldf <- cbind(lgdp, lfd, lco2)
head(ldf, 6)
tail(ldf, 6)

autoplot(ldf, facets = TRUE) +
  ggtitle("Time Series Plots for lgdp, lfd, and lco2") +
  xlab("Time") 

#### Since the variables are I(1), inclusion of the intercept in the VAR model is only appropriate when at least one of the variables has a deterministic trend. Given that lgp, and lfd
### appear to have deterministic trends, the intercept must be included in the VAR model to account for these systematic changes over time.

### We will arbitrarily choose our maximum lag to select to be 4
# Select initial VAR order --> VAR (2)
VARselect(ldf, lag.max = 4, type = "const") #const means include intercept

#Estimate VAR(2)
var2 <- VAR(ldf, p = 2, type = "const")

#Portmanteau test for serial correlation in residuals
serialtest <- serial.test(var2, lags.pt = 8, type = "PT.asymptotic")
print(serialtest)

# p-value > α, we fail to reject the null and conclude that the errors are vector white noise.

# In order to ensure that VAR of order 2 is the most parsimonious model, we test it against VAR of order 1.
# Test the null hypothesis that the VAR order is 1 against the alternative that it is 2.

LU <- logLik(var2)
print(LU)

var1 <- VAR(ldf, p = 1, type = "const")
LR <- logLik(var1)
print(LR)

LRcalc <- 2*(LU - LR)
print(LRcalc)

#Obtain chi-squared value at 5% level of significance with 9 df
qchisq(0.95, 9)
# LRcalc = 30.11783 > 16.919, we reject H0 that the VAR order is 1 in favour of H1 that the VAR order is 2.

### Now, we know that we must estimate a VECM model of order 1 if there is cointegration.

### Now, we will perform our cointegration test using the VECM model

# As lco2, lgdp and lfd all appear to have non-zero means, we should include an intercept in the cointegrating equation. 
# Since two out of the three variables (lgdp, and lfd) also appear to have a linear deterministic trend, we should have a linear trend in the VECM.
# Thus, we will either fit a case 3 VECM model (with non-zero mean, linear trend, and non-zero mean cointegrating relation)  or case 4 VECM model (with non-zero mean, linear trend, and trend cointegrating relation)  

# We will estimate the case 3 VECM model first
fit_vecm <- VECM(ldf, lag = 1, include = "const", estim = "ML", LRinclude = "none")
summary(fit_vecm)
# AIC =  -922.5832, BIC -893.8723

# Case 4 VECM model
fit_vecm_1 <- VECM(ldf, lag = 1, include = "const", estim = "ML", LRinclude = "trend")
summary(fit_vecm_1)
# AIC = -926.8192 BIC = -898.1083 
## It is clearly seen the case 3 model is a better fit with lower AIC and BIC values.

# Now, we perform the cointegration test (Trace and Maximum Eigenvalue Tests)
coint_eigen <- rank.test(fit_vecm, type = "eigen", cval = 0.05)
coint_trace <- rank.test(fit_vecm, type = "trace", cval = 0.05)

#Print the rank selected by each test 
coint_eigen 
coint_trace
# Both conclude that the cointegrating rank is 1.

#Print the test-statistics and associated p-values 
summary(coint_eigen) 
summary(coint_trace)

#Estimate the coefficients of VECM model with r = 1
fit_vecm <- VECM(ldf, lag = 1, r = 1, include = "const", estim = "ML", LRinclude = "none")

#Print estimation output 
summary(fit_vecm)

#Extract cointegrating vector 
fit_vecm$model.specific$coint

###########################
# Granger causality
###########################

# Extract dataframe that also contains ECT

vecm.df <- fit_vecm$model

# Combine vecm.df with dldf 

vecm.df2 <- cbind(vecm.df, dldf)
colnames(vecm.df2)

# update header names 

colnames(vecm.df2) <- c("lgdp", "lfd", "lco2", "ECT", 
                          "Intercept", "dlgdp.L1", "dlfd.L1", 
                          "dlco2.L1", "dlco2", "dlgdp", 
                          "dlfd")

# Estimate UR eqn by OLS
library(dynlm)

UR.fit_dlco2 <- dynlm(dlco2 ~ ECT + L(dlgdp, 1) + L(dlfd, 1) + L(dlco2, 1), 
                data = vecm.df2)
summary(UR.fit_dlco2)
# ECT is significant, convergence to the long-run equilibrium will be at a relatively slow speed 

# Perform wald coeff test 

coef(UR.fit_dlco2)
coef(UR.fit_dlco2)[3]

library(aod)

"Does lgdp GC lco2?"
wald.test(b = coef(UR.fit_dlco2), Sigma = vcov(UR.fit_dlco2), Terms = 3)
# H0: gamma 1,1 = 0
# H1: gamma 1,1 not = 0
# do not reject null - lgdp does not granger cause lco2 (p-value = 0.75)

"Does lfd GC lco2?"
wald.test(b = coef(UR.fit_dlco2), Sigma = vcov(UR.fit_dlco2), Terms = 4)
#reject null - lfd granger causes lco2 (p-value = 0.0098)
# H0: gamma 1,2 = 0
# H1: gamma 1,2 not = 0

UR.fit_dlgdp <- dynlm(dlgdp ~ ECT + L(dlgdp, 1) + L(dlfd, 1) + L(dlco2, 1), 
                data = vecm.df2)
summary(UR.fit_dlgdp)
# ECT is significant, convergence to the long-run equilibrium will be at a relatively slow speed 
coef(UR.fit_dlgdp)
coef(UR.fit_dlgdp)[5]

UR.fit_dlfd <- dynlm(dlfd ~ ECT + L(dlgdp, 1) + L(dlfd, 1) + L(dlco2, 1), 
                data = vecm.df2)
summary(UR.fit_dlfd)
# ECT is significant, convergence to the long-run equilibrium will be at a relatively slow speed
coef(UR.fit_dlfd)
coef(UR.fit_dlfd)[5]

"Does lco2 GC lgdp?"
wald.test(b = coef(UR.fit_dlgdp), Sigma = vcov(UR.fit_dlgdp), Terms = 5)
# H0: gamma 2,1 = 0
# H1: gamma 2,1 not = 0
# fail to reject null at 5% level of significance - lco2 does not granger cause lgdp (p-value = 0.059)

"Does lco2 GC lfd?"

wald.test(b = coef(UR.fit_dlfd), Sigma = vcov(UR.fit_dlfd), Terms = 5)
# H0: gamma 3,1 = 0
# H1: gamma 3,1 not = 0
#reject null - lc02 granger causes lfd (p-value = 0.029)
# there is bi-directional granger causality

"Does lgdp GC lfd?"
wald.test(b = coef(UR.fit_dlfd), Sigma = vcov(UR.fit_dlfd), Terms = 3)
# H0: gamma 3,2 = 0
# H1: gamma 3,2 not = 0
# reject null at 1% level of significance - lgdp does granger cause lfd (p-value = 0.047)

"Does lfd GC lgdp?"
wald.test(b = coef(UR.fit_dlgdp), Sigma = vcov(UR.fit_dlgdp), Terms = 4)
# H0: gamma 2,3 = 0
# H1: gamma 2,3 not = 0
# reject null at 5% level of significance - lfd does granger cause lgdp (p-value = 0.011)

########################################
# Orthogonal Impulse Response Analysis
########################################

set.seed(123)

#Calculate impulse response functions of a shock to GDP on carbon emissions
irf_VECM1 <- irf(fit_vecm, impulse = "lgdp", response = "lco2")

# Extract the response and CIs
response_lco2 <- irf_VECM1$irf[["lgdp"]][, "lco2"]
ci_lower <- irf_VECM1$Lower[["lgdp"]][, "lco2"]
ci_upper <- irf_VECM1$Upper[["lgdp"]][, "lco2"]

# Plot the Orthogonal Impulse Response
plot(response_lco2, type = "l", col = "blue", ylim = range(c(ci_lower, ci_upper)),
     main = "Orthogonal Impulse Response from dlgdp", ylab = "dlco2", xlab = "95% Bootstrap CI, 100 Runs")

# Add confidence intervals
lines(ci_lower, col = "red", lty = 2)
lines(ci_upper, col = "red", lty = 2)

#Look at several lower CIs to tentatively check whether they differ from zero
print(ci_lower)

set.seed(123)

#Calculate impulse response functions of a shock to financial development on carbon emissions
irf_VECM2 <- irf(fit_vecm, impulse = "lfd", response = "lco2")

# Extract the response
response_lco2_1 <- irf_VECM2$irf[["lfd"]][, "lco2"]
ci_lower1 <- irf_VECM2$Lower[["lfd"]][, "lco2"]
ci_upper1 <- irf_VECM2$Upper[["lfd"]][, "lco2"]

# Plot the Orthogonal Impulse Response
plot(response_lco2_1, type = "l", col = "blue", ylim = range(c(ci_lower1, ci_upper1 + 0.010)),
     main = "Orthogonal Impulse Response from dlfd", ylab = "dlco2", xlab = "95% Bootstrap CI, 100 Runs")

# Add confidence intervals
lines(ci_lower1, col = "red", lty = 2)
lines(ci_upper1, col = "red", lty = 2)

#Look at several lower CIs to tentatively check whether they differ from zero
print(ci_lower1)
