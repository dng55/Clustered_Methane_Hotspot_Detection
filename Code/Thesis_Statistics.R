# library(reticulate)
# use_python("/Library/Frameworks/Python.framework/Versions/3.9/bin/python3.9")
# source_python("/Users/darianng/Documents/MSc_Geography/MSc Thesis/Code/pickle_reader.py")
# data <- read_pickle_file("/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Compilation/Young_Annual/2022_growing_season/compiled_data.p")
rm(list=ls())
library(reticulate)
pd <- import("pandas")
data <- pd$read_pickle("/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Compilation/Young_Annual/2022_growing_season/compiled_data.p")

# Consolidating data into a dataframe
NDVI <- unlist(data$collapsed_landsat$NDVI)
NDWI <- unlist(data$collapsed_landsat$NDWI)
NDMI <- unlist(data$collapsed_landsat$NDMI)
temp <- unlist(data$collapsed_landsat$temp)
month.f <- as.factor(data$collapsed_landsat$month)
CH4 <- unlist(data$collapsed_ffp)
df<-data.frame(NDVI,NDWI,NDMI,temp,month.f,CH4)
# Introducing Month as a factor


# Boxplot mean CH4
tapply(df$CH4, df$month.f, mean) # simple means
boxplot(CH4 ~ month.f, data=df, col='pink', 
        main="CH4 by month")

# Scatter plot matrix
pairs(~ NDVI + NDWI + NDMI + temp + month.f, data=df, main="Scatterplot Matrix for CH4 data",pch = 20)
pairs(~ CH4+NDVI + NDWI + NDMI + temp + month.f, data=df, main="Scatterplot Matrix for CH4 data",pch = 20)


# Linear Model
model <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp+ month.f*NDVI+month.f*NDWI+month.f*NDMI+month.f*temp, data=df)
model2 <- lm(CH4 ~+NDVI+NDWI+NDMI+temp, data=df)
summary(model)
anova(model)
Anova(model, type="III")
df$yhat <- fitted(model)
df$resid <- resid(model)
p_val <- 1-pf(summary(model)$fstatistic[1],summary(model)$fstatistic[2],summary(model)$fstatistic[3])
p_val2 <- 1-pf(summary(model2)$fstatistic[1],summary(model2)$fstatistic[2],summary(model2)$fstatistic[3])

# Compare b/n model and model2
par(mfrow=c(2,1), cex=1.0, mai=c(0.8, 1.3, 0.4, 0.8))
# Model
plot(yhat ~ CH4, data=df, xlab="Observed CH4", 
     ylab="Predicted CH4", main="Predicted vs. observed CH4 (with month.f)",
     pch=19)
abline(a=0,b=1,col="red") 
text(80,20,paste("R-squared:",round(summary(model)$r.squared,4)))
# Model 2
plot(df$CH4,fitted(model2), xlab="Observed CH4", 
     ylab="Predicted CH4", main="Predicted vs. observed CH4",
     pch=19)
abline(a=0,b=1,col="red") 
text(80,20,paste("R-squared:",round(summary(model2)$r.squared,4)))


#  Assumptions met for Analysis 3?
#####################################################################
# save the yhats and residuals.
df$yhat <- fitted(model)
df$resid <- resid(model)
cbind(df$CH4, df$month, df$NDVI,
      df$yhat, df$resid) # show yhat's (emmeans) & residuals.
# Diagnostic plots -- Assumptions met?
par(mfrow=c(2,2), cex=1.0, mai=c(1.0,1.0,0.6,0.6))
# Fitted line plot. Predicted y versus observed y.
plot(df$yhat ~ df$CH4, xlab="y", 
     ylab="yhat", main="Predicted vs Actual y",
     pch=19)
abline(a=0,b=1,col="red") 
plot(df$resid ~ df$yhat, main="Residual Plot")  
abline(a=0,b=0, col="red")  # residual plot
qqnorm(df$resid)     # normality plot
qqline(df$resid, col=2)
hist(df$resid, density=10, main="Residuals Distribution",
     col="green", border="black") 
par(mfrow=c(1,1), cex=1.0, mai=c(1.0,1.0,1.0,1.0) )

# Do a normality test
shapiro.test(df$resid) # Shapiro-Wilk normality test

require(car) # This must be loaded in advance!
Anova(model, type=('III')) # get the Type III SS tests


# Save Model Results
df_coeff <- coefficients(model)
savedir = "/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Model/"
setwd(savedir)
write.csv(df, file = "Young2022_fitted.csv", row.names = TRUE)
write.csv(df_coeff, file = "Young2022_coeff.csv", row.names = TRUE)


# Fitted line plot: predicted y versus the observed y.
plot(yhat ~ CH4, data=df, xlab="Observed CH4", 
     ylab="Predicted CH4", main="Predicted vs. observed CH4",
     pch=19)
abline(a=0,b=1,col="red") 
