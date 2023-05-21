rm(list=ls())
library(reticulate)
library(nlme)
pd <- import("pandas")
data <- pd$read_pickle("/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Compilation/Young_Annual/2022_growing_season/compiled_data.p")

############################ PREP WORK #############################
resid_plot <- function (model,df){
  par(mfrow=c(2,2), cex=1.0, mai=c(1.0,1.0,0.6,0.6))
  # Fitted line plot. Predicted y versus observed y.
  plot(df$CH4, fitted(model), xlab="y", 
       ylab="yhat", main="Predicted vs Actual y",
       pch=19)
  abline(a=0,b=1,col="red") 
  plot(fitted(model),resid(model), main="Residual Plot")  
  abline(a=0,b=0, col="red")  # residual plot
  qqnorm(resid(model))     # normality plot
  qqline(resid(model), col=2)
  hist(resid(model), density=10, main="Residuals Distribution",
       col="green", border="black") 
  par(mfrow=c(1,1), cex=1.0, mai=c(1.0,1.0,1.0,1.0) )
  mtext(deparse(model$call[2]),side = 3,line = 4,outer = FALSE)
  }

# Consolidating data into a dataframe
NDVI <- unlist(data$collapsed_landsat$NDVI)
NDWI <- unlist(data$collapsed_landsat$NDWI)
NDMI <- unlist(data$collapsed_landsat$NDMI)
temp <- unlist(data$collapsed_landsat$temp)
# Introducing Month as a factor
month.f <- as.factor(data$collapsed_landsat$month)
CH4 <- unlist(data$collapsed_ffp)
df<-data.frame(NDVI,NDWI,NDMI,temp,month.f,CH4)

# Test transformations [Can't log(CH4)]
pairs(~ CH4+NDVI + NDWI + NDMI + temp + month.f, data=df, main="Scatterplot Matrix for CH4 data",pch = 20)
# Compare original to transformed
par(mfrow=c(2,1), cex=1.0, mai=c(0.8, 1.3, 0.4, 0.8))
plot(df$CH4,df$NDVI,pch=20)
plot(df$CH4,(df$NDVI)^3,pch=8)

# Transforming
df$sqNDVI <- df$NDVI^2
df$cubNDVI <- df$NDVI^3
df$expNDWI <- exp(df$NDWI)
df$invNDWI <- 1/(df$NDWI)
df$sqNDWI <- (df$NDWI)^2
df$sqNDMI <- df$NDMI^2
df$cubNDMI <- df$NDMI^3

############################ LM MODEL TESTING #############################

# Base model, no month.f
model0 <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp+sqNDWI, data=df)
shapiro.test(resid(model0))
resid_plot(model0,df)
anova(model0)

# [4] interactions
model1 <- lme(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
              +month.f*NDVI
              +month.f*NDWI
               +month.f*NDMI
               +month.f*temp
             , data=df)
shapiro.test(resid(model1))
resid_plot(model1,df)

# All transformed
model2 <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp +
               sqNDVI+cubNDVI+expNDWI+invNDWI+sqNDWI+sqNDMI+cubNDMI+
               month.f*NDVI+month.f*NDWI+month.f*NDMI+month.f*temp+
               month.f*sqNDVI+month.f*cubNDVI+month.f*expNDWI+month.f*invNDWI+
               month.f*sqNDWI+month.f*sqNDMI+month.f*cubNDMI
             , data=df)
shapiro.test(resid(model2))
resid_plot(model2,df)
anova(model2)

# [1] Significant significant transforms
model3 <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
              +month.f*NDVI+month.f*NDWI+month.f*NDMI+month.f*temp
              , data=df)
shapiro.test(resid(model3))
anova(model3)

############################ LME MODEL TESTING #############################

lme.model <- lme(CH4 ~ NDVI+NDWI+NDMI+temp
                   +month.f*NDVI+month.f*NDWI+month.f*NDMI+month.f*temp,
                   data=df,
                   random=~1|month.f, method="ML")
shapiro.test(resid(lme.model))
anova(lme.model)
resid_plot(lme.model,df)

############################ MLR Stepwise #############################
# IN-OUT
library(MASS) # Need to load the package MASS first!
model.null <- lm(CH4 ~ 1,data = df) # no x variables, the null model
model.full <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
                 +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                 , data=df) # all x variables


model.step <- stepAIC(model.full,
                      scope = list(upper = ~month.f+NDVI+NDWI+NDMI+temp
                                   +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp,
                                   lower = ~month.f+NDVI+NDWI+NDMI+temp),trace = TRUE)

# SUBSETS
library(leaps) # The leaps package must be loaded first!
CH4.subsets <- regsubsets(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
                          +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                          , data=df,nbest=3)
A<- summary(CH4.subsets) # make sense of the outputs.
results<-data.frame(A$which,A$rsq,A$adjr2,A$bic) #save some of them
names(results) 
# Show the nbest=3 for 1 x variable, best 3 for 2 x variables, etc.
results
C<-order(-results$A.adjr2) # get the order of results by increasing bic
results[C,]  # show the results by bic

# Try OLSRR - regsubset only works with quantitative vars. OLSRR works with categorical vars
library(olsrr)
model.olsrr <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
                  +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp, data=df)
ols_step_best_subset(model.olsrr)

# model3 <- lm(CH4 ~ month.f*NDVI+month.f*NDWI, data = df)
# summary(model3)

############################ MLR Stepwise-Found models #############################
# By backward AIC
model.backward <- lm(CH4 ~ month.f + NDVI + NDWI + month.f*NDVI + month.f*NDWI,data=df)
# By forward AIC
model.forward <- lm(CH4 ~ NDMI + temp,data=df)
# By backward AIC, lower limited by All Fixed Terms
model.forward_allFixed <- lm(CH4 ~ month.f + NDVI + NDWI + NDMI + temp + month.f*NDVI + month.f*NDWI,data=df)

shapiro.test(resid(model.good))
dev.new()
resid_plot(model.good,df)
anova(model.good)
summary(model.good)


# Save Model Results
df_coeff <- coefficients(model.good)
savedir = "/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Model/"
setwd(savedir)
write.csv(df, file = "Young2022_fitted_bestAIC_limAllFixed.csv", row.names = TRUE)
write.csv(df_coeff, file = "Young2022_coeff_bestAIC_limAllFixed.csv", row.names = TRUE)