rm(list=ls())
library(reticulate)
library(nlme)
pd <- import("pandas")

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
############################ PREP WORK #############################

# Loading data
site_name = 'Young'
dataPath = '/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Compilation/'
data <- pd$read_pickle(paste(dataPath,site_name,"_Annual/21to23_growing_season/",site_name,"_compiled_data.p",sep=""))

# Consolidating data into a dataframe
NDVI <- unlist(data$collapsed_landsat$NDVI)
NDWI <- unlist(data$collapsed_landsat$NDWI)
NDMI <- unlist(data$collapsed_landsat$NDMI)
temp <- unlist(data$collapsed_landsat$temp)
# Introducing Month and Year as factors
month.f <- as.factor(data$collapsed_landsat$month)
year.f <- as.factor(data$collapsed_landsat$year)
CH4 <- unlist(data$collapsed_ffp)
df<-data.frame(NDVI,NDWI,NDMI,temp,month.f,CH4)
df<-data.frame(NDVI,NDWI,NDMI,temp,month.f,year.f,CH4)

############ Exclude data from model ########################
data$collapsed_landsat$period
exclude_date <- c(20210925)
keep_idx <- !(1:length(data$collapsed_landsat$NDVI) %in% which(data$collapsed_landsat$period %in% exclude_date))
# Store data not in these dates
NDVI <- unlist(data$collapsed_landsat$NDVI[keep_idx])
NDWI <- unlist(data$collapsed_landsat$NDWI[keep_idx])
NDMI <- unlist(data$collapsed_landsat$NDMI[keep_idx])
temp <- unlist(data$collapsed_landsat$temp[keep_idx])
month.f <- as.factor(data$collapsed_landsat$month[keep_idx])
year.f <- as.factor(data$collapsed_landsat$year[keep_idx])
CH4 <- unlist(data$collapsed_ffp[keep_idx])
removed_L8 <- rep(exclude_date,times=length(CH4))
df_reduced<-data.frame(NDVI,NDWI,NDMI,temp,month.f,year.f,CH4,removed_L8)

# CHECK HERE IF df IS CORRECT
model_df = df_reduced
model.full <- lm(CH4 ~ month.f+year.f+NDVI+NDWI+NDMI+temp
                 +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                 +year.f*NDVI+ year.f*NDWI+ year.f*NDMI+ year.f*temp
                 , data=model_df) # all x variables
model.step <- stepAIC(model.full,
                      scope = list(upper = ~month.f+year.f+NDVI+NDWI+NDMI+temp
                                   +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                                   +year.f*NDVI+ year.f*NDWI+ year.f*NDMI+ year.f*temp,
                                   lower = ~month.f+year.f+NDVI+NDWI+NDMI+temp),trace = TRUE)
model.good <- lm(formula=formula(model.step),data=model_df)
model.good <- lm(CH4 ~ month.f+year.f+NDVI+NDWI+NDMI+temp+
                   month.f*NDWI + month.f*NDMI + year.f*temp,data=model_df)

shapiro.test(resid(model.good))
resid_plot(model.good,model_df)
summary(model.good)

# SAVE
save_model <-model.good
save_df <- model_df
df_coeff <- coefficients(save_model)
save_df$fitted <- fitted(save_model)

savedir = "/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Model/"
setwd(savedir)
write.csv(save_df, file = "Young2factor_fitted_trained_20210925.csv", row.names = TRUE)
write.csv(df_coeff, file = "Young2factor_coeff_trained_20210925.csv", row.names = TRUE)
print('Saved')
############ ^^ Exclude data from model ^^ ########################


###### Test factor --> significantly different means ######
tapply(df$CH4, df$month.f, mean) # simple means
par(mfrow=c(2,1), cex=1, mai=c(0.8,0.8,0.8,0.8))
boxplot(CH4 ~ month.f,data=df, col='pink',main="CH4 by Month")
tapply(df$CH4, df$year.f, mean)
boxplot(CH4 ~ year.f, data=df, col='pink',main="CH4 by Year")
par(mfrow=c(1,1))
###### Test factor --> significantly different means ######

pairs(~ CH4+NDVI + NDWI + NDMI + temp + month.f+year.f, data=df, main="Scatterplot Matrix for CH4 data",pch = 20)

# Multicollinearity
predictors <- df[, c("NDVI", "NDWI", "NDMI", "temp")]
cor_matrix <- cor(predictors)
vif <- sapply(predictors, function(x) vif(lm(x ~ ., data = predictors)))


########################## INCLUDE YEAR.F #######################
model.full <- lm(CH4 ~ month.f+year.f+NDVI+NDWI+NDMI+temp
              +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
              +year.f*NDVI+ year.f*NDWI+ year.f*NDMI+ year.f*temp
              , data=df) # all x variables

model.step <- stepAIC(model.full,
                      scope = list(upper = ~month.f+year.f+NDVI+NDWI+NDMI+temp
                                   +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                                   +year.f*NDVI+ year.f*NDWI+ year.f*NDMI+ year.f*temp,
                                   lower = ~month.f+year.f+NDVI+NDWI+NDMI+temp),trace = TRUE)

shapiro.test(resid(model.full))
resid_plot(model.full,df_reduced)
summary(model.full)


############################ MLR Stepwise #############################
# IN-OUT
library(MASS) # Need to load the package MASS first!
model.null <- lm(CH4 ~ 1,data = df) # no x variables, the null model
model.full <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
                 +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                 , data=df) # all x variables


shapiro.test(resid(model.good))
dev.new()
resid_plot(model.good,df)
anova(model.good)
summary(model.good)

# # Try OLSRR - regsubset only works with quantitative vars. OLSRR works with categorical vars
# library(olsrr)
# model.olsrr <- lm(CH4 ~ month.f+NDVI+NDWI+NDMI+temp
#                   +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp, data=df)
# ols_step_best_subset(model.olsrr)
# 
# # model3 <- lm(CH4 ~ month.f*NDVI+month.f*NDWI, data = df)
# # summary(model3)




#################### Mixed Model, repeated in time, TEMPORAL AUTOCORRELATION ###################
library(nlme)
y.LME0<-gls(CH4~month.f+year.f+NDVI+NDWI+NDMI+temp+month.f*NDVI+month.f*NDWI+month.f*NDMI+month.f*temp
            +year.f*NDVI+year.f*NDWI+year.f*NDMI+year.f*temp,data=df,method="REML") 

y.LME0<-gls(CH4~month.f+year.f+NDVI+NDWI+NDMI+temp+month.f*NDWI+month.f*NDMI
            +year.f*temp,data=df,method="REML") 

corStruct <- corCAR1(form = ~ 1 | month.f/year.f)
# Fit the GLS model
model <- gls(response ~ NDVI + NDWI + NDMI + temp, 
             data = data,
             correlation = corStruct)


model.step <- stepAIC(y.LME2,
                      scope = list(upper = ~month.f+year.f+NDVI+NDWI+NDMI+temp
                                   +month.f*NDVI+ month.f*NDWI+ month.f*NDMI+ month.f*temp
                                   +year.f*NDVI+ year.f*NDWI+ year.f*NDMI+ year.f*temp,
                                   lower = ~1)
                      ,direction='both',trace = TRUE)

anova(y.LME0)
summary( y.LME0) 
# graph the y versus yhat for each month
plot(y.LME0,CH4~fitted(.)|month.f, abline=c(0,1))
# graph the y versus yhat for each year
plot(y.LME0,CH4~fitted(.)|year.f, abline=c(0,1))

## MIXED MODEL 2: One level not hierarchical.  However, CAR1 within Person
## and assuming this is the same for all persons. Note that AR1 could
## be used here, since the measures are 2 years apart (regular in time).
y.LME2<-update(y.LME0, correlation = corCAR1(form = ~ NDVI+NDWI+NDMI | month.f))
y.LME2  # get the estimated coefficients for the LME model
anova(y.LME2)
summary(y.LME2)

df$yhat.LME2.0<-fitted( y.LME2,level=0)  # estimated population averaged yhat
df$resid.LME2.0<-resid( y.LME2,level=0)  # estimated residuals
df$stdresid.LME2.0<-resid( y.LME2,level=0, type="n")  # so-called normalized residuals

shapiro.test(resid( y.LME2,level=0))
# y.LME2 has better logLik and AIC than y.LME3. No need to further complicate.
y.LME3<-update(y.LME0, correlation = corCompSymm(form = ~ 1 | month.f))  
shapiro.test(resid( y.LME3,level=0))

par(mfrow=c(2,2),mai=c(0.6,0.6,0.6,0.6),cex=0.55)
plot(df$yhat.LME2.0,df$stdresid.LME2.0, main="LME Model 2 Standardized Residuals Plot",
     xlab="yhat", ylab="residual")
plot(df$CH4,df$yhat.LME2.0, main="LME Model 2, Fitted line plot",
     ylab="yhat", xlab="y")
qqnorm(df$stdresid.LME2.0, main="LME Model 2 Standardized residuals, Normality plot")
abline(0, 1,lwd=2,lty=2,col='red')
hist(df$stdresid.LME2.0, breaks =8 , density=10,col="green", border="black",
     main="LME Model 2 , Error Distribution")
par(mfrow=c(1,1),mai=c(1.0,1.0,1.0,1.0),cex=1.0)


# SAVE
save_model <- y.LME2
save_df <- df
df_coeff <- coefficients(save_model)
save_df$fitted <- fitted(save_model,level=0)

savedir = "/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Model/"
setwd(savedir)
write.csv(save_df, file = "Young2factor_fitted_autocorr_bestAIC.csv", row.names = TRUE)
write.csv(df_coeff, file = "Young2factor_coeff_autocorr_bestAIC.csv", row.names = TRUE)
print('Saved')




############################ MLR Stepwise-Found models #############################

# By backward AIC
model.backward <- lm(CH4 ~ month.f + NDVI + NDWI + month.f*NDVI + month.f*NDWI,data=df)
lme.model.backward <- lme(CH4 ~ NDVI + NDWI + month.f*NDVI + month.f*NDWI,data=df, random=~1|month.f,method='ML')

# By forward AIC
model.forward <- lm(CH4 ~ NDMI + temp,data=df)

# By backward AIC, lower limited by All Fixed Terms
model.backward_allFixed <- lm(CH4 ~ month.f + NDVI + NDWI + NDMI + temp + month.f*NDVI + month.f*NDWI,data=df)
lme.model.backward_allFixed <- lme(CH4 ~ NDVI + NDWI + NDMI + temp + month.f*NDVI + month.f*NDWI,data=df, random=~1|month.f,method='ML')

# By backward AIC, must include NDVI,NDWI,NDMI only
model.backward_allFixednotemp <- lm(CH4 ~ month.f + NDVI + NDWI + NDMI + month.f*NDVI + month.f*NDWI,data=df)
lme.model.backward_allFixednotemp <- lme(CH4 ~ NDVI + NDWI + NDMI + month.f*NDVI + month.f*NDWI,data=df, random=~1|month.f,method='ML')

#2 factor: backward AIC (best)
model.backward.2factor <- lm(CH4 ~ month.f + year.f + NDVI + NDWI+ NDMI+ temp 
                             + month.f*NDVI+ month.f*NDWI + year.f*temp,data=df)

#2 factor (reduced): backward AIC (best)
model.backward.2factor.training <- lm(CH4 ~ month.f + year.f + NDVI + NDWI + month.f:NDVI + month.f:NDWI
                                     ,data=df_reduced)

model.good <-model.backward.2factor
shapiro.test(resid(model.good))
dev.new()
resid_plot(model.good,df)
anova(model.good)
summary(model.good)
anova(model.good,model.full)
anova(lme.model.full, lme.model.backward, lme.model.backward_allFixed, lme.model.backward_allFixednotemp)

# Save Model Results
save_df <- df
save_model <-model.full
# Make sure model matches its df
if (length(save_df$CH4) == length(fitted(save_model))){
  df_coeff <- coefficients(save_model)
  # MAKE SURE THIS IS CORRECT df!!
  save_df$fitted <- fitted(save_model)
  savedir = "/Users/darianng/Documents/MSc_Geography/MSc Thesis/Data/Model/"
  setwd(savedir)
  write.csv(save_df, file = "Young2factor_fitted_full.csv", row.names = TRUE)
  write.csv(df_coeff, file = "Young2factor_coeff_full.csv", row.names = TRUE)
  print('Saved')
}