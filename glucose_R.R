#########################################
#Project Title: Glycemic control in critically ill patients with and without diabetes
#Date created: 8/1/2022
#Date updated: 12/2/2022
#########################################
library("mgcv")
library("dplyr")
library("tidymv")
library("oddsratio")
library("ggplot2")
library("tidyverse")
library("BB")
library("pROC")
library("visreg")
#########################################
# 1. Whole cohort
df<- read_csv("df_uniquepatient_morethan2_nonan.csv")
df$diabetes <- ifelse(df$diabetes==TRUE,"DM","Non DM")
df$diabetes <- factor(df$diabetes)
df$operative <- factor(df$operative)
df$ventday1<- factor(df$ventday1)
df$inotropevasopressor <- factor(df$inotropevasopressor)
df$BMI_cat <- factor(df$BMI_cat)

# 1.1 TWA mean
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean, by=diabetes) + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_twamean",
          by = "diabetes", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, glucose_twamean < 240)
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Time-weighted average glucose (mg/dL)", y="Probability of hospital mortality") +
  scale_x_continuous(breaks = seq(50, 240, by = 10)) 
concurvity(model_all_twamean)

# 1.2 Min
model_all_min<-gam(hospitaldischargestatus ~ s(glucose_min, by=diabetes) + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                   data = df,
                   family = binomial, method="REML")
v<-visreg(model_all_min, xvar = "glucose_min",
          by = "diabetes", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_min > 20) & (glucose_min < 90))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Minimum glucose (mg/dL)", y="Probability of hospital mortality")+
  scale_x_continuous(breaks = seq(20, 90, by = 10))

# 1.3 CV
model_all_cv<-gam(hospitaldischargestatus ~ s(glucose_cv, by=diabetes)+ s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                  data = df,
                  family = binomial, method="REML")
v<-visreg(model_all_cv, xvar = "glucose_cv",
          by = "diabetes", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_cv > 10) & (glucose_cv < 50))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Coefficient of variation (%)", y="Probability of hospital mortality")

#########################################
# 2. Non-DM
df_nonDM<- read_csv("df_uniquepatient_nonDMonly_nonan.csv")
df_nonDM$operative <- factor(df_nonDM$operative)
df_nonDM$BMI_cat <- factor(df_nonDM$BMI_cat)
df_nonDM$ventday1<- factor(df_nonDM$ventday1)
df_nonDM$inotropevasopressor <- factor(df_nonDM$inotropevasopressor)

# 2.1 TWA mean
# Crude
model_nonDM_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean),
                         data = df_nonDM,
                         family = binomial,
                         method='REML')
plot(model_nonDM_twamean,shade=T,scale=0,
     xlab="glucose_twamean", ylab="f(glucose_twamean)")

# Adjusted
model_nonDM_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean) + operative  + s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                         data = df_nonDM,
                         family = binomial)

v_nonDM<-visreg(model_nonDM_twamean, xvar = "glucose_twamean", data= df_nonDM, trans = plogis,  plot=FALSE)
v2_nonDM <- subset(v_nonDM, glucose_twamean < 200)
plot(v2_nonDM, rug=FALSE, overlay=TRUE, gg=TRUE)  + 
  labs(x = "Time-weighted average glucose (mg/dL)", y="Mortality", title="Time-weighted average level  and mortality in Non-DM patients") 

#Applied separately on crude and adjusted models
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(100,60))
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(100,80))
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(100,120))
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(100,160))
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(100,180))
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(120,180))
or_gam(data = df_nonDM, model = model_nonDM_twamean, pred = "glucose_twamean",values = c(150,180))

# 2.2 Min
# Crude
model_nonDM_min<-gam(hospitaldischargestatus ~ s(glucose_min) , data = df_nonDM,
                     family = binomial, method='REML')

# Adjusted
model_nonDM_min<-gam(hospitaldischargestatus ~ s(glucose_min) + operative + s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                     data = df_nonDM,
                     family = binomial, method='REML')
v_nonDM<-visreg(model_nonDM_min, xvar = "glucose_min", by = "operative", data= df_nonDM, trans = plogis,  plot=FALSE)
v2_nonDM <- subset(v_nonDM, (glucose_min > 20) & (glucose_min < 80))
plot(v2_nonDM, rug=FALSE, overlay=TRUE, gg=TRUE)  + 
  labs(x = "Minimum glucose (mg/dL)", y="Mortality", title="Minimum glucose and mortality in Non-DM patients")

#Applied separately on crude and adjusted models
or_gam(data = df_nonDM, model = model_nonDM_min, pred = "glucose_min",values = c(80,70))
or_gam(data = df_nonDM, model = model_nonDM_min, pred = "glucose_min",values = c(80,60))
or_gam(data = df_nonDM, model = model_nonDM_min, pred = "glucose_min",values = c(80,50))
or_gam(data = df_nonDM, model = model_nonDM_min, pred = "glucose_min",values = c(80,40))
or_gam(data = df_nonDM, model = model_nonDM_min, pred = "glucose_min",values = c(80,30))

# 2.3 CV
#Crude
model_nonDM_cv<-gam(hospitaldischargestatus ~ s(glucose_cv), data = df_nonDM,
                    family = binomial, method='REML')

#Adjusted
model_nonDM_cv<-gam(hospitaldischargestatus ~ s(glucose_cv)  +operative+ s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                    data = df_nonDM,
                    family = binomial, method='REML')

#Applied separately on crude and adjusted models
or_gam(data = df_nonDM, model = model_nonDM_cv, pred = "glucose_cv",values = c(10,20))
or_gam(data = df_nonDM, model = model_nonDM_cv, pred = "glucose_cv",values = c(10,30))
or_gam(data = df_nonDM, model = model_nonDM_cv, pred = "glucose_cv",values = c(10,40))
or_gam(data = df_nonDM, model = model_nonDM_cv, pred = "glucose_cv",values = c(10,50))

#########################################
# 3. DM
df_DM<- read_csv("df_uniquepatient_DMonly_nonan.csv")
df_DM$operative <- factor(df_DM$operative)
df_DM$BMI_cat <- factor(df_DM$BMI_cat)
df_DM$ventday1<- factor(df_DM$ventday1)
df_DM$inotropevasopressor <- factor(df_DM$inotropevasopressor)

# 3.1 TWA mean
#Crude
model_DM_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean)  , data = df_DM,
                      family = binomial, method='REML')
#Adjusted
model_DM_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean)  + operative+ s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                      data = df_DM,
                      family = binomial, method='REML')

#Applied separately on crude and adjusted models
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(100,60))
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(100,80))
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(100,120))
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(100,160))
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(100,180))
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(120,180))
or_gam(data = df_DM, model = model_DM_twamean, pred = "glucose_twamean",values = c(150,180))

# 3.2 Min
# Crude
model_DM_min<-gam(hospitaldischargestatus ~ s(glucose_min) , data = df_DM,
                  family = binomial, method='REML')
# Adjusted
model_DM_min<-gam(hospitaldischargestatus ~ s(glucose_min) + operative + s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                  data = df_DM,
                  family = binomial, method='REML')

#Applied separately on crude and adjusted models
or_gam(data = df_DM, model = model_DM_min, pred = "glucose_min",values = c(80,70))
or_gam(data = df_DM, model = model_DM_min, pred = "glucose_min",values = c(80,60))
or_gam(data = df_DM, model = model_DM_min, pred = "glucose_min",values = c(80,50))
or_gam(data = df_DM, model = model_DM_min, pred = "glucose_min",values = c(80,40))
or_gam(data = df_DM, model = model_DM_min, pred = "glucose_min",values = c(80,30))

# 3.3 CV
# Crude
model_DM_cv<-gam(hospitaldischargestatus ~ s(glucose_cv) , data = df_DM,
                 family = binomial, method='REML')
# Adjusted
model_DM_cv<-gam(hospitaldischargestatus ~ s(glucose_cv) + operative  + s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                 data = df_DM,
                 family = binomial, method='REML')

#Applied separately on crude and adjusted models
or_gam(data = df_DM, model = model_DM_cv, pred = "glucose_cv",values = c(10,20))
or_gam(data = df_DM, model = model_DM_cv, pred = "glucose_cv",values = c(10,30))
or_gam(data = df_DM, model = model_DM_cv, pred = "glucose_cv",values = c(10,40))
or_gam(data = df_DM, model = model_DM_cv, pred = "glucose_cv",values = c(10,50))

#################################
#Finding glucose range
# 1 Non-DM
df_train_nonDM<- read_csv("df_uniquepatient_nonDMonly_nonan.csv")
df_train_nonDM$operative <- factor(df_train_nonDM$operative)
df_train_nonDM$BMI_cat <- factor(df_train_nonDM$BMI_cat)
df_train_nonDM$ventday1<- factor(df_train_nonDM$ventday1)
df_train_nonDM$inotropevasopressor <- factor(df_train_nonDM$inotropevasopressor)
model_nonDM_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean) + operative+ s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor , 
                         data = df_train_nonDM,
                         family = binomial, method='REML')
plot(model_nonDM_twamean, select=1, ylim=c(-1,3), xlim=c(50,240), seWithMean = TRUE, rug = FALSE, shade = TRUE, shade.col = "grey", xaxt='n',
     xlab="Time-weighted average glucose (mg/dL)", ylab="f(TWA glucose)")
abline(h=0,lty=2,lwd=0.5)
abline(v=80,lty=2,lwd=0.5)
abline(v=122,lty=2,lwd=0.5)
axis(1, at = seq(50, 240, by = 10), las=2)

# 2 DM
df_train_DM<- read_csv("df_uniquepatient_DMonly_nonan.csv")
df_train_DM$hospitaldischargestatus <- factor(df_train_DM$hospitaldischargestatus)
df_train_DM$admission_type <- factor(df_train_DM$admission_type)
df_train_DM$operative <- factor(df_train_DM$operative)
df_train_DM$BMI_cat <- factor(df_train_DM$BMI_cat)
df_train_DM$ventday1<- factor(df_train_DM$ventday1)
df_train_DM$inotropevasopressor <- factor(df_train_DM$inotropevasopressor)

model_DM_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean)  + operative+ s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor
                      , data = df_train_DM,
                      family = binomial, method='REML')
plot(model_DM_twamean, xlab="Time-weighted average glucose (mg/dL)", ylab="f(TWA glucose)",
     select=1, ylim=c(-1,3), xlim=c(50,240), seWithMean = TRUE, rug = FALSE, shade = TRUE, shade.col = "grey", xaxt='n')
abline(h=0,lty=2,lwd=0.5)
abline(v=85,lty=2,lwd=0.5)
abline(v=148,lty=2,lwd=0.5)
axis(1, at = seq(50, 240, by = 10), las=2)
#################################
#Function calculating AUROC, De Long's test
# 1. Non-DM
df_train_nonDM<- read_csv("X_train_nonDM.csv")
df_test_nonDM<- read_csv("X_test_nonDM.csv")
# 1.1 Crude
glucose_twamean_80120<-cut(df_train_nonDM$glucose_twamean,c(0,80,120,Inf))
glucose_twamean_80120 <- factor(glucose_twamean_80120 , levels=c("(80,120]","(0,80]","(120,Inf]"))
model_80120 <- glm(hospitaldischargestatus ~ glucose_twamean_80120,data=df_train_nonDM,
                   family="binomial")
glucose_twamean_80110<-cut(df_train_nonDM$glucose_twamean,c(0,80,110,Inf))
glucose_twamean_80110 <- factor(glucose_twamean_80110 , levels=c("(80,110]","(0,80]","(110,Inf]"))
model_80110 <- glm(hospitaldischargestatus ~ glucose_twamean_80110,data=df_train_nonDM,
                   family="binomial")
glucose_twamean_180<-cut(df_train_nonDM$glucose_twamean,c(0,180,Inf))
model_180 <- glm(hospitaldischargestatus ~ glucose_twamean_180,data=df_train_nonDM,
                 family="binomial")
glucose_twamean_180200<-cut(df_train_nonDM$glucose_twamean,c(0,180,200,Inf))
glucose_twamean_180200 <- factor(glucose_twamean_180200 , levels=c("(180,200]","(0,180]","(200,Inf]"))
model_180200 <- glm(hospitaldischargestatus ~ glucose_twamean_180200,data=df_train_nonDM,
                    family="binomial")

# 1.2 Adjusted
glucose_twamean_80120<-cut(df_train_nonDM$glucose_twamean,c(0,80,120,Inf))
glucose_twamean_80120 <- factor(glucose_twamean_80120 , levels=c("(80,120]","(0,80]","(120,Inf]"))
df_train_nonDM$BMI_cat <- factor(df_train_nonDM$BMI_cat , levels=c("18.5-25","<18.5","25-<30","30-<35", "35-<40", ">=40"))
model_80120 <- glm(hospitaldischargestatus ~ glucose_twamean_80120+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor 
                   ,data=df_train_nonDM,
                   family="binomial")
summary(model_80120)
exp(model_80120$coefficients)
exp(confint.default(model_80120))

glucose_twamean_80110<-cut(df_train_nonDM$glucose_twamean,c(0,80,110,Inf))
glucose_twamean_80110 <- factor(glucose_twamean_80110 , levels=c("(80,110]","(0,80]","(110,Inf]"))
model_80110 <- glm(hospitaldischargestatus ~ glucose_twamean_80110+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor 
                   ,data=df_train_nonDM,
                   family="binomial")

glucose_twamean_180<-cut(df_train_nonDM$glucose_twamean,c(0,180,Inf))
model_180 <- glm(hospitaldischargestatus ~ glucose_twamean_180+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor 
                 ,data=df_train_nonDM,
                 family="binomial")

glucose_twamean_180200<-cut(df_train_nonDM$glucose_twamean,c(0,180,200,Inf))
glucose_twamean_180200 <- factor(glucose_twamean_180200 , levels=c("(180,200]","(0,180]","(200,Inf]"))
model_180200 <- glm(hospitaldischargestatus ~ glucose_twamean_180200+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor 
                    ,data=df_train_nonDM,
                    family="binomial")

#1.3 Predictions
glucose_twamean_80120<-cut(df_test_nonDM$glucose_twamean,c(0,80,120,Inf))
logit_p<-predict.glm(model_80120 , newdata=df_test_nonDM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_nonDM.80120<-roc(df_test_nonDM$hospitaldischargestatus,p0)
roc_nonDM.80120
ci(roc_nonDM.80120)

glucose_twamean_80110<-cut(df_test_nonDM$glucose_twamean,c(0,80,110,Inf))
logit_p<-predict.glm(model_80110, newdata=df_test_nonDM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_nonDM.80110<-roc(df_test_nonDM$hospitaldischargestatus,p0)
roc_nonDM.80110
ci(roc_nonDM.80110)

glucose_twamean_180<-cut(df_test_nonDM$glucose_twamean,c(0,180,Inf))
logit_p<-predict.glm(model_180, newdata=df_test_nonDM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_nonDM.180<-roc(df_test_nonDM$hospitaldischargestatus,p0)
roc_nonDM.180
ci(roc_nonDM.180)

glucose_twamean_180200<-cut(df_test_nonDM$glucose_twamean,c(0,180,200,Inf))
logit_p<-predict.glm(model_180200, newdata=df_test_nonDM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_nonDM.180200<-roc(df_test_nonDM$hospitaldischargestatus,p0)
roc_nonDM.180200
ci(roc_nonDM.180200)

# 1.4 Delong's test
roc.test(roc_nonDM.80120,roc_nonDM.80110)
roc.test(roc_nonDM.80120,roc_nonDM.180)
roc.test(roc_nonDM.80120,roc_nonDM.180200)

################################
## 2. DM
df_train_DM<- read_csv("X_train_DM.csv")
df_test_DM<- read_csv("X_test_DM.csv")
# 2.1 Crude
glucose_twamean_cut<-cut(df_train_DM$glucose_twamean,c(0,90,150,Inf))
glucose_twamean_cut <- factor(glucose_twamean_cut , levels=c("(90,150]","(0,90]","(150,Inf]"))
df_train_DM$BMI_cat <- factor(df_train_DM$BMI_cat , levels=c("18.5-25","<18.5","25-<30","30-<35", "35-<40", ">=40"))
glmmodel_DM <- glm(hospitaldischargestatus ~ glucose_twamean_cut,data=df_train_DM, family="binomial")
glucose_twamean_80110<-cut(df_train_DM$glucose_twamean,c(0,80,110,Inf))
model_80110 <- glm(hospitaldischargestatus ~ glucose_twamean_80110,data=df_train_DM,family="binomial")
glucose_twamean_180<-cut(df_train_DM$glucose_twamean,c(0,180,Inf))
model_180 <- glm(hospitaldischargestatus ~ glucose_twamean_180,data=df_train_DM, family="binomial")
glucose_twamean_180200<-cut(df_train_DM$glucose_twamean,c(0,180, 200,Inf))
model_180200 <- glm(hospitaldischargestatus ~ glucose_twamean_180,data=df_train_DM,family="binomial")

#2.2 Adjusted
glucose_twamean_cut<-cut(df_train_DM$glucose_twamean,c(0,90,150,Inf))
glucose_twamean_cut <- factor(glucose_twamean_cut , levels=c("(90,150]","(0,90]","(150,Inf]"))
df_train_DM$BMI_cat <- factor(df_train_DM$BMI_cat , levels=c("18.5-25","<18.5","25-<30","30-<35", "35-<40", ">=40"))
glmmodel_DM <- glm(hospitaldischargestatus ~ glucose_twamean_cut+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor 
                   ,data=df_train_DM,
                   family="binomial")
summary(glmmodel_DM)
exp(glmmodel_DM$coefficients)
exp(confint.default(glmmodel_DM))

glucose_twamean_80110<-cut(df_train_DM$glucose_twamean,c(0,80,110,Inf))
glucose_twamean_80110 <- factor(glucose_twamean_80110 , levels=c("(80,110]","(0,80]","(110,Inf]"))
model_80110 <- glm(hospitaldischargestatus ~ glucose_twamean_80110 + operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor
                   ,data=df_train_DM,
                   family="binomial")

glucose_twamean_180<-cut(df_train_DM$glucose_twamean,c(0,180,Inf))
model_180 <- glm(hospitaldischargestatus ~ glucose_twamean_180+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor
                 ,data=df_train_DM,
                 family="binomial")

glucose_twamean_180200<-cut(df_train_DM$glucose_twamean,c(0,180, 200,Inf))
glucose_twamean_180200 <- factor(glucose_twamean_180200 , levels=c("(80,200]","(0,180]","(200,Inf]"))
model_180200 <- glm(hospitaldischargestatus ~ glucose_twamean_180+ operative+ age + apachescore + BMI_cat + ventday1 + inotropevasopressor
                    ,data=df_train_DM,
                    family="binomial")

glucose_twamean_cut<-cut(df_test_DM$glucose_twamean,c(0,90,150,Inf))
glucose_twamean_80110<-cut(df_test_DM$glucose_twamean,c(0,80,110,Inf))
glucose_twamean_180<-cut(df_test_DM$glucose_twamean,c(0,180,Inf))
glucose_twamean_180200<-cut(df_test_DM$glucose_twamean,c(0,180, 200,Inf))

#2.3 Predictions
logit_p<-predict.glm(glmmodel_DM, newdata=df_test_DM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_DM.gam<-roc(df_test_DM$hospitaldischargestatus,p0)
roc_DM.gam
ci(roc_DM.gam)

logit_p<-predict.glm(model_80110, newdata=df_test_DM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_DM.80110<-roc(df_test_DM$hospitaldischargestatus,p0)
roc_DM.80110
ci(roc_DM.80110)

logit_p<-predict.glm(model_180, newdata=df_test_DM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_DM.180<-roc(df_test_DM$hospitaldischargestatus,p0)
roc_DM.180
ci(roc_DM.180)

logit_p<-predict.glm(model_180200, newdata=df_test_DM)
p0<-exp(logit_p)/(1+exp(logit_p))
roc_DM.180200<-roc(df_test_DM$hospitaldischargestatus,p0)
roc_DM.180200
ci(roc_DM.180200)

#2.4 De Long's test
roc.test(roc_DM.gam,roc_DM.80110)
roc.test(roc_DM.gam,roc_DM.180)
roc.test(roc_DM.gam,roc_DM.180200)

#################################
#Subgroup analysis: medical vs. surgical
df<- read_csv("df_uniquepatient_morethan2_nonan.csv")
df$diabetes <- factor(df$diabetes)
df$operative <- ifelse(df$operative==1,"Surgical","Medical")
df$operative <- factor(df$operative)
df$ventday1<- factor(df$ventday1)
df$inotropevasopressor <- factor(df$inotropevasopressor)
df$BMI_cat <- factor(df$BMI_cat)
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean, by=operative) + diabetes + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_twamean",
          by = "operative", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, glucose_twamean < 240)
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Time-weighted average glucose (mg/dL)", y="Probability of hospital mortality") +
  scale_x_continuous(breaks = seq(50, 240, by = 10))

#Min
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_min, by=operative) + diabetes + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_min",
          by = "operative", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, (glucose_min > 20) & (glucose_min < 80))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Minimum glucose (mg/dL)", y="Probability of hospital mortality")+
  scale_x_continuous(breaks = seq(20, 80, by = 10))

#CV
model_all_cv<-gam(hospitaldischargestatus ~ s(glucose_cv, by=operative)+ diabetes+ s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                  data = df,
                  family = binomial, method="REML")
v<-visreg(model_all_cv, xvar = "glucose_cv",
          by = "operative", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_cv > 10) & (glucose_cv < 50))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Coefficient of variation (%)", y="Probability of hospital mortality")
##########################################
#Subgroup analysis: trauma vs. non_trauma
df<- read_csv("df_uniquepatient_morethan2_nonan.csv")
df$diabetes <- factor(df$diabetes)
df$trauma <- ifelse(df$trauma==1,"Trauma","Non-trauma")
df$trauma <- factor(df$trauma)
df$operative<- factor(df$operative)
df$ventday1<- factor(df$ventday1)
df$inotropevasopressor <- factor(df$inotropevasopressor)
df$BMI_cat <- factor(df$BMI_cat)
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean, by=trauma) + diabetes + s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_twamean",
          by = "trauma", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, glucose_twamean < 240)
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Time-weighted average glucose (mg/dL)", y="Probability of hospital mortality") +
  scale_x_continuous(breaks = seq(50, 240, by = 10))

#Min
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_min, by=trauma) + diabetes + s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_min",
          by = "trauma", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, (glucose_min > 20) & (glucose_min < 80))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Minimum glucose (mg/dL)", y="Probability of hospital mortality")+
  scale_x_continuous(breaks = seq(20, 80, by = 10))

#CV
model_all_cv<-gam(hospitaldischargestatus ~ s(glucose_cv, by=trauma)+ diabetes+ s(age) + s(apachescore) + BMI_cat + ventday1 + inotropevasopressor ,
                  data = df,
                  family = binomial, method="REML")
v<-visreg(model_all_cv, xvar = "glucose_cv",
          by = "trauma", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_cv > 10) & (glucose_cv < 50))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Coefficient of variation (%)", y="Probability of hospital mortality")

###################################################
#Subgroup analysis: steroid vs. no_steroid
df<- read_csv("df_uniquepatient_morethan2_nonan.csv")
df$diabetes <- factor(df$diabetes)
df$steroid_use <- ifelse(df$steroid_use==1,"Steroid","No steroid")
df$steroid_use <- factor(df$steroid_use)
df$operative <- factor(df$operative)
df$ventday1<- factor(df$ventday1)
df$inotropevasopressor <- factor(df$inotropevasopressor)
df$BMI_cat <- factor(df$BMI_cat)
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean, by=steroid_use) + diabetes + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_twamean",
          by = "steroid_use", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, glucose_twamean < 240)
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Time-weighted average glucose (mg/dL)", y="Probability of hospital mortality") +
  scale_x_continuous(breaks = seq(50, 240, by = 10))

#Min
model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_min, by=steroid_use) + diabetes + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_min",
          by = "operative", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, (glucose_min > 20) & (glucose_min < 80))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Minimum glucose (mg/dL)", y="Probability of hospital mortality")+
  scale_x_continuous(breaks = seq(20, 80, by = 10))

#CV
model_all_cv<-gam(hospitaldischargestatus ~ s(glucose_cv, by=steroid_use)+ diabetes+ s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                  data = df,
                  family = binomial, method="REML")
v<-visreg(model_all_cv, xvar = "glucose_cv",
          by = "operative", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_cv > 10) & (glucose_cv < 50))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Coefficient of variation (%)", y="Probability of hospital mortality")


#################################
#Sensitivity analysis: IDDM vs. NIDDM
df_DMstatusonly<-read_csv("df_DMstatusonly_nonan.csv")
df_DMstatusonly$DM_status <- ifelse(df_DMstatusonly$DM_status=="NIDDM_diet","NIDDM on diet",
                                    (ifelse(df_DMstatusonly$DM_status=="NIDDM_OHA", "NIDDM on OHA",
                                            (ifelse(df_DMstatusonly$DM_status=="No_DM", "Non DM",
                                                    (ifelse(df_DMstatusonly$DM_status=="IDDM", "IDDM", "")))))))
df_DMstatusonly$DM_status <- factor(df_DMstatusonly$DM_status)
df_DMstatusonly$operative <- factor(df_DMstatusonly$operative)
df_DMstatusonly$ventday1<- factor(df_DMstatusonly$ventday1)
df_DMstatusonly$inotropevasopressor <- factor(df_DMstatusonly$inotropevasopressor)
df_DMstatusonly$BMI_cat <- factor(df_DMstatusonly$BMI_cat)

DMstatus_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean, by=DM_status) +  s(age) + s(apachescore) + operative+ ventday1 + inotropevasopressor,
                      data = df_DMstatusonly,
                      family = binomial, method="REML")
v<-visreg(DMstatus_twamean, xvar = "glucose_twamean",
          by = "DM_status", data= df_DMstatusonly, trans = plogis,  plot=FALSE)
v2 <- subset(v, glucose_twamean < 240)
plot(v2, rug=FALSE, overlay=TRUE, gg=TRUE) +theme(legend.title = element_blank())+
  labs(x = "Time-weighted average glucose (mg/dL)", y="Probability of hospital mortality")

DMstatus_min<-gam(hospitaldischargestatus ~ s(glucose_min, by=DM_status) + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor,
                  data = df_DMstatusonly,
                  family = binomial, method="REML")
v<-visreg(DMstatus_min, xvar = "glucose_min",
          by = "DM_status", data= df_DMstatusonly, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_min > 20) & (glucose_min < 80))
plot(v2, rug=FALSE, overlay=TRUE, gg=TRUE) +theme(legend.title = element_blank())+
  labs(x = "Minimum glucose (mg/dL)", y="Probability of hospital mortality")+
  scale_x_continuous(breaks = seq(20, 80, by = 10))

DMstatus_cv<-gam(hospitaldischargestatus ~ s(glucose_cv, by=DM_status) + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor,
                 data = df_DMstatusonly,
                 family = binomial, method="REML")
v<-visreg(DMstatus_cv, xvar = "glucose_cv",
          by = "DM_status", data= df_DMstatusonly, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_cv > 10) & (glucose_cv < 50))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Coefficient of variation (%)", y="Probability of hospital mortality")
#################################
#Sensitivity analysis: all patients (including LOS<2 days)
df<- read_csv("df_uniquepatient_all_nonan.csv")
df$diabetes <- ifelse(df$diabetes==TRUE,"DM","Non DM")
df$diabetes <- factor(df$diabetes)
df$operative <- factor(df$operative)
df$ventday1<- factor(df$ventday1)
df$inotropevasopressor <- factor(df$inotropevasopressor)
df$BMI_cat <- factor(df$BMI_cat)

model_all_twamean<-gam(hospitaldischargestatus ~ s(glucose_twamean, by=diabetes) + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                       data = df,
                       family = binomial, method="REML")
v<-visreg(model_all_twamean, xvar = "glucose_twamean",
          by = "diabetes", data= df,scale="response",  plot=FALSE)
v2 <- subset(v, (glucose_twamean > 50) & (glucose_twamean < 240))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Time-weighted average glucose (mg/dL)", y="Probability of hospital mortality") +
  scale_x_continuous(breaks = seq(50, 240, by = 10))

concurvity(model_all_twamean)

model_all_min<-gam(hospitaldischargestatus ~ s(glucose_min, by=diabetes) + s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                   data = df,
                   family = binomial, method="REML")
v<-visreg(model_all_min, xvar = "glucose_min",
          by = "diabetes", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_min > 20) & (glucose_min < 90))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Minimum glucose (mg/dL)", y="Probability of hospital mortality")+
  scale_x_continuous(breaks = seq(20, 90, by = 10))

model_all_cv<-gam(hospitaldischargestatus ~ s(glucose_cv, by=diabetes)+ s(age) + s(apachescore) + operative+ BMI_cat + ventday1 + inotropevasopressor ,
                  data = df,
                  family = binomial, method="REML")
v<-visreg(model_all_cv, xvar = "glucose_cv",
          by = "diabetes", data= df, trans = plogis,  plot=FALSE)
v2 <- subset(v, (glucose_cv > 10) & (glucose_cv < 50))
plot(v2, overlay=TRUE, rug=FALSE, gg=TRUE)  +theme(legend.title = element_blank())+
  labs(x = "Coefficient of variation (%)", y="Probability of hospital mortality")