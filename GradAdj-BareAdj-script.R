## This scripts contains all the descriptive analyses on bare adjectives, and some basic inferential stats, but no modelling.

library(tidyverse)
library(multidplyr) # To speed things up
library(lme4)

options(tibble.width=Inf)


# Set up cluster for mutlidplyr:
cluster <- new_cluster(parallel::detectCores()-1)
cluster_library(cluster, "tidyverse")


# Color palettes for graphs:
TwoColorPalette = c(rgb(.5,0,.7),rgb(.9,.65,0))
FourColorPalette = c("#F5793A","#A95AA1","#85C0F9","#0F2080")


# Load and clean the data:
source("GradAdj-load-data.R")

## Define a few useful functions ##

# Standard error function for graphs
se <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))}

# Find the first non-NA value, or return NA if there is only NA:
FindValue=function(x)x[which.max(!is.na(x))]

custom_qqplot <- function(x,col="red"){qqnorm(x);qqline(x,col=col)}

# Define a prettier ScaleType factor for graph labels:
AdjData <- AdjData %>% mutate(
  Scale=factor(ScaleType2,
               levels=c("Relative","MinMax","Min","Max"),
               labels=c("Unbound","Double bound","Lower-bound","Upper-bound"))
)

####################################################
# Analysis 1: Fit sigmoids and analyse parameters
#     Check effect of negation on threshold
####################################################

# Sigmoid fitting function:
SigmoidParamFit <- function(.data){
  scales=unique(.data$ScaleType)
  ScaleNum=as.numeric(factor(.data$ScaleType))
  # Init:
  par=c(3.2,rep(.05,length(scales)),rep(50,length(scales)))
  names(par)=c("sigmaSigmoid",paste0(scales,"kSigmoid"),paste0(scales,"muSigmoid"))
  negloglik <- function(par){
    Pred <- 100/(1+exp(-par[1+ScaleNum]*(100*.data$Degree-par[1+length(scales)+ScaleNum])))
    -sum(case_when(
      .data$PolValue==0 ~ pnorm(0,Pred,exp(par[1]),log.p = T),
      .data$PolValue==100 ~ pnorm(100,Pred,exp(par[1]),lower.tail = F,log.p = T),
      T ~ dnorm(.data$PolValue,Pred,exp(par[1]),log = T)
    ))
  }
  summarize(.data,data.frame(t(optim(par,negloglik,control = list(parscale=c(25,rep(.03,length(scales)),rep(50,length(scales))),ndeps=1e-4))$par)))
}


# Censored Sigmoid fitting function:
CensoredSigmoidParamFit <- function(.data){
  scales=unique(.data$ScaleType)
  ScaleNum=as.numeric(factor(.data$ScaleType))
  # Init:
  par=c(3.2,rep(.05,length(scales)),rep(50,length(scales)))
  names(par)=c("sigmaCenSig",paste0(scales,"kCenSig"),paste0(scales,"muCenSig"))
  negloglik <- function(par){
    Pred <- case_when(.data$Degree==0~0,
                      .data$Degree==1~100,
                      T~100/(1+exp(-par[1+ScaleNum]*(100*.data$Degree-par[1+length(scales)+ScaleNum]))))
    -sum(case_when(
      .data$PolValue==0 ~ pnorm(0,Pred,exp(par[1]),log.p = T),
      .data$PolValue==100 ~ pnorm(100,Pred,exp(par[1]),lower.tail = F,log.p = T),
      T ~ dnorm(.data$PolValue,Pred,exp(par[1]),log = T)
    ))
  }
  summarize(.data,data.frame(t(optim(par,negloglik,control = list(parscale=c(25,rep(.03,length(scales)),rep(50,length(scales))),ndeps=1e-4))$par)))
}


cluster_copy(cluster, c("SigmoidParamFit","CensoredSigmoidParamFit"))

sigmoidFits <- AdjData %>%
  group_by(Results.index,ScaleType,Polarity) %>%
  partition(cluster) %>%
  do(SigmoidParamFit(.)) %>%
  collect() %>%
  ungroup()

cenSigFits <- AdjData %>%
  group_by(Results.index,ScaleType,Polarity) %>%
  partition(cluster) %>%
  do(CensoredSigmoidParamFit(.)) %>%
  collect() %>%
  ungroup()

AdjDataTMP <- AdjData %>%
  select(-ends_with("Sigmoid"),-ends_with("CenSig")) %>% # just for safety in case the code is run in wrong order
  # Merge information about sigmoids:
  left_join(sigmoidFits,by = c("Results.index", "ScaleType", "Polarity")) %>%
  left_join(cenSigFits,by = c("Results.index", "ScaleType", "Polarity")) %>%
  rowwise() %>%
  mutate(
    kSigmoid=FindValue(c_across(ends_with("kSigmoid"))),
    muSigmoid=FindValue(c_across(ends_with("muSigmoid"))),
    kCenSig=FindValue(c_across(ends_with("kCenSig"))),
    muCenSig=FindValue(c_across(ends_with("muCenSig")))
  ) %>%
  ungroup() %>%
  select(-starts_with("Relative"),-starts_with("Min",ignore.case = F),-starts_with("Max",ignore.case = F)) %>%
  # Find best sigmoid that best describes each block of data (censored or not)
  mutate(PredSigmoid = 100/(1+exp(-kSigmoid*(100*Degree-muSigmoid))),
         PredCenSig = case_when(Degree==0~0,
                                Degree==1~100,
                                T~100/(1+exp(-kCenSig*(100*Degree-muCenSig)))),
         bestPred=if_else(sigmaCenSig<sigmaSigmoid,PredCenSig,PredSigmoid),
         priorMean=(1-p1-p0)*100*(a/(a+b))+100*p1)

# t-test on threshold midpoint:
TTestData <- AdjDataTMP %>%
  mutate(bestMu=if_else(sigmaCenSig<sigmaSigmoid,muCenSig,muSigmoid)) %>%
  filter(bestMu> -50&bestMu< 150) %>% # Remove about 7% data on each side (outliers, usually flat distributions).
  group_by(Results.index,ScaleType) %>%
  summarize(muAff=first(bestMu[Polarity=="Aff"]),
            muNeg=first(bestMu[Polarity=="Neg"]))

t.test(x=TTestData$muAff,y=TTestData$muNeg,paired=T)


# t-test on threshold steepness:
TTestData <- AdjDataTMP %>%
  mutate(bestk=if_else(sigmaCenSig<sigmaSigmoid,kCenSig,kSigmoid)) %>%
  filter(bestk>0 & bestk<.5) %>% # Remove about 2.5% data on each side.
  group_by(Results.index,ScaleType) %>%
  summarize(kAff=first(bestk[Polarity=="Aff"]),
            kNeg=first(bestk[Polarity=="Neg"]))

t.test(x=TTestData$kAff,y=TTestData$kNeg,paired=T)


## Since negation has no effect, refit sigmoid on pooled data ##
rm(AdjDataTMP)
sigmoidFits <- AdjData %>%
  group_by(Results.index,ScaleType) %>%
  partition(cluster) %>%
  do(SigmoidParamFit(.)) %>%
  collect() %>%
  ungroup()
cenSigFits <- AdjData %>%
  group_by(Results.index,ScaleType) %>%
  partition(cluster) %>%
  do(CensoredSigmoidParamFit(.)) %>%
  collect() %>%
  ungroup()

AdjData <- AdjData %>%
  select(-ends_with("Sigmoid"),-ends_with("CenSig")) %>% # just for safety in case the code is run in wrong order
  # Merge information about sigmoids:
  left_join(sigmoidFits,by = c("Results.index", "ScaleType")) %>%
  left_join(cenSigFits,by = c("Results.index", "ScaleType")) %>%
  rowwise() %>%
  mutate(
    kSigmoid=FindValue(c_across(ends_with("kSigmoid"))),
    muSigmoid=FindValue(c_across(ends_with("muSigmoid"))),
    kCenSig=FindValue(c_across(ends_with("kCenSig"))),
    muCenSig=FindValue(c_across(ends_with("muCenSig")))
  ) %>%
  ungroup() %>%
  select(-starts_with("Relative"),-starts_with("Min",ignore.case = F),-starts_with("Max",ignore.case = F)) %>%
  # Find best sigmoid that best describes each block of data (censored or not)
  mutate(
    bestCensored = (sigmaCenSig<sigmaSigmoid),
    bestPred=if_else(bestCensored,
                     case_when(Degree==0~0,
                               Degree==1~100,
                               T~100/(1+exp(-kCenSig*(100*Degree-muCenSig)))),
                     100/(1+exp(-kSigmoid*(100*Degree-muSigmoid)))),
    priorMean=(1-p1-p0)*100*(a/(a+b))+100*p1
    )


#######################################
# Define a measure of "absoluteness":
#    Slopes near extreme degrees 
#######################################

FindSlope <- function(values,degrees,n=3){
  tmp <- tibble(degrees,values) %>%
    group_by(degrees) %>%
    summarize(val = mean(values)/100)
  Slope1 = as.numeric(coefficients(lm(val~degrees,data=head(tmp,n)))[2])
  Slope2 = as.numeric(coefficients(lm(val~degrees,data=tail(tmp,n)))[2])
  return(tibble(Slope1,Slope2))
}

# Example with one participant for the slides:
AdjData %>%
  filter(Results.index==106) %>%
  mutate(NormDegree=(Degree-minDeg)/(maxDeg-minDeg)) %>%
  group_by(NormDegree,ScaleType,Scale) %>%
  summarize(Response = mean(PolValue)) %>%
  ungroup() %>%
  ggplot(aes(x=NormDegree,y=Response,col=Scale)) +
  facet_wrap(.~Scale) +
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_color_manual(name="Distribution",values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Normalized Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))

AdjData %>%
  filter(Results.index==106) %>%
  mutate(NormDegree=((Degree-minDeg)/(maxDeg-minDeg))) %>%
  group_by(ScaleType) %>%
  summarize(FindSlope(PolValue,NormDegree))

# Apply function to whole dataset:
slope_data <- AdjData %>%
  mutate(NormDegree=((Degree-minDeg)/(maxDeg-minDeg))) %>%
  group_by(Results.index,ScaleType,ScaleType2,Scale,Adj,p0,p1,a,b) %>%
  summarize(FindSlope(PolValue,NormDegree)) %>%
  ungroup()
              


# Plot the slopes for each type of cc:
pdf(file="GraphSlopes.pdf", width=8,height=5)
slope_data %>%
  group_by(Scale) %>%
  summarize(x=mean(Slope1),y=mean(Slope2),se_x=se(Slope1),se_y=se(Slope2)) %>%
  mutate(xmin=(x-se_x),xmax=(x+se_x),
         ymin=(y-se_y),ymax=(y+se_y)) %>%
  # summarize(x=mean(ihs(Slope1)),y=mean(ihs(Slope2)),se_x=se(ihs(Slope1)),se_y=se(ihs(Slope2))) %>%
  # mutate(xmin=sinh(x-se_x),xmax=sinh(x+se_x),
  #        ymin=sinh(y-se_y),ymax=sinh(y+se_y),
  #        x=sinh(x),y=sinh(y)) %>%
  ggplot(aes(x=x,y=y,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,col=Scale))+
  geom_point(size=1) +
  geom_errorbar(aes(width=.5*(xmax-xmin)))+
  geom_errorbarh(aes(height=.5*(ymax-ymin)))+
  xlab("Bottom slope")+
  ylab("Top slope")+
  #coord_cartesian(xlim=c(0,20),ylim=c(0,20))+
  theme_bw()+
  scale_color_manual(values=FourColorPalette,name="Distribution")
dev.off()

# Graph by ScaleType (distinguishing max1 an max2)
pdf(file="GraphSlopesFull.pdf", width=8,height=5)
slope_data %>%
  #mutate(Slope1=ihs(Slope1),Slope2=ihs(Slope2)) %>%
  group_by(ScaleType) %>%
  summarize(x=mean(Slope1),y=mean(Slope2),se_x=se(Slope1),se_y=se(Slope2)) %>%
  mutate(xmin=(x-se_x),xmax=(x+se_x),
         ymin=(y-se_y),ymax=(y+se_y)) %>%
  mutate(ScaleType = factor(ScaleType,
                            levels=c("Relative1","Relative2","MinMax1","MinMax2","Min1","Min2","Max1","Max2"),
                            labels=c("Unbound 1","Unbound 2","Double bound 1","Double bound 2","Lower-bound 1","Lower-bound 2","Upper-bound 1","Upper-bound 2")
  )) %>%
  ggplot(aes(x=x,y=y,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,col=ScaleType,shape=ScaleType))+
  geom_point(size=3) +
  geom_errorbar(aes(width=.5*(xmax-xmin)))+
  geom_errorbarh(aes(height=.5*(ymax-ymin)))+
  xlab("Bottom slope")+
  ylab("Top slope")+
  #coord_cartesian(xlim=c(0,20),ylim=c(0,20))+
  theme_bw()+
  scale_color_manual(values=rep(FourColorPalette,each=2))+
  scale_shape_manual(values=rep(c(15,18),4))
dev.off()


# Slopes are not normally distributed, mainly because positive outliers
slope_data %>% pull(Slope1) %>% sort %>% custom_qqplot
slope_data %>% pull(Slope2) %>% sort %>% custom_qqplot

# Inverse hyperbolic sine (IHS) function
ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}

# IHS transformation does an ok job, though still not perfect.
slope_data %>% pull(Slope1)  %>% ihs%>% sort %>% custom_qqplot
slope_data %>% pull(Slope2) %>% ihs %>% sort %>% custom_qqplot

model_slope1 <- (lmer(ihs(Slope1)~Scale+(1|Results.index),data=slope_data))
model_slope2 <- (lmer(ihs(Slope2)~Scale+(1|Results.index),data=slope_data))
summary(model_slope1)
summary(model_slope1)[[10]] |>
  as_tibble() %>%
  mutate(p_value = 2*pnorm(abs(`t value`),lower.tail = F),Scale=levels(AdjData$Scale)) %>%
  select(Scale,everything())
summary(model_slope2)
summary(model_slope2)[[10]] |>
  as_tibble() %>%
  mutate(p_value = 2*pnorm(abs(`t value`),lower.tail = F),Scale=levels(AdjData$Scale)) %>%
  select(Scale,everything())

plot(residuals(model_slope1)~fitted(model_slope1));abline(0,0,col="red")
plot(residuals(model_slope2)~fitted(model_slope2));abline(0,0,col="red")



###################################################
# Effect of probability mass on slopes at boundary
###################################################


# Effect of probability mass at 0 and 1:
pdf(file="BotSlope-by-p0.pdf", width=6,height=5)
slope_data %>%
  filter(p0>0) %>%
  ggplot(aes(x=p0,y=ihs(Slope1),col=Scale,fill=Scale,group=Scale))+
  facet_grid(.~Scale,scale="free_x")+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_manual(values = FourColorPalette[2:3],aesthetics = c("color","fill"),guide="none")+
  scale_x_continuous(name="Probability mass at 0",labels = scales::label_percent(accuracy=1))+
  scale_y_continuous(name="Bottom slope (IHS-transformed)")
dev.off()


pdf(file="TopSlope-by-p1.pdf", width=6,height=5)
slope_data %>%
  filter(p1>0) %>%
  ggplot(aes(x=p1,y=ihs(Slope2),col=ScaleType,fill=ScaleType,group=ScaleType))+
  facet_grid(.~Scale,scale="free_x")+
  geom_point()+
  geom_smooth(method = "lm")+
  theme_bw()+
  scale_color_manual(values = FourColorPalette[c(2,4)],aesthetics = c("color","fill"),guide="none")+
  scale_x_continuous(name="Probability mass at 1",labels = scales::label_percent(accuracy=1))+
  scale_y_continuous(name="Top slope (IHS-transformed)")
dev.off()


model_slope1_p0 <- update(model_slope1, .~.+scale(p0))
anova(model_slope1_p0,model_slope1)
summary(model_slope1_p0)


model_slope2_p1 <- update(model_slope2, .~.+scale(p1))
anova(model_slope2_p1,model_slope2)
summary(model_slope2_p1)


## In response to some comments, check differences between adjectives: ##
slope_data  %>%
  mutate(Adj=factor(Adj,levels=0:7,labels=c("roagly","vibble","drok","scrop","plard","hif","tepable","plawic"))) %>%
  ggplot(aes(x=Adj,y=ihs(Slope1),color=Adj,group=paste(Scale,Adj))) +
  facet_wrap(.~Scale) + 
  geom_boxplot()+
  theme_bw()
slope_data  %>%
  mutate(Adj=factor(Adj,levels=0:7,labels=c("roagly","vibble","drok","scrop","plard","hif","tepable","plawic"))) %>%
  ggplot(aes(x=Adj,y=ihs(Slope2),color=Adj,group=paste(Scale,Adj))) +
  facet_wrap(.~Scale) + 
  geom_boxplot()+
  theme_bw()

adj_model_1 <- update(model_slope1, .~.+factor(Adj))
anova(adj_model_1,model_slope1)

adj_model_2 <- update(model_slope2, .~.+factor(Adj))
anova(adj_model_2,model_slope2)


## Regarding the ambiguity of double-bound adjectives
double_bound_lm <- slope_data %>%
  filter(ScaleType2=="MinMax") %>%
  lm(ihs(Slope2)~ihs(Slope1),data=.)

summary(double_bound_lm)

###############
# Fit beta-cdf 
###############

# Attempt to fit beta curves instead of simple sigmoids, which could possibly be an improvement on the RH-R model as well.

# Function to fit a beta cdf to participants responses
beta_fit_FUN <- function(degrees,response,scale=T){
  # We need to scale the responses, otherwise weird things happen on the part of scale that isn't used (which results in very low a/b values)
  if(scale){response = (response-min(response))/(max(response)-min(response))}
  target <- function(x){
    -sum(tobitLogLik(response,pbeta(degrees,exp(x[1]),exp(x[2])),exp(x[3])))
  }
  out <- exp(optim(par=c(a=1,b=1,sigma=.2),target)$par)
  return(tibble(a=out[1],b=out[2],sigma=out[3]))
}

cluster_copy(cluster,"beta_fit_FUN")

BetaFit <- AdjData %>%
  group_by(Results.index,ScaleType) %>%
  partition(cluster) %>%
  summarize(beta_fit_FUN(Degree,PolValue/100)) %>%
  collect()%>%
  ungroup()

summary(BetaFit)


BetaFit %>%
  mutate(a=log(a),b=log(b)) %>%
  group_by(ScaleType) %>%
  summarize(x=exp(mean(a)),y=exp(mean(b)),se_x=exp(se(a)),se_y=exp(se(b))) %>%
  mutate(xmin=(x/se_x),xmax=(x*se_x),
         ymin=(y/se_y),ymax=(y*se_y)) %>%
  mutate(ScaleType = factor(ScaleType,
                            levels=c("Relative1","Relative2","MinMax1","MinMax2","Min1","Min2","Max1","Max2"),
                            labels=c("Unbound 1","Unbound 2","Double bound 1","Double bound 2","Lower-bound 1","Lower-bound 2","Upper-bound 1","Upper-bound 2")
  )) %>%
  ggplot(aes(x=x,y=y,xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,col=ScaleType,shape=ScaleType))+
  geom_point(size=3) +
  geom_errorbar(aes(width=.5*(xmax-xmin)))+
  geom_errorbarh(aes(height=.5*(ymax-ymin)))+
  xlab("a")+
  ylab("b")+
  #coord_cartesian(xlim=c(0,20),ylim=c(0,20))+
  theme_bw()+
  scale_color_manual(values=rep(FourColorPalette,each=2))+
  scale_shape_manual(values=rep(c(15,18),4))




bet_plot_data <- BetaFit %>%
  rename(resp_a = a, resp_b = b) %>%
  left_join(
    AdjData %>% group_by(Results.index,ScaleType) %>% summarise(a=a[1],b=b[1])
  )
bet_plot_data %>%
  ggplot(aes(x=log(a),y=log(resp_a),col=ScaleType))+
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype=2)+
  scale_color_manual(values=rep(FourColorPalette,each=2))
  #coord_cartesian(ylim=c(0,50))

bet_plot_data %>%
  ggplot(aes(x=log(b),y=log(resp_b),col=ScaleType))+
  geom_point()+
  geom_abline(intercept=0,slope=1,linetype=2)+
  scale_color_manual(values=rep(FourColorPalette,each=2))


bet_plot_data %>%
  ggplot(aes(x=log(a),y=log(b),col=ScaleType,shape=ScaleType))+
  geom_point()+
  scale_color_manual(values=rep(FourColorPalette,each=2))+
  scale_shape_manual(values=rep(c(16,18),4))


BetaFit %>%
  group_by(ScaleType) %>%
  summarize(mean_log_a = mean(log(a)),mean_log_b=mean(log(b)),se_a=se(log(a)),se_b=se(log(b))) %>%
  ggplot(aes(x=mean_log_a,y=mean_log_b,
             xmin=mean_log_a-se_a,xmax=mean_log_a+se_a,
             ymin=mean_log_b-se_b,ymax=mean_log_b+se_b,
             col=ScaleType,shape=ScaleType))+
  geom_point(size=2) +
  geom_errorbar(width=.05)+
  geom_errorbarh(height=.05)+
  scale_color_manual(values=rep(FourColorPalette,each=2))+
  scale_shape_manual(values=rep(c(15,18),4))+
  theme_bw()


BetaFit %>%
  filter(ScaleType%in%c("MinMax1","MinMax2")) %>%
  ggplot(aes(x=log(a),y=log(b),col=ScaleType))+
  geom_point()


