# This script contains descriptive plots of the data
# For plots of computed slopes at scale ends, see the main script.

library(tidyverse)
options(tibble.width=Inf)

# Color palettes for graphs:
TwoColorPalette = c(rgb(.5,0,.7),rgb(.9,.65,0))
FourColorPalette = c("#F5793A","#A95AA1","#85C0F9","#0F2080")

# Standard error function for graphs
se <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))}

# Load and clean the data:
source("GradAdj-load-data.R")


############################################################
# Comparing plots by raw degree, normed degree, normed rank
#                       Binned plots
############################################################

### BY BARE DEGREE ###

# Bare adjective data, by binned degree and Scale
AdjData %>%
  mutate(BinnedDegree = cut(Degree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(0,seq(0.05,1,.1),1))),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  group_by(Scale,BinnedDegree) %>%
  summarize(Value=mean(PolValue),Error=se(PolValue)) %>%
  ungroup() %>%
  mutate(BinnedDegree=as.numeric(as.character(BinnedDegree))) %>%
  ggplot(aes(x=BinnedDegree,y=Value,ymin=Value-Error,ymax=Value+Error,color=Scale,fill=Scale))+
  facet_wrap(.~Scale)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_line()+
  geom_ribbon(alpha=.2)+
  theme_bw() +
  scale_color_manual(values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))


### BY NORMED DEGREE ###

# Bare adjective data, by binned normed degree and Scale
AdjData %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         #BinnedNormedDegree = cut(NormedDegree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(seq(0,1,.1),.99))),
         BinnedNormedDegree = cut(NormedDegree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(0,seq(0.05,1,.1),1))),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
         ) %>%
  # Pooling together extremities where there isn't enough data:
  mutate(BinnedNormedDegree = case_when(
    ScaleType2=="Min"&BinnedNormedDegree%in%c(.95,1) ~ 1,
    ScaleType2=="Max"&BinnedNormedDegree%in%c(.05,.15) ~ .1,
    T ~ as.numeric(as.character(BinnedNormedDegree))
  )) %>%
  group_by(Scale,BinnedNormedDegree) %>%
  summarize(Value=mean(PolValue),Error=se(PolValue)) %>%
  ungroup() %>%
  mutate(BinnedNormedDegree=as.numeric(as.character(BinnedNormedDegree))) %>%
  ggplot(aes(x=BinnedNormedDegree,y=Value,ymin=Value-Error,ymax=Value+Error,color=Scale,fill=Scale))+
  facet_wrap(.~Scale)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black",inherit.aes=F)+
  geom_line()+
  geom_ribbon(alpha=.2)+
  theme_bw() +
  scale_color_manual(values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Normalized degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))


### BY RANK ###

# Bare adjective data, by binned normalized rank and Scale
AdjData %>%
  mutate(Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  group_by(Results.index,ScaleType) %>%
  mutate(ComputedRank = as.numeric(factor(Rank)), # Get ranks of the degrees, ignoring repeats
         ScaledRank = (ComputedRank-min(ComputedRank))/(max(ComputedRank)-min(ComputedRank)),
         BinnedScaledRank = cut(ScaledRank,breaks=sort(c(seq(-0.1,1,.1))),labels=sort(c(0,seq(0.05,1,.1)))),
  ) %>%
  ungroup() %>%
  group_by(Scale,BinnedScaledRank) %>%
  summarize(Value=mean(PolValue),Error=se(PolValue)) %>%
  ungroup() %>%
  mutate(BinnedScaledRank=as.numeric(as.character(BinnedScaledRank))) %>%
  ggplot(aes(x=BinnedScaledRank,y=Value,ymin=Value-Error,ymax=Value+Error,color=Scale,fill=Scale))+
  facet_wrap(.~Scale)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_line()+
  geom_ribbon(alpha=.2)+
  theme_bw() +
  scale_color_manual(values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Normalized Rank",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))


############################################################
# Comparing plots by raw degree, normed degree, normed rank
#              Raw data with smooth geom
############################################################

### BY RAW DEGREE ###

# Bare adjective data, by raw degree and Scale with smoothed GAM
AdjData %>%
  mutate(Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  ggplot(aes(x=Degree,y=PolValue,color=Scale,fill=Scale))+
  facet_wrap(.~Scale)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam')+
  theme_bw() +
  scale_color_manual(name="Distribution",values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))

### BY NORMED DEGREE ###

# Bare adjective data, by normed degree and Scale with smoothed GAM
AdjData %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue,color=Scale,fill=Scale))+
  facet_wrap(.~Scale)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam')+
  theme_bw() +
  scale_color_manual(name="Distribution",values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Normalized Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))


### BY RANK ###

# Bare adjective data, by normalized rank and Scale with smoothed GAM
AdjData %>%
  mutate(Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  group_by(Results.index,ScaleType) %>%
  mutate(ComputedRank = as.numeric(factor(Rank)), # Get ranks of the degrees, ignoring repeats
         ScaledRank = (ComputedRank-min(ComputedRank))/(max(ComputedRank)-min(ComputedRank))
  ) %>%
  ungroup() %>%
  ggplot(aes(x=ScaledRank,y=PolValue,color=Scale,fill=Scale))+
  facet_wrap(.~Scale)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam')+
  theme_bw() +
  scale_color_manual(name="Distribution",values=FourColorPalette,aesthetics = c("colour","fill"))+
  scale_x_continuous(name="Normalized Rank",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))



##################
# Effect of p0/p1
##################

# Not enough data to see a clear pattern for each value, so we bin.

# Min adjectives:
AdjData %>%
  filter(ScaleType2=="Min") %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  mutate(BinnedP0=cut(p0,
                      breaks=c(0,0.2,0.35,0.50,0.8),
                      labels=c("[0.1 , 0.2]","[0.25 , 0.35]","[0.4 , 0.5]","[0.55 , 0.7]"))
         ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue, group=BinnedP0))+
  facet_wrap(.~BinnedP0)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam',col=FourColorPalette[3],fill=FourColorPalette[3])+
  theme_bw() +
  scale_x_continuous(name="Normalized Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))


# Max adjectives:
AdjData %>%
  filter(ScaleType2=="Max") %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  mutate(BinnedP1=cut(p1,
                      breaks=c(0,0.2,0.35,0.50,0.8),
                      labels=c("[0.1 , 0.2]","[0.25 , 0.35]","[0.4 , 0.5]","[0.55 , 0.7]"))
  ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue, group=BinnedP1))+
  facet_wrap(.~BinnedP1)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam',col=FourColorPalette[4],fill=FourColorPalette[4])+
  theme_bw() +
  scale_x_continuous(name="Normalized Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))


# MinMax adjectives:
AdjData %>%
  filter(ScaleType2=="MinMax") %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  mutate(BinnedP0=cut(p0,
                      breaks=c(0,0.2,0.35,0.50,0.8),
                      labels=c("[0.1 , 0.2]","[0.25 , 0.35]","[0.4 , 0.5]","[0.55 , 0.7]")),
         BinnedP1=cut(p1,
                      breaks=c(0,0.2,0.35,0.50,0.8),
                      labels=c("[0.1 , 0.2]","[0.25 , 0.35]","[0.4 , 0.5]","[0.55 , 0.7]")),
         BinnedP01=cut(round(p0+p1,2),
                       breaks=c(c(0,.35,.45,.65))),
         p0=BinnedP0,p1=BinnedP1
  ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue, group=paste(BinnedP0,BinnedP1)))+
  facet_grid(p1~p0,labeller=label_both)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam',col=FourColorPalette[2],fill=FourColorPalette[2])+
  theme_bw() +
  scale_x_continuous(name="Normalized Degree",labels=scales::percent)+
  scale_y_continuous(name="Acceptability",labels=scales::label_percent(scale=1))




##################################
# Difference between type 1 and 2
##################################

# Min adjectives:
AdjData %>%
  filter(ScaleType2=="Min") %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  mutate(BinnedP0=cut(p0,
                      breaks=c(0,0.2,0.35,0.50,0.8),
                      labels=c("[0.1 , 0.2]","[0.25 , 0.35]","[0.4 , 0.5]","[0.55 , 0.7]"))
  ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue, group=BinnedP0))+
  facet_grid(ScaleType~BinnedP0)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam')+
  theme_bw() +
  xlab("Normalized degree")+ylab("Acceptability")

AdjData %>%
  filter(ScaleType2=="Min") %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  mutate(BinnedP0=cut(p0,
                      breaks=c(0,0.35,0.8),
                      labels=c("Low p0","High p0"))
  ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue, group=BinnedP0))+
  facet_grid(ScaleType~BinnedP0)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam')+
  theme_bw() +
  xlab("Normalized degree")+ylab("Acceptability")

AdjData %>%
  filter(ScaleType2=="Max") %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative'))
  ) %>%
  mutate(BinnedP1=cut(p1,
                      breaks=c(0,0.35,0.8),
                      labels=c("Low p1","High p1"))
  ) %>%
  ggplot(aes(x=NormedDegree,y=PolValue, group=BinnedP1))+
  facet_grid(ScaleType~BinnedP1)+
  geom_abline(slope=100,intercept=0,linetype=2,color="black")+
  geom_point(alpha=.1,pch=16)+
  geom_smooth(method = 'gam',color="red")+
  theme_bw() +
  xlab("Normalized degree")+ylab("Acceptability")



##########################
# Graphs for intensifiers
##########################


# Graph by Scale and adverb (a bit too noisy)
# Pre-process errors:
tmp <- data %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         BinnedNormedDegree = cut(NormedDegree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(0,seq(0.05,1,.1),1))),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Adverb = as.character(Adverb),
         Adverb = case_when(
           Adverb=="aBitAdjAtAll" & Polarity == "Aff" ~ "aBitAdj",
           Adverb=="aBitAdjAtAll" & Polarity == "Neg" ~ "atAllAdj",
           T ~ Adverb
         ),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative')),
         Adverb=factor(Adverb,levels=c("Adj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj","atAllAdj"))
  ) %>%
  mutate(BinnedNormedDegree = case_when(
    ScaleType2=="Min"&BinnedNormedDegree%in%c(.95,1) ~ 1,
    ScaleType2=="Max"&BinnedNormedDegree%in%c(.05,.15) ~ .1,
    T ~ as.numeric(as.character(BinnedNormedDegree))
  )) %>%
  group_by(Adverb,BinnedNormedDegree,Scale) %>%
  summarize(Error=se(PolValue)) %>%
  pivot_wider(id_cols = c(BinnedNormedDegree,Scale),names_from = Adverb,values_from=Error) %>%
  rename(AdjErr=Adj) %>%
  pivot_longer(cols=ends_with("Adj"),names_to = "Adverb",values_to = "Error") %>%
  mutate(Adverb=factor(Adverb,levels=c("atAllAdj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj")))

data %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         BinnedNormedDegree = cut(NormedDegree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(0,seq(0.05,1,.1),1))),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Adverb = as.character(Adverb),
         Adverb = case_when(
           Adverb=="aBitAdjAtAll" & Polarity == "Aff" ~ "aBitAdj",
           Adverb=="aBitAdjAtAll" & Polarity == "Neg" ~ "atAllAdj",
           T ~ Adverb
         ),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative')),
         Adverb=factor(Adverb,levels=c("Adj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj","atAllAdj"))
  ) %>%
  mutate(BinnedNormedDegree = case_when(
    ScaleType2=="Min"&BinnedNormedDegree%in%c(.95,1) ~ 1,
    ScaleType2=="Max"&BinnedNormedDegree%in%c(.05,.15) ~ .1,
    T ~ as.numeric(as.character(BinnedNormedDegree))
  )) %>%
  group_by(Adverb,BinnedNormedDegree,Scale) %>%
  summarize(Value=mean(PolValue)) %>% #,Error=se(PolValue)
  pivot_wider(id_cols = c(BinnedNormedDegree,Scale),names_from = Adverb,values_from=Value) %>%
  rename(AdjRef=Adj) %>%
  pivot_longer(cols=ends_with("Adj"),names_to = "Adverb",values_to = "Value") %>%
  mutate(Adverb=factor(Adverb,levels=c("atAllAdj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj"))) %>%
  left_join(tmp) %>%
  mutate(Adverb=factor(Adverb,
                       levels=c("atAllAdj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj"),
                       labels=c("at all","a bit","quite","very","extremely","absolutely")),
  ) %>%
  ggplot(aes(x=AdjRef,y=Value,ymin=Value-Error,ymax=Value+Error,group=Scale,color=Scale,fill=Scale)) +
  facet_wrap(.~Adverb)+
  geom_line()+
  geom_ribbon(alpha=.2)+
  geom_ribbon(aes(ymin=AdjRef-AdjErr,ymax=AdjRef+AdjErr),alpha=.2,color="transparent",fill="black")+
  geom_line(aes(y=AdjRef),color="black",linetype=2,alpha=.5)+
  theme_bw() +
  scale_color_manual(values=FourColorPalette,aesthetics = c("fill","color"))



# Graph by Adverb only
tmp <- data %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         BinnedNormedDegree = cut(NormedDegree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(0,seq(0.05,1,.1),1))),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Adverb = as.character(Adverb),
         Adverb = case_when(
           Adverb=="aBitAdjAtAll" & Polarity == "Aff" ~ "aBitAdj",
           Adverb=="aBitAdjAtAll" & Polarity == "Neg" ~ "atAllAdj",
           T ~ Adverb
         ),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative')),
         Adverb=factor(Adverb,levels=c("Adj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj","atAllAdj"))
  ) %>%
  group_by(Adverb,BinnedNormedDegree) %>%
  summarize(Error=se(PolValue)) %>%
  pivot_wider(id_cols = c(BinnedNormedDegree),names_from = Adverb,values_from=Error) %>%
  rename(AdjErr=Adj) %>%
  pivot_longer(cols=ends_with("Adj"),names_to = "Adverb",values_to = "Error")

data %>%
  mutate(NormedDegree = (Degree-minDeg)/(maxDeg-minDeg),
         BinnedNormedDegree = cut(NormedDegree,breaks=sort(c(seq(-0.1,1,.1),.999)),labels=sort(c(0,seq(0.05,1,.1),1))),
         Scale=factor(ScaleType2,
                      levels=c("Relative","MinMax","Min","Max"),
                      labels=c("Unbound","Double bound","Lower-bound","Upper-bound")),
         Adverb = as.character(Adverb),
         Adverb = case_when(
           Adverb=="aBitAdjAtAll" & Polarity == "Aff" ~ "aBitAdj",
           Adverb=="aBitAdjAtAll" & Polarity == "Neg" ~ "atAllAdj",
           T ~ Adverb
         ),
         Polarity=factor(Polarity,levels=c("Aff","Neg"),labels=c("Affirmative",'Negative')),
         Adverb=factor(Adverb,levels=c("Adj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj","atAllAdj"))
  ) %>%
  group_by(Adverb,BinnedNormedDegree) %>%
  summarize(Value=mean(PolValue)) %>% #,Error=se(PolValue)
  pivot_wider(id_cols = c(BinnedNormedDegree),names_from = Adverb,values_from=Value) %>%
  rename(AdjRef=Adj) %>%
  pivot_longer(cols=ends_with("Adj"),names_to = "Adverb",values_to = "Value") %>%
  left_join(tmp) %>%
  mutate(Adverb=factor(Adverb,
                       levels=c("atAllAdj","aBitAdj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj"),
                       labels=c("at all","a bit","quite","very","extremely","absolutely")),
         ) %>%
  ggplot(aes(x=AdjRef,y=Value,ymin=Value-Error,ymax=Value+Error)) +
  facet_wrap(.~Adverb)+
  geom_line(color=TwoColorPalette[1])+
  geom_ribbon(alpha=.2,color=TwoColorPalette[1],fill=TwoColorPalette[1])+
  geom_ribbon(aes(ymin=AdjRef-AdjErr,ymax=AdjRef+AdjErr),alpha=.2,color=TwoColorPalette[2],fill=TwoColorPalette[2])+
  geom_line(aes(y=AdjRef),color=TwoColorPalette[2],linetype=2)+
  #xlim(0,100)+ylim(0,100)+
  xlab("Rating for bare adjective")+
  ylab("Rating for modified adjective")+
  theme_bw()

