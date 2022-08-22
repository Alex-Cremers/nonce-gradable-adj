##################################
# Load, Clean and Process the data
##################################


# Load raw anonymized results:
# This scripts loads the data, and that's about it. It's meant to be called by the other scripts, not run on its own.

subjectInfo <- read_csv("GradAdj-SubjectInfo.csv",show_col_types = FALSE)
data <- read_csv("GradAdj-AnonData.csv",show_col_types = FALSE)


# define Polarity-controlled response and Flag cases where participants missed the negation
tmp <- data %>%
  mutate(PolValue=if_else(Polarity=="Neg",100L-value,value)) %>%
  group_by(Results.index,ScaleType,Polarity) %>%
  summarize(Coef=lm(PolValue~Degree)$coefficients[[2]])
Threshold = mean(tmp$Coef) - sd(tmp$Coef)
print(paste0("Threshold for polarity errors: ",round(Threshold/100,2)))

data <- data %>%
  left_join(tmp) %>% 
  mutate(PolarityError = Coef < Threshold) %>%
  mutate(PolValue=if_else(Polarity=="Neg",100L-value,value)) %>%
  ungroup()

rm(tmp)

subjectInfo <- left_join(subjectInfo,data %>% group_by(Results.index) %>%
                           summarize(PolErrors=sum(PolarityError)/10))
# Excluded participants: (too fast, too many duplicates, or missed polarity on more than half of the blocks)
excluded <- subjectInfo %>%
  filter(duplicates>.5|medianRT<1000|PolErrors>=8) %>%
  pull(Results.index)

print(paste0("Participants excluded from analyses: ",length(excluded)))

print("Polarity errors by polarity:")
data %>%
  group_by(Polarity) %>%
  summarize(`Error rate` = mean(PolarityError))

# Filter out excluded participants and clean up data:
# We skip the response time filter because we want to keep either all or none of the trials in a block
# (otherwise this messes up with the degree distribution when computing the Klein model)
N <- nrow(data)
data <- data %>%
  filter(!Results.index %in% excluded) %>%
  filter(!PolarityError) %>%
  select(-PolarityError)
print(paste0("Proportion of data thrown away: ",round(1-nrow(data)/N,3)))


# We now focus on the bare adjectives:
AdjData <- data %>%
  filter(Adverb=="Adj") %>%
  select(-Adverb)


#### Prepare data for Stan ####

StanData <- AdjData %>%
  mutate(y=PolValue/100,
         Adj=as.numeric(factor(Adj)),
         Results.index=as.numeric(factor(Results.index)),
         NormedDegree = (Degree-minDeg)/(maxDeg-minDeg))%>%
  mutate(
    ScaleType2 = substr(ScaleType,1,nchar(ScaleType)-1),
    Scale=factor(ScaleType2,
                 levels=c("Relative","MinMax","Min","Max"),
                 labels=c("Unbound","Double bound","Lower-bound","Upper-bound"))
  ) %>%
  group_by(Results.index,ScaleType) %>%
  mutate(NormedRank=as.numeric(factor(rank(Degree,ties.method="min")))) %>%
  ungroup()


StanAdverbData <- data %>%
  filter(Adverb!="aBitAdjAtAll") %>%
  mutate(y=PolValue/100,
         Adj=as.numeric(factor(Adj)),
         Adverb=as.numeric(factor(Adverb,levels=c("Adj","quiteAdj","veryAdj","extremelyAdj","absolutelyAdj"))),
         Results.index=as.numeric(factor(Results.index)),
         NormedDegree = (Degree-minDeg)/(maxDeg-minDeg))%>%
  mutate(
    ScaleType2 = substr(ScaleType,1,nchar(ScaleType)-1),
    Scale=factor(ScaleType2,
                 levels=c("Relative","MinMax","Min","Max"),
                 labels=c("Unbound","Double bound","Lower-bound","Upper-bound"))
  )%>%
  group_by(Results.index,ScaleType) %>%
  mutate(NormedRank=as.numeric(factor(rank(Degree,ties.method="min")))) %>%
  ungroup()



