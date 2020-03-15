# Analysis
# Dani Gargya, daniela@gargya.de
# March 2020


# calculating jaccard turnover ----
# Load packages ----
library(tidyverse) # (Contains loads of useful functions)
library(ggplot2) # (Package for making nice graphs)
library(vegan)

# load data ----
avatar_data <- read.csv("data/avatar_fauna.csv")


# bio <- read.csv("data/bio.csv")
bio_turnover <- bio %>% 
  select(SAMPLE_DESC, STUDY_ID_PLOT, STUDY_ID, PLOT, YEAR, GENUS_SPECIES, ID_SPECIES, sum.allrawdata.ABUNDANCE) %>% 
  group_by(STUDY_ID_PLOT) %>% 
  filter(YEAR %in% c(max(YEAR), min(YEAR))) %>% 
  group_by(STUDY_ID_PLOT, YEAR, GENUS_SPECIES) %>% 
  summarise(Abundance = sum(sum.allrawdata.ABUNDANCE))
  
bio_t_matrix <- bio_turnover %>% 
  group_by(STUDY_ID_PLOT) %>% 
  spread(GENUS_SPECIES, Abundance)







# from biotime website ----
# example with study ID 10
bio_10 <- bio %>%
  filter(STUDY_ID == 10) #%>% 
  #rename(Year = YEAR,
         #SampleID = SAMPLE_DESC,
         #Species = GENUS_SPECIES,
         #Abundance = sum.allrawdata.ABUNDANCE,
         #resamps = NUMBER_OF_SAMPLES)

# Maria Dornelas 09.06.2015
rarefysamples<-function(Year, SampleID, Species, Abundance, resamps) {
  #######################################################################
  # takes as input a  Year, SampleID, Species, Abundance and number of resamples
  # which should be in dataframe so that elements match
  # calculates turnover:
  # 1) between each year and the first year 
  # 2) between pairs of adjacent years 
  # 3) between each year and the lasy year of the time series
  # for the rarefied pooled samples
  ###########################################################################
  
  rareftab<-data.frame(array(NA,dim=c(0,3)))
  # getting vector with number of samples per year
  nsamples<-c()
  for(y in unique(Year)){
    nsamples<-c(nsamples, length(unique(SampleID[Year==y])))
  }
  t<-1
  minsample<-min(nsamples)
  for(repeats in 1:resamps){
    raref<-data.frame(array(NA,dim=c(1,3)))
    names(raref)<-c("Year","Species","Abundance")
    for(y in unique(Year)){
      #getting samples for this year
      samps<-unique(SampleID[Year==y])
      # re-sampling samples to equalize number of samples
      sam<-as.character(sample(samps,minsample,replace=T))
      # getting data that belongs to bootstraped samples
      rarefyear<-data.frame(SampleID[which(SampleID %in% sam & Year == y)],Species[which(SampleID %in% sam &Year
                                                                                         ==y)],Abundance[which(SampleID %in% sam & Year == y)])
      names(rarefyear)<-c("SampleID", "Species", "Abundance")
      # calculating pooled abundances of each species to store
      spabun<-tapply(as.numeric(rarefyear[,3]),as.character(rarefyear[,2]),sum)
      spar<-data.frame(rep(y, length(spabun)),names(spabun),spabun, row.names=NULL)
      names(spar)<-c("Year","Species","Abundance")
      raref<-rbind(raref,spar)
    }
    # calculating year by species table of abundance
    rareftab<-rbind(rareftab,cbind(rep(repeats,dim(raref)[1]),raref))
  }
  return (rareftab)
}

# run rarefaction code ----
TSrf<-list()
IDs<-unique(bio_10$STUDY_ID)

for(i in 1:length(IDs)){
  data<-bio_10[bio_10$STUDY_ID==IDs[i],]
  TSrf[[i]]<-rarefysamples(data$YEAR, data$SAMPLE_DESC, data$GENUS_SPECIES, data$sum.allrawdata.ABUNDANCE, 1) 
}
names(TSrf)<-IDs

rf<-do.call(rbind, TSrf)
rf<-data.frame(rf, ID=rep(names(TSrf), times=unlist(lapply(TSrf, nrow))))
rf<-rf[!is.na(rf$Year),-1]

#### prepare the rarefied output for diversity metric code
t1<-with(rf, tapply(Abundance, list(Year, Species, ID), sum))
t2<-unique(rf[, c("ID", "Year")])

#### produces a list of matrices for each study - in this case is only a single dataset
dsList<-list()

for(i in 1:dim(t1)[3]){
  id<-dimnames(t1)[[3]][i]
  a<-subset(t2, ID==id)$Year
  b<-t1[dimnames(t1)[[1]] %in% as.character(a),,i]
  dsList[[i]]<-data.frame(Year = rownames(b), b)
}

names(dsList) <- dimnames(t1)[[3]]

#### replacing NA with zero
for(i in 1:(length(dsList))){
  dsList[[i]][is.na(dsList[[i]])]<-0
}

# investigate biodiversity trends ----
#### Functions to calculate alpha and beta diversity metrics

####calculates species richness and total abundance
getAlpha<-function(x){ 
  data.frame(S=apply(x>0, 1, sum),
             N=apply(x, 1, sum))
}

####calculates Jaccard and Morisita-Horn similarity indices
getBeta<-function(x, origin=1, consec=F){
  size<-dim(x)[1]-1
  ifelse(consec==T, id<-cbind(1:size, 1:size+1), id<-cbind(1:dim(x)[1], origin))
  output<-data.frame(Jaccard=as.matrix(1-vegdist(x, method="jaccard"))[id],
                     MorHorn=as.matrix(1-vegdist(x, method="horn"))[id])
  ifelse(consec==T, rownames(output)<-rownames(x)[-1], rownames(output)<-rownames(x))
  output
}


#### get the alpha diversity metrics
alphaX<-lapply(dsList, function(x)getAlpha(x[,-1]))
alpha<-do.call(rbind,alphaX)
alpha$Year<-substring(rownames(alpha), nchar(rownames(alpha))-3)

#### get the beta diversity metrics
betaX<-lapply(dsList, function(x)getBeta(x[,-1]))
beta<-do.call(rbind, lapply(betaX, function(x)x[-1,]))
beta$Year<-substring(rownames(beta), nchar(rownames(beta))-3)

# function clean theme
themeNoGridAxes <- function() {
  theme_bw()+
    theme(axis.text=element_text(size=16,color="black"),
          axis.title=element_text(size=16,face="bold"),
          legend.position="none",     
          plot.title=element_text(size=18,face="bold", hjust=0.5), 
          plot.background=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(color="black")
    )
}

#### plot species richness temporal trend with simple regression line
fitSR<-lm(alpha$S~alpha$Year)

srPlot<-ggplot(alpha, aes(x=Year, y=S))+geom_line(size=1, colour="gray30")+
  geom_point(aes(x=Year, y=S), colour="darkgreen", size=5)+
  scale_x_continuous(breaks=c(1993, 1995, 1997, 1999, 2001, 2003, 2005))+
  geom_abline(intercept=fitSR$coef[1], slope=fitSR$coef[2], colour="gray10", size=1, linetype=3)+
  themeNoGridAxes()+ylab("S")+xlab("Year")
srPlot


#### plot the Jaccard similarity temporal trend with a simple regression line
fitJ<-lm(beta$Jaccard~beta$Year)

jPlot<-ggplot(beta, aes(x=Year, y=Jaccard))+geom_line(size=1, colour="gray30")+
  geom_point(aes(x=Year, y=Jaccard), colour="dodgerblue", size=5)+
  #scale_x_continuous(breaks=c(1992, 1995, 1996, 1999, 2001, 2003, 2005))+
  geom_abline(intercept=fitJ$coef[1], slope=fitJ$coef[2], colour="gray10", size=1, linetype=3)+
  themeNoGridAxes()+ylab("Jaccard Similarity")+xlab("Year")+ylim(0, 1)
jPlot
