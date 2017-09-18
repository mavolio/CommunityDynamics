library(tidyverse)
library(clibrary)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)


#work
#sim<-read.csv('C:\\Users\\megha\\Dropbox\\SESYNC\\SESYNC_RACs\\R Files/SimCom_Sept.csv')%>%
sim<-read.csv('SimCom_Sept.csv')%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)


#####CALCULATING DIVERSITY METRICS WITHIN A TIME STEP FOR EACH REPLICATE AND THEN AVERAGING LATER

#1) function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

# 2) function to calculate EQ evenness from Smith and Wilson 1996
#' @x the vector of abundances of each species
#' if all abundances are equal it returns a 1
E_q<-function(x){
  x1<-x[x!=0]
  if (length(x1)==1) {
    return(NA)
  }
  if (abs(max(x1) - min(x1)) < .Machine$double.eps^0.5) {##bad idea to test for zero, so this is basically doing the same thing testing for a very small number
    return(1)
  }
  r<-rank(x1, ties.method = "average")
  r_scale<-r/max(r)
  x_log<-log(x1)
  fit<-lm(r_scale~x_log)
  b<-fit$coefficients[[2]]
  2/pi*atan(b)
}

#3 I would like this to be functions to look at reordering
####New appraoch to Rank Shifts
###ranks - taking into account that all speices are not always present.
###Give all species with zerio abundace the S+1 rank for that year.
###includes species that are not present in year X but appear in year X+1 or are present in year X and disappear in year X+1


#NOTE this is where I tried to make a loop to rank for each dataset.

##add ranks 
sim_rank_pres<-sim%>%
  filter(abundance!=0)%>%
  tbl_df()%>%
  group_by(id, time, site)%>%
  mutate(rank=rank(-abundance, ties.method = "average"))%>%
  tbl_df()

#adding zeros
sim_addzero <- sim %>%
  group_by(id) %>%
  nest() %>%
  mutate(spread_df = purrr::map(data, ~spread(., key=species, value=abundance, fill=0) %>%
                                  gather(key=species, value=abundance, -site, -time, -id2))) %>%
  unnest(spread_df)

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
sim_zeros<-sim_addzero%>%
  filter(abundance==0)
##get species richness for each year
sim_S<-group_by(sim, id, time, site, id2)%>%
  summarize(S=S(abundance))
##merge together make zero abundances rank S+1
sim_zero_rank<-merge(sim_zeros, sim_S, by=c("id","time","site", "id2"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
sim_rank<-rbind(sim_rank_pres, sim_zero_rank)

##calculating re-ordering

reordering=data.frame(id=c(), time=c(), MRSc=c())#expeiment year is year of timestep2

spc_id<-unique(sim_rank$id2)
Rprof("loop1-profile.log")
initialtime <- proc.time()[3]
for (i in 1:length(spc_id)){
  subset<-sim_rank%>%
    filter(id2==spc_id[i])
  id2<-spc_id[i]
  #now get all timestep within an experiment
  timestep<-sort(unique(subset$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(time==timestep[i])
    
    subset_t2<-subset%>%
      filter(time==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id2"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    metrics<-data.frame(id2=id2, time=timestep[i+1], MRSc=MRSc)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

Rprof(NULL)
elapsedTime <- proc.time()[3] - initialtime
message(sprintf("Elapsed time: %s", elapsedTime))

