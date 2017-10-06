library(tidyverse)
library(codyn)
library(vegan)
library(Kendall)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)


#work
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
#' @param sim foo dataset with columns for time, plot, species, and abundance, (and optionally an id column for grouping)
add_ranks_for_non_present_and_present_species <- function(sim) {
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
}

##calculating re-ordering
calculate_reordering <- function(sim_rank) {
  reordering=data.frame(id=c(), time=c(), MRSc=c())#expeiment year is year of timestep2
  
  spc_id<-unique(sim_rank$id2)
  
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
  return(reordering)
}

# TODO define functions throughout whole file
# TODO rename data frame argument to df
# TODO add arguments to name the time, site, species, and abundance columns
# TODO refactor to allow users to use arbitrary column names (remove dplyr?)
# TODO refactor to optimize speed of function
# TODO add new function to codyn
# TODO document the functions fully
# TODO add tests to codyn for the new functions
mean_rank_shift_corrected <- function(df) {
  sim_rank <- add_ranks_for_non_present_and_present_species(df)
  result <- calculate_reordering(sim_rank)
  return(result)
}

ranked <- mean_rank_shift_corrected(sim)

#sim_reorder<-reordering%>%
#  separate(id2, c("id","site"), sep="::")%>%
#  group_by(id, time)%>%
#  summarise(MRSc=mean(MRSc))

#4/5)Bray-Curtis Mean Change & Bray-Curtis Dissimilarity 

#####Calculating Bray-Curtis both comparing the mean community change between consequtive time steps and the change in dispersion between two time steps.
##Doing this for all years of an experiment at one time point because want to ensure all points are in the same space.
###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years


#make a new dataframe with just the label;
id_u<-unique(sim$id)

#makes an empty dataframe
bray_curtis=data.frame(id=c(), time=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 

#Calculating bc mean change and dispersion
for(i in 1:length(id_u)) {
  
  #subsets out each dataset
  subset=sim%>%
    filter(id==id_u[i])%>%
    select(id, time, species, abundance, site)
  
  #get years
  timestep<-sort(unique(subset$time))
  
  #transpose data
  species=subset%>%
    spread(species, abundance, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,4:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$time, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))
  
  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(id=id_u[i],
                           time=timestep[2:length(timestep)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(id=id_u[i],
                   time=species$time,
                   site=species$site,
                   dist=disp$distances)%>%
    tbl_df%>%
    group_by(id, time)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(id=id_u[i],
                      time=timestep[2:length(timestep)],
                      dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("id","time"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}
sim_bray_curtis<-bray_curtis

#6) curve comparision
####Looking at the shape of the curve - cc
###this compares the areas of difference between two curves that are consequtive time steps for a plot.

d_output=data.frame(id=c(), time=c(), site=c(), Darea=c())#expeiment year is year of timestep2

com<-unique(sim$id)

for (i in 1:length(com)){
  subset<-sim%>%
    filter(id==com[i])
  
  ranks<-subset%>%
    filter(abundance!=0)%>%
    group_by(time, site)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()
  
  com_id2<-com[i]
  
  timestep<-sort(unique(ranks$time))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-ranks%>%
      filter(time==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(site)%>%
      unique()
    
    subset_t2<-ranks%>%
      filter(time==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(site)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="site")
    #dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
    ##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="site")
    
    #dropping plots with only 1 species in any of the two years    
    drop<-subset_t12_2%>%
      group_by(time, site)%>%
      mutate(numplots=length(site))%>%
      ungroup()%>%
      group_by(site)%>%
      mutate(min=min(numplots))%>%
      select(site, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="site")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(time, site)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
      group_by(site) %>%
      do({
        y <- unique(.$time)###assumption this is a length 2 list
        df1 <- filter(., time==y[[1]])
        df2 <- filter(., time==y[[2]])
        sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
        sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
        r <- sort(unique(c(0, df1$relrank, df2$relrank)))
        h <- abs(sf1(r) - sf2(r))
        w <- c(diff(r), 0)
        data.frame(Dstar=sum(w*h))#do has to output a dataframe
      })
    
    d_output1=data.frame(id=com_id2, time=timestep[i+1], site=result$site, Dstar=result$Dstar)#expeiment year is year of timestep2
    
    d_output<-rbind(d_output, d_output1)
  }
}

sim_dstar<-d_output%>% 
  group_by(id, time)%>%
  summarise(Dstar=mean(Dstar))

########comparing control versus treatment plots

#
###need a contorl and treatment column

data<-read.csv("subsetexperiments.csv")


data2<-data%>%
  mutate(trt=ifelse(site_code=="KNZ"&treatment=="N1P0", "C", ifelse(site_code=="NIN"&treatment=="1NF","C", ifelse(site_code=="SEV"&treatment=="C","C","T"))))

reordering_ct=data.frame(site_project_comm=c(), treatment=c(), calendar_year=c(), MRSc_diff=c(), spdiffc=c())

explist<-unique(data2$id)

for (i in 1:length(explist)){
  ##get zero abundances to be filled in for all species.
  ##this works the first time only
  subset<-data2%>%
    filter(id==explist[i])%>%
    spread(genus_species, relcov, fill=0)
  
  spc<-explist[i]
  
  ##make wide and get averages of each species by treatment
  wide<-subset%>%
    gather(genus_species, relcov, 12:ncol(subset))%>%
    group_by(site_project_comm, calendar_year, treatment, trt,genus_species)%>%
    summarize(relcov=mean(relcov))
  
  ##add ranks dropping zeros
  rank_pres<-wide%>%
    filter(relcov!=0)%>%
    tbl_df()%>%
    group_by(site_project_comm, calendar_year, treatment, trt)%>%
    mutate(rank=rank(-relcov, ties.method = "average"))%>%
    tbl_df()
  
  ###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
  ##pull out zeros
  zeros<-wide%>%
    filter(relcov==0)
  ##get species richness for each year
  rich<-group_by(wide, site_project_comm, calendar_year, treatment)%>%
    summarize(S=S(relcov))
  ##merge together make zero abundances rank S+1
  zero_rank<-merge(zeros, rich, by=c("site_project_comm","calendar_year", "treatment"))%>%
    mutate(rank=S+1)%>%
    select(-S)%>%
    tbl_df()
  ##combine all
  rank<-rbind(rank_pres, zero_rank)
  
  timestep<-sort(unique(rank$calendar_year)) 
  for(i in 1:(length(timestep))){
    
    time<-rank%>%
      filter(calendar_year==timestep[i])
    
    time_id<-timestep[i]
    
    #fitler out control plots
    control<-time%>%
      filter(trt=="C")
    
    treat_list<-unique(subset(time, trt!="C")$treatment)
    
    for (i in 1:length(treat_list)){
      treat<-time%>%
        filter(treatment==treat_list[i])
      
      treat_id<-treat_list[i]
      
      subset_ct<-merge(control, treat, by=c("site_project_comm", "calendar_year","genus_species"), all=T)%>%
        filter(relcov.x!=0|relcov.y!=0)
      
      MRSc_diff<-mean(abs(subset_ct$rank.x-subset_ct$rank.y))/nrow(subset_ct)
      
      spdiff<-subset_ct%>%
        filter(relcov.x==0|relcov.y==0)
      
      spdiffc<-nrow(spdiff)/nrow(subset_ct)
      
      metrics<-data.frame(site_project_comm=spc, treatment=treat_id, calendar_year=time_id, MRSc_diff=MRSc_diff, spdiffc=spdiffc)#spc_id
      ##calculate differences for these year comparison and rbind to what I want.
      
      reordering_ct=rbind(metrics, reordering_ct)  
    }
  }
}
