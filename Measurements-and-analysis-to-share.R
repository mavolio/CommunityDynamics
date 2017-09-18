library(tidyr)
library(dplyr)
library(codyn)
library(vegan)
library(Kendall)
library(ggplot2)
library(gridExtra)
library(reldist)
library(grid)
library(gtable)
library(purrr)

sim<-read.csv("~/Documents/SESYNC/SESYNC_RACs/R Files/SimCom_June.csv")%>%
  mutate(time=as.numeric(iteration),
         id2=paste(id, site, sep="::"))%>%
  select(-X, -sample, -iteration)
  
codyndat<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/relative cover_nceas and converge_12012015_cleaned.csv")%>%
  gather(species, abundance, sp1:sp99)%>%
  filter(site_code!="MISS")

codyndat_info<-read.csv("~/Dropbox/CoDyn/R Files/11_06_2015_v7/siteinfo_key.csv")%>%
  filter(site_project_comm!="")

###CLEANING CODYN DATASET
#restrict to species that are present in an experiment
splist<-codyndat%>%
  group_by(site_code, project_name, community_type, species)%>%
  summarize(present=sum(abundance))%>%
  filter(present!=0)%>%
  select(-present)

#merge back and will drop species that do not exist in a dataset
codyndat_clean<-merge(codyndat, splist, by=c("site_code","project_name","community_type","species"))%>%
  select(-X, -sitesubplot, -site_code, -project_name, -community_type)%>%
  mutate(id=paste(site_project_comm, plot_id, sep="::"))

#####CALCULATING DIVERSITY METRICS WITHIN A TIME STEP FOR EACH REPLICATE AND THEN AVERAGING LATER

#function to calculate richness
#' @x the vector of abundances of each species
S<-function(x){
  x1<-x[x!=0]
  length(x1)
}

#function to calculate EQ evenness from Smith and Wilson 1996
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


#####CALCULATING DIVERSITY METRICS ACROSS CONSECUTIVE TIME STEPS
#gains and losses
codyndat_loss<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="disappearance")
codyndat_gain<-turnover(df=codyndat_clean, time.var="experiment_year", species.var="species", abundance.var="abundance", replicate.var="id", metric="appearance")
codyndat_gains_loss<-merge(codyndat_gain, codyndat_loss, by=c("experiment_year","id"))%>%
  separate(id, c("site_project_comm", "plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarize(gain=mean(appearance),
            loss=mean(disappearance))

####New appraoch to Rank Shifts
###ranks - taking into account that all speices are not always present.
###Give all species with zerio abundace the S+1 rank for that year.
###includes species that are not present in year X but appear in year X+1 or are present in year X and disappear in year X+1

#####if the datset does not have zero abundances for species that are not present year to year do this first.
#adding zeros
sim_addzero <- sim %>%
  group_by(id) %>%
  nest() %>%
  mutate(spread_df = purrr::map(data, ~spread(., key=species, value=abundance, fill=0) %>%
                                  gather(key=species, value=abundance, -site, -time, -id2))) %>%
  unnest(spread_df)

###if dataset has zeros start here.

##add ranks dropping zeros
codyndat_rank_pres<-codyndat_clean%>%
  filter(abundance!=0)%>%
  tbl_df()%>%
  group_by(site_project_comm, experiment_year, plot_id)%>%
  mutate(rank=rank(-abundance, ties.method = "average"))%>%
  tbl_df()

###make zero abundant species have the rank S+1 (the size of the species pool plus 1)
##pull out zeros
codyndat_zeros<-codyndat_clean%>%
  filter(abundance==0)
##get species richness for each year
codyndat_S<-group_by(codyndat_clean, site_project_comm, experiment_year, plot_id)%>%
  summarize(S=S(abundance))
##merge together make zero abundances rank S+1
codyndat_zero_rank<-merge(codyndat_zeros, codyndat_S, by=c("site_project_comm","experiment_year","plot_id"))%>%
  mutate(rank=S+1)%>%
  select(-S)%>%
  tbl_df()
##combine all
codyndat_rank<-rbind(codyndat_rank_pres, codyndat_zero_rank)

##calculate reordering between time steps 3 ways, rank correlations, mean rank shifts not corrected, and mean ranks shifts corrected for the size of the speceis pool

reordering=data.frame(id=c(), experiment_year=c(), MRSc=c())#expeiment year is year of timestep2

spc_id<-unique(codyndat_rank$id)
  
for (i in 1:length(spc_id)){
  subset<-codyndat_rank%>%
    filter(id==spc_id[i])
  id<-spc_id[i]
  
  splist<-subset%>%
    select(species)%>%
    unique()
  sppool<-length(splist$species)
  
#now get all timestep within an experiment
  timestep<-sort(unique(subset$experiment_year))    

  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-subset%>%
      filter(experiment_year==timestep[i])
    
    subset_t2<-subset%>%
      filter(experiment_year==timestep[i+1])
    
    subset_t12<-merge(subset_t1, subset_t2, by=c("species","id"), all=T)%>%
      filter(abundance.x!=0|abundance.y!=0)
    
    MRSc<-mean(abs(subset_t12$rank.x-subset_t12$rank.y))/nrow(subset_t12)
    
    metrics<-data.frame(id=id, experiment_year=timestep[i+1], MRSc=MRSc)#spc_id
    ##calculate differences for these year comparison and rbind to what I want.
    
    reordering=rbind(metrics, reordering)  
  }
}

codyndat_reorder<-reordering%>%
  separate(id, c("site_project_comm","plot_id"), sep="::")%>%
  group_by(site_project_comm, experiment_year)%>%
  summarise(MRSc=mean(MRSc))



#####Calculating Bray-Curtis both comparing the mean community change between consequtive time steps and the change in dispersion between two time steps.
##Doing this for all years of an experiment at one time point because want to ensure all points are in the same space.
###first, get bray curtis dissimilarity values for each all years within each experiment between all combinations of plots
###second, get distance of each plot to its year centroid 
###third: mean_change is the distance the centroids of consequtive years
####fourth: dispersion_diff is the average dispersion of plots within a treatment to treatment centriod then compared between consequtive years

#Codyn dataset first

#make a new dataframe with just the label;
site_project_comm_u<-unique(codyndat_clean$site_project_comm)

#makes an empty dataframe
bray_curtis=data.frame(site_project_comm=c(), experiment_year=c(), bc_mean_change=c(), bc_dispersion_diff=c()) 

##calculating bray-curtis mean change and disperison differecnes
for(i in 1:length(site_project_comm_u)) {
  
  #subsets out each dataset
  subset=codyndat_clean%>%
    filter(site_project_comm==site_project_comm_u[i])%>%
    select(site_project_comm, experiment_year, species, abundance, plot_id)
  
  #get years
  experiment_years<-sort(unique(subset$experiment_year))
  
  #transpose data
  species=subset%>%
    spread(species, abundance, fill=0)
  
  #calculate bray-curtis dissimilarities
  bc=vegdist(species[,4:ncol(species)], method="bray")
  
  #calculate distances of each plot to year centroid (i.e., dispersion)
  disp=betadisper(bc, species$experiment_year, type="centroid")
  
  #getting distances between centroids over years; these centroids are in BC space, so that's why this uses euclidean distances
  cent_dist=as.matrix(vegdist(disp$centroids, method="euclidean"))

  ##extracting only the comparisions we want year x to year x=1.
  ###(experiment_year is year x+1
  cent_dist_yrs=data.frame(site_project_comm=site_project_comm_u[i],
                           experiment_year=experiment_years[2:length(experiment_years)],
                           mean_change=diag(cent_dist[2:nrow(cent_dist),1:(ncol(cent_dist)-1)]))
  
  #collecting and labeling distances to centroid from betadisper to get a measure of dispersion and then take the mean for a year
  disp2=data.frame(site_project_comm=site_project_comm_u[i],
                  experiment_year=species$experiment_year,
                  plot_id=species$plot_id,
                  dist=disp$distances)%>%
    tbl_df%>%
    group_by(site_project_comm, experiment_year)%>%
    summarize(dispersion=mean(dist))
  
  ##subtract consequtive years subtracts year x+1 - x. So if it is positive there was greater dispersion in year x+1 and if negative less dispersion in year x+1
  disp_yrs=data.frame(site_project_comm=site_project_comm_u[i],
                           experiment_year=experiment_years[2:length(experiment_years)],
                           dispersion_diff=diff(disp2$dispersion))
  
  #merge together change in mean and dispersion data
  distances<-merge(cent_dist_yrs, disp_yrs, by=c("site_project_comm","experiment_year"))
  
  #pasting dispersions into the dataframe made for this analysis
  bray_curtis=rbind(distances, bray_curtis)  
}

codyndat_braycurtis<-bray_curtis


####Looking at the shape of the curve - cc
###this compares the areas of difference between two curves that are consequtive time steps for a plot.

#codyn dat first
d_output=data.frame(site_project_comm=c(), experiment_year=c(), plot_id=c(), Darea=c())#expeiment year is year of timestep2

spc<-unique(codyndat_clean$site_project_comm)

for (i in 1:length(spc)){
  subset<-codyndat_clean%>%
    filter(site_project_comm==spc[i])
 
  ranks<-subset%>%
    filter(abundance!=0)%>%
    group_by(experiment_year, plot_id)%>%
    mutate(rank=rank(-abundance, ties.method="average"),
           maxrank=max(rank),
           relrank=rank/maxrank)%>%
    arrange(-abundance)%>%
    mutate(cumabund=cumsum(abundance))%>%
    ungroup()

  spc_id2<-spc[i]
  
  timestep<-sort(unique(ranks$experiment_year))    
  
  for(i in 1:(length(timestep)-1)) {#minus 1 will keep me in year bounds NOT WORKING
    subset_t1<-ranks%>%
      filter(experiment_year==timestep[i])
    
    plots_t1<-subset_t1%>%
      select(plot_id)%>%
      unique()
    
    subset_t2<-ranks%>%
      filter(experiment_year==timestep[i+1])
    
    plots_t2<-subset_t2%>%
      select(plot_id)%>%
      unique()
    
    plots_bothyrs<-merge(plots_t1, plots_t2, by="plot_id")
#dataset of two years    
    subset_t12<-rbind(subset_t1, subset_t2)
    
##dropping plots that were not measured both years
    subset_t12_2<-merge(plots_bothyrs, subset_t12, by="plot_id")
    
#dropping plots with only 1 species in any of the two years    
    drop<-subset_t12_2%>%
      group_by(experiment_year, plot_id)%>%
      mutate(numplots=length(plot_id))%>%
      ungroup()%>%
      group_by(plot_id)%>%
      mutate(min=min(numplots))%>%
      select(plot_id, min)%>%
      unique()
    
    subset_t12_3<-merge(subset_t12_2, drop, by="plot_id")%>%
      filter(min!=1)%>%
      ungroup()%>%
      group_by(experiment_year, plot_id)%>%
      arrange(rank)%>%
      ungroup()
    
    result <- subset_t12_3 %>%
    group_by(plot_id) %>%
    do({
      y <- unique(.$experiment_year)###assumption this is a length 2 list
      df1 <- filter(., experiment_year==y[[1]])
      df2 <- filter(., experiment_year==y[[2]])
      sf1 <- stepfun(df1$relrank, c(0, df1$cumabund))
      sf2 <- stepfun(df2$relrank, c(0, df2$cumabund))
      r <- sort(unique(c(0, df1$relrank, df2$relrank)))
      h <- abs(sf1(r) - sf2(r))
      w <- c(diff(r), 0)
      data.frame(Dstar=sum(w*h))#do has to output a dataframe
  })

    d_output1=data.frame(site_project_comm=spc_id2, experiment_year=timestep[i+1], plot_id=result$plot_id, Dstar=result$Dstar)#expeiment year is year of timestep2
    
    d_output<-rbind(d_output, d_output1)
  }
}
codyndat_dstar<-d_output%>% 
  group_by(site_project_comm, experiment_year)%>%
  summarise(Dstar=mean(Dstar))