#####################################################################################
## Example code and accompanying dataset for estimating risk to migratory birds
## using a three-component framework: multispecies migratory connectivity as a 
## measure of EXPOSURE, conservation assessment scores as a measure of VULNERABILTIY,
## and temperature and land-cover change as a measure of HAZARD. The product of these
## three components represents risk of population declines for the migratory bird
## species considered in the analysis. 
#
## As an example, the 5 species included here are blackpoll warbler, grasshopper sparrow, 
## ovenbird, prothonotary warbler and tree swallow. For the paper, a total of 112 species 
## were included but data sharing agreements with individual data holders preclude our 
## ability to share a complete dataset publicly. Thus, results produced here will differ 
## from the full analysis. Regardless, this code can be used as a template for calculating 
## risk and its three components for any set of species for which data are available. 
## Finally, note that this script is the second in the workflow; first, species-specific 
## connectivity is estimated following the MC_ex script using available movement data for 
## a given species.
#
### DIRECTIONS to run script: set desired working directory in LINE 73 and run ###
# note: may need to install necessary libraries if not already installed
#
# Inputs:
# 1. Folder of species-specific connectivity results (Connectivity_ExSpGroup), which are produced
# from the first script in this workflow (MC_ex.R). Note that 5 species results are provided for
# this example, but due to data restrictions, only raw data for bkpwar are provided as an example 
# to run the first script. Thus, only the bkpwar file (_MigConnect_Results.csv) can be reproduced
# from the MC_ex code. The other 4 species results files are shared here with permission.
# For more details on data sharing, please see the Data Availability statement. 
# Note that outputs from the full analysis of 112 species are
# also available on GitHub and Zenodo (see Data Availability statement for links).
#
# 2. Folder of summarized vulnerability, and temperature and land-cover change per MCR 
# (VulHaz_data) used to calculate vulnerability and hazard per MCR pair, respectively.
#
# Outputs: 
# 1. Data frame of exposure (i.e. multispecies migratory connectivity) estimates and 
# uncertainty estimates (if applicable; only for MCR pairs with >1 species) per MCR pair (MSMC_est.csv).
# If only 1 species for a given MCR pair, exposure is the mean of the connectivity parameter for the given species
# and (bootstrapped) CIs are listed as NA. Note that in the paper, for connections with only 1 species, we report the
# 95% CI of the connectivity parameter estimated for the given species. See Methods for details.
#
# 2. Example figure (MSMC_WeightMeansCI.png) of exposure and uncertainty estimates per MCR pair for those
# pairs with > 1 species (i.e. weighted means and 80% and 95% bootstrapped CIs)
#
# 3. Data frame of risk calculations (Risk_components.csv): exposure, hazard, vulnerability, and 
# risk, both unscaled and scaled (normalized) values, for each MCR pair
#
# 4. Model summary from example post-hoc linear regression estimating relative importance of the
# three componenets of risk (exposure, vulnerability, hazard) in explaining variation in risk
#
# Outputs are saved within an output folder created within the working directory. Note that the script
# iterates through folders within the directory that have connectivity results for different species groups
# as desired, if separate exposure/risk outputs are desired for different focal species groups.
#
## Code developed by SP Saunders with assistance from coauthors DeLuca, Meehan,
## Taylor, and Michel.
#
## Analysis pertains to the manuscript entitled "Multispecies migratory connectivity
## reveals hemispheric risk to birds" by S. Saunders, W. DeLuca, B. Bateman, J. Deppe,
## J. Grand, E. Knight, T. Meehan, N. Michel, N. Seavy, M. Smith, L. Taylor, C. Witko,
## and C. Wilsey.
###################################################################################

#load libraries
library(dplyr) 
library(tidyr) 
library(ggplot2)
library(purrr) 
library(tidytext)
library(scales) 
library(jagsUI) 

setwd("Z:/5_MigratoryBirdInitiative/Migratory Connectivity/Connectivity files_June2024/") 
#working directory where data folders (Connectivity_data_ExSp and VulHaz_data) are located

dirs <- list.dirs()
dirs <- dirs[grepl("Connectivity_*", dirs)] 
# The above lists all folders with this pattern of path name (Connectivity_) in working directory. Script is designed to iterate through
# folders of species groups to estimate multispecies connectivity (exposure) per species group of interest
# e.g. waterbirds, landbirds; or separately per taxa if multitaxa data available, i.e. separate folders for
# mammals and birds. Risk is then estimated per species group.
# In this example, only 1 folder with 5 example bird species connectivity results is provided due to
# data persmissions that limit sharing of species-specific information.

for (direct in dirs) { #loop through folders to produce separate multispecies connectivity and risk estimates per group
  dir_code <- direct
  
# create output folder for results within species group folders
dir.create(path = paste0(dir_code,"/Output/"), recursive = TRUE)

#######################################################################################################################################
## Step 1. Calculating multispecies connectivity (i.e. weighted means of species-specific connectivities) and uncertainty per MCR pair
## NOTE: 5 example species considered here: blackpoll warbler, grasshopper sparrow, ovenbird, prothonotary warbler, tree swallow
######################################################################################################################################

#pull species-specific result files for each group (i.e. folders within working directory)
conprops <- data.frame()
lf <- list.files(dir_code,pattern = "*_MigConnect_Results.csv")
for (i in 1:length(lf)){
  temp <- read.csv(paste0(dir_code,"/",lf[i]))
  conprops <- rbind(conprops, data.frame(Species=substr(lf[i],1,6), Mean=temp$mean, SDev=temp$sdev, Lower95=temp$low95,
                                         Upper95=temp$up95,Br_MCR=temp$brdg, NB_MCR=temp$nonbrdg
  ))                                
}

#add grouping column to indicate MCR pair
conprops.allspp <- conprops %>% unite('Group', Br_MCR:NB_MCR, remove=FALSE)

#Calculate weighted means and weighted variances (and CI of weighted means) per mcr pairs. This is an iterative procedure to estimate mean and 
#variance from a set of measurements `x` (species-specific mean connectivity estimates) made with known measurement errors `sigma.i` (species-
#specific standard deviations fo connectivity estimates). See Methods for more details.

varmean.wt <- function(Group.Name, data){
  # limit to rows for that group
  df<- data %>% filter(Group == Group.Name)
  # Initialize.
  v <- var(df$x)
  n <- length(df$x)
  S <- mean(df$sigma.i^2) # Avoids recomputing this constant within the loop
  for (i in 1:400) { #will iterate 400 times to update (based on prelim runs; convergence reached w/in 100-200 iterations)
    w <- 1/(v + df$sigma.i^2); w <- w / sum(w) #compute weights
    m <- sum(df$x * w) #update the estimated (weighted) mean
    v.0 <- max(0, sum((df$x-m)^2 * w) * n / (n-1) - S) #bias-adjusted weighted variance
    v <- v.0
  }
  
  # put results in dataframe
  df <- data.frame("Group" = Group.Name, "WeightVar" = v, "WeightMean" = m)
  
  # format Group Col
  df$Group <- as.character(df$Group)
  return(df)
}

#rename Mean and SDev for use of function
conprops.allspp <- conprops.allspp %>% 
  mutate(x = Mean) %>% mutate(sigma.i = SDev) %>%
  as.data.frame()

# unique group names need to loop through for calcs
Group.Names <- unique(conprops.allspp$Group)

# take each group name and use the function above using purrr package
weightmeanvars.allspp <- Group.Names %>% map_df(~ varmean.wt(Group.Name = .x, data = conprops.allspp))

#conduct parametric bootstrap to get 95% CI and 80% CI of weighted mean for each mcr pair
#first merge conprops with weightmeanvars by group
conprops.allspp.merge <- merge(conprops.allspp, weightmeanvars.allspp, by = "Group")

#var.wt function to use (same as above but not within group loop)
var.wt <- function(x, sigma.i) {
  # Initialize
  v <- var(x)
  n <- length(x)
  S <- mean(sigma.i^2) 
  for (i in 1:400) {  
    w <- 1/(v + sigma.i^2); w <- w / sum(w)
    m <- sum(x * w)
    v.0 <- max(0, sum((x-m)^2 * w) * n / (n-1) - S) 
    v <- v.0
  }
  list(Variance=v, Mean=m, Iterations=i)
}

#function to run bootstrap
boot.ci <- function(Group.Name, data){
  # limit to rows for that group
  df<- data %>% filter(Group == Group.Name)
  X <- replicate(10000, {  #set to 10k draws based on prelim runs producing virtually identical CIs to at least the thousandth decimal place 
    y <- rnorm(nrow(df), mean(df$WeightMean), sqrt(mean(df$WeightVar))) #take means of WeightMean and WeightVar cols since identical per group
    epsilon <- rnorm(nrow(df), 0, df$sigma.i) 
    x <- y + epsilon #apply random measurement error to weighted means with the given SDs
    var.wt(x, df$sigma.i)
  })
  X <- matrix(unlist(X), 3)
  lowci <- quantile(X[2,], 0.025) #95% CI of weighted mean
  upci <- quantile(X[2,], 0.975) #95% CI
  lowci80 <- quantile(X[2,], 0.10) #80% CI of weighted mean
  upci80 <- quantile(X[2,], 0.90)  #80% CI
  
  # put results in dataframe
  df <- data.frame("Group" = Group.Name, "Low95CI.WeightMean" = lowci, "Up95CI.WeightMean" = upci, "Low80CI.WeightMean" = lowci80, "Up80CI.WeightMean" = upci80)
  
  # format Group Col
  df$Group <- as.character(df$Group)
  return(df)
}

# unique group names to loop through for calcs that are groups with >1 species 
names.red <- conprops.allspp.merge %>% group_by(Group) %>% count() %>% filter(n > 1)
Group.Names <- unique(names.red$Group) #keep only groups with >1 species (otherwise use nonweighted mean and CI)

# run bootstrap on all remaining mcr pairs using purrr to get 95% CI and 80% CI of weighted mean. NOTE est run time: ~15 mins
weightmeancis.allspp <- Group.Names %>% map_df(~ boot.ci(Group.Name = .x, data = conprops.allspp.merge))

#merge weightmeanvars with weightmeancis by group, adding NA when group has only 1 entry 
weightdf.allspp <- merge(weightmeanvars.allspp, weightmeancis.allspp, by =  "Group", all.x = TRUE) 

#calculate standard (unweighted) means when weighted mean is NA (only 1 species)
multimeans.allspp <- conprops.allspp %>% group_by(Br_MCR, NB_MCR) %>% summarise(Multimean = mean(Mean)) %>% as.data.frame()

#create matching group column to merge with weightdf
multimeans.allspp <- multimeans.allspp %>% unite('Group', Br_MCR:NB_MCR, remove=FALSE) 

#merge unweighted and weighted dfs
mergedf.allspp <- merge(multimeans.allspp, weightdf.allspp, by =  "Group", all.x = TRUE) #merge dfs of orginal calcs plus weighted calcs
mergedf.allspp <- mergedf.allspp %>% mutate(AllMean = coalesce(WeightMean,Multimean)) #create single col for all means (either weight or unweight)

#clean up df
finaldf.new <- mergedf.allspp %>%
  select(Group, Br_MCR, NB_MCR, AllMean, Multimean, WeightMean,  
         Low80CI.WeightMean, Up80CI.WeightMean, Low95CI.WeightMean, Up95CI.WeightMean) %>%
  arrange(Br_MCR, NB_MCR)

###Example plot of exposure estimates for MCR pairs with >1 species #################################################################
#example caterpillar plot of weighted means and 80% and 95% CIs for MCR pairs with > 1 species (not including MCR pairs with 1 species)
finaldf.new.wm <- finaldf.new %>% filter(!is.na(WeightMean))
#plot
ex_plot <- ggplot(data = finaldf.new.wm, aes(x = tidytext::reorder_within(NB_MCR,-WeightMean,Br_MCR))) + 
  geom_errorbar(aes(ymin=Low95CI.WeightMean, ymax=Up95CI.WeightMean),size=0.6,width=0)+ 
  geom_errorbar(aes(ymin=Low80CI.WeightMean, ymax=Up80CI.WeightMean),size=1.6,width=0)+ 
  geom_point(aes(y=WeightMean),size=2, col="black")+ 
  tidytext::scale_x_reordered()+
  coord_flip()+
  labs(y="Weighted mean proportion",x="Nonbreeding MCR")+
  ggtitle("Exposure and uncertainty for connections with >1 species")+
  theme_bw() +
  theme(axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        axis.title.y=element_text(size = 11), 
        axis.title.x=element_text(size = 11),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour="black"),
        panel.background = element_rect(fill = 'white', colour = 'white'),
        plot.margin = margin(0.5,1,0.5,0.5,"cm")) +
  facet_wrap(~Br_MCR, ncol = 6, scales = "free_y")

png(file=paste0(dir_code,"/Output/MSMC_WeightMeanCI.png"), width=12, height=5, units="in", res=300)
print(ex_plot)
dev.off()

#organize and save final df of multispecies connectivity estimates per mcr combos [EXPOSURE] #####################################
finaldf.save <- finaldf.new %>% 
  select(Group, Br_MCR, NB_MCR, AllMean,  
         Low80CI.WeightMean, Up80CI.WeightMean, Low95CI.WeightMean, Up95CI.WeightMean) %>%
  rename(Exposure = AllMean) %>%
  arrange(Br_MCR, NB_MCR)
write.csv(finaldf.save, file=paste0(dir_code,"/Output/MSMC_est.csv"), row.names = FALSE) #save in Output folder in WD
#NOTE: the AllMean field represents Exposure, or the average connectivity per MCR pair (either weighted mean if >1 species or the unweighted mean if 1 species)
# 80% and 95% CIs are also saved; note that if CIs are NA, the given connection represents only 1 species (i.e. bootstrapped CIs of weighted mean not calculated)

######################################################################################################
## Step 2. Calculate average vulnerability per MCR pair
# Vulnerability is defined as PIF ACAD CCS max scores per species. See Methods for more details.
######################################################################################################
#read in vulnerability scores dataset from ACAD [for 5 example species]
vulspp <- read.csv("VulHaz_data/Vul_Group_ExSpecies.csv") #file located in working directory

#summarize by average vulnerability per MCR pair [note this will differ from full analysis using 112 species]
vulsppavg <- vulspp %>% group_by(Group) %>% summarise(MeanCCSmax=mean(ccs.max)) %>% as.data.frame()

############################################################################################
## Step 3. Combine temp change and land cover vulnerability per MCR pair to yield hazard
# See Methods for more details on IPCC temp change layer and ESRI Clark Labs LC change layer
#############################################################################################
#temp change
tempchange_mcr_df <- read.csv(file="VulHaz_data/TempChange_MCR.csv", header = TRUE) #file located in working directory

#calculate average temp change per MCR pair (group)
tempchange_groups <- tempchange_mcr_df %>%
  rename(Br_MCR=MCR_Code) %>%
  mutate(NB_MCR=Br_MCR) %>%
  as.data.frame()
tempchange_brmcrs <- tempchange_groups %>%
  select(Br_MCR, meantempchange_exex)
tempchange_nbmcrs <- tempchange_groups %>%
  select(NB_MCR, meantempchange_exex)
finaldf.tc <- finaldf.new %>%
  select(Group, Br_MCR, NB_MCR) %>%
  left_join(., tempchange_brmcrs, by="Br_MCR") %>%
  rename(Br_tc=meantempchange_exex) %>%
  left_join(., tempchange_nbmcrs, by="NB_MCR") %>%
  rename(NB_tc=meantempchange_exex) %>%
  mutate(BNB_meantc=(Br_tc+NB_tc)/2) %>%
  as.data.frame()

#lc change
lcvul_mcr_df <- read.csv(file="VulHaz_data/LandCoverVul2050_MCR.csv", header = TRUE) #file in working directory

#calculate average LC change per MCR pair
lulcvul_groups <- lcvul_mcr_df %>%
  rename(Br_MCR=MCR_Code) %>%
  mutate(NB_MCR=Br_MCR) %>%
  as.data.frame()
lcvul_brmcrs <- lulcvul_groups %>%
  select(Br_MCR, mean_lcv)
lcvul_nbmcrs <- lulcvul_groups %>%
  select(NB_MCR, mean_lcv)
finaldf.lc <- finaldf.new %>%
  select(Group, Br_MCR, NB_MCR) %>%
  left_join(., lcvul_brmcrs, by="Br_MCR") %>%
  rename(Br_lccvul=mean_lcv) %>%
  left_join(., lcvul_nbmcrs, by="NB_MCR") %>%
  rename(NB_lcvul=mean_lcv) %>%
  mutate(BNB_meanlcvul=(Br_lccvul+NB_lcvul)/2) %>%
  as.data.frame()

#first normalize (using min max scaling) both metrics so they are able to be summed equivalently
#pull cols from respective dfs:
finaldf.tc.red <- finaldf.tc %>% 
  mutate(BNB_meantc.norm = (BNB_meantc-min(BNB_meantc, na.rm=TRUE))/(max(BNB_meantc, na.rm=TRUE)-min(BNB_meantc, na.rm=TRUE))) %>%
  select(Group, Br_MCR, NB_MCR, BNB_meantc, BNB_meantc.norm)

finaldf.lc.red <- finaldf.lc %>% 
  mutate(BNB_meanlc.norm = (BNB_meanlcvul-min(BNB_meanlcvul, na.rm=TRUE))/(max(BNB_meanlcvul, na.rm=TRUE)-min(BNB_meanlcvul, na.rm=TRUE))) %>%
  select(Group, Br_MCR, NB_MCR, BNB_meanlcvul, BNB_meanlc.norm)

#merge dfs and create one hazard metric as summed normalized metrics [note this will differ from full analysis using full 921 MCR pairs]
finaldf.haz <- left_join(finaldf.tc.red, finaldf.lc.red) %>%
  mutate(BNB_MeanHazard=BNB_meantc.norm+BNB_meanlc.norm)

#########################################################################
## Multiplying expsoure, vulnerability and hazard to yield risk 
########################################################################

#exposure, vulnerability and hazard by group 
hazdf.risk <- finaldf.haz %>% select(Group, Br_MCR, NB_MCR, BNB_MeanHazard)
finaldf.vulavg <- left_join(finaldf.new,vulsppavg)
convuldf.risk <- finaldf.vulavg %>% select(Group, Br_MCR, NB_MCR, AllMean, MeanCCSmax)
finaldf.risk <- left_join(convuldf.risk, hazdf.risk)

#rescale each metric between 1 and 100 before multiplying; create risk metric as product of rescaled metrics [note this will differ from full analysis results]
finaldf.risk <- finaldf.risk %>%
  mutate(connect.rescale = scales::rescale(AllMean, to = c(1, 100)),
         vul.rescale = scales::rescale(MeanCCSmax, to = c(1, 100)),
         haz.rescale = scales::rescale(BNB_MeanHazard, to = c(1, 100)),
         mbrisk = connect.rescale*vul.rescale*haz.rescale,
         mbrisk.norm = (mbrisk-min(mbrisk,na.rm=TRUE))/(max(mbrisk,na.rm=TRUE)-min(mbrisk,na.rm=TRUE))) #normalize risk between 0 and 1 for interpretability

finaldf.risk.save <- finaldf.risk %>%
  rename(Exposure=AllMean, Vulnerability=MeanCCSmax, Hazard=BNB_MeanHazard,Risk=mbrisk,Norm_Risk=mbrisk.norm,
       Norm_Exposure=connect.rescale, Norm_Vulnerability=vul.rescale, Norm_Hazard=haz.rescale) #rename fields for interpretability

#save final risk calculations per MCR pair
write.csv(finaldf.risk.save, file=paste0(dir_code,"/Output/Risk_components.csv"), row.names = FALSE) #save in Output file in WD

#############################################################################################################################################
# Example post-hoc linear regression to evaluate drivers of risk
# Note that results differ from full analysis, but exposure is still the strongest driver of variation in risk (i.e. largest effect size)
# for this set of 5 example species, followed by vulnerability then hazard
##############################################################################################################################################
# prep data for jags
finaldf.risk.jags <- finaldf.risk %>%
  select(connect.rescale,vul.rescale,haz.rescale,mbrisk.norm) %>%
  mutate(connect.scale=(connect.rescale-mean(connect.rescale,na.rm = TRUE))/sd(connect.rescale,na.rm = TRUE), #standardize predictors for regression
         vul.scale=(vul.rescale-mean(vul.rescale,na.rm = TRUE))/sd(vul.rescale,na.rm = TRUE),
         haz.scale=(haz.rescale-mean(haz.rescale,na.rm = TRUE))/sd(haz.rescale,na.rm = TRUE)) %>%
  select(connect.scale,vul.scale,haz.scale,mbrisk.norm) %>%
  as.data.frame()
jagsdata_s1 <- with(finaldf.risk.jags, list(mbrisk.norm = mbrisk.norm, connect.scale = connect.scale,vul.scale=vul.scale, haz.scale=haz.scale,N = length(mbrisk.norm)))

# write model file
modfile <- tempfile()
writeLines("
model{
  # Likelihood:
  for (i in 1:N){
    mbrisk.norm[i] ~ dnorm(mu[i], tau) # tau is precision (1 / variance)
    mu[i] <- alpha + beta1 * connect.scale[i] + beta2 * vul.scale[i] + beta3 * haz.scale[i]
  }
  # Priors:
  alpha ~ dnorm(0, 0.01) # intercept
  beta1 ~ dnorm(0, 0.01) # slope
  beta2 ~ dnorm(0, 0.01) # slope
  beta3 ~ dnorm(0, 0.01) # slope
  sigma ~ dunif(0, 100) # standard deviation
  tau <- 1 / (sigma * sigma) 
}
", con = modfile)

init_values <- function(){
  list(alpha = rnorm(1), beta1 = rnorm(1), beta2 = rnorm(1), beta3 = rnorm(1), sigma = runif(1))
}
params <- c("alpha", "beta1","beta2","beta3", "sigma")

fit_lm1 <- jags(data = jagsdata_s1, inits = init_values, parameters.to.save = params, model.file = modfile,
                n.chains = 3, n.iter = 12000, n.burnin = 2000, n.thin = 10)

#save output (beta1=exposure; beta2=vulnerability; beta3=hazard)
summary.save <- as.data.frame(fit_lm1$summary)
write.csv(summary.save, file=paste0(dir_code,"/Output/Risk_regression_modsum.csv")) #write to output folder 

} #End loop through species group folders

## END SCRIPT