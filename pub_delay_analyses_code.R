library('data.table')
library('vcd')
library('lme4')
library('MASS')
library('optimx')
library('scales')
library('ggpubr')
library('ggplot2')
library('DHARMa')
library('AER')
library('multcomp')
library('RColorBrewer')
library('ggridges')
library('ggeffects')
library('emmeans')
library('MuMIn')


##### read in ce_data_pubdelay_update_corr.csv
masterall <- fread(choose.files())

##### replace NAs for end date with start date if only has start date (i.e. study conducted within one year)
masterall[is.na(stdate2)&!is.na(stdate1),stdate2:=stdate1]

##### note studies where there are no dates available - we will exclude these (made up ~3% of all studies in database)
masterall[is.na(stdate2)&is.na(stdate1),nodate:="ND"]

##### calculate difference between years (publication year - year data collection ended (stdate2))
masterall$datediff = masterall$pubdate-masterall$stdate2

##### remove reviews and studies without dates.
masterall_trim <- unique(masterall[review==0&!is.na(stdate1)&!is.na(stdate2)&!is.na(datediff),])


##### create dataset for main model (testing predictors of publication delay including publication source, publication date, and Conservation Evidence synopsis (i.e., conservation subject))
masterall_trim_mainmod <- unique(masterall_trim[,list(pageid,before,controlled,randomised,datediff,syn,int,pubdate,peer_review)])

######## test and run poisson model first
mod1<-glm(datediff~pubdate+syn+peer_review,family=poisson(link="log"),data=masterall_trim_mainmod)

summary(mod1)
#plot(mod1)
#mod1$coefficients
#mod1$deviance/mod1$df.residual
# 
# ####### check proportion of 0's in the data for zero inflation
# data.tab<-table(masterall_trim_mainmod$datediff==0)
# data.tab/sum(data.tab)
# #proportion of 0's expected from a Poisson distribution
# mu <- mean(masterall_trim_mainmod$datediff)
# cnts <- rpois(1000, mu)
# data.tab <- table(cnts == 0)
# data.tab/sum(data.tab)
# 
# ##### some zero inflation. Expected and observed zeros are slightly different.
# 
# 
# ###### check for overdispersion
# mod1$deviance/mod1$df.resid
# ###### or alternatively,mod1 via Pearson's residuals
# Resid <- resid(mod1, type = "pearson")
# sum(Resid^2) / (nrow(masterall_trim_mainmod) - length(coef(mod1)))
# 
# testDispersion(simulateResiduals(mod1,  refit=T))
# testZeroInflation(simulateResiduals(mod1, refit=T))
# 
# dispersiontest(mod1)
# 
# #### darma tests
# 
# dat.sim <- simulateResiduals(mod1, n=1000)
# plot(dat.sim)
# testOutliers(dat.sim,type='bootstrap',nBoot = 1000)
# 
# ###### there seems to be an issue with overdispersion with theta much greater than 1.
# 
# 
# ##########Pearson's χ2 residuals - explores whether there are any significant patterns remaining in the residuals
# dat.resid <- sum(resid(mod1, type = "pearson")^2)
# 1 - pchisq(dat.resid, mod1$df.resid)
# 
# 
# #Deviance (G2) - similar to the χ2 test above, yet uses deviance
# 1-pchisq(mod1$deviance, mod1$df.resid)
# 
# #The DHARMa package also has a routine for running a Kologorov-Smirnov test test to explore overall uniformity of the residuals as a goodness-of-fit test on the scaled residuals.
# testUniformity(dat.sim)
# 
# 
# ###fitted vs observed data inspection
# tempdata <- data.frame(obs=masterall_trim_mainmod$datediff, fitted=fitted(mod1))
# 
# ggplot(tempdata, aes(y=fitted, x=obs)) + geom_point()
# 
# #### quasi-R2
# 1-(mod1$deviance/mod1$null)
# 
# ####### overall poisson model seems poor and suffers from overdispersion
# 
# 
# ####### try poisson without big outliers (i.e., with delay > 20 years)
# masterall_trim_mainmod_noout <- unique(masterall_trim[datediff<=20,list(pageid,before,controlled,randomised,datediff,syn,int,pubdate,peer_review)])
# 
# mod1<-glm(datediff~pubdate+syn+peer_review,family=poisson(link="log"),data=masterall_trim_mainmod_noout)
# summary(mod1)
# #plot(mod1)
# #mod1$coefficients
# #mod1$deviance/mod1$df.residual
# 
# #proportion of 0's in the data
# data.tab<-table(masterall_trim_mainmod_noout$datediff==0)
# data.tab/sum(data.tab)
# #proportion of 0's expected from a Poisson distribution
# mu <- mean(masterall_trim_mainmod_noout$datediff)
# cnts <- rpois(1000, mu)
# data.tab <- table(cnts == 0)
# data.tab/sum(data.tab)
# 
# ##### some zero inflation. Expected and observed zeros are slighlty different.
# 
# ######overdispersion
# mod1$deviance/mod1$df.resid
# #Or alternatively,mod1 via Pearson's residuals
# Resid <- resid(mod1, type = "pearson")
# sum(Resid^2) / (nrow(masterall_trim_mainmod_noout) - length(coef(mod1)))
# 
# testDispersion(simulateResiduals(mod1,  refit=T))
# testZeroInflation(simulateResiduals(mod1, refit=T))
# 
# dispersiontest(mod1)
# 
# #### darma tests
# 
# dat.sim <- simulateResiduals(mod1, n=1000)
# plot(dat.sim)
# testOutliers(dat.sim,type='bootstrap',nBoot = 1000)
# 
# #### still issues with overdispersion
# 
# ##########Pearson's χ2 residuals - explores whether there are any significant patterns remaining in the residuals
# dat.resid <- sum(resid(mod1, type = "pearson")^2)
# 1 - pchisq(dat.resid, mod1$df.resid)
# 
# #Deviance (G2) - similar to the χ2 test above, yet uses deviance
# 1-pchisq(mod1$deviance, mod1$df.resid)
# 
# #The DHARMa package also has a routine for running a Kologorov-Smirnov test test to explore overall uniformity of the residuals as a goodness-of-fit test on the scaled residuals.
# testUniformity(dat.sim)
# 
# ###fitted vs observed data inspection
# tempdata <- data.frame(obs=masterall_trim_mainmod_noout$datediff, fitted=fitted(mod1))
# 
# ggplot(tempdata, aes(y=fitted, x=obs)) + geom_point()
# 
# #### quasi-R2
# 1-(mod1$deviance/mod1$null)
# 
# 
########## poisson had overdispersion, even with removed outliers, try quasi-Poisson

########
masterall_trim_mainmodqp <- unique(masterall_trim[,list(pageid,before,controlled,randomised,datediff,syn,int,pubdate,peer_review)])

######## create factor levels in order to aid model interpretation and multiple comparisons - compare against category with lowest mean delay
######## data exploration suggests Bee Conservation has lowest mean delay, as does non-peer-reviewed (peer_review==0)

masterall_trim_mainmodqp$syn = factor(masterall_trim_mainmodqp$syn,levels=c("Bee Conservation",unique(masterall_trim_mainmodqp[syn!="Bee Conservation",syn])))
masterall_trim_mainmodqp$peer_review = factor(masterall_trim_mainmodqp$peer_review, levels=c(0,1))
mod1 <- glm(datediff~pubdate+syn+peer_review, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp)
summary(mod1)

######quasi-AIC model selection
#deviance calc
dfun <- function(object) {
  with(object,sum((weights * residuals^2)[weights > 0])/df.residual)
}
#reuses AIC from poisson family estimation
x.quasipoisson <- function(...) {
  res <- quasipoisson(...)
  res$aic <- poisson(...)$aic
  res
}
gmod1 <- update(mod1,family="x.quasipoisson",
       na.action=na.fail)
(gg <- dredge(gmod1,rank="QAIC", chat=dfun(gmod1)))


#####################################################################################################

####### journals multiple comparisons based on quasi-Poisson

emm_pr <- emmeans(mod1, ~peer_review, type="response")
confint(emm_pr)
test(emm_pr)
contrast(emm_pr, method="pairwise")


#####################################################################################################
##############publication date###########################################
#########################################################################
ggemmeans(mod1, ~pubdate)

########################################
#######sensitivity analysis#############
########################################
mod1980 <- glm(datediff~pubdate+syn+peer_review, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp[pubdate>=1980])
mod1990 <- glm(datediff~pubdate+syn+peer_review, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp[pubdate>=1990])
mod2000 <- glm(datediff~pubdate+syn+peer_review, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp[pubdate>=2000])
mod2010 <- glm(datediff~pubdate+syn+peer_review, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp[pubdate>=2010])

summary(mod1980)
ggemmeans(mod1980, ~pubdate)

summary(mod1990)
ggemmeans(mod1990, ~pubdate)

summary(mod2000)
ggemmeans(mod2000, ~pubdate)

summary(mod2010)
ggemmeans(mod2010, ~pubdate)


####################################################################
#############synopsis###############################################
####################################################################
emm_syn <- emmeans(mod1, ~syn, type="response")
confint(emm_syn)
contrast(emm_syn, method="pairwise")


#######################################################################################################
###### taxonomic model ##################################################################
####################################################################################
###### read in ampdat.csv, birddat.csv, and mamdat.csv (contains IUCN Red List names and statuses to match to our dataset on delays)
ampdat <- fread(choose.files())
birddat <- fread(choose.files())
mamdat <- fread(choose.files())

masterall_trim_taxa_all <- unique(masterall_trim[,list(pageid,before,controlled,randomised,int,syn,stdate1,stdate2,datediff,binom,class,peer_review,pubdate)])

ampdat$scientificName<-tolower(ampdat$scientificName)
namesamp <- unique(masterall_trim_taxa_all[syn=="Amphibian Conservation"&class=="Amphibia"&!is.na(masterall_trim_taxa_all$binom),binom])
for(i in 1:length(namesamp)){masterall_trim_taxa_all[binom==namesamp[i],iucncat := unique(ampdat[scientificName==namesamp[i],redlistCategory])]}

unique(masterall_trim_taxa_all[syn=="Amphibian Conservation"&class=="Amphibia"&!is.na(masterall_trim_taxa_all$binom),iucncat])


mamdat$scientificName<-tolower(mamdat$scientificName)
namesmam <- unique(masterall_trim_taxa_all[(syn=="Bat Conservation"|syn=="Primate Conservation"|syn=="Terrestrial Mammal Conservation")&class=="Mammalia"&!is.na(masterall_trim_taxa_all$binom),binom])
for(i in 1:length(namesmam)){
  if(length(unique(mamdat[scientificName==namesmam[i],redlistCategory]))>0){
    masterall_trim_taxa_all[binom==namesmam[i],iucncat := unique(mamdat[scientificName==namesmam[i],redlistCategory])]
  }
}
unique(masterall_trim_taxa_all[(syn=="Bat Conservation"|syn=="Primate Conservation"|syn=="Terrestrial Mammal Conservation")&class=="Mammalia"&!is.na(masterall_trim_taxa_all$binom),iucncat])
unique(masterall_trim_taxa_all[is.na(iucncat)&(syn=="Bat Conservation"|syn=="Primate Conservation"|syn=="Terrestrial Mammal Conservation")&class=="Mammalia"&!is.na(masterall_trim_taxa_all$binom),binom])
setdiff(namesmam,mamdat$scientificName)
#### some mammal species cannot be matched. These tend to be domesticated or hybrid animals.

birddat$scientificName<-tolower(birddat$scientificName)
namesbird <- unique(masterall_trim_taxa_all[syn=="Bird Conservation"&class=="Aves"&!is.na(masterall_trim_taxa_all$binom),binom])
for(i in 1:length(namesbird)){
  masterall_trim_taxa_all[binom==namesbird[i],iucncat := unique(birddat[scientificName==namesbird[i],redlistCategory])]
}
unique(masterall_trim_taxa_all[syn=="Bird Conservation"&class=="Aves"&!is.na(masterall_trim_taxa_all$binom),iucncat])

############# only include studies on amphibians, birds, or mammals from relevant taxa-specific synopses 
masterall_trim_taxa_allclean <- masterall_trim_taxa_all[(syn=="Bird Conservation"|syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"|syn=="Amphibian Conservation")&!is.na(iucncat)&(class=="Aves"|class=="Mammalia"|class=="Amphibia")&!is.na(datediff)]

############# shorten IUCN categories to acronyms
masterall_trim_taxa_allclean[iucncat=="Data Deficient",iucncat2:="DD"]
masterall_trim_taxa_allclean[iucncat=="Least Concern",iucncat2:="LC"]
masterall_trim_taxa_allclean[iucncat=="Near Threatened",iucncat2:="NT"]
masterall_trim_taxa_allclean[iucncat=="Vulnerable",iucncat2:="VU"]
masterall_trim_taxa_allclean[iucncat=="Endangered",iucncat2:="EN"]
masterall_trim_taxa_allclean[iucncat=="Critically Endangered",iucncat2:="CR"]
masterall_trim_taxa_allclean[iucncat=="Extinct in the Wild",iucncat2:="EW"]

masterall_trim_taxa_allclean$iucncat2<-factor(masterall_trim_taxa_allclean$iucncat2,levels=c("DD","LC","NT","VU","EN","CR","EW"))
#masterall_trim_taxa_allclean[is.na(iucncat2)&!is.na(iucncat)]

###too few studies on EW and DD
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="EW",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="DD",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="NT",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="VU",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="EN",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="CR",list(pageid,datediff,before,controlled,randomised,pubdate)]))

##################### remove EW and DD category studies
masterall_trim_taxa_allclean <- unique(masterall_trim_taxa_allclean[iucncat2!="EW" & iucncat2!="DD",list(iucncat,iucncat2,syn,int,class,binom,pageid,datediff,before,controlled,randomised,peer_review,pubdate)])

######### only consider most threatened species in each study
uniqstuds <- unique(masterall_trim_taxa_allclean[,list(syn,int,pageid,datediff,before,controlled,randomised,class,peer_review,pubdate)])
masterall_trim_taxa_allcleanmost <- masterall_trim_taxa_allclean

for(i in 1:nrow(uniqstuds)){
  tokeep<-max(as.numeric(masterall_trim_taxa_allcleanmost[class==uniqstuds$class[i]&syn==uniqstuds$syn[i]&int==uniqstuds$int[i]&pageid==uniqstuds$pageid[i]&datediff==uniqstuds$datediff[i]& before==uniqstuds$before[i]&controlled==uniqstuds$controlled[i]&randomised==uniqstuds$randomised[i],iucncat2]))
  masterall_trim_taxa_allclean[class==uniqstuds$class[i]&syn==uniqstuds$syn[i]&int==uniqstuds$int[i]&pageid==uniqstuds$pageid[i]&datediff==uniqstuds$datediff[i]& before==uniqstuds$before[i]&controlled==uniqstuds$controlled[i]&randomised==uniqstuds$randomised[i],iucncat2:=tokeep]
}

masterall_trim_taxa_allcleanmostuniq <- masterall_trim_taxa_allcleanmost[,list(int,syn,pageid,datediff,before,controlled,randomised,iucncat2,class,peer_review,pubdate)]

######### order different factors
masterall_trim_taxa_allcleanmostuniq$syn = factor(masterall_trim_taxa_allcleanmostuniq$syn,levels=c("Bee Conservation",unique(masterall_trim_taxa_allcleanmostuniq[syn!="Bee Conservation",syn])))
masterall_trim_taxa_allcleanmostuniq$peer_review = factor(masterall_trim_taxa_allcleanmostuniq$peer_review,levels=sort(unique(masterall_trim_taxa_allcleanmostuniq$peer_review)))
masterall_trim_taxa_allcleanmostuniq$iucncat2 = factor(masterall_trim_taxa_allcleanmostuniq$iucncat2,levels=c("LC","NT","VU","EN","CR","EW"))

###### run a GLM for taxa including other variables from original model
mod2 <- glm(datediff~pubdate+syn+peer_review+iucncat2, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq)
summary(mod2)
#plot(mod1)

gmod2 <- update(mod2,family="x.quasipoisson",
                na.action=na.fail)
(gg <- dredge(gmod2,rank="QAIC", chat=dfun(gmod2)))

############# run mulitple comparisons for IUCN categories and output table
emm_iucn <- emmeans(mod2, ~iucncat2, type="response")
confint(emm_iucn)
contrast(emm_iucn,method="pairwise")


############# now run taxa specific analyses

###amphibian
amp_only <- masterall_trim_taxa_allcleanmostuniq[syn=="Amphibian Conservation"]
mod2a <- glm(datediff~pubdate+peer_review+iucncat2, family=quasipoisson(link="log"),data=amp_only)
summary(mod2a)

gmod2a <- update(mod2a,family="x.quasipoisson",
                na.action=na.fail)
(gg <- dredge(gmod2a,rank="QAIC", chat=dfun(gmod2a)))

#mod2a <- glm(datediff~pubdate+peer_review, family=quasipoisson(link="log"),data=amp_only)
#summary(mod2a)

#best model doesn't include iucn category so for reference just use original model to show comparisons
emm_iucn_amp <- emmeans(mod2a, ~iucncat2, type="response")
confint(emm_iucn_amp)
contrast(emm_iucn_amp,method="pairwise")


###bird
bird_only <- masterall_trim_taxa_allcleanmostuniq[syn=="Bird Conservation"]
mod2b <- glm(datediff~pubdate+peer_review+iucncat2, family=quasipoisson(link="log"),data=bird_only)
summary(mod2b)

gmod2b <- update(mod2b,family="x.quasipoisson",
                 na.action=na.fail)
(gg <- dredge(gmod2b,rank="QAIC", chat=dfun(gmod2b)))

mod2b <- glm(datediff~iucncat2, family=quasipoisson(link="log"),data=bird_only)
summary(mod2b)

emm_iucn_bird <- emmeans(mod2b, ~iucncat2, type="response")
confint(emm_iucn_bird)
contrast(emm_iucn_bird,method="pairwise")


###mammal
mam_only <- masterall_trim_taxa_allcleanmostuniq[syn=="Terrestrial Mammal Conservation"|syn=="Primate Conservation"|syn=="Bat Conservation"]
mod2m <- glm(datediff~pubdate+peer_review+iucncat2+syn, family=quasipoisson(link="log"),data=mam_only)
summary(mod2m)

gmod2m <- update(mod2m,family="x.quasipoisson",
                 na.action=na.fail)
(gg <- dredge(gmod2m,rank="QAIC", chat=dfun(gmod2m)))

mod2m <- glm(datediff~pubdate+peer_review+iucncat2+syn, family=quasipoisson(link="log"),data=mam_only)
summary(mod2m)

emm_iucn_mam <- emmeans(mod2m, ~iucncat2, type="response")
confint(emm_iucn_mam)
contrast(emm_iucn_mam,method="pairwise")
