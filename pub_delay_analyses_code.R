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

##### read in ce_data_pubdelay.csv
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
masterall_trim_mainmod <- unique(masterall_trim[,list(pageid,before,controlled,randomised,datediff,syn,int,pubdate,pubtype)])

######## test and run poisson model first
mod1<-glm(datediff~pubdate+syn+pubtype,family=poisson(link="log"),data=masterall_trim_mainmod)

summary(mod1)
#plot(mod1)
#mod1$coefficients
#mod1$deviance/mod1$df.residual

####### check proportion of 0's in the data for zero inflation
data.tab<-table(masterall_trim_mainmod$datediff==0)
data.tab/sum(data.tab)
#proportion of 0's expected from a Poisson distribution
mu <- mean(masterall_trim_mainmod$datediff)
cnts <- rpois(1000, mu)
data.tab <- table(cnts == 0)
data.tab/sum(data.tab)

##### no need to worry about zero inflation. Expected and observed zeros are similar.


###### check for overdispersion
mod1$deviance/mod1$df.resid
###### or alternatively,mod1 via Pearson's residuals
Resid <- resid(mod1, type = "pearson")
sum(Resid^2) / (nrow(masterall_trim_mainmod) - length(coef(mod1)))

testDispersion(simulateResiduals(mod1,  refit=T))
testZeroInflation(simulateResiduals(mod1, refit=T))

dispersiontest(mod1)

#### darma tests

dat.sim <- simulateResiduals(mod1, n=1000)
plot(dat.sim)
testOutliers(dat.sim,type='bootstrap',nBoot = 1000)

###### there seems to be an issue with overdispersion with theta much greater than 1.


##########Pearson's χ2 residuals - explores whether there are any significant patterns remaining in the residuals
dat.resid <- sum(resid(mod1, type = "pearson")^2)
1 - pchisq(dat.resid, mod1$df.resid)


#Deviance (G2) - similar to the χ2 test above, yet uses deviance
1-pchisq(mod1$deviance, mod1$df.resid)

#The DHARMa package also has a routine for running a Kologorov-Smirnov test test to explore overall uniformity of the residuals as a goodness-of-fit test on the scaled residuals.
testUniformity(dat.sim)


###fitted vs observed data inspection
tempdata <- data.frame(obs=masterall_trim_mainmod$datediff, fitted=fitted(mod1))

ggplot(tempdata, aes(y=fitted, x=obs)) + geom_point()

#### quasi-R2
1-(mod1$deviance/mod1$null)

####### overall poisson model seems poor and suffers from overdispersion


####### try poisson without big outliers (i.e., with delay > 20 years)
masterall_trim_mainmod_noout <- unique(masterall_trim[datediff<=20,list(pageid,before,controlled,randomised,datediff,syn,int,pubdate,pubtype)])

mod1<-glm(datediff~pubdate+syn+pubtype,family=poisson(link="log"),data=masterall_trim_mainmod_noout)
summary(mod1)
#plot(mod1)
#mod1$coefficients
#mod1$deviance/mod1$df.residual

#proportion of 0's in the data
data.tab<-table(masterall_trim_mainmod_noout$datediff==0)
data.tab/sum(data.tab)
#proportion of 0's expected from a Poisson distribution
mu <- mean(masterall_trim_mainmod_noout$datediff)
cnts <- rpois(1000, mu)
data.tab <- table(cnts == 0)
data.tab/sum(data.tab)

#####no need to worry about zeros. Expected and observed zeros are similar.

######overdispersion
mod1$deviance/mod1$df.resid
#Or alternatively,mod1 via Pearson's residuals
Resid <- resid(mod1, type = "pearson")
sum(Resid^2) / (nrow(masterall_trim_mainmod_noout) - length(coef(mod1)))

testDispersion(simulateResiduals(mod1,  refit=T))
testZeroInflation(simulateResiduals(mod1, refit=T))

dispersiontest(mod1)

#### darma tests

dat.sim <- simulateResiduals(mod1, n=1000)
plot(dat.sim)
testOutliers(dat.sim,type='bootstrap',nBoot = 1000)

#### still issues with overdispersion

##########Pearson's χ2 residuals - explores whether there are any significant patterns remaining in the residuals
dat.resid <- sum(resid(mod1, type = "pearson")^2)
1 - pchisq(dat.resid, mod1$df.resid)

#Deviance (G2) - similar to the χ2 test above, yet uses deviance
1-pchisq(mod1$deviance, mod1$df.resid)

#The DHARMa package also has a routine for running a Kologorov-Smirnov test test to explore overall uniformity of the residuals as a goodness-of-fit test on the scaled residuals.
testUniformity(dat.sim)

###fitted vs observed data inspection
tempdata <- data.frame(obs=masterall_trim_mainmod_noout$datediff, fitted=fitted(mod1))

ggplot(tempdata, aes(y=fitted, x=obs)) + geom_point()

#### quasi-R2
1-(mod1$deviance/mod1$null)


########## poisson had overdispersion, even with removed outliers, try quasi-Poisson

######## 
masterall_trim_mainmodqp <- unique(masterall_trim[,list(pageid,before,controlled,randomised,datediff,syn,int,pubdate,pubtype)])

######## create factor levels in order to aid model interpretation and multiple comparisons - compare against category with lowest mean delay
######## data exploration suggests Management of Captive Animals has lowest mean delay
######## pubtype of non_journal_lit also has lowest mean delay
masterall_trim_mainmodqp$syn = factor(masterall_trim_mainmodqp$syn,levels=c("Management of Captive Animals",unique(masterall_trim_mainmodqp[syn!="Management of Captive Animals",syn])))
masterall_trim_mainmodqp$pubtype = factor(masterall_trim_mainmodqp$pubtype,levels=c("non_journal_lit",unique(masterall_trim_mainmodqp[pubtype!="non_journal_lit",pubtype])))
mod1 <- glm(datediff~pubdate+syn+pubtype, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp)
summary(mod1)
#plot(mod1)

#write.csv(summary(mod1)$coefficients,"modresults.csv")


#####################################################################################################

####### journals multiple comparisons based on quasi-Poisson

####### run and output table
tmctestpubtype <- glht(mod1, linfct = mcp(pubtype = "Tukey"))
#summary(tmctestpubtype)
#confint(tmctestpubtype)
sumtmctestpubtype <- summary(tmctestpubtype)
tabtmcpubtype <- data.table(comp=names(sumtmctestpubtype$test$coefficients),est=round(sumtmctestpubtype$test$coefficients,digits=3),sigma=round(sumtmctestpubtype$test$sigma,digits=3),tstat=round(sumtmctestpubtype$test$tstat,digits=3),pval=round(sumtmctestpubtype$test$pvalues,digits=3))
#write.csv(tabtmcpubtype,"pubtypetukeytests.csv")

##### make categories for x axis of publication delay to aid plotting - combining delays 11-20, and >20 years into separate categories.
masterall_trim_mainmodqp[datediff<=10, datedifffac := as.character(datediff)]
masterall_trim_mainmodqp[datediff>10 & datediff<=20, datedifffac := "11-20"]
masterall_trim_mainmodqp[datediff>20, datedifffac := "20+"]

####### plot findings for journals
####### obtain maxvalues for each category for plotting
mxvals2 <- sapply(unique(masterall_trim_mainmodqp$pubtype),function(y){return(max(table(masterall_trim_mainmodqp[pubtype==y,datedifffac])))})

####### add these max values to relevant pubtype category in dataset to aid plotting
for(i in 1:length(unique(masterall_trim_mainmodqp$pubtype))){masterall_trim_mainmodqp[pubtype==unique(masterall_trim_mainmodqp$pubtype)[i],maxdatediffjr:=mxvals2[i]]}

###### run simple glm with only pubtype (and without intercept) to obtain predicted values from model for plotting
mod1pubtype <- glm(datediff~pubtype-1, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp)
preds=predict(mod1pubtype,data.frame(pubtype=gsub("pubtype","",names(mod1pubtype$coefficients))),type="response",se.fit=TRUE)
modelvalspt <- data.table(pubtype=gsub("pubtype","",names(mod1pubtype$coefficients)),mns=preds$fit,upr=preds$fit+qnorm(0.975)*preds$se.fit,lwr=preds$fit-qnorm(0.975)*preds$se.fit)

###### add these predicted values to relevant pubtype category in dataset to aid plotting
for(i in 1:nrow(modelvalspt)){
  masterall_trim_mainmodqp[pubtype==modelvalspt$pubtype[i],modelmnpt:=modelvalspt$mns[i]]
  masterall_trim_mainmodqp[pubtype==modelvalspt$pubtype[i],modelmnuprpt:=modelvalspt$upr[i]]
  masterall_trim_mainmodqp[pubtype==modelvalspt$pubtype[i],modelmnlwrpt:=modelvalspt$lwr[i]]
}

###### Add new labels to categories for plotting
masterall_trim_mainmodqp[pubtype=="non_journal_lit",pubtypelabel:="Non-journal literature"]
masterall_trim_mainmodqp[pubtype=="journal",pubtypelabel:="Recognised journals"]
masterall_trim_mainmodqp[pubtype=="journal_notrec",pubtypelabel:="Unrecognised journals"]

##### order publication type labels to aid plotting
masterall_trim_mainmodqp$pubtypelabel <- factor(masterall_trim_mainmodqp$pubtypelabel, levels=c("Recognised journals","Unrecognised journals","Non-journal literature"))

masterall_trim_mainmodqp$datedifffac<-factor(masterall_trim_mainmodqp$datedifffac,levels=as.character(c(sort(as.numeric(unique(masterall_trim_mainmodqp[datediff<=10,datediff]))),"11-20","20+")))


##### create publication delay plot by journal
ggplot(masterall_trim_mainmodqp)+ 
  geom_density(aes(x=datedifffac,y=..count..,group=pubtypelabel),fill="grey50",stat="binline",binwidth=1,show.legend=FALSE)+
  scale_x_discrete(name="Publication delay (years)")+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies",breaks=pretty_breaks(n=4))+
  geom_segment(aes(x=modelmnpt+1,xend=modelmnpt+1,y=0,yend=maxdatediffjr,group=pubtypelabel),colour="red",linetype="solid",lwd=1.05)+
  geom_segment(aes(x=modelmnuprpt+1,xend=modelmnuprpt+1,y=0,yend=maxdatediffjr,group=pubtypelabel),colour="red",linetype="dashed",lwd=0.95)+
  geom_segment(aes(x=modelmnlwrpt+1,xend=modelmnlwrpt+1,y=0,yend=maxdatediffjr,group=pubtypelabel),colour="red",linetype="dashed",lwd=0.95)+
  facet_wrap(.~pubtypelabel,ncol=1,scales="free_y")+
  theme(aspect.ratio=0.33,strip.text = element_text(size=19,margin=margin(5,5,5,5,unit="pt")),axis.title=element_text(size=21),axis.text = element_text(size=18),strip.background = element_rect(fill="grey90"))

#ggsave("delayplotsforpaperbypubtype.svg",height=45,width=45,units="cm",device="svg")


#####################################################################################################

###### synopses multiple comparisons based on quasi-Poisson

####### run and output table
tmctestsyn <- glht(mod1, linfct = mcp(syn = "Tukey"))
#summary(tmctestsyn)
#confint(tmctestsyn)
sumtmctestsyn <- summary(tmctestsyn)
tabtmcsyn <- data.table(comp=names(sumtmctestsyn$test$coefficients),est=round(sumtmctestsyn$test$coefficients,digits=3),sigma=round(sumtmctestsyn$test$sigma,digits=3),tstat=round(sumtmctestsyn$test$tstat,digits=3),pval=round(sumtmctestsyn$test$pvalues,digits=3))
#write.csv(tabtmcsyn,"synopsistukeytests.csv")

####### create two separate datasets (one by synopsis to create separate plots, another containing all synopses to create overall plot) then merge with labels to denote difference 
masterall_trim_mainmodqp_all <- masterall_trim_mainmodqp
masterall_trim_mainmodqp_all$syn <- "All synopses"
masterall_syn_plot <- rbind(masterall_trim_mainmodqp_all,masterall_trim_mainmodqp)

###### create factor levels for synopses in alphabetical order
masterall_syn_plot$syn<-factor(masterall_syn_plot$syn,levels=sort(unique(masterall_syn_plot$syn)))

###### createfind and add summary statistics for each synopsis to dataset
for(i in 1:length(unique(masterall_syn_plot$syn))){
  masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],meddatediff:=as.integer(median(masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],datediff]))]
  masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],iqr_upr:=quantile(masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],datediff],probs=c(0.75))]
  masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],iqr_lwr:=quantile(masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],datediff],probs=c(0.25))]
  masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],meandatediff:=mean(masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],datediff])]
  masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],st_err:=sd(masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],datediff])/sqrt(length(masterall_syn_plot[syn==unique(masterall_syn_plot$syn)[i],datediff]))]
}

###### for purposes of plotting, summarise delays in bins - check these are still here  
unique(masterall_syn_plot$datediffac)

##### separate datasets into one for plot of whole database and plots by separate synopses
masterall_syn_plotnotall <- masterall_syn_plot[syn!="All synopses"]
masterall_syn_plotall <- masterall_syn_plot[syn=="All synopses"]
##### make plot of whole database appear at one point on x axis to aid plotting
masterall_syn_plotall[,syn:=as.numeric(syn)]
masterall_syn_plotall[,syn:=0]

##### run simple GLM with intercept only to find overall mean delay for whole database plot
mod1synall <- glm(datediff~1, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp)
meanallsyn <- exp(coef(mod1synall))
meanallsynupr <- exp(coef(mod1synall))+(qnorm(0.975)*unique(predict(mod1synall,type="response",se.fit=TRUE)[[2]]))
meanallsynlwr <- exp(coef(mod1synall))-(qnorm(0.975)*unique(predict(mod1synall,type="response",se.fit=TRUE)[[2]]))
meanallsyndat <- data.table(meanallsyn,meanallsynupr,meanallsynlwr)

##### plot the delay for whole database
plot1a<-ggplot() + 
  geom_density(data=masterall_syn_plotall,aes(x=datedifffac,y=..count..,group=syn),stat="binline",binwidth=1,show.legend=FALSE,fill="grey50")+
  scale_x_discrete(name="Publication delay (years)")+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies")+ggtitle("All synopses")+
  geom_segment(data=meanallsyndat,aes(x=meanallsyn+1,xend=meanallsyn+1,y=0,yend=1615),colour="red",linetype="solid",lwd=1.05)+
  geom_segment(data=meanallsyndat,aes(x=meanallsynupr+1,xend=meanallsynupr+1,y=0,yend=1615),colour="red",linetype="dashed",lwd=0.95)+
  geom_segment(data=meanallsyndat,aes(x=meanallsynlwr+1,xend=meanallsynlwr+1,y=0,yend=1615),colour="red",linetype="dashed",lwd=0.95)+
  theme(plot.title = element_text(size=25),strip.text = element_text(size=19),axis.title=element_text(size=23),axis.text=element_text(size=21),panel.grid.major = element_line(size=0.25,colour="grey70"))

##### find max delay and add data to dataset by synopsis
mxvals1 <- sapply(unique(masterall_syn_plotnotall$syn),function(y){  return(max(table(masterall_syn_plotnotall[syn==y,datedifffac])))})
for(i in 1:length(unique(mxvals1))){masterall_syn_plotnotall[syn==unique(masterall_syn_plotnotall$syn)[i],maxddatediff:=mxvals1[i]]}

##### run simple GLM for synopsis only without intercept to get mean delay for each synopsis for plotting
mod1syn <- glm(datediff~syn-1, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp)
preds=predict(mod1syn,data.frame(syn=gsub("syn","",names(mod1syn$coefficients))),type="response",se.fit=TRUE)
modelvalssyn <- data.table(syn=gsub("syn","",names(mod1syn$coefficients)),
                           mns=preds$fit,
                           upr=preds$fit+qnorm(0.975)*preds$se.fit,
                           lwr=preds$fit-qnorm(0.975)*preds$se.fit)

##### add predictions by synopsis to dataset for plotting
for(i in 1:nrow(modelvalssyn)){
  masterall_syn_plotnotall[syn==modelvalssyn$syn[i],modelmn:=modelvalssyn$mns[i]]
  masterall_syn_plotnotall[syn==modelvalssyn$syn[i],modelmnupr:=modelvalssyn$upr[i]]
  masterall_syn_plotnotall[syn==modelvalssyn$syn[i],modelmnlwr:=modelvalssyn$lwr[i]]
}

##### plot publication delay by separate synopses (+1 needed to adjust mean lines on x axis relative to histograms)
plot1b<-ggplot(masterall_syn_plotnotall) + 
  geom_density(aes(x=datedifffac,y=..count..,group=syn),fill="grey90",stat="binline",binwidth=1,show.legend=FALSE)+
  scale_x_discrete(name="Publication delay (years)")+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies",breaks=pretty_breaks(n=4))+
  geom_segment(aes(x=modelmn+1,xend=modelmn+1,y=0,yend=maxddatediff,group=syn),colour="red",linetype="solid",lwd=1.05)+
  geom_segment(aes(x=modelmnupr+1,xend=modelmnupr+1,y=0,yend=maxddatediff,group=syn),colour="red",linetype="dashed",lwd=0.95)+
  geom_segment(aes(x=modelmnlwr+1,xend=modelmnlwr+1,y=0,yend=maxddatediff,group=syn),colour="red",linetype="dashed",lwd=0.95)+
  facet_wrap(.~syn,ncol=3,scales="free_y")+
  theme(plot.title=element_text(size=25),strip.background = element_rect(fill="grey70"),strip.text = element_text(size=21),axis.title=element_text(size=23),axis.text.y=element_text(size=20),axis.text.x=element_text(size=17),panel.grid.major = element_line(size=0.25,colour="grey70"))

##### arrange overall plot with plot by separate synopses
ggarrange(plot1a,plot1b, labels=c("A","B"), nrow=2,ncol=1, heights=c(1,2), font.label=list(size=25))

#ggsave("delayplotsforpaper.svg",height=45,width=60,units="cm",device="svg")




#####################################################################################

###### obtain predicted values from simple GLM for publication date only for plotting
mod1time <- glm(datediff~pubdate, family=quasipoisson(link="log"),data=masterall_trim_mainmodqp)
preds=predict(mod1time,data.frame(pubdate=seq(1912,2020,1)),type="response",se.fit=TRUE)

##### add these to data table for plotting
modelvalstime <- data.table(pubdate=seq(1912,2020,1),mns=preds$fit,upr=preds$fit+qnorm(0.975)*preds$se.fit,lwr=preds$fit-qnorm(0.975)*preds$se.fit)

##### create plot of publication delay over time  - restrict to publication delay <=20 years
ggplot() + geom_hex(data=masterall_trim_mainmodqp,aes(y=datediff,x=pubdate),colour="black",bins=21)+
  # geom_ribbon(data=pred,aes(x=pubdate,ymin=lwr,ymax=upr),fill="grey",alpha=0.5)+
  geom_line(data=modelvalstime,aes(x=pubdate,y=mns),colour="red",size=1)+
  geom_line(data=modelvalstime,aes(x=pubdate,y=lwr),colour="red",linetype=2,size=1)+
  geom_line(data=modelvalstime,aes(x=pubdate,y=upr),colour="red",linetype=2,size=1)+
  scale_fill_gradient(low = "white", high = "black",breaks=seq(50,450,100))+
  theme_classic()+scale_y_continuous(name="Publication delay (years)",
                                     breaks=pretty_breaks(n=10),limits=c(-0.5,20.5),expand = c(0.01, 0.01))+
  scale_x_continuous(name="Publication year",pretty_breaks(n=5))+
  guides(fill=guide_colourbar(title="Number of studies"))+
  theme(aspect.ratio=1,axis.title=element_text(size=23),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20),legend.text = element_text(size=20),legend.title=element_text(size=20))

#ggsave("pubdelay_overtimeplots.svg",height=45,width=45,units="cm",device="svg")

##### create plot of publication delay over time - use all data including delay >20 years
ggplot() + geom_hex(data=masterall_trim_mainmodqp,aes(y=datediff,x=pubdate),colour="black",bins=101)+
  # geom_ribbon(data=pred,aes(x=pubdate,ymin=lwr,ymax=upr),fill="grey",alpha=0.5)+
  geom_line(data=modelvalstime,aes(x=pubdate,y=mns),colour="red",size=1)+
  geom_line(data=modelvalstime,aes(x=pubdate,y=lwr),colour="red",linetype=2,size=1)+
  geom_line(data=modelvalstime,aes(x=pubdate,y=upr),colour="red",linetype=2,size=1)+
  scale_fill_gradient(low = "white", high = "black",breaks=seq(50,450,100))+
  theme_classic()+scale_y_continuous(name="Publication delay (years)",
                                     breaks=pretty_breaks(n=10),limits=c(-0.5,100.5),expand = c(0.01, 0.01))+
  scale_x_continuous(name="Publication year",pretty_breaks(n=5))+
  guides(fill=guide_colourbar(title="Number of studies"))+
  theme(aspect.ratio=1,axis.title=element_text(size=23),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20),legend.text = element_text(size=20),legend.title=element_text(size=20))

#ggsave("pubdelay_overtimeplots_alldata.svg",height=45,width=45,units="cm",device="svg")


#######################################################################################################
###### taxanomic model

###### read in ampdat.csv, birddat.csv, and mamdat.csv (contains IUCN Red List names and statuses to match to our dataset on delays)
ampdat <- fread(choose.files())
birddat <- fread(choose.files())
mamdat <- fread(choose.files())

masterall_trim_taxa_all <- unique(masterall_trim[,list(pageid,lat,long,before,controlled,randomised,int,syn,stdate1,stdate2,datediff,binom,class,pubtype,pubdate)])

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

masterall_trim_taxa_allclean$iucncat2<-factor(masterall_trim_taxa_allclean$iucncat2,levels=c("DD", "LC","NT","VU","EN","CR","EW"))
#masterall_trim_taxa_allclean[is.na(iucncat2)&!is.na(iucncat)]

###too few studies on EW and DD
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="EW",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="DD",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="NT",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="VU",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="EN",list(pageid,datediff,before,controlled,randomised,pubdate)]))
nrow(unique(masterall_trim_taxa_allclean[iucncat2=="CR",list(pageid,datediff,before,controlled,randomised,pubdate)]))

##################### remove EW and DD category studies
masterall_trim_taxa_allclean <- unique(masterall_trim_taxa_allclean[iucncat2!="EW" & iucncat2!="DD",list(iucncat,iucncat2,syn,int,class,binom,pageid,datediff,before,controlled,randomised,pubtype,pubdate)])

######### only consider most threatened species in each study
uniqstuds <- unique(masterall_trim_taxa_allclean[,list(syn,int,pageid,datediff,before,controlled,randomised,class,pubtype,pubdate)])
masterall_trim_taxa_allcleanmost <- masterall_trim_taxa_allclean

for(i in 1:nrow(uniqstuds)){
  tokeep<-max(as.numeric(masterall_trim_taxa_allcleanmost[class==uniqstuds$class[i]&syn==uniqstuds$syn[i]&int==uniqstuds$int[i]&pageid==uniqstuds$pageid[i]&datediff==uniqstuds$datediff[i]& before==uniqstuds$before[i]&controlled==uniqstuds$controlled[i]&randomised==uniqstuds$randomised[i],iucncat2]))
  masterall_trim_taxa_allclean[class==uniqstuds$class[i]&syn==uniqstuds$syn[i]&int==uniqstuds$int[i]&pageid==uniqstuds$pageid[i]&datediff==uniqstuds$datediff[i]& before==uniqstuds$before[i]&controlled==uniqstuds$controlled[i]&randomised==uniqstuds$randomised[i],iucncat2:=tokeep]
}

masterall_trim_taxa_allcleanmostuniq <- masterall_trim_taxa_allcleanmost[,list(int,syn,pageid,datediff,before,controlled,randomised,iucncat2,class,pubtype,pubdate)]

######### order different factors
masterall_trim_taxa_allcleanmostuniq$syn = factor(masterall_trim_taxa_allcleanmostuniq$syn,levels=c("Management of Captive Animals",unique(masterall_trim_taxa_allcleanmostuniq[syn!="Management of Captive Animals",syn])))
masterall_trim_taxa_allcleanmostuniq$pubtype = factor(masterall_trim_taxa_allcleanmostuniq$pubtype,levels=sort(unique(masterall_trim_taxa_allcleanmostuniq$pubtype)))
masterall_trim_taxa_allcleanmostuniq$iucncat2 = factor(masterall_trim_taxa_allcleanmostuniq$iucncat2,levels=c("LC","NT","VU","EN","CR","EW"))

###### run a GLM for taxa including other variables from original model
mod2 <- glm(datediff~pubdate+syn+pubtype+iucncat2, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq)
summary(mod2)
#plot(mod1)

############# run mulitple comparisons for IUCN categories and output table
tmctestiucn <- glht(mod2, linfct = mcp(iucncat2 = "Tukey"))
#summary(tmctestiucn)
#confint(tmctestiucn)
sumtmctestiucn <- summary(tmctestiucn)
tabtmciucn <- data.table(comp=names(sumtmctestiucn$test$coefficients),est=round(sumtmctestiucn$test$coefficients,digits=3),sigma=round(sumtmctestiucn$test$sigma,digits=3),tstat=round(sumtmctestiucn$test$tstat,digits=3),pval=round(sumtmctestiucn$test$pvalues,digits=3))
#write.csv(tabtmciucn,"iucntukeytests.csv")


############# now run taxa specific analyses

###amphibian
mod2a <- glm(datediff~pubdate+pubtype+iucncat2, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq[syn=="Amphibian Conservation"])
summary(mod2)
#plot(mod1a)

tmctestiucna <- glht(mod2a, linfct = mcp(iucncat2 = "Tukey"))
#summary(tmctestiucna)
#confint(tmctestiucna)
sumtmctestiucna <- summary(tmctestiucna)
tabtmciucna <- data.table(comp=names(sumtmctestiucna$test$coefficients),est=round(sumtmctestiucna$test$coefficients,digits=3),sigma=round(sumtmctestiucna$test$sigma,digits=3),tstat=round(sumtmctestiucna$test$tstat,digits=3),pval=round(sumtmctestiucna$test$pvalues,digits=3))
#write.csv(tabtmciucna,"iucntukeytestsamp.csv")

###bird
mod2b <- glm(datediff~pubdate+pubtype+iucncat2, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq[syn=="Bird Conservation"])
summary(mod2b)
#plot(mod1b)

tmctestiucnb <- glht(mod2b, linfct = mcp(iucncat2 = "Tukey"))
#summary(tmctestiucnb)
#confint(tmctestiucnb)
sumtmctestiucnb <- summary(tmctestiucnb)
tabtmciucnb <- data.table(comp=names(sumtmctestiucnb$test$coefficients),est=round(sumtmctestiucnb$test$coefficients,digits=3),sigma=round(sumtmctestiucnb$test$sigma,digits=3),tstat=round(sumtmctestiucnb$test$tstat,digits=3),pval=round(sumtmctestiucnb$test$pvalues,digits=3))
#write.csv(tabtmciucnb,"iucntukeytestsbird.csv")

###mammal
mod2m <- glm(datediff~pubdate+pubtype+iucncat2, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq[syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"])
summary(mod2m)
#plot(mod1m)

tmctestiucnm <- glht(mod2m, linfct = mcp(iucncat2 = "Tukey"))
#summary(tmctestiucnm)
#confint(tmctestiucnm)
sumtmctestiucnm <- summary(tmctestiucnm)
tabtmciucn <- data.table(comp=names(sumtmctestiucnm$test$coefficients),est=round(sumtmctestiucnm$test$coefficients,digits=3),sigma=round(sumtmctestiucnm$test$sigma,digits=3),tstat=round(sumtmctestiucnm$test$tstat,digits=3),pval=round(sumtmctestiucnm$test$pvalues,digits=3))
#write.csv(tabtmciucn,"iucntukeytestsmam.csv")

##################### now plot different taxanomic analyses

############# make sure histogram categories are created as in main analyses
masterall_trim_taxa_allcleanmostuniq[datediff<=10, datedifffac := as.character(datediff)]
masterall_trim_taxa_allcleanmostuniq[datediff>10 & datediff<=20, datedifffac := "11-20"]
masterall_trim_taxa_allcleanmostuniq[datediff>20, datedifffac := "20+"]
masterall_trim_taxa_allcleanmostuniq$datedifffac<-factor(masterall_trim_taxa_allcleanmostuniq$datedifffac,as.character(c(sort(as.numeric(unique(masterall_trim_taxa_allcleanmostuniq[datediff<=10,datediff]))),"11-20","20+")))

############ run simple GLM for IUCN category (without intercept) to get mean values for separate taxa plots
mod2aa <- glm(datediff~iucncat2-1, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq[syn=="Amphibian Conservation"])
preds=predict(mod2aa,data.frame(iucncat2=gsub("iucncat2","",names(mod2aa$coefficients))),type="response",se.fit=TRUE)
modelvalsaa <- data.table(iucncat2=gsub("iucncat2","",names(mod2aa$coefficients)),mns=preds$fit,upr=preds$fit+qnorm(0.975)*preds$se.fit,lwr=preds$fit-qnorm(0.975)*preds$se.fit)

############## add these to dataset
for(i in 1:nrow(modelvalsaa)){
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsaa$iucncat2[i]&syn=="Amphibian Conservation",modelmn:=modelvalsaa$mns[i]]
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsaa$iucncat2[i]&syn=="Amphibian Conservation",modelmnupr:=modelvalsaa$upr[i]]
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsaa$iucncat2[i]&syn=="Amphibian Conservation",modelmnlwr:=modelvalsaa$lwr[i]]
}

######## repeat for birds
mod2bb <- glm(datediff~iucncat2-1, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq[syn=="Bird Conservation"])
preds=predict(mod2bb,data.frame(iucncat2=gsub("iucncat2","",names(mod2bb$coefficients))),type="response",se.fit=TRUE)
modelvalsbb <- data.table(iucncat2=gsub("iucncat2","",names(mod2bb$coefficients)),
                          mns=preds$fit,
                          upr=preds$fit+qnorm(0.975)*preds$se.fit,
                          lwr=preds$fit-qnorm(0.975)*preds$se.fit)

for(i in 1:nrow(modelvalsbb)){
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsbb$iucncat2[i]&syn=="Bird Conservation",modelmn:=modelvalsbb$mns[i]]
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsbb$iucncat2[i]&syn=="Bird Conservation",modelmnupr:=modelvalsbb$upr[i]]
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsbb$iucncat2[i]&syn=="Bird Conservation",modelmnlwr:=modelvalsbb$lwr[i]]
}

######## repeat for mammals
mod2mm <- glm(datediff~iucncat2-1, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniq[syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"])
preds=predict(mod2mm,data.frame(iucncat2=gsub("iucncat2","",names(mod2mm$coefficients))),type="response",se.fit=TRUE)
modelvalsmm <- data.table(iucncat2=gsub("iucncat2","",names(mod2mm$coefficients)),
                          mns=preds$fit,
                          upr=preds$fit+qnorm(0.975)*preds$se.fit,
                          lwr=preds$fit-qnorm(0.975)*preds$se.fit)

for(i in 1:nrow(modelvalsmm)){
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsmm$iucncat2[i]&(syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"),modelmn:=modelvalsmm$mns[i]]
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsmm$iucncat2[i]&(syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"),modelmnupr:=modelvalsmm$upr[i]]
  masterall_trim_taxa_allcleanmostuniq[iucncat2==modelvalsmm$iucncat2[i]&(syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"),modelmnlwr:=modelvalsmm$lwr[i]]
}

############ plot data for amphibians
plotamp<-ggplot(masterall_trim_taxa_allcleanmostuniq[syn=="Amphibian Conservation"]) + 
  geom_density(aes(x=datedifffac,y=..count..,group=iucncat2,fill=iucncat2),stat="binline",binwidth=1,show.legend=FALSE)+
  scale_x_discrete(name="Publication delay (years)",limits=levels(masterall_trim_taxa_allcleanmostuniq$datedifffac))+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies",breaks=pretty_breaks(n=4))+
  geom_vline(aes(xintercept=modelmn+1,group=iucncat2),colour="black",linetype="solid",lwd=1.05)+
  geom_vline(aes(xintercept=modelmnupr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  geom_vline(aes(xintercept=modelmnlwr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  facet_wrap(.~iucncat2,scales="fixed",ncol=1,strip.position = "right")+
  theme(strip.text = element_text(size=18),axis.title=element_text(size=18),
        legend.position = "none",plot.title=element_text(size=20))+  
  scale_fill_manual(values=c((brewer.pal(n=5, name="OrRd"))[1:5]))+ggtitle("Amphibians")

############ plot data for birds
plotbird<-ggplot(masterall_trim_taxa_allcleanmostuniq[syn=="Bird Conservation"]) + 
  geom_density(aes(x=datedifffac,y=..count..,group=iucncat2,fill=iucncat2),stat="binline",binwidth=1,show.legend=FALSE)+
  scale_x_discrete(name="Publication delay (years)")+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies",breaks=pretty_breaks(n=4))+
  geom_vline(aes(xintercept=modelmn+1,group=iucncat2),colour="black",linetype="solid",lwd=1.05)+
  geom_vline(aes(xintercept=modelmnupr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  geom_vline(aes(xintercept=modelmnlwr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  facet_wrap(.~iucncat2,scales="fixed",ncol=1,strip.position = "right")+
  theme(strip.text = element_text(size=18),axis.title=element_text(size=18),
        legend.position = "none",plot.title=element_text(size=20))+  
  scale_fill_manual(values=c((brewer.pal(n=5, name="OrRd"))[1:5]))+ggtitle("Birds")

############ plot data for mammals
plotmam<-ggplot(masterall_trim_taxa_allcleanmostuniq[syn=="Terrestrial Mammal Conservation"|syn=="Bat Conservation"|syn=="Primate Conservation"]) + 
  geom_density(aes(x=datedifffac,y=..count..,group=iucncat2,fill=iucncat2),stat="binline",binwidth=1,show.legend=FALSE)+
  scale_x_discrete(name="Publication delay (years)")+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies",breaks=pretty_breaks(n=4))+
  geom_vline(aes(xintercept=modelmn+1,group=iucncat2),colour="black",linetype="solid",lwd=1.05)+
  geom_vline(aes(xintercept=modelmnupr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  geom_vline(aes(xintercept=modelmnlwr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  facet_wrap(.~iucncat2,scales="fixed",ncol=1,strip.position = "right")+
  theme(strip.text = element_text(size=18),axis.title=element_text(size=18),
        legend.position = "none",plot.title=element_text(size=20))+  
  scale_fill_manual(values=c((brewer.pal(n=5, name="OrRd"))[1:5]))+ggtitle("Mammals")

###### combine plots
ggarrange(plotamp,plotbird,plotmam,nrow=1)
#ggsave("pubdelay_iucnplots.svg",height=45,width=45,units="cm",device="svg")




#### do the same but plotting all taxa together as a whole
masterall_trim_taxa_allcleanmostuniqall <- masterall_trim_taxa_allcleanmostuniq
mod2d <- glm(datediff~iucncat2-1, family=quasipoisson(link="log"),data=masterall_trim_taxa_allcleanmostuniqall)
preds=predict(mod2d,data.frame(iucncat2=gsub("iucncat2","",names(mod2d$coefficients))),type="response",se.fit=TRUE)
modelvalsd <- data.table(iucncat2=gsub("iucncat2","",names(mod2d$coefficients)),mns=preds$fit,upr=preds$fit+qnorm(0.975)*preds$se.fit,lwr=preds$fit-qnorm(0.975)*preds$se.fit)

for(i in 1:nrow(modelvalsd)){
  masterall_trim_taxa_allcleanmostuniqall[iucncat2==modelvalsd$iucncat2[i],modelmn:=modelvalsd$mns[i]]
  masterall_trim_taxa_allcleanmostuniqall[iucncat2==modelvalsd$iucncat2[i],modelmnupr:=modelvalsd$upr[i]]
  masterall_trim_taxa_allcleanmostuniqall[iucncat2==modelvalsd$iucncat2[i],modelmnlwr:=modelvalsd$lwr[i]]
}

ggplot(masterall_trim_taxa_allcleanmostuniqall) + 
  geom_density(aes(x=datedifffac,y=..count..,group=iucncat2,fill=iucncat2),stat="binline",binwidth=1,show.legend=FALSE)+
  scale_x_discrete(name="Publication delay (years)")+theme_ridges(center_axis_labels = TRUE)+
  scale_y_continuous(name="Number of studies",breaks=pretty_breaks(n=4))+
  geom_vline(aes(xintercept=modelmn+1,group=iucncat2),colour="black",linetype="solid",lwd=1.05)+
  geom_vline(aes(xintercept=modelmnupr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  geom_vline(aes(xintercept=modelmnlwr+1,group=iucncat2),colour="black",linetype="dashed",lwd=0.95)+
  facet_wrap(.~iucncat2,scales="fixed",ncol=1,strip.position = "right")+
  theme(strip.text = element_text(size=18),axis.title=element_text(size=18),
        legend.position = "none",plot.title=element_text(size=20))+  
  scale_fill_manual(values=c((brewer.pal(n=5, name="OrRd"))[1:5]))+ggtitle("Amphibians, Birds, and Mammals")

#ggsave("pubdelay_iucnplots_all.svg",height=45,width=45,units="cm",device="svg")
