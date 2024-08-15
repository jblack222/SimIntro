##################################################################################
##################################################################################
## SimIntro.R
##
## R-code to calculate heterozygosity from simulations using the raw genetic data.
##
## Adapted from microsatellite code by Jakub Stoklosa (School of Mathematics and Statistics, The University of New South Wales, NSW 2052, Australia).
## 
## This program is meant to be used for non-commercial purposes only.
##################################################################################
##################################################################################

##############################################################################

############# NOTE FROM JOHN ###########
# When going from vcf -> genind, vcfR will remove loci with 100% missing data
# Critical to ensure that any sites are located and removed for ALL datasets

rm(list=ls());

setwd("~/path/to/simulated/introduction");
source("~/path/to/sims1_prog.R");
library(psych);
library(Rlab);
library(poppr);
library(vcfR); 

pops4genalex <- read.delim("popmap.txt", he_pop1r = F, sep = "\t")

all_vcf<-read.vcfR("all_samples.vcf")
all_genind<-vcfR2genind(all_vcf)
all_genalex<-genind2genalex(all_genind, pop = pops4genalex[,2], overwrite = TRUE, filename = "all_samples.genalex.csv")

dataALL_0<-as.matrix(read.csv("SBB.all.genalex.csv",blank.lines.skip=TRUE,colClasses=NA, header = FALSE))

dataALL_1<-(dataALL_0[-1,-1])
dataALL_2<-(dataALL_1[-1,-1])
dataALL_3<-(dataALL_2[-1,])

dataALL_F <- cbind(pops4genalex, dataALL_3)

data_pop1 <- subset(dataALL_F, V2 == 'Melbourne')
data_pop1 <- data_pop1[,-1:-2]
data_pop2 <- subset(dataALL_F, V2 == 'Sydney')
data_pop2 <- data_pop2[,-1:-2]
data_pop3 <- subset(dataALL_F, V2 == 'Darwin')
data_pop3 <- data_pop3[,-1:-2]

N_pop1 <- nrow(data_pop1)
N_pop2 <- nrow(data_pop2)
N_pop3 <- nrow(data_pop3)
# No. of individuals in each population.

n<-ncol(dataALL_3)/2      
# No. of Locus for each population

data_pop1<-as.matrix(data_pop1)
data_pop1<-matrix(as.numeric(c(data_pop1)),nrow=N_pop1,ncol=2*n,byrow=F);
data_pop2<-as.matrix(data_pop2)
data_pop2<-matrix(as.numeric(c(data_pop2)),nrow=N_pop2,ncol=2*n,byrow=F);
data_pop3<-as.matrix(data_pop3)
data_pop3<-matrix(as.numeric(c(data_pop3)),nrow=N_pop3,ncol=2*n,byrow=F);

data_pop1[data_pop1==0]=2 
data_pop2[data_pop2==0]=2 
data_pop3[data_pop3==0]=2 

## CHOOSE YOUR POPULATION TO COMBINE
data1<-data_pop1
data2<-data_pop2

## CHOOSE HOW MANY POPULATIONS TO MIX  
no.data.sets<-2;

## CHOOSE HOW MANY SAMPLES FROM EACH POPULATION
no.ind.sample<-c(20,4);

## LIST ALL DATASETS TO BE USED eg data1, data2, data99
all.data<-list(data1,data2); 


new.raw.dat<-c();

for(i in 1:no.data.sets)
  {
  data.row<-c(1:nrow(all.data[[i]]));
  rand.chosen<-sample(data.row,no.ind.sample[i],replace=T);
  new.chosen<-all.data[[i]][rand.chosen,];
  new.raw.dat<-rbind(new.raw.dat,new.chosen);
  }  
  
## SET YOUR GENERATION TIME AND GROWTH RATE
T.sim.pan<-2; 
# No. of the panmictic generation time.
growth.rate<-2.0;  
# Growth rate for each t-th gen. time.


N.sim.vec<-c(nrow(new.raw.dat));

for(i in 1:T.sim.pan){N.sim.vec<-c(N.sim.vec,round(N.sim.vec[i]*growth.rate));}

N.sim.vec;

data_pan<-simulation_pan(new.raw.dat,N.sim.vec,n,2,T.sim.pan)[[1]];

## NAME YOUR PRE-MIXING FILE
write.csv(data_pan,"20pop_1.4pop_2.pan.csv"); 
# this outputs the mic data after panmixia and pop growth as a csv


hete.conf_big<-c();
Na.conf_big<-c();

## SET YOUR MAXIMUM POPULATION SIZE
pop.size_vec<-c(300);  
# Change the pop size here.
lnN<-length(pop.size_vec);

## SET PARAMETERS BELOW
for(l in 1:lnN)
  {
  N.sim<-pop.size_vec[l];  
# Population size at each generation.
  n.sim<-7333;  
# No. of locus for each individual at each generation.
  p.sim<-2;  
# No. of alleles per locus for each individual at each generation.
  T.sim<-50;  
# No. of generations.
  Rep.sim<-5; 
# No. of repeated sims. for each gen.
  
  start<-Sys.time();
  data.sims4<-het_simulation3(data_pan,N.sim,n.sim,p.sim,T.sim,Rep.sim);
  end<-Sys.time(); end-start;
  print("Currently generating sims. for a pop. size of:")
  print(pop.size_vec[l+1]);
  
  hete.conf.info<-data.sims4$het.conf.info;
  Na.conf.info<-data.sims4$Na.conf.info;
  
  hete.conf_big<-cbind(hete.conf_big,hete.conf.info[,1]);
  Na.conf_big<-cbind(Na.conf_big,Na.conf.info[,1]);
  }

par(mfrow=c(1,1),las=1);

matplot(1:T.sim,hete.conf_big,xlab="Gen. time",ylab="",type="l",main="H_e",
        col=c(1:lnN),lty=c(1:lnN),lwd=2,cex.lab=1.5,font.lab=1,cex.axis=1.2,cex.main=1.5,las=1);

legend("bottomleft",legend=pop.size_vec,lwd=2,col=c(1:lnN),lty=c(1:lnN),merge=TRUE,bty="n",cex=1.5);

x11(); 
# This makes two plots in the R-console, they are stacked on one another. 

par(mfrow=c(1,1),las=1);

matplot(1:T.sim,Na.conf_big,xlab="Gen. time",ylab="",type="l",main="Na",
        col=c(1:length(pop.size_vec)),lty=c(1:length(pop.size_vec)),
        lwd=2,cex.lab=1.5,font.lab=1,cex.axis=1.2,cex.main=1.5,las=1);

write.csv(Na.conf_big,"Na_file.pop_1.pop_2.csv"); 
# this outputs the NA to a csv for the various pop sizes
write.csv(hete.conf_big,"He_file.pop_1.pop_2.csv"); 
# this outputs HE to a csv for the various pop sizes
