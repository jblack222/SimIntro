################################################################################
################################################################################
## Proj2_sims1_prog.r
##
## Some functions for running simulations to calculate heterozygosity 
## from simulations using the raw genetic data.
##
## Authors: Jakub Stoklosa (School of Mathematics and Statistics, 
## The University of New South Wales, NSW 2052, Australia).
##
## Please report any problems/suggestions to Jakub Stoklosa at:
## j.stoklosa@unsw.unimelb.edu.au
##
## This program is meant to be used for non-commercial purposes only.
################################################################################
################################################################################

## Simulation functions.

het_function<-function(data)
  {
  n_o<-nrow(data);
  H_o<-sum(data[,1]!=data[,2])/n_o;
  allele.no<-min(data):max(data);
  match.vec<-c();
  
  for(i in 1:length(allele.no))
    {
    match.vec1<-c(((data[,1]==allele.no[i])+0)+((data[,2]==allele.no[i])+0));
    match.vec<-cbind(match.vec,match.vec1);
    }
  
  q<-apply(match.vec,2,sum)/(n_o*2);
  q.sq<-q^2;
  sum_q.sq<-sum(q^2);
  H_e<-1-sum_q.sq;
  H_e_corr<-H_e*((2*n_o)/(2*n_o-1));
  Na<-sum(q>0);
  c(n_o,H_o,H_e,H_e_corr,Na);
  }

het_simulation1<-function(data1,N.sim,n.sim,p.sim,T.sim)
  {
  data.track<-list();
  het.track<-list();
  het.mean.track<-c();
  
  data.old<-data1;
  for(t in 1:T.sim)
    {
    data.new<-matrix(NA,nrow=N.sim,ncol=n.sim*p.sim);
    for(i in 1:N.sim)
      {
      for(j in 1:n.sim)
        {
        k<-(p.sim*j)-1;
        sample.locus<-c(data.old[,c(k:(k+1))])[!is.na(c(data.old[,c(k:(k+1))]))];
        data.new[i,k:(k+1)]<-sample(sample.locus,size=p.sim,replace=T);
        }
      }
    
    het.mat<-c();
    for(j in 1:n.sim)
      {
      k<-(p.sim*j)-1;
      het.mat<-rbind(het.mat,het_function(data.new[,c(k:(k+1))]));
      }
    colnames(het.mat)<-c("n_o","H_o","H_e","H_e_corr","Na");
    het.mean.track<-rbind(het.mean.track,(apply(het.mat,2,mean)[-1]));
    
    if(t==10 | t==20 | t==50 | t==100 | t==250 | t==500 | t==1000){data.track<-c(data.track,list(data.new));}
    if(t==10 | t==20 | t==50 | t==100 | t==250 | t==500 | t==1000){het.track<-c(het.track,list(het.mat));}
    
    data.old<-data.new;
    }
  names(data.track)[[1]]<-"10th Generation";   names(het.track)[[1]]<-"10th Generation";
  names(data.track)[[2]]<-"20th Generation";   names(het.track)[[2]]<-"20th Generation";
  names(data.track)[[3]]<-"50th Generation";   names(het.track)[[3]]<-"50th Generation"; 
  names(data.track)[[4]]<-"100th Generation";  names(het.track)[[4]]<-"100th Generation";
  names(data.track)[[5]]<-"250th Generation";  names(het.track)[[5]]<-"250th Generation";
  names(data.track)[[6]]<-"500th Generation";  names(het.track)[[6]]<-"500th Generation";
  names(data.track)[[7]]<-"1000th Generation"; names(het.track)[[7]]<-"1000th Generation";
  
  list(data=data.track,het=het.track,het.mean.track=het.mean.track);
  }

het_simulation2<-function(data1,N.sim,n.sim,p.sim,T.sim,Rep.sim)
  {
  data.old<-list();
  het.mean.track<-list();
  het.conf.info<-c();
  
  for(u in 1:Rep.sim){data.old<-c(data.old,list(data1))}
  
  for(t in 1:T.sim)
    {
    het.mean.track1<-c();
    data.new.track<-list();
    
    for(u in 1:Rep.sim)
      {
      data.new<-matrix(NA,nrow=N.sim,ncol=n.sim*p.sim);
      for(i in 1:N.sim)
        {
        for(j in 1:n.sim)
          {
          k<-(p.sim*j)-1;
          sample.locus<-c(data.old[[u]][,c(k:(k+1))])[!is.na(c(data.old[[u]][,c(k:(k+1))]))];
          data.new[i,k:(k+1)]<-sample(sample.locus,size=p.sim,replace=T);
          }
        }
      
      het.mat<-c();
      for(j in 1:n.sim)
        {
        k<-(p.sim*j)-1;
        het.mat<-rbind(het.mat,het_function(data.new[,c(k:(k+1))]));
        }
      colnames(het.mat)<-c("n_o","H_o","H_e","H_e_corr","Na");
      het.mean.track1<-rbind(het.mean.track1,(apply(het.mat,2,mean)[-1]));  
      data.new.track<-c(data.new.track,list(data.new));
      }
    
    data.old<-data.new.track;
    het.conf.info1<-c(mean(het.mean.track1[,2]),mean(het.mean.track1[,2])-1.96*sd(het.mean.track1[,2]),mean(het.mean.track1[,2])+1.96*sd(het.mean.track1[,2]));
    het.conf.info<-rbind(het.conf.info,het.conf.info1);
    }
  
  list(het.conf.info=het.conf.info);
  }

het_simulation3<-function(data1,N.sim,n.sim,p.sim,T.sim,Rep.sim)
  {
  data.old<-list();
  het.mean.track<-list();
  het.conf.info<-c();
  Na.conf.info<-c();
  
  for(u in 1:Rep.sim){data.old<-c(data.old,list(data1))}
  
  for(t in 1:T.sim)
    {
    het.mean.track1<-c();
    data.new.track<-list();
    
    for(u in 1:Rep.sim)
      {
      data.new<-matrix(NA,nrow=N.sim,ncol=n.sim*p.sim);
      for(i in 1:N.sim)
        {
        for(j in 1:n.sim)
          {
          k<-(p.sim*j)-1;
          sample.locus<-c(data.old[[u]][,c(k:(k+1))])[!is.na(c(data.old[[u]][,c(k:(k+1))]))];
          data.new[i,k:(k+1)]<-sample(sample.locus,size=p.sim,replace=T);
          }
        }
    
      het.mat<-c();
      for(j in 1:n.sim)
        {
        k<-(p.sim*j)-1;
        het.mat<-rbind(het.mat,het_function(data.new[,c(k:(k+1))]));
        }
      colnames(het.mat)<-c("n_o","H_o","H_e","H_e_corr","Na");
      het.mean.track1<-rbind(het.mean.track1,(apply(het.mat,2,mean)[-1]));  
      data.new.track<-c(data.new.track,list(data.new));
      }

    data.old<-data.new.track;
    het.conf.info1<-c(mean(het.mean.track1[,3]),mean(het.mean.track1[,3])-1.96*sd(het.mean.track1[,3]),mean(het.mean.track1[,3])+1.96*sd(het.mean.track1[,3]));
    het.conf.info<-rbind(het.conf.info,het.conf.info1);
    
    Na.conf.info1<-c(mean(het.mean.track1[,4]),mean(het.mean.track1[,4])-1.96*sd(het.mean.track1[,4]),mean(het.mean.track1[,4])+1.96*sd(het.mean.track1[,4]));
    Na.conf.info<-rbind(Na.conf.info,Na.conf.info1);
    }

  list(het.conf.info=het.conf.info,Na.conf.info=Na.conf.info);
  }

## Same function as above but only reports output for one specified time gen.

simulation_pan<-function(data1,N.sim,n.sim,p.sim,T.sim)
  {
  data.track<-list();
  data.old<-data1;
  
  for(t in 1:T.sim)
    {
    data.new<-matrix(NA,nrow=N.sim[t+1],ncol=(n.sim*p.sim));
    
    for(i in 1:N.sim[t+1])
      {
      for(j in 1:n.sim)
        {
        k<-(p.sim*j)-1;
        sample.locus<-c(data.old[,c(k:(k+1))])[!is.na(c(data.old[,c(k:(k+1))]))];
        data.new[i,k:(k+1)]<-sample(sample.locus,size=p.sim,replace=T);
        }
      }
    if(t==T.sim){data.track<-c(data.track,list(data.new));}
    data.old<-data.new;
    }
  data.track;
  }
