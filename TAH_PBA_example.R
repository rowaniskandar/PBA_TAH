# PBA for TAH model
# developed by rowan iskandar
# modified 18012022
# this file is for both par and non-par version
# use TAH_functions.R


currpath <- dirname(rstudioapi::callFun("getActiveDocumentContext")$path)  
setwd(currpath)
library(rjags) #install https://sourceforge.net/projects/mcmc-jags/ first
library(R2OpenBUGS)
library(R2jags)
library(dplyr)
library(reshape2)
#library(ggplot2)
library(stringr)
#library(INLA, lib.loc = "C:/R_libraries")
library(flexsurv)
library(zip)
library(plyr)
library(BCEA)
#library(ggplot2)
library(nloptr) #optimization
library(doParallel)  
library(future)
library(future.apply)
library(purrr)
library(truncnorm)
# library(gmailr)
library(support)
# library(openxlsx)
# library(beepr)
# library(RColorBrewer)
# library(scales)

#source("TAH_functions_new_model_v7.R") #load 

################################################################################
#functions
########################################################################################
#pbox functions given different minimal data
########################################################################################

pbox_mean_inv <- function(u,type,min,max,mean){
  crit <- (max-mean)/(max-min)
  #print(crit)
  if (type=="upper"){
    #need to double check
    if(u==1){
      x=max
      #x=mean
    }
    else if(crit<u && u<1){
      x=(u*max-max+mean)/u
    }
    else if(u==crit){
      x=min
    }
    else if(u<crit){
      x=min
    }
  }
  else{
    if (u==1){
      x=max
    }
    else if (0<u && u< crit){
      x=(u*min-mean)/(u-1)
    }
    else if(u<=0){
      x=(u*min-mean)/(u-1) #not sure, min <= x < mean will yield u=0
    }
    else if(crit <=u && u<1){ #not sure about the condition
      x=max
    }
  }
  return(x)
}

#inverse for minmax
#newer version may 19 2021
pbox_minmax_inv <-function(u,type,min,max){
  if (type=="upper"){
    if(u<1){
      x=min
    }
    else{x=max}
  }
  else{
    if(u==0){
      x=min}
    else{
      x=max
    }
  }
  return (x)
}

#inverse for median
pbox_median_inv <-function(u,type,min,max,median){
  if (type=="upper"){ 
    if (u<0.5){
      x=min
    }
    else if(u >= 0.5 && u<1){
      x=median
    }
    else{ #different from ver 1
      x=max
    }
  }
  else{
    if (u<0.5){
      x=median
    }
    else if(0.5 <= u){
      x=max
    }
    else{
      x=min
    }
  }
  return (x)
}

pbox_std_inv <- function(u,type, min, max, mean, std){
  F1 <- (std^2)/((mean-min)^2+std^2)
  #print(F1)
  F2 <- (max-mean)^2/((max-mean)^2 + std^2)
  #print(F2)
  if (type=="upper"){
    if (u < F1 ){
      x=min
    }
    else if ( F1 <=u  && u < F2 ){
      x=mean-std*sqrt((1-u)/u)
    }
    else if (u >= F2){ 
      num <- (max-mean)*(mean-min)-std^2
      denom <- u*(max-min)+mean-max
      x <- max - num/denom
    }
  }
  else {
    if (u <= F1 ){
      num <- (max-mean)*(mean-min)-std^2
      denom <- max-mean-u*(max-min)
      x<- min + num/denom
    }
    else if ( F1 < u  && u <= F2 ){
      x=mean+std*sqrt(u/(1-u))
    }
    else if (u> F2) {
      x=max
    }
  }
  return(x)
}

pbox_median_std_inv <-function(u,type,min,max, median, mean, std){
  F1 <- (std^2)/((mean-min)^2+std^2)
  F2 <- (max-mean)^2/((max-mean)^2 + std^2)
  if (type=="upper"){ #
    if(u < F1){
      x=min
    }
    else if (F1 <=u && u< 0.5){
      x=mean-std*sqrt((1-u)/u)
    }
    else if (0.5 <= u && u< F2){
      x=median
    }
    else if ( F2 <= u && u <1){
      num <- (max-mean)*(mean-min)-std^2
      denom <- u*(max-min)+mean-max
      xtemp <- max - num/denom
      if (xtemp < median){
        x=median
      }
      else{x<-xtemp}
    }
    else if(u==1){
      x=max
    }
  }
  else{
    if(u<F1){
      num <- (max-mean)*(mean-min)-std^2
      denom <- max-mean-u*(max-min)
      xtemp<- min + num/denom
      if (xtemp >median){
        x<-median
      } else {
        x <-xtemp
      }
    }
    else if (F1 <= u && u <0.5){
      x<-median
    }
    else if(0.5 <= u && u <F2){
      x=mean+std*sqrt(u/(1-u))
    }
    else if (F2<=u){
      x=max
    }
  }
  return(x)
}

########################################################################################
#TAH CEA model for the case study
########################################################################################

markov.model <- function(tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync,
                         short.term.mortality.crmt, long.term.mortality.crmt,
                         short.term.mortality.sync, long.term.mortality.sync,
                         tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke,
                         p.dis.stroke.crmt, p.dis.stroke.sync,
                         TH.pre=Time.horiz.pre,TH.post=Time.horiz.post){
  
  MM.pre.crmt <- array(0,dim= c(3, 3, TH.pre-1))
  MM.pre.crmt[1,1,]<- (1 - tp.alive2death_HTx.crmt)
  MM.pre.crmt[1,3,]<- tp.alive2death_HTx.crmt*c(short.term.mortality.crmt,rep(long.term.mortality.crmt, TH.pre-4))
  MM.pre.crmt[1,2,]<- tp.alive2death_HTx.crmt*(1-c(short.term.mortality.crmt,rep(long.term.mortality.crmt, TH.pre-4)))
  MM.pre.crmt[2,2,]<- 1
  MM.pre.crmt[3,3,]<- 1
  
  
  MM.pre.sync <- array(0,dim= c(3, 3, TH.pre-1))
  MM.pre.sync[1,1,]<- (1 - tp.alive2death_HTx.sync)
  MM.pre.sync[1,3,]<- tp.alive2death_HTx.sync*c(short.term.mortality.sync,rep(long.term.mortality.sync, TH.pre-4))
  MM.pre.sync[1,2,]<- tp.alive2death_HTx.sync*(1-c(short.term.mortality.sync,rep(long.term.mortality.sync, TH.pre-4)))
  MM.pre.sync[2,2,]<- 1
  MM.pre.sync[3,3,]<- 1
  
  trace.pre <- array(0, dim = c(TH.pre, 3,  2))
  trace.pre[1, 1, ] <- cohort
  
  
  for(i in 2:(TH.pre)){
    trace.pre[i, , 1] <- trace.pre[i-1,, 1] %*% MM.pre.crmt[,,i-1]
    trace.pre[i, , 2] <- trace.pre[i-1,, 2] %*% MM.pre.sync[,,i-1]
  }
  
  
  #Post-transplantation 
  
  #calculate probability of transplantation with novel TAH and Syncardia
  p.HTx.crmt<-vector()
  for(i in 1:TH.pre-1){p.HTx.crmt[i]<-MM.pre.crmt[1,2,i]*prod(MM.pre.crmt[1,1,1:i-1])}
  p.HTx.crmt<-sum(p.HTx.crmt)
  
  p.HTx.sync<-vector()
  for(i in 1:TH.pre-1){p.HTx.sync[i]<-MM.pre.sync[1,2,i]*prod(MM.pre.sync[1,1,1:i-1])}
  p.HTx.sync<-sum(p.HTx.sync)
  
  
  #post transplantation transition probability matrix (state 1: alive free from stroke, state 2: alive with disabling stroke, state 3: dead)
  MM.post <- array(0,dim= c(3, 3, TH.post-1))
  MM.post[1,1,]<- (1 - tp.HTx2death.stroke.free)
  MM.post[1,3,]<- tp.HTx2death.stroke.free
  MM.post[2,2,]<- 1-tp.HTx2death.dis.stroke
  MM.post[2,3,]<- tp.HTx2death.dis.stroke
  MM.post[3,3,]<- 1
  
  
  
  trace.post <- array(0, dim = c(TH.post, 3,  2))
  trace.post[1, , 1] <- c(p.HTx.crmt*(1-p.dis.stroke.crmt),p.HTx.crmt*p.dis.stroke.crmt,1-p.HTx.crmt)
  trace.post[1, , 2] <- c(p.HTx.sync*(1-p.dis.stroke.sync),p.HTx.sync*p.dis.stroke.sync,1-p.HTx.sync)
  
  
  for(i in 2:(TH.post)){
    trace.post[i, , 1] <- trace.post[i-1,, 1] %*% MM.post[,,i-1]
    trace.post[i, , 2] <- trace.post[i-1,, 2] %*% MM.post[,,i-1]
  }
  
  trace<-list(trace.pre=trace.pre,trace.post=trace.post,p.HTx.crmt=p.HTx.crmt,p.HTx.sync=p.HTx.sync)
  return(trace)
  
}


costs <- function(tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync,
                  short.term.mortality.crmt, long.term.mortality.crmt,
                  short.term.mortality.sync, long.term.mortality.sync,
                  tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke,
                  p.dis.stroke.crmt, p.dis.stroke.sync,
                  # pre-transplantation
                  LOS.dis.crmt , LOS.dis.sync,
                  p.surgr , p.neurd , p.majin , p.majdm,
                  p.dis,
                  c.norm.ward,c.neurd, c.majin, c.majdm,c.surgr,
                  c.m1,
                  # post-transplantation
                  c.HTx, c_pc.post,
                  #
                  TH.pre,TH.post, N, dc)
{
  
  #calculate pre and post transplantations markov traces
  trace <- markov.model(tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync,
                        short.term.mortality.crmt, long.term.mortality.crmt,
                        short.term.mortality.sync, long.term.mortality.sync,
                        tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke,
                        p.dis.stroke.crmt, p.dis.stroke.sync,
                        TH.pre,TH.post)
  
  #pre-transplantation costs
  # calculate LOS days per cycle for those discharged and those remaining in the hospital
  # cost of LOS (patients either discharged or remaining in the hospital)
  LOS_pc.dis.crmt<-c(rep(30,LOS.dis.crmt%/%30),LOS.dis.crmt%%30) #number of days per cycle that discharged patients spend in the hospital before being sent home
  LOS_pc.dis.sync<-c(rep(30,LOS.dis.sync%/%30),LOS.dis.sync%%30)
  
  
  #create vector of monthly LOS costs
  LOS_pc<-matrix(0, nrow=TH.pre,ncol=2)
  LOS_pc[1:length(LOS_pc.dis.crmt),1]<-LOS_pc.dis.crmt
  LOS_pc[1:length(LOS_pc.dis.sync),2]<-LOS_pc.dis.sync
  c.LOS_pc<-(LOS_pc%*%diag(p.dis)+
               rep(30*(1-p.dis),
                   each=nrow(LOS_pc)))*c.norm.ward
  
  # Cost of adverse events (assumed they all happen within 1 year - no discount)
  
  
  p.advev  <- rbind(p.surgr,p.neurd,p.majin,p.majdm)
  c.advev        <- apply(N*p.advev*c(c.surgr,c.neurd,c.majin,c.majdm),2,sum)
  
  
  
  # Multiplying the costs of  LOS  for the proportions of patients in the alive with support state 
  
  disc_cx.pre<-1/(1+dc)^(0:(TH.pre-1))
  
  c_pc.pre<-array(NA, dim = 2)
  c_pc.pre[1]<-sum((trace$trace.pre[,1,1]*((c.LOS_pc)[,1])*disc_cx.pre))  #novel TAH
  c_pc.pre[2]<-sum((trace$trace.pre[,1,2]*((c.LOS_pc)[,2])*disc_cx.pre))  #syncardia
  
  # add treatment and adverse events costs
  cost.pre<-c_pc.pre+c.m1+c.advev
  
  #post-transplantation costs
  disc_cx.post<-1/(1+dc)^(5:(TH.post+5-1))
  
  cost.post<-vector()
  cost.post[1]<-c.HTx*trace$p.HTx.crmt+sum((trace$trace.post[,1,1]*(c_pc.post))*disc_cx.post)
  cost.post[2]<-c.HTx*trace$p.HTx.sync+sum((trace$trace.post[,1,2]*(c_pc.post))*disc_cx.post)
  
  
  cost<-cost.pre+cost.post
  
  return(cost)
}

effects <- function(tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync,
                    short.term.mortality.crmt, long.term.mortality.crmt,
                    short.term.mortality.sync, long.term.mortality.sync,
                    tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke,
                    p.dis.stroke.crmt, p.dis.stroke.sync,
                    # pre-transplantation
                    LOS.dis.crmt , LOS.dis.sync,
                    p.surgr , p.neurd , p.majin ,p.majdm,
                    e.alive.dis , e.alive.hos ,    
                    e.dis.surgr , e.dis.neurd , e.dis.majin , e.dis.majdm, 
                    # post-transplantation
                    e.post.dis.stroke, e.post.stroke.free,
                    #
                    TH.pre,TH.post, N, de)
{
  
  trace <- markov.model(tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync,
                        short.term.mortality.crmt, long.term.mortality.crmt,
                        short.term.mortality.sync, long.term.mortality.sync,
                        tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke,
                        p.dis.stroke.crmt, p.dis.stroke.sync,
                        TH.pre,TH.post)
  #pre-transplantation utilities
  LOS_pc.dis.crmt<-c(rep(30,LOS.dis.crmt%/%30),LOS.dis.crmt%%30) #number of days per cycle that discharged patients spend in the hospital before being sent home
  LOS_pc.dis.sync<-c(rep(30,LOS.dis.sync%/%30),LOS.dis.sync%%30) 
  
  LOS_pc<-matrix(0, nrow=TH.pre,ncol=2)
  LOS_pc[1:length(LOS_pc.dis.crmt),1]<-LOS_pc.dis.crmt
  LOS_pc[1:length(LOS_pc.dis.sync),2]<-LOS_pc.dis.sync
  
  e.alive<-matrix(0, nrow=TH.pre,ncol=2)
  
  
  e.alive[,1]<-     LOS_pc[,1]   *e.alive.hos[1]/(12*30)+                # Average QALYs based on number of days in and out hospital - novel TAH
    (30-LOS_pc[,1])  *e.alive.dis[1]/(12*30)
  
  e.alive[,2]<-     LOS_pc[,2]   *e.alive.hos[2]/(12*30)+                # Average QALYs based on number of days in and out hospital - Syncardia
    (30-LOS_pc[,2])  *e.alive.dis[2]/(12*30)
  
  
  # disutilities of adverse events
  e.advev_pc<-matrix(0, nrow=TH.pre,ncol=2)
  
  p.advev  <- rbind(p.surgr,p.neurd,p.majin,p.majdm)
  e.advev        <- apply(N*p.advev*c(e.dis.surgr , e.dis.neurd , e.dis.majin , e.dis.majdm),2,sum)
  
  
  disc_fx.pre<-1/(1+de)^(0:(TH.pre-1))
  
  e_pc.pre<-array(NA, dim = 2)
  
  e_pc.pre[1]<-sum((trace$trace.pre[,1,1]*(e.alive[,1]))*disc_fx.pre)
  e_pc.pre[2]<-sum((trace$trace.pre[,1,2]*(e.alive[,1]))*disc_fx.pre)
  
  #post-transplantation utilities
  disc_fx.post<-1/(1+de)^(5:(TH.post+5-1))
  
  e_pc.post<-array(NA, dim = 2)
  
  e_pc.post[1]<-sum((trace$trace.post[,1,1]*(e.post.stroke.free)*disc_fx.post))
  e_pc.post[2]<-sum((trace$trace.post[,1,2]*(e.post.dis.stroke )*disc_fx.post))
  
  e_pc<-e_pc.pre-e.advev+e_pc.post
  
  return(e_pc)
}

########################################################################################
#functions to do parallelization
########################################################################################
#function to run the model given a set of parameters
TAH.model.single.run_par_inmb<-function(list_params,
                                        objective.index, #change sign objective function: 1: positive, 2: negative
                                        trt.index, #which strategy: 1: novel TAH, 2:syncardia
                                        tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync,
                                        long.term.mortality.crmt, #short.term.mortality.crmt
                                        short.term.mortality.sync, long.term.mortality.sync,
                                        tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke,
                                        p.dis.stroke.crmt, p.dis.stroke.sync,
                                        # pre-transplantation
                                        LOS.dis.crmt , LOS.dis.sync,
                                        p.surgr , p.neurd , p.majin ,p.majdm,
                                        p.dis,
                                        c.norm.ward,c.neurd, c.majin, c.majdm,c.surgr,
                                        c.m1,
                                        # post-transplantation
                                        c.HTx,c_pc.post,
                                        e.alive.dis , e.alive.hos ,    
                                        e.dis.surgr , e.dis.neurd , e.dis.majin , e.dis.majdm, 
                                        # post-transplantation
                                        e.post.dis.stroke, e.post.stroke.free,
                                        TH.pre,TH.post, N, dc, de,
                                        WTP){
  
  short.term.mortality.crmt <-  c(list_params[1],list_params[2],list_params[3])
  p.majdm <- c(list_params[4],p.majdm)
  
  res.costs<-costs(tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
                   short.term.mortality.crmt=short.term.mortality.crmt, long.term.mortality.crmt=long.term.mortality.crmt,
                   short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
                   tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
                   p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
                   # pre-transplantation
                   LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
                   p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
                   p.dis=p.dis,
                   c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
                   c.m1=c.m1,
                   # post-transplantation
                   c.HTx=c.HTx,c_pc.post=c_pc.post,
                   TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc                   
  )
  
  #cost_sim[s,]<-res.costs
  
  res.eff<-effects(tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
                   short.term.mortality.crmt=short.term.mortality.crmt, long.term.mortality.crmt=long.term.mortality.crmt,
                   short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
                   tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
                   p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
                   # pre-transplantation
                   LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
                   p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
                   e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,    
                   e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm, 
                   # post-transplantation
                   e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
                   TH.pre=TH.pre,TH.post=TH.post, N=N, de=de                   
  )  
  #NB <- trt.effect-trt.cost/WTP
  NMB <- res.costs - res.eff*WTP
  INMB <- NMB[2]-NMB[1]
  if(objective.index==1){
    return(INMB)
  }
  else{
    return(-INMB)
  }
}

#function that will be parallelized
#inmb + Bound Optimization BY Quadratic Approximation optim version (faster) - parallel version
fun_pba_par<- function(list.index,list.params, #list of parameters that are varied - PSA
                       #list.params.other, # list ofparams - not varied
                       TH.pre,TH.post, N, dc, de,
                       WTP, 
                       index_slice, pbox.int, index.array){
  # trt.index, #which strategy: 1: novel TAH, 2:syncardia
  # #outcome.index=1, #which outcome: 1:cost, 2: effect
  #########################################
  shape <-list.params[[1]]
  scale.crmt <-list.params[[2]]
  scale.sync <-list.params[[3]]
  short.term.mortality.sync <-list.params[[4]]
  long.term.mortality.crmt <- list.params[[5]]
  long.term.mortality.sync <- list.params[[6]]
  p.surgr <-list.params[[7]]
  p.neurd <-list.params[[8]]
  p.majin <-list.params[[9]]
  p.majdm <-list.params[[10]]
  p.dis.stroke.crmt <-list.params[[11]]
  p.dis.stroke.sync <-list.params[[12]]
  LOS.dis.crmt <-list.params[[13]]
  LOS.dis.sync <-list.params[[14]]
  p.dis <- list.params[[15]]
  c.surgr <-list.params[[16]]
  c.neurd <-list.params[[17]]
  c.majin <-list.params[[18]]
  c.majdm <-list.params[[19]]
  c.m1 <-list.params[[20]]
  c_pc.post <-list.params[[21]]
  e.alive.dis <-list.params[[22]]
  e.alive.hos <-list.params[[23]]
  e.post.dis.stroke <-list.params[[24]]
  e.post.stroke.free <-list.params[[25]]
  tp.HTx2death.stroke.free <-list.params[[26]]
  tp.HTx2death.dis.stroke <-list.params[[27]]
  
  tp.alive2death_HTx.crmt<-make.surv.weibull(Time.horiz.pre,shape = shape, scale=scale.crmt,1)
  tp.alive2death_HTx.sync<-make.surv.weibull(Time.horiz.pre,shape = shape, scale=scale.sync,1)
  
  #########################################
  dim_slice <- 1 # length(index_slice) #size of set K (equation 12) or the number of possible partitions
  result_list <- c(rep(0,(dim_slice*3)))
  dim(result_list) <- c(dim_slice,3)
  #pbox.int <- pbox.int[ , , ,trt.index]
  #pbox.int <- pbox.int[ , , ,trt.index]
  num.pbox.param <- dim(pbox.int)[1] #get the number of pbox params
  for (i in 1:dim_slice){
    k <- list.index[[5]]
    lbound <- matrix(rep(0,num.pbox.param*2),nrow=num.pbox.param,ncol=2)
    ubound <- matrix(rep(0,num.pbox.param*2),nrow=num.pbox.param,ncol=2)
    
    for (j in 1:num.pbox.param){
      for (jj in 1:2){
        lbound[j,jj] <- pbox.int[j,index.array[k,j],1,jj]
        ubound[j,jj] <- pbox.int[j,index.array[k,j],2,jj]
      }
    }
    
    #using auglag
    objective.index=2
    x.0 <-(lbound[,1]+ubound[,1])/2
    #temporary fix due to lbound is higher than ubound (vice versa) - not yet fixed as of 15 May 2021
    bound.compare <- (lbound[,1] < ubound[,1])
    lbound_2 <- lbound
    ubound_2 <- ubound
    
    for (i.bound in 1:length(ubound[,1])){
      if(bound.compare[i.bound]==FALSE){
        lbound_2[i.bound,1] <- ubound[i.bound,1]
        ubound_2[i.bound,1] <- lbound[i.bound,1]
      }
    }
    ubound <- ubound_2
    lbound <- lbound_2
    
    optim.max <- bobyqa(lbound[,1],TAH.model.single.run_par_inmb,lbound[,1],ubound[,1],
                        nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000),
                        objective.index=objective.index, #change sign objective function: 1: positive, 2: negative
                        trt.index=1, #which strategy: 1: novel TAH, 2:syncardia
                        tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
                        #short.term.mortality.crmt=short.term.mortality.crmt,
                        long.term.mortality.crmt=long.term.mortality.crmt,
                        short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
                        tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
                        p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
                        # pre-transplantation
                        LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
                        p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
                        p.dis=p.dis,
                        c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
                        c.m1=c.m1,
                        # post-transplantation
                        c.HTx=c.HTx,c_pc.post=c_pc.post,
                        e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,
                        e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm,
                        # post-transplantation
                        e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
                        TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc, de=de,
                        WTP=WTP)
    
    
    #min
    objective.index=1
    optim.min <- bobyqa(lbound[,1],TAH.model.single.run_par_inmb,lbound[,1],ubound[,1],
                        nl.info = FALSE, control=list(xtol_rel=1e-9, maxeval=1000),
                        objective.index=objective.index, #change sign objective function: 1: positive, 2: negative
                        trt.index=1, #which strategy: 1: novel TAH, 2:syncardia
                        tp.alive2death_HTx.crmt=tp.alive2death_HTx.crmt, tp.alive2death_HTx.sync=tp.alive2death_HTx.sync,
                        #short.term.mortality.crmt=short.term.mortality.crmt,
                        long.term.mortality.crmt=long.term.mortality.crmt,
                        short.term.mortality.sync=short.term.mortality.sync, long.term.mortality.sync=long.term.mortality.sync,
                        tp.HTx2death.stroke.free=tp.HTx2death.stroke.free, tp.HTx2death.dis.stroke=tp.HTx2death.dis.stroke,
                        p.dis.stroke.crmt=p.dis.stroke.crmt, p.dis.stroke.sync=p.dis.stroke.sync,
                        # pre-transplantation
                        LOS.dis.crmt=LOS.dis.crmt , LOS.dis.sync=LOS.dis.sync,
                        p.surgr=p.surgr , p.neurd=p.neurd , p.majin=p.majin ,p.majdm=p.majdm,
                        p.dis=p.dis,
                        c.norm.ward=c.norm.ward,c.neurd=c.neurd, c.majin=c.majin, c.majdm=c.majdm,c.surgr=c.surgr,
                        c.m1=c.m1,
                        # post-transplantation
                        c.HTx=c.HTx,c_pc.post=c_pc.post,
                        e.alive.dis=e.alive.dis , e.alive.hos=e.alive.hos ,
                        e.dis.surgr=e.dis.surgr , e.dis.neurd=e.dis.neurd , e.dis.majin=e.dis.majin , e.dis.majdm=e.dis.majdm,
                        # post-transplantation
                        e.post.dis.stroke=e.post.dis.stroke, e.post.stroke.free=e.post.stroke.free,
                        TH.pre=TH.pre,TH.post=TH.post, N=N, dc=dc, de=de,
                        WTP=WTP)
    
    #calculate weights (assume independence)
    w=1
    for (j in 1:num.pbox.param){
      #w1 <- pbox.int[3,index.array[k,1],3]
      #w2 <- pbox.int[5,index.array[k,2],3]
      w <- w*pbox.int[j,index.array[k,j],3,1]
    }
    #optim.results[k,3] <- w1*w2
    result_list <- c(optim.min$value,optim.max$value, w)
  }
  return(result_list)
}

################################################################################
#data
################################################################################
c.TAH.dev     <-c(200000,100000)  

Time.horiz.pre  <-35             #time horizon for the pre-trasnplant section (months)
Time.horiz.post <-300          #time horizon for the post-trasnplant section (months)

nsims       <-10000
cohort      <-1                  #cohort size
disc.fx     <-(1+0.035)^(1/12)-1 # (monthly) discount rate for effect
disc.cx     <-(1+0.035)^(1/12)-1 # (monthly) discount rate for costs

pat.year<- 50                    # number of expected patients * year


dir.op.costs <-50000             # direct operating costs (cost of selling goods)
ind.op.costs <-1000000           # indirect operating costs
mean.mktshr  <-4.5               # parameters of the distribution linking effectiveness with mkt share - mean
sd.mktshr    <-2              # parameters of the distribution linking effectiveness with mkt share - sd

t.rel.ind    <-10                      # time during which the decision problem is relevant - INDUSTRY
disc.evi.ind <-1/(1+0.035)^0:t.rel.ind # vector of discount factors

t.rel.pay    <-10                      # time during which the decision problem is relevant - PAYER
disc.evi.pay <-1/(1+0.035)^0:t.rel.pay # vector of discount factors


k            <-seq(0,60000,100)    # k vector for sensitivity analysis

#loda data from prior.bugs model
load("bugs.prior.model.Rdata")
#plot(prior.model)

prior.model.data<-prior.model$BUGSoutput$sims.matrix
prior.model.data<-as.data.frame(prior.model.data)

#Take samples
sat.seed= 123
samples<-sample(1:nrow(prior.model.data),nsims,replace=FALSE)
prior.model.data<-prior.model.data[samples,]

#transition probabilities from heart transplant to death
#load(paste(path,"tp.HTx2death.Rdata",sep=""))
load("tp.HTx2death.Rdata")

tp.HTx2death.stroke.free   <- tp.HTx2death[1:Time.horiz.post-1]                                   
tp.HTx2death.dis.stroke    <- tp.HTx2death[1:Time.horiz.post-1]  #mortality from transplanted with or without stroke assumed equal (to avoid double counting)

#costs
c.HTx         <-20545                                       # cost of heart transplant procedure (sharpless et al.)
c.norm.ward   <-251                                         #cost of day in normal ward after TAH implant (tariffa lungo degenza per malattie cardiocircolatorie)

#HRQpL
e.post.stroke.free         <-0.83/12                                                                         #monthly HRQoL Patients with HTx                    

#disutilities - assumed already taken into account in alive with support utility values to avoid double counting
e.dis.surgr   <-0                #disutility from bleeding/surgical repair event
e.dis.neurd   <-0                #disutility from neurological disfunction event
e.dis.majin   <-0                #disutility major infection               event
e.dis.majdm   <-0                #disutility from infection                event

k.bench  <-30000

cost_sim<-array(NA,c(nsims,2))
effs_sim<-array(NA,c(nsims,2))

basecase.data<-array(NA,c(nsims,13))
colnames(basecase.data)<-c("effs.crmt","effs.sync","cost.crmt","cost.sync","NMB",
                           "delta.effs", "delta.cost","approved","mrkt.share.sim",
                           "pop.revs","pop.cost","go.dec","nogo.dec")

################################################################################
##### preparing data for functions ######
################################################################################
# Parameter	Description
# e.post.dis.stroke	Permanent disutility following a disabling stroke while on support
# p.dis.stroke	Proability of disabling (moderate or severe) stroke (conditional on having a stroke) while on support
# pi.death.1.1.	Probability of death in the first month after implant - conditional on having an event (death or transplant) - novel TAH
# pi.death.2.1.	Probability of death in the second month after implant - conditional on having an event (death or transplant) - novel TAH
# pi.death.3.1.	Probability of death in the third month after implant - conditional on having an event (death or transplant) - novel TAH
# 
# pi.se.crmt.4.	Probability of having a major device malfunction
# scale.crmt	location parameter of the survival curve before transplantation - competing events: death or transplantation - NOVEL TAH
# shape	shape parameter of the survival curve before transplantation - competing events: death or transplantation - NOVEL TAH
# e.post.stroke.free	Quality of life after transplant - without disabling stroke
# 

#transition probabilities from alive with support to death or transplant
shape<-prior.model.data$shape[]
scale.crmt<-prior.model.data$scale.crmt[]
scale.sync<-prior.model.data$scale.sync[]

tp.alive2death_HTx.crmt<-make.surv.weibull(Time.horiz.pre,shape = shape, scale=scale.crmt,1)
tp.alive2death_HTx.sync<-make.surv.weibull(Time.horiz.pre,shape = shape, scale=scale.sync,1)

short.term.mortality.crmt<-prior.model$BUGSoutput$sims.list$pi.death[samples,1:3,1]  #<- 3 PBA parameters
short.term.mortality.sync<-prior.model$BUGSoutput$sims.list$pi.death[samples,1:3,2]

long.term.mortality.crmt<-prior.model$BUGSoutput$sims.list$pi.death[samples,4,1]
long.term.mortality.sync<-prior.model$BUGSoutput$sims.list$pi.death[samples,4,2]

#Probabilities of AE  for novel TAH [1] and Syncardia [2]
p.surgr <-cbind(prior.model.data$`pi.se.crmt[1]`[],prior.model.data$`pi.se.sync[1]`[])  #bleeding/surgical repair
p.neurd <-cbind(prior.model.data$`pi.se.crmt[2]`[],prior.model.data$`pi.se.sync[2]`[])  #stroke
p.majin <-cbind(prior.model.data$`pi.se.crmt[3]`[],prior.model.data$`pi.se.sync[3]`[])  #major infection
p.majdm <-cbind(prior.model.data$`pi.se.crmt[4]`[],prior.model.data$`pi.se.sync[4]`[])  #major device malfunction #<- one PBA parameter

p.dis.stroke.crmt<-p.dis.stroke.sync<-prior.model.data$p.dis.stroke[]

#resurce use and costs
LOS.dis.crmt   <-prior.model.data$LOS.dis.crmt[]                                          #LOS of discharged patients novel TAH
LOS.dis.sync   <-LOS.dis.crmt                                                              #LOS of discharged patients Syncardia (assumed equal)
p.dis          <-cbind(prior.model.data$`p.dis[1]`[],prior.model.data$`p.dis[2]`[])          #probability of being discharged for novel TAH [1,] and Syncardia [2,]
b          <-cbind(prior.model.data$`p.dis[1]`[],prior.model.data$`p.dis[2]`[])          #probability of being discharged for novel TAH [1,] and Syncardia [2,]

c.surgr        <-prior.model.data$`cost[1]`[]                                             #cost of bleeding/surgical repair event
c.neurd        <-prior.model.data$`cost[2]`[]                                             #cost of neurological disfunction
c.majin        <-prior.model.data$`cost[3]`[]                                             #cost of major infection
c.majdm        <-prior.model.data$`cost[4]`[]                                             #cost of infection
c.m1           <-prior.model$BUGSoutput$sims.list$c.m1[,]                                  #cost of TAH implant + 1st month
# c_pc.post[,1:7] <-prior.model.data$BUGSoutput$sims.list$c_pc.post[,1:7]
# c_pc.post[8:24]<- prior.model.data$BUGSoutput$sims.list$c_pc.post[1,7]

#c.m1           <-c(prior.model.data$`c.m1[1]`,prior.model.data$`c.m1[2]`)                                  #cost of TAH implant + 1st month
c.m1           <-cbind(prior.model.data$`c.m1[1]`,prior.model.data$`c.m1[2]`)  

#c_pc.post      <-rep(0,Time.horiz.post)                                                    #monthly cost after HTx (Sharpless et al.)
c_pc.post      <-matrix(rep(0,nsims*Time.horiz.post), nrow=nsims, ncol=Time.horiz.post)     #monthly cost after HTx (Sharpless et al.)

c_pc.post[,1:7] <-c(prior.model.data$`c_pc.post[1]`[],prior.model.data$`c_pc.post[2]`[],prior.model.data$`c_pc.post[3]`[],
                    prior.model.data$`c_pc.post[4]`[],prior.model.data$`c_pc.post[5]`[],prior.model.data$`c_pc.post[6]`[],
                    prior.model.data$`c_pc.post[7]`[])

c_pc.post[,8:24]<- prior.model.data$`c_pc.post[7]`[]

#HRQpL
e.alive.dis         <-cbind(prior.model.data$`e.alive.dis[1]`[],prior.model.data$`e.alive.dis[2]`[])/12   #monthly HRQoL Patients with TAH discharged at home - Novel TAH [1] and Syncardia [2]
e.alive.hos         <-e.alive.dis                                                                     #monthly HRQoL Patients with TAH at hospital        - Novel TAH [1] and Syncardia [2]
e.post.dis.stroke   <-prior.model.data$e.post.dis.stroke[]
e.post.stroke.free   <-prior.model.data$e.post.stroke.free[]

################################################################################
#start of PBA
##############################################################################
#1. uncertainty representation
num.pbox.params <- 4 #only 4 parameters that are treated as PBA parameters
#default: summary statistics calculated from PSA data

#bound scenarios (variable: bound):
#"extreme": [0,1]
#"PSA": [min, max]
#"user" supplied

#summary statistics supplied, choose one of the following (variable: statistics)
#minmax: min, max
#"mean": min, max, mean
#"std": min, max, mean, std
#"median", min, max, median
#"median_std": min, max, mean, std, median

#looping over number of pbox intervals
pbox.int.loop <- c(10)
statistics.loop <-c("mean")

#number of intervals for inmb
min.x <- -200000 #min(df.nmb.upper,df.nmb.lower)
max.x <- 250000#max(df.nmb.upper,df.nmb.lower)
n.x <-100000
steps <- (max.x-min.x)/n.x
steps.vec <- c(rep(steps,n.x))
x.values <- seq(min.x,max.x,length.out=n.x)

#preparing matrices for results by combinations of pbox intervals and summary statistics
#columns: number of x, bounds: lower upper, statistics, pbox int
cdf.all <- array(c(rep(0,n.x*2*length(statistics.loop),length(pbox.int.loop))), dim=c(n.x,2,length(statistics.loop),length(pbox.int.loop)))
#start of loop over pbox intervals (pbox.int.loop) and over which statistics (statistics.loop)
for (i.int in 1:length(pbox.int.loop)){
  for (i.stat in 1:length(statistics.loop)){
    
    bound <- "extreme"
    statistics <- statistics.loop[i.stat]
    num.pbox.int <- pbox.int.loop[i.int] #number of p-box discretization for UP
    #combination.med.std <- "median_std"
    #################################################
    if(bound=="extreme"){
      short.term.mortality.crmt.1.min<-0
      short.term.mortality.crmt.2.min<-0
      short.term.mortality.crmt.3.min<-0
      p.majdm.crmt.min <- 0
      
      short.term.mortality.crmt.1.max<- 1
      short.term.mortality.crmt.2.max<- 1
      short.term.mortality.crmt.3.max<- 1
      p.majdm.crmt.max <- 1
    } else if(bound=="PSA"){
      short.term.mortality.crmt.1.min<-  min(short.term.mortality.crmt[,1])
      short.term.mortality.crmt.2.min<-  min(short.term.mortality.crmt[,2])
      short.term.mortality.crmt.3.min<- min(short.term.mortality.crmt[,3])
      p.majdm.crmt.min <- min(p.majdm[,1])
      #maximum
      short.term.mortality.crmt.1.max<- max(short.term.mortality.crmt[,1])
      short.term.mortality.crmt.2.max<- max(short.term.mortality.crmt[,2])
      short.term.mortality.crmt.3.max<- max(short.term.mortality.crmt[,3])
      p.majdm.crmt.max <- max(p.majdm[,1])
    } else{
      #minimum
      short.term.mortality.crmt.1.min<-  0 #min(short.term.mortality.crmt[,1])
      short.term.mortality.crmt.2.min<- 0  #min(short.term.mortality.crmt[,2])
      short.term.mortality.crmt.3.min<-0 #min(short.term.mortality.crmt[,3])
      p.majdm.crmt.min <- 0 #min(p.majdm[,1])
      #maximum
      short.term.mortality.crmt.1.max<-1 #max(short.term.mortality.crmt[,1])
      short.term.mortality.crmt.2.max<- 1#max(short.term.mortality.crmt[,2])
      short.term.mortality.crmt.3.max<- 1 #max(short.term.mortality.crmt[,3])
      p.majdm.crmt.max <-0.5 #max(p.majdm[,1])
    }
    #median
    short.term.mortality.crmt.1.median<-median(short.term.mortality.crmt[,1])
    short.term.mortality.crmt.2.median<-median(short.term.mortality.crmt[,2])
    short.term.mortality.crmt.3.median<-median(short.term.mortality.crmt[,3])
    p.majdm.crmt.median <- median(p.majdm[,1])
    #mean
    short.term.mortality.crmt.1.mean <-mean(short.term.mortality.crmt[,1])
    short.term.mortality.crmt.2.mean <-mean(short.term.mortality.crmt[,2])
    short.term.mortality.crmt.3.mean <-mean(short.term.mortality.crmt[,3])
    p.majdm.crmt.mean <- mean(p.majdm[,1])
    #standard deviation
    short.term.mortality.crmt.1.std<-sd(short.term.mortality.crmt[,1])
    short.term.mortality.crmt.2.std <-sd(short.term.mortality.crmt[,2])
    short.term.mortality.crmt.3.std <-sd(short.term.mortality.crmt[,3])
    p.majdm.crmt.std <- sd(p.majdm[,1])
    
    ##################################################################
    #steps S1 and S2
    
    #params.table: table of parameters which are represented by pboxs and the data available for each parameter
    #rows: PBA parameter
    #columns: min, max, median, mean, std (available summary statistics)
    #3rd dim: treatment #this is an extra dimension for any further stratification of the cohort. It is not necessary for p-box implementation and model-specific.
    #params.table <- matrix(rep(0,num.pbox.params*5),nrow=num.pbox.params,ncol=5)
    params.table <- array(rep(0,num.pbox.params*5),dim=c(num.pbox.params, 5, 2))
    mean.std <- 0.2 #assumption: std as a fraction of mean as no data is available
    #crmt params for the PBA params
    
    params.table[1, ,1] = c(short.term.mortality.crmt.1.min,short.term.mortality.crmt.1.max,short.term.mortality.crmt.1.median,short.term.mortality.crmt.1.mean, short.term.mortality.crmt.1.std )
    params.table[2, ,1] = c(short.term.mortality.crmt.2.min,short.term.mortality.crmt.2.max,short.term.mortality.crmt.2.median,short.term.mortality.crmt.2.mean, short.term.mortality.crmt.2.std )
    params.table[3, ,1]= c(short.term.mortality.crmt.3.min,short.term.mortality.crmt.3.max,short.term.mortality.crmt.3.median,short.term.mortality.crmt.3.mean, short.term.mortality.crmt.3.std )
    params.table[4, ,1]= c(p.majdm.crmt.min,p.majdm.crmt.max,p.majdm.crmt.median,p.majdm.crmt.mean,p.majdm.crmt.std)
    
    #steps S3 and S4
    #inverse sampling
    #discretization over the cdf values (u) (NOTE: this is range of p-box)
    #partitioning u domain
    #columns of inv.int matrix: (1) lower bound - interval (2) upper bound - interval (3) width of interval
    inv.int <- matrix(rep(0,(num.pbox.int-1)*3),nrow=num.pbox.int-1,ncol=3) 
    dummy <- seq(0,1,length.out=num.pbox.int)
    for (i in 1:(num.pbox.int-1)){
      inv.int[i,1]=dummy[i]
      inv.int[i,2]=dummy[i+1]
      inv.int[i,3]=dummy[i+1]-dummy[i]
    }
    #sampling over x values, using inverse pbox functions, applied to the specified u values from above 
    #dimensions: (1) variable, (2) intervals, (3) lower or upper bounds or weight
    #pbox.int <-array(rep(0,num.pbox.params*(num.pbox.int-1)*3), dim=c(num.pbox.params, num.pbox.int-1, 3))
    pbox.int <-array(rep(0,num.pbox.params*(num.pbox.int-1)*3*2), dim=c(num.pbox.params, num.pbox.int-1, 3,2))
    for (j in 1:num.pbox.params){ #repeat over model parameters subjected to pbox
      for (i in 1:(num.pbox.int-1)){
        for (jj in 1:1){
          if(statistics=="minmax"){
            pbox.int[j,i,1,jj]=pbox_minmax_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj])
            pbox.int[j,i,2,jj]=pbox_minmax_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj])
          } else if (statistics=="median"){
            pbox.int[j,i,1,jj]=pbox_median_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj])
            pbox.int[j,i,2,jj]=pbox_median_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj])
          } else if (statistics=="mean"){
            pbox.int[j,i,1,jj]=pbox_mean_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj])
            pbox.int[j,i,2,jj]=pbox_mean_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj])
          } else if (statistics=="std"){
            pbox.int[j,i,1,jj]=pbox_std_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj],params.table[j,5,jj])
            pbox.int[j,i,2,jj]=pbox_std_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,4,jj],params.table[j,5,jj])
          } else if (statistics=="median_std"){
            pbox.int[j,i,1,jj]=pbox_median_std_inv(inv.int[i,1],"upper",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj],params.table[j,4,jj],params.table[j,5,jj])
            pbox.int[j,i,2,jj]=pbox_median_std_inv(inv.int[i,2],"lower",params.table[j,1,jj],params.table[j,2,jj],params.table[j,3,jj],params.table[j,4,jj],params.table[j,5,jj])
          }  
          pbox.int[j,i,3,jj]=inv.int[i,3]
        }
      }
    }
    ################################################################################
    #2. uncertainty propagation
    ################################################################################
    #create multi-indices (for each p-box parameter)
    num.param.pba <- num.pbox.params #number of model parameters with p-box
    n.indices <- (num.pbox.int-1)**num.param.pba #number of unique indices - for each parameter, there will be num.pbox.int index
    input_data <- 1:n.indices

    #steps S5
    vec <- c(1:(num.pbox.int-1))
    lst <- lapply(numeric(num.param.pba), function(x) vec)
    index.array <- as.matrix(expand.grid(lst))
    #print(dim(index.array)[1])
    index.array.last <- rep(0,dim(index.array)[1])
    index.array <- cbind(index.array,index.array.last)
    #print(dim(index.array)[1])
    #index for par version
    list.index <- list()
    index_in_index.array <- seq(1,n.indices)
    for (zz in 1:n.indices){
      list.index[[zz]] <- list(index.array[zz,1],index.array[zz,2],index.array[zz,3],index.array[zz,4],index_in_index.array[zz]) 
    }
    ##############################################################################
    #par version
    list_test.pbox <- list()
    # list_test <- list()
    #fixed PSA parameters at their average values
    list_test.pbox <- list(mean(shape),
                           mean(scale.crmt),
                           mean(scale.sync),
                           #short.term.mortality.crmt[i,], 3 PBA parameters
                           apply(short.term.mortality.sync,2 ,mean),
                           mean(long.term.mortality.crmt),
                           mean(long.term.mortality.sync),
                           apply(p.surgr,2 ,mean),
                           apply(p.neurd,2 ,mean),
                           apply(p.majin,2 ,mean),
                           mean(p.majdm[,2]), #p.majdm for novTAH is PBA param
                           mean(p.dis.stroke.crmt),
                           mean(p.dis.stroke.sync),
                           mean(LOS.dis.crmt),
                           mean(LOS.dis.sync),
                           apply(p.dis,2 ,mean),
                           mean(c.surgr),
                           mean(c.neurd),
                           mean(c.majin),
                           mean(c.majdm),
                           apply(c.m1,2,mean),
                           #mean(c.m1),
                           apply(c_pc.post,2 ,mean),
                           apply(e.alive.dis,2 ,mean),
                           apply(e.alive.hos,2 ,mean),
                           mean(e.post.dis.stroke),
                           mean(e.post.stroke.free),
                           tp.HTx2death.stroke.free,
                           tp.HTx2death.dis.stroke
    )
    
    
    params_pdf.pbox <- list_test.pbox
    # params_pdf.MC <- list_test.MC
    ##############################################################################
    #running PBA
    ################################################################################
    #steps S6
    #approach 1 with parallel over K
    cl <- makeClusterPSOCK(
      # Public IP number of EC2 instance
      workers=availableCores(),
      dryrun = FALSE,
      connectTimeout = 120,
      verbose=TRUE,
      outfile=NULL,
      rshlogfile=TRUE
    )
    plan(list(tweak(cluster, workers = cl), multiprocess))
    #running the optimization
    result_list_outer.sp.pbox <- future.apply::future_lapply(list.index, fun_pba_par ,params_pdf.pbox,
                                                             TH.pre=Time.horiz.pre,TH.post=Time.horiz.post, N=cohort, dc=disc.cx, de=disc.fx,
                                                             WTP=k.bench,
                                                             index_slice=input_data,pbox.int=pbox.int,index.array=index.array)

    ################################################################################
    #processing PBA outputs 
    
    #parallel over K
    nmb <- sapply (result_list_outer.sp.pbox, function (x) {length (x) <- 3; return (x)})
    nmb.t <- t(nmb)
    #INMB approach (as output) 16 May 2021
    df.nmb.lower <- data.frame(nmb.t[,1])
    df.nmb.upper <- data.frame(abs(nmb.t[,2]))
    df.weight <- data.frame(nmb.t[,3]) #weight
    
    #step S7
    #calculate empirical cdf of INMB
    for (j in 1:n.x){
      #cdf.x.lower.sp[j]<- sum(df.x$V2[df.x$V1<=x.values[j]]) #note reverse, see Ferson 2015
      cdf.x.lower.sp.pbox[j]<- sum(df.weight[df.nmb.upper<=x.values[j]]) #note reverse, see Ferson 2015
      cdf.x.upper.sp.pbox[j]<- sum(df.weight[df.nmb.lower<=x.values[j]])
    }
    
    x.mean.lower.sp.df.pbox <- as.data.frame(cbind(x.values,cdf.x.lower.sp.pbox))
    colnames(x.mean.lower.sp.df.pbox) <- c("x","y")
    
    x.mean.upper.sp.df.pbox <- as.data.frame(cbind(x.values,cdf.x.upper.sp.pbox))
    colnames(x.mean.upper.sp.df.pbox) <- c("x","y")
    #collecitng the result for the current loop
    cdf.all[,1,i.stat,i.int] <-cdf.x.lower.sp.pbox
    cdf.all[,2,i.stat,i.int] <-cdf.x.upper.sp.pbox
    
  }
}
#get the results for post-processing
pbox1.upper.std <- as.data.frame(cbind(x.values,cdf.all[,2,1,1]))
pbox1.lower.std <- as.data.frame(cbind(x.values,cdf.all[,1,1,1]))

colnames(pbox1.upper.std) <- c("V1","V2")
colnames(pbox1.lower.std) <- c("V1","V2")

#calculate the minimum (maximum) of INMB by averaging UBF (LBF) 
mean.upper.std <-x.values[1:n.x-1]%*%diff(pbox1.upper.std$V2)
mean.lower.std <-x.values[1:n.x-1]%*%diff(pbox1.lower.std$V2)

########################################################################
########################################################################
#hurwicz
n.alpha <- 1001
alphas <- c(seq(0,1,(n.alpha-1)^(-1)))
mean.std <- c(rep(0,n.alpha))
mean.upper.std.s <- mean.upper.std[1,1]
mean.lower.std.s <- mean.lower.std[1,1]

for(a in 1:n.alpha){
  mean.std[a] <- alphas[a]*mean.upper.std.s+(1-alphas[a])*mean.lower.std.s
}

########################################################################
#end of code
########################################################################