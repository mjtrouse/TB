#Rcodes for replication clock plasmid/lung-spleen Mtb dissemination

#dependencies for ODE mods
library(FME)
library(graphics)

#dependencies for simulation mods
library(GillespieSSA2)
library(ggplot2)
library(scales)
library(dplyr)
library(adaptivetau)

#Data
load("~/tb_repo_data.RData")

#plotting fits
#change dataset and solution DF
{
#spleen plots
plot(dataset.both[,1], log10(dataset.both[,4]),
     xlab="day", ylab="free spleen", col = 3, pch = 22, ylim=c(0,6))#free spleen
points(averages.log$day, averages.log$mean.log.Fs, 
       col="blue", pch=17)
lines(mytimes0,log10(solmodfit.m04[,4]), col=2, lwd = 2)

plot(dataset.both[,1], log10(dataset.both[,5]),
     xlab="day", ylab="plasmid spleen", col = 3, pch = 22, ylim=c(0,6))  #plasmid spleen
points(averages.by.day$day, averages.log$mean.log.Ps, 
       col="blue", pch=17)
lines(mytimes0,log10(solmodfit.m04[,5]), col=2, lwd = 2)

#lung plots
plot(dataset.both[[1]], log10(dataset.both[[2]]),
     xlab="day", ylab="free lung", col = 3, pch = 22, ylim=c(0,8))#free lung
points(averages.log$day, averages.log$mean.log.Fl, 
       col="blue", pch=17)
lines(mytimes0,log10(solmodfit.m04i[,2]), col=2, lwd = 2)

plot(dataset.both[,1], log10(dataset.both[,3]),
     xlab="day", ylab="plasmid lung", col = 3, pch = 22, ylim=c(0,8))  #plasmid lung
points(averages.log$day, averages.log$mean.log.Pl, 
       col="blue", pch=17)
lines(mytimes0,log10(solmodfit.m04[,3]), col=2, lwd = 2)
}

#################################################################
##Independent spleen infection model fitting starting day 13
#functions for independent model
{
#residual function
resids=function(fit, fix, dataset, modelname,  modelcols, datacols, init.vector)
{
  timepoints=dataset[[1]]
  timename=names(dataset)[1]
  if(is.null(datacols)) datacols=2:(dim(dataset)[2]) # default to all x1...xn in the data
  if (is.null(init.vector)) init.vector=dataset[1, datacols] # take init conds from the data if necessary
  
  trans<-function(x) log10(x) # takes the log10 of each value to fit over scale rather than specific points. 
  soln=ode.solve(fit=fit,fix=fix,modelname=modelname, init.vector=init.vector, timepoints=timepoints, timename=timename)
  
  resids=as.vector(as.matrix(trans(soln[,modelcols])-trans(dataset[,datacols])))
  resids=na.omit(resids) # omits NA values in the data
  resids=resids[!is.infinite(resids)] # removes infinite residuals
  #resids<-resids[-1]   # removing first fake data point
}

#ODE solve framework for independent model

ode.solve=function(fit, fix, modelname, init.vector, timepoints, timename){ 
  allpars <- c(fit, fix)
  init.conds <- c()
  for(i in 1:length(init.vector)){
    init.conds <- c(init.conds, allpars[init.vector[[i]]])
  }
  soln=as.data.frame(lsoda(y=init.conds, times=timepoints, func=modelname, parms=allpars))
  soln$allp<-soln[,3]  #count all plasmid bearing bacteria
  soln$total<-soln[,2] + soln[,3]  #count total bacteria
  soln$perc<-soln$allp/soln$total #count percent plasmid
  #soln$nonq<-soln[,2] + soln[,3] #count non quiescent bacteria
  if (is.null(timename)) timename="Hours" 
  names(soln)[1]=timename
  return(soln)
}

#replication function 
bacrep <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- 0
      dPs.dt <- 0
    }
    if (time >13 & time < 26 ){
      dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs
      dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps
    }
    if (time >26 & time < 69 ){
      dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs
      dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps
    }
    if (time >69){
      dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs
      dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps
    }   
    return(list(c(dFs.dt, dPs.dt)))})
}
}

#indepedent spleen infection model
{
mytimes<- seq(13,111,.1) #spleen data from day 13-111
init.vector <- c("Fs", "Ps")
fix <- c( s1=0.18, Fs = 21, Ps = 7)  #averages of day 13
fit <- c(del1= 1.09, del2= .04, del3=0.18, rep1= 1.7, rep2= .07, rep3= .21)

modelfitq0<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrep,dataset=(dataset.independent),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfitq0<-ode.solve(coef(modelfitq0),fix,bacrep,init.vector,mytimes,"day")
}
################################################################
#Quiescence model fitting spleen data only from day 13

#functions for quiescence mods
{
#residual function same

#ode solve includes counts for number of quiescent bacteria in output
ode.solveq=function(fit, fix, modelname, init.vector, timepoints, timename){ 
  allpars <- c(fit, fix)
  init.conds <- c()
  for(i in 1:length(init.vector)){
    init.conds <- c(init.conds, allpars[init.vector[[i]]])
  }
  soln=as.data.frame(lsoda(y=init.conds, times=timepoints, func=modelname, parms=allpars))
  soln$allp<-soln[,3] + soln[,5]  #count all plasmid bearing bacteria
  soln$total<-soln[,2] + soln[,3] + soln[,4] + soln[,5] #count total bacteria
  soln$perc<-soln[,6]/soln[,7] #count percent plasmid
  soln$nonq<-soln[,2] + soln[,3] #count non quiescent bacteria
  if (is.null(timename)) timename="Hours" 
  names(soln)[1]=timename
  return(soln)
}

#replication function 
bacrepq <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- 0
      dPs.dt <- 0
      dFsq.dt<-0
      dPsq.dt<-0
    }
    if (time >13 & time < 26 ){
      dFs.dt <- rep1*Fs + s1*rep1*Ps - q1*rep1*Fs-del1*Fs
      dPs.dt <- rep1*Ps - s1*rep1*Ps - q1*rep1*Ps-del1*Ps
      dFsq.dt<- q1*rep1*Fs
      dPsq.dt<- q1*rep1*Ps
    }
    if (time >26 & time < 69 ){
      dFs.dt <- rep2*Fs + s1*rep2*Ps - q1*rep2*Fs- del2*Fs
      dPs.dt <- rep2*Ps - s1*rep2*Ps - q1*rep2*Ps- del2*Ps
      dFsq.dt<-q1*rep2*Fs
      dPsq.dt<-q1*rep2*Ps
      
    }
    if (time >69){
      dFs.dt <- rep3*Fs + s1*rep3*Ps - q1*rep3*Fs-del3*Fs
      dPs.dt <- rep3*Ps - s1*rep3*Ps -q1*rep3*Ps-del3*Ps
      dFsq.dt<- q1*rep3*Fs
      dPsq.dt<- q1*rep3*Ps
    }   
    return(list(c(dFs.dt, dPs.dt, dFsq.dt, dPsq.dt)))})} 
}

#quiescence model 
{
mytimes<- seq(13,111,.1) #spleen data from day 13-111
init.vector <- c("Fs", "Ps", "Fsq", "Psq")
fix <- c(Fs = 21, Ps = 7, Fsq=0, Psq=0, s1=0.18, q1=0)  #day 13 avg, change q1 for rate
#fit <- c(del1= 1.09, del2= .04, del3=0.18, rep1= 1.7, rep2= .07, rep3= .21) #calculated estimates
#using independent model fit rates
fit<-c(del1=0.99, del2=0.04, del3=0.24, rep1=1.40, rep2=0.07, rep3=0.27)

modelfitq0.0<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepq,dataset=(dataset.independent),modelcols=c(6,7),datacols=c(3,4),init.vector=init.vector)
solmodfitq0.0<-ode.solveq(coef(modelfitq0.0),fix,bacrepq,init.vector,mytimes,"day")
}
############################################################################
#Influx models fitting spleen data only from day 13

#influx from day 13 functions
{
#resids and odesolve same as independent spleen model

#constant influx rate
bacrepisc <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- 0
      dPs.dt <- 0
    }
    if (time >13 & time < 26 ){
      dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs+If
      dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps+Ip
    }
    if (time >26 & time < 69 ){
      dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs+If
      dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps+Ip
    }
    if (time >69){
      dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs+If
      dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps+Ip
    }   
    return(list(c(dFs.dt, dPs.dt)))})
}

#two influx rates
bacrepis2 <- function(time, state, parms){
    with(as.list(c(state,parms)),{
      if (time <= 13) {
        dFs.dt <- 0
        dPs.dt <- 0
      }
      if (time >13 & time < 26 ){
        dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs+If
        dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps+Ip
      }
      if (time >26 & time < 69 ){
        dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs+If2
        dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps+Ip2
      }
      if (time >69){
        dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs+If2
        dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps+Ip2
      }   
      return(list(c(dFs.dt, dPs.dt)))})
  }

#three influx rates
bacrepis3 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- 0
      dPs.dt <- 0
    }
    if (time >13 & time < 26 ){
      dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs+If
      dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps+Ip
    }
    if (time >26 & time < 69 ){
      dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs+If2
      dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps+Ip2
    }
    if (time >69){
      dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs+If3
      dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps+Ip3
    }   
    return(list(c(dFs.dt, dPs.dt)))})
}
}

#same initial vector and times for all influx models from day 13
mytimes<- seq(13,111,.1) #spleen data from day 13-111
init.vector <- c("Fs", "Ps")

#constant influx rate model from day 13
{
fix<-c(Fs = 21,Ps = 7, s1=0.18,If=2.9680743 , Ip=0.5100717 )

fit <- c(del1= 1.16228066 , del2= 0.01194348 , del3=0.18870354 ,
         rep1= 1.56614220 , rep2= 0.03806798 , rep3= 0.21992302 )

modelfit.isc<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepisc,dataset=(dataset.independent),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.isc<-ode.solve(coef(modelfit.isc),fix,bacrepisc,init.vector,mytimes,"day")
}
#two influx rate model from day 13
{
fix<-c(Fs = 21,Ps = 7, s1=0.18,If=20.0479587   , Ip=1.4893622 , 
       If2=122.4857412    , Ip2=0.3538574 )
fit <- c( del1 = 1.13370367 , del2= 0.05126102, del3=0.20285011, rep1= 1.47610914, 
          rep2=0.07022283, rep3=0.23505741)

modelfit.inf2<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepis2,dataset=(dataset.independent),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.inf2<-ode.solve(coef(modelfit.inf2),fix,bacrepis2,init.vector,mytimes,"day")
}

#three influx rate model from day 13
{
fix<-c(Fs = 21,Ps = 7, s1=0.18,If=18.8754690  , Ip=1.5676295 , 
       If2=142.3095933    , Ip2=0.2309642   ,
       If3= 50.7927315  , Ip3=0.4893566)

fit <- c(del1 = 1.13370367   , del2= 0.05126102   , del3=0.20285011    , 
         rep1= 1.47610914   , rep2=0.07022283  , rep3= 0.23505741)


modelfit.inf3<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepis3,dataset=(dataset.independent),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.inf3<-ode.solve(coef(modelfit.inf3),fix,bacrepis3,init.vector,mytimes,"day")
}

#####################################################################
#influx models fitting to spleen data only from day 0

#functions for day 0 influx mods
{
#resids and odesolve same and independent spleen 

#constant influx from day 0
bacrepisc0 <- function(time, state, parms){
      with(as.list(c(state,parms)),{
        if (time <= 13) {
          dFs.dt <- If
          dPs.dt <- Ip
        }
        if (time >13 & time < 26 ){
          dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs+If
          dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps+Ip
        }
        if (time >26 & time < 69 ){
          dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs+If
          dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps+Ip
        }
        if (time >69){
          dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs+If
          dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps+Ip
        }   
        return(list(c(dFs.dt, dPs.dt)))})
}

#2 influx rates from day 0
bacrepis20 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- If
      dPs.dt <- Ip
    }
    if (time >13 & time < 26 ){
      dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs+If
      dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps+Ip
    }
    if (time >26 & time < 69 ){
      dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs+If2
      dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps+Ip2
    }
    if (time >69){
      dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs+If2
      dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps+Ip2
    }   
    return(list(c(dFs.dt, dPs.dt)))})

}

#3 rate model only when starting from day 13

#4 influx rates from day 0
bacrepis40 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- If0
      dPs.dt <- Ip0
    }
    if (time >13 & time < 26 ){
      dFs.dt <- rep1*Fs + s1*rep1*Ps - del1*Fs+If
      dPs.dt <- rep1*Ps - s1*rep1*Ps - del1*Ps+Ip
    }
    if (time >26 & time < 69 ){
      dFs.dt <- rep2*Fs + s1*rep2*Ps - del2*Fs+If2
      dPs.dt <- rep2*Ps - s1*rep2*Ps -del2*Ps+Ip2
    }
    if (time >69){
      dFs.dt <- rep3*Fs + s1*rep3*Ps - del3*Fs+If3
      dPs.dt <- rep3*Ps - s1*rep3*Ps -del3*Ps+Ip3
    }   
    return(list(c(dFs.dt, dPs.dt)))})
  
}
}

#conditions same for all influx mods from day 0
mytimes0<-seq(0,111,0.1)
init.vector <- c("Fs", "Ps")

#constant influx from day 0 model
{
fix<-c(Fs = 0,Ps = 0, s1=0.18,If=1.98141027 , Ip=0.99750098)
fit <- c(del1= 1.17332857    , del2= 0.01379888    , del3=0.19065362,   
         rep1= 1.56707412    , rep2= 0.04085744    , rep3= 0.21860093)

modelfit.isc0<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepisc0,dataset=(dataset.indepedent.d0),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.isc0<-ode.solve(coef(modelfit.isc0),fix,bacrepisc0,init.vector,mytimes0,"day")
}

#2 influx rates from day 0 model
{
#fitted influx then fitted rates
fix<-c(Fs = 0,Ps = 0, s1=0.18, If=3.5620404  , Ip=1.5001247   , 
       If2=11.2663606    , Ip2=0.4717342    )
fit <- c( del1 = 1.13370367    , del2= 0.05126102    , del3=0.20285011    , 
          rep1= 1.47610914    , rep2=0.07022283    , rep3= 0.23505741)

modelfit.is20<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepis20,dataset=(dataset.indpendent.d0),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.is20<-ode.solve(coef(modelfit.is20),fix,bacrepis20,init.vector,mytimes0,"day")
}

#4 influx rates from day 0 model
{
fix<-c(Fs = 0,Ps = 0, s1=0.18,If0=2.0038628  , Ip0=0.5049933, If=5.0083361 , Ip=0.9894203  , 
       If2=19.4401025, Ip2=0.3495042,
       If3=329.2480605, Ip3=0.7846695)

fit <- c(del1 = 1.13370367   , del2= 0.05126102   , del3=0.20285011    , 
         rep1= 1.47610914   , rep2=0.07022283  , rep3= 0.23505741 )


modelfit.is40<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepis40,dataset=(dataset.independent.d0),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.is40<-ode.solve(coef(modelfit.is40),fix,bacrepis40,init.vector,mytimes0,"day")
}

####################################################################
#influx only- 4 rates from day 0

#functions for influx only models
{
#resids and odesolve same as independent model

#ode function
bacrepio4 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFs.dt <- If0
      dPs.dt <- Ip0
    }
    if (time >13 & time < 26 ){
      dFs.dt <- If
      dPs.dt <- Ip
    }
    if (time >26 & time < 69 ){
      dFs.dt <- If2
      dPs.dt <- Ip2
    }
    if (time >69){
      dFs.dt <- If3
      dPs.dt <- Ip3
    }   
    return(list(c(dFs.dt, dPs.dt)))})
}
}

#starting conditions(same as influx model)
mytimes0<-seq(0,111,0.1)
init.vector <- c("Fs", "Ps")

#influx only- 4 rates from day 0 model
{
fix<-c(Fs = 0,Ps = 0)
fit <- c(If0=2, Ip0=0.5, If=20 , Ip=3 , 
         If2=500, Ip2=1,
         If3=1200, Ip3=2)


modelfit.io4<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepio4,dataset=(dataset.independent.d0),modelcols=c(2,3),datacols=c(2,3),init.vector=init.vector)
solmodfit.io4<-ode.solve(coef(modelfit.io4),fix,bacrepio4,init.vector,mytimes0,"day")
}

######################################################################
#migration models fitting lung and spleen data from day 0
#spleen replication/death starting at day 13

#migration model functions
{
  #resids same as indepedent model
  
  #need a new odesolve with lung adn spleen counts
ode.solve2=function(fit, fix, modelname, init.vector, timepoints, timename){ 
    allpars <- c(fit, fix)
    init.conds <- c()
    for(i in 1:length(init.vector)){
      init.conds <- c(init.conds, allpars[init.vector[[i]]])
    }
    soln=as.data.frame(lsoda(y=init.conds, times=timepoints, func=modelname, parms=allpars))
    soln$tl<-soln[,2] + soln[,3]  #count total bacteria lung
    soln$ts<-soln[,4]+soln[,5]  #count total bacteria spleen
    soln$ppl<-soln$Pl/soln$tl #percent plasmid lung
    soln$pps<-soln$Ps/soln$ts #percent plasmid spleen
    if (is.null(timename)) timename="Days" 
    names(soln)[1]=timename
    return(soln)
  }

#constant migration rate all time periods
bacrepm01 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFl.dt<- rep1l*Fl +s1*rep1l*Pl-del1l*Fl-mf*Fl
      dPl.dt<- rep1l*Pl - s1*rep1l*Pl-del1l*Pl-mp*Pl
      dFs.dt <- mf*Fl
      dPs.dt <- mp*Pl
    }
    if (time >13 & time < 26 ){
      dFl.dt<- rep2l*Fl +s1*rep2l*Pl-del2l*Fl-mf*Fl
      dPl.dt<- rep2l*Pl - s1*rep2l*Pl-del2l*Pl-mp*Pl
      dFs.dt <- rep1s*Fs + s1*rep1s*Ps - del1s*Fs+mf*Fl
      dPs.dt <- rep1s*Ps - s1*rep1s*Ps - del1s*Ps+mp*Pl
    }
    if (time >26 & time < 69 ){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp*Pl
      dFs.dt <- rep2s*Fs + s1*rep2s*Ps - del2s*Fs+mf*Fl
      dPs.dt <- rep2s*Ps - s1*rep2s*Ps -del2s*Ps+mp*Pl
    }
    if (time >69){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp*Pl
      dFs.dt <- rep3s*Fs + s1*rep3s*Ps - del3s*Fs+mf*Fl
      dPs.dt <- rep3s*Ps - s1*rep3s*Ps -del3s*Ps+mp*Pl
    }   
    return(list(c(dFl.dt, dPl.dt, dFs.dt, dPs.dt)))})
}

#2 migration rate function
bacrepm02 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFl.dt<- rep1l*Fl +s1*rep1l*Pl-del1l*Fl-mf*Fl
      dPl.dt<- rep1l*Pl - s1*rep1l*Pl-del1l*Pl-mp*Pl
      dFs.dt <- mf*Fl
      dPs.dt <- mp*Pl
    }
    if (time >13 & time < 26 ){
      dFl.dt<- rep2l*Fl +s1*rep2l*Pl-del2l*Fl-mf*Fl
      dPl.dt<- rep2l*Pl - s1*rep2l*Pl-del2l*Pl-mp*Pl
      dFs.dt <- rep1s*Fs + s1*rep1s*Ps - del1s*Fs+mf*Fl
      dPs.dt <- rep1s*Ps - s1*rep1s*Ps - del1s*Ps+mp*Pl
    }
    if (time >26 & time < 69 ){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf2*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp2*Pl
      dFs.dt <- rep2s*Fs + s1*rep2s*Ps - del2s*Fs+mf2*Fl
      dPs.dt <- rep2s*Ps - s1*rep2s*Ps -del2s*Ps+mp2*Pl
    }
    if (time >69){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf2*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp2*Pl
      dFs.dt <- rep3s*Fs + s1*rep3s*Ps - del3s*Fs+mf2*Fl
      dPs.dt <- rep3s*Ps - s1*rep3s*Ps -del3s*Ps+mp2*Pl
    }   
    return(list(c(dFl.dt, dPl.dt, dFs.dt, dPs.dt)))})
}

#3 migration rate function
bacrepm03.1 <- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFl.dt<- rep1l*Fl +s1*rep1l*Pl-del1l*Fl-mf0*Fl
      dPl.dt<- rep1l*Pl - s1*rep1l*Pl-del1l*Pl-mp0*Pl
      dFs.dt <- mf0*Fl
      dPs.dt <- mp0*Pl
    }
    if (time >13 & time < 26 ){
      dFl.dt<- rep2l*Fl +s1*rep2l*Pl-del2l*Fl-mf*Fl
      dPl.dt<- rep2l*Pl - s1*rep2l*Pl-del2l*Pl-mp*Pl
      dFs.dt <- rep1s*Fs + s1*rep1s*Ps - del1s*Fs+mf*Fl
      dPs.dt <- rep1s*Ps - s1*rep1s*Ps - del1s*Ps+mp*Pl
    }
    if (time >26 & time < 69 ){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf2*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp2*Pl
      dFs.dt <- rep2s*Fs + s1*rep2s*Ps - del2s*Fs+mf2*Fl
      dPs.dt <- rep2s*Ps - s1*rep2s*Ps -del2s*Ps+mp2*Pl
    }
    if (time >69){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf2*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp2*Pl
      dFs.dt <- rep3s*Fs + s1*rep3s*Ps - del3s*Fs+mf2*Fl
      dPs.dt <- rep3s*Ps - s1*rep3s*Ps -del3s*Ps+mp2*Pl
    }   
    return(list(c(dFl.dt, dPl.dt, dFs.dt, dPs.dt)))})
}

#4 migration rate function
bacrepm04<- function(time, state, parms){
  with(as.list(c(state,parms)),{
    if (time <= 13) {
      dFl.dt<- rep1l*Fl +s1*rep1l*Pl-del1l*Fl-mf0*Fl
      dPl.dt<- rep1l*Pl - s1*rep1l*Pl-del1l*Pl-mp0*Pl
      dFs.dt <- mf0*Fl
      dPs.dt <- mp0*Pl
    }
    if (time >13 & time < 26 ){
      dFl.dt<- rep2l*Fl +s1*rep2l*Pl-del2l*Fl-mf*Fl
      dPl.dt<- rep2l*Pl - s1*rep2l*Pl-del2l*Pl-mp*Pl
      dFs.dt <- rep1s*Fs + s1*rep1s*Ps - del1s*Fs+mf*Fl
      dPs.dt <- rep1s*Ps - s1*rep1s*Ps - del1s*Ps+mp*Pl
    }
    if (time >26 & time < 69 ){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf2*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp2*Pl
      dFs.dt <- rep2s*Fs + s1*rep2s*Ps - del2s*Fs+mf2*Fl
      dPs.dt <- rep2s*Ps - s1*rep2s*Ps -del2s*Ps+mp2*Pl
    }
    if (time >69){
      dFl.dt<- rep3l*Fl +s1*rep3l*Pl-del2l*Fl-mf3*Fl
      dPl.dt<- rep3l*Pl - s1*rep3l*Pl-del2l*Pl-mp3*Pl
      dFs.dt <- rep3s*Fs + s1*rep3s*Ps - del3s*Fs+mf3*Fl
      dPs.dt <- rep3s*Ps - s1*rep3s*Ps -del3s*Ps+mp3*Pl
    }   
    return(list(c(dFl.dt, dPl.dt, dFs.dt, dPs.dt)))})
}

}

#conditions
init.vector <- c("Fl", "Pl", "Fs", "Ps")

#model fits
#migration rates fitted then fixed for all 
{
#one migration rate model fit
fix <- c(Fl=67, Pl=200, Fs = 0, Ps = 0, s1=0.18,
         mf=0.002533678 , mp=0.001367941
)
fit <- c(  del1s= 0.4546721    , del2s= 0.5525795   , del3s=0.7886913    , 
           rep1s= 0.3980166    , rep2s= 0.0000100 , rep3s= 0.7733606  ,
           del1l= 0.4395634    , del2l= 0.1132162     , 
           rep1l= 0.8201337    , rep2l= 0.2660779     , rep3l= 0.1379681  
)


modelfit.m01<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepm01,dataset=(dataset.both),
                     modelcols=c(2,3,4,5),datacols=c(2,3,4,5),init.vector=init.vector, lower=c(rep(0,length(fit))))
solmodfit.m01<-ode.solve2(coef(modelfit.m01),fix,bacrepm01,init.vector,mytimes0,"hour")

#two migration rate fitting
fix <- c(Fl=67, Pl=200, Fs = 0, Ps = 0, s1=0.18,
         mf=0.0019745061    , mp=0.0020974351   ,
         mf2=0.0042508071   , mp2=0.0006152705 
)

fit <- c( del1l= 0.63000217   , del2l= 0.05037954    , 
          rep1l= 0.95230550    , rep2l= 0.35689430   , rep3l= 0.05779962  ,
          del1s= 0.48446582     , del2s=0.15057246    , del3s=1.02407188     , 
          rep1s= 0.12509146     , rep2s= 0.00008624  , rep3s= 0.93896781 )


modelfit.m02<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepm02,dataset=(dataset.both),
                     modelcols=c(2,3,4,5),datacols=c(2,3,4,5),init.vector=init.vector, lower=c(rep(0,length(fit))))
solmodfit.m02<-ode.solve2(coef(modelfit.m02),fix,bacrepm02,init.vector,mytimes0,"day")


#3 migration rate model fit
fix <- c(Fl=67, Pl=200, Fs = 0, Ps = 0, s1=0.18,
         mf0=0.0002958677 , mp0=0.0052754090 ,
         mf=0.0046876933     , mp=0.0014301969    ,
         mf2=0.0061780689    , mp2=0.0006395014 
)

fit <- c(  del1l= 0.63009877    , del2l= 0.07413043     , 
           rep1l= 0.95172435     , rep2l= 0.37280870    , rep3l= 0.08240063   ,
           del1s= 0.49390455      , del2s=0.15661676     , del3s=0.93418873      , 
           rep1s= 0.12290471      , rep2s= 0.00000006   , rep3s= 0.89412269 )

modelfit.m03.1<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepm03.1,dataset=(dataset.both),
                       modelcols=c(2,3,4,5),datacols=c(2,3,4,5),init.vector=init.vector, lower=c(rep(0,length(fit))))
solmodfit.m03.1<-ode.solve2(coef(modelfit.m03.1),fix,bacrepm03.1,init.vector,mytimes0,"day")

#4 migration rate model fit
fix <- c(Fl=67, Pl=200, Fs = 0, Ps = 0, s1=0.18,
         mf0=0.0003108439  , mp0=0.0048276895  ,
         mf=0.0056706353 , mp=0.0012732144     ,
         mf2=0.0042046571, mp2=0.0006395144 ,
         mf3= 0.0047763864  , mp3= 0.0009930757 
)

fit <- c( del1l= 0.62610837 , del2l= 0.09893231, 
          rep1l= 0.95260463, rep2l= 0.39975506 , rep3l= 0.10630213  ,
          del1s= 0.50461278, del2s=0.13031000 , del3s=0.68215110   , 
          rep1s= 0.09288570  , rep2s= 0.00000007, rep3s= 0.64439897  
)

modelfit.m04<-modFit(method = "Marq", f=resids,p=fit,fix=fix,modelname=bacrepm04,dataset=(dataset.both),
                     modelcols=c(2,3,4,5),datacols=c(2,3,4,5),init.vector=init.vector, lower=c(rep(0,length(fit))))
solmodfit.m04<-ode.solve2(coef(modelfit.m04),fix,bacrepm04,init.vector,mytimes0,"day")
}

#######################################################################
#######################################################################
#######################################################################
##Stochastic simulation models
#comparing GillespieSSA and adaptivetau
#Ran using params from 3-migration-rate model of dissemination
#simulated plasmid-bearing and plasmid-free bacteria counts
#compared sims to ODE model

########################################################################
#Gillespie

#initial conditions
init.fp<-c(Fl=67 , Pl=200, Fs=0, Ps=0)

pars.fp<-c(rL1=m03par$rep1l,rL2=m03par$rep2l,rL3=m03par$rep3l,
           dL1=m03par$del1l,dL2=m03par$del2l,
           rS1=m03par$rep1s, rS2=m03par$rep2s, rS3=m03par$rep3s,
           dS1=m03par$del1s, dS2=m03par$del2s, dS3=m03par$del3s,
           mf1=0.000295868, mp1=0.005275409,  
           mf2=0.004687693, mp2=0.001430197, 
           mf3=0.006178069, mp3=0.000639501,
           s=0.18, #seg probability
           t1=13, t2=26, t3=69, t4=111) 

tmax<-111

#defining a function to run one time course of Gillespie
onerun.m3fp <- function (params=params, tmax=tmax){
  reactions<-list(
    reaction("(time < t1 ? rL1 :
                time < t2 ? rL2 : rL3) * Fl", c(Fl = +1), name="lung_rep1_free"),
    reaction("(time < t1 ? rL1 :
                time < t2 ? rL2 : rL3) * Pl", c(Pl = +1), name="lung_rep1_plas"),
    
    reaction("(time < t1 ? rL1 :
                time < t2 ? rL2 : rL3) *s* Pl", c(Pl = -1, Fl=+1), name="lung_plas_seg"),
    
    reaction("(time < t1 ? dL1 : dL2) * Fl", c(Fl= -1), name="lung_death_free"), 
    reaction("(time < t1 ? dL1 : dL2) * Pl", c(Pl= -1), name="lung_death_plas"),
    
    reaction("(time < t1 ? mf1 : 
               time < t2? mf2: mf3) * Fl"  , c(Fl = -1, Fs = +1), name="mig_free" ),
    reaction("(time < t1 ? mp1 : 
               time < t2? mp2: mp3) * Pl"  , c(Pl = -1, Ps = +1), name="mig_plas" ),
    
    
    reaction("(time < t1 ? rS1 : 
               time < t2 ? rS2 : rS3) * Fs" , c(Fs= +1), name="spleen_rep_free" ),
    reaction("(time < t1 ? rS1 : 
               time < t2 ? rS2 : rS3) * Ps" , c(Ps= +1), name="spleen_rep_plas" ),
    
    reaction("(time < t1 ? rS1 : 
               time < t2 ? rS2 : rS3) *s* Ps" , c(Ps= -1, Pl=+1), name="spleen_plas_seg" ),
    
    reaction("(time < t1 ? dS1 :
               time < t2 ? dS2 : dS3) * Fs" , c(Fs= -1), name="spleen_death_free" ),
    reaction("(time < t1 ? dS1 :
               time < t2 ? dS2 : dS3) * Ps" , c(Ps= -1), name="spleen_death_plas" )
  )
  
  sim = ssa(initial_state = init.fp,
            reactions = reactions,
            params = pars.fp,
            tmax,
            method=ssa_etl(tau=0.3),
            census_interval = 0.1,
            verbose=F)
  sim
} #onerun function with gillespie ETL 

#get solutions for one run of Gillespie
r<-(onerun.m3fp(pars.fp, tmax))

#Bootstrap Gillespie ETL
{
nboot<-200
set.seed(412)
tmax<-111

allruns<-c()
for (n in 1:nboot) {
  one<-onerun.m3fp(pars.fp,tmax)
  run<-data.frame(time=unlist(one$time),  #unlist output from ssa 
                  Fl=unlist(one$state[,1]), 
                  Pl=unlist(one$state[,2]),
                  Fs=unlist(one$state[,3]), 
                  Ps=unlist(one$state[,4]),
                  sim.number=rep(n,length(one$time)))
  ifelse(n==1,allruns<-run, allruns<-rbind(allruns,run))
}

#add total bacteria counts
allruns$Tl<-allruns$Fl+allruns$Pl
allruns$Ts<-allruns$Fs+allruns$Ps

#Get averages for plotting
#average total lung count 
avg.Tl<-allruns %>%
  group_by(sim.number) %>%
  group_by(time=floor(time)) %>%
  summarise(avg=mean(Tl))

#average total spleen count
avg.Ts<-allruns.m3fp2 %>%
  group_by(sim.number) %>%
  group_by(time=floor(time)) %>%
  summarise(avg=mean(Ts))

#average plasmid lung 
avg.Pl<-allruns %>%
  group_by(sim.number) %>%
  group_by(time=floor(time)) %>%
  summarise(avg=mean(Pl))

#average plasmid spleen 
avg.Ps<-allruns %>%
  group_by(sim.number) %>%
  group_by(time=floor(time)) %>%
  summarise(avg=mean(Ps))
}

##########################################################################
#Adaptivetau simulations

#adaptive tau function
onerunav <- function (params=params, tmax=tmax){
  rateFun<- function(x, params, t) {
    return( c(rL(params,t) * x[c("Fl","Pl")],        #lung rep           
              dL(params,t) * x[c("Fl", "Pl")],       #lung death
              params$s     *rL(params,t)* x[c("Pl")],             #segregation lung
              mf(params,t) * x["Fl"],                #free mig
              mp(params,t) * x["Pl"],                #plasmid mig
              rS(params,t) * x[c("Fs", "Ps")],          #spleen rep
              dS(params,t) * x[c("Fs", "Ps")],          #spleen death 
              params$s     *rS(params,t)* x[c("Ps")])        #segregation spleen
            
    )
  }#rateFun
  
  transitions= list(c(Fl = +1),          # trans 1:bacteria grow Lung 
                    c(Pl= +1),           
                    c(Fl = -1), 
                    c(Pl = -1),         
                    c(Fl = +1, Pl = -1),         # trans 5:segregation lung
                    c(Fl = -1, Fs = +1),         
                    c(Pl = -1, Ps = +1),         
                    c(Fs = +1), 
                    c(Ps = +1),         
                    c(Fs = -1),
                    c(Ps = -1),         
                    c(Fs = +1, Ps = -1)          # trans 12: segregation spleen
  ) 
  
  
  simResults = ssa.adaptivetau(init.values = init.fp,
                               transitions, rateFun,
                               params = params,
                               tf=tmax)
  simResults
} #onerun adaptive function

#run bootstraps with adaptivetau
{ 
  nboot<-200
  set.seed(412)
  tmax<-111
  
  init<-c(Fl=67, Pl=200, Fs=0, Ps=0)
  
  allruns.ad<-c()
  for (n in 1:nboot) {
    one<-as.data.frame(onerunav(params=as.list(pars.fp),tmax=tmax))
    one$sim.number<-rep(n,length(one$time))
    ifelse(n==1,allruns.ad32<-one, allruns.ad<-rbind(allruns.ad,one))
  }
  
  #add totals
  #add total counts
  allruns.ad$Tl<-allruns.ad$Fl+allruns.ad$Pl
  allruns.ad$Ts<-allruns.ad$Fs+allruns.ad$Ps
  
  #get bootstrap averages for plotting
  #average total lung count for fp4
  avg.Tl.ad<-allruns.a %>%
    group_by(sim.number) %>%
    group_by(time=floor(time)) %>%
    summarise(avg=mean(Tl))
  
  #average total spleen count for fp4
  avg.Ts.ad<-allruns.ad %>%
    group_by(sim.number) %>%
    group_by(time=floor(time)) %>%
    summarise(avg=mean(Ts))
  
  #average plasmid lung for fp4
  avg.Pl.ad<-allruns.ad %>%
    group_by(sim.number) %>%
    group_by(time=floor(time)) %>%
    summarise(avg=mean(Pl))
  
  #average plasmid spleen for fp4
  avg.Ps.ad<-allruns.ad %>%
    group_by(sim.number) %>%
    group_by(time=floor(time)) %>%
    summarise(avg=mean(Ps))
  
}

#########################################################################
#Run ODE model with same params as simulations
{
#time dependent params
{
  t1<-13
  t2<-26
  t3<-69
  t4<-111
  
  rL<-function(params, t) {
    return (ifelse(t<t1,params$rL1,
                   ifelse(t<t2,params$rL2,params$rL3)))
  }
  
  dL<-function(params, t) {
    return (ifelse(t<t1,params$dL1,params$dL2))
  }
  
  rS<-function(params, t) {
    return (ifelse(t<t1,0,
                   ifelse(t<t2,params$rS1,
                          ifelse(t<t3,params$rS2, params$rS3))))
  }
  
  #no replication or death in spleen during first time period
  
  dS<-function(params, t) {
    return (ifelse(t<t1,0,
                   ifelse(t<t2, params$dS1,
                          ifelse(t<t3,params$dS2,params$dS3))))
  }
  
  mf<-function(params, t) {
    return (ifelse(t<t1, params$mf1, 
                   ifelse(t<t2, params$mf2,params$mf3)))
  } 
  mp<-function(params, t) {
    return (ifelse(t<t1, params$mp1, 
                   ifelse(t<t2, params$mp2,params$mp3)))
  } 
  
} # parameter functions (don't need for gillespie set up but need for ODE)

tmax<-111
times<-seq(0,tmax,1)

#ODE model function
ODEmodel.fp3 <- function(t, state, parms){
  with(as.list(c(t,state,parms)),
       {
         dLf = (rL(parms,t)-dL(parms,t) - mf(parms,t))*Fl +s*rL(parms,t)*Pl
         dLp = (rL(parms, t)-dL(parms,t) - mp(parms,t))*Pl-s*rL(parms,t)*Pl
         dSf = (rS(parms,t)-dS(parms,t))*Fs + mf(parms,t)*Fl + s*rS(parms,t)*Ps
         dSp = (rS(parms,t)-dS(parms,t))*Ps + mp(parms,t)*Pl -s*rS(parms,t)*Ps
         return(list(c(dLf,dLp, dSf, dSp)))
       }
  )
}   

#running ODE model(initial params from literature)
soln.fp3 <- as.data.frame(ode(y=c(Fl=67, Pl=200, Fs=0, Ps=0), 
                              times=times, 
                              func=ODEmodel.fp3, 
                              method='lsoda',
                              parms=as.list(pars.fp), atol = 1e-17, rtol = 1e-12)
)


#make total counts from ODE results
soln.fp3$Tl<-soln.fp3$Pl+soln.fp3$Fl
soln.fp3$Ts<-soln.fp3$Ps+soln.fp3$Fs
}


#######################################################################
#Plotting simulations vs. ODE model

#Gillespie vs. ODE plots
{
  #lung total
  ggplot()+
    scale_x_continuous(breaks = c(1,13,26, 111))+
    scale_y_continuous(label=math_format())+
    theme_bw()+
    theme(legend.position=c(0.8,0.2), legend.background = element_blank(),
          legend.key=element_rect(fill=NA),
          legend.text = element_text(size=13),
          axis.title = element_text(size=13), axis.text=element_text(size=13))+
    geom_line(data=allruns, aes(time, log10(Tl), group=sim.number, col="a"))+
    geom_line(data=avg.Tl, aes(time, log10(avg), col="b"))+
    geom_point(data=dataset.both, aes(day, log10(Tl), col="c"))+
    geom_line(data=soln.fp3, aes(time, log10(Tl), col="d"), lty="dashed")+
    geom_segment(aes(x=c(-1,11,24,67,109),y=averages.log$mean.log.Tl, 
                     xend=c(3,15,28,71,113),yend=averages.log$mean.log.Tl), color="red", lwd=1)+
    labs(x="Days after infection", y="Mtb Total (lung)", title="Plasmid/Free Sims (Starting CFU=267)")+
    scale_color_manual(name="", 
                       values=c("a"="gray", "b"="navy", "c"="black", "d"="darkgreen"),
                       guide=guide_legend(override.aes = 
                                            list(shape=c(NA,NA,16, NA), 
                                                 linetype=c("solid", "solid", "blank", "dashed"))),
                       labels=c("simulations(n=200)", "sim avg", "raw data", "ODE mod"))  
  
  #lung plasmid
  ggplot()+
    scale_x_continuous(breaks = c(1,13,26, 111))+
    scale_y_continuous(label=math_format())+
    theme_bw()+
    theme(legend.position=c(0.8,0.2), legend.background = element_blank(),
          legend.key=element_rect(fill=NA),
          legend.text = element_text(size=13),
          axis.title = element_text(size=13), axis.text=element_text(size=13))+
    geom_line(data=allruns, aes(time, log10(Pl), group=sim.number, col="a"))+
    geom_line(data=avg.Pl, aes(time, log10(avg), col="b"))+
    geom_point(data=dataset.both, aes(day, log10(Pl), col="c"))+
    geom_line(data=soln.fp3, aes(time, log10(Pl), col="d"), lty="dashed")+
    geom_segment(aes(x=c(-1,11,24,67,109),y=averages.log$mean.log.Pl, 
                     xend=c(3,15,28,71,113),yend=averages.log$mean.log.Pl), color="red", lwd=1)+
    labs(x="Days after infection", y="Mtb Plasmid(lung)", title="Plasmid/Free Sims(Starting CFU=267)")+
    scale_color_manual(name="", 
                       values=c("a"="gray", "b"="navy", "c"="black", "d"="darkgreen"),
                       guide=guide_legend(override.aes = 
                                            list(shape=c(NA,NA,16, NA), 
                                                 linetype=c("solid", "solid", "blank", "dashed"))),
                       labels=c("simulations(n=200)", "sim avg", "raw data", "ODE mod")) 
  
  #spleen total
  ggplot()+
    scale_x_continuous(breaks = c(1,13,26, 111))+
    scale_y_continuous(label=math_format())+
    theme_bw()+
    theme(legend.position=c(0.8,0.2), legend.background = element_blank(),
          legend.key=element_rect(fill=NA),
          legend.text = element_text(size=13),
          axis.title = element_text(size=13), axis.text=element_text(size=13))+
    geom_line(data=allruns, aes(time, log10(Ts), group=sim.number, col="a"))+
    geom_line(data=avg.Ts, aes(time, log10(avg), col="b"))+
    geom_point(data=dataset.both, aes(day, log10(Ts), col="c"))+
    geom_line(data=soln.fp3, aes(time, log10(Ts), col="d"), lty="dashed")+
    geom_segment(aes(x=c(-1,11,24,67,109),y=averages.log$mean.log.Ts, 
                     xend=c(3,15,28,71,113),yend=averages.log$mean.log.Ts), color="red", lwd=1)+
    labs(x="Days after infection", y="Mtb Total(spleen)", title="Plasmid/Free Sims(Starting CFU=267)")+
    scale_color_manual(name="", 
                       values=c("a"="gray", "b"="navy", "c"="black", "d"="darkgreen"),
                       guide=guide_legend(override.aes = 
                                            list(shape=c(NA,NA,16, NA), 
                                                 linetype=c("solid", "solid", "blank", "dashed"))),
                       labels=c("simulations(n=200)", "sim avg", "raw data", "ODE mod")) 
  
  #spleen plasmid
  ggplot()+
    scale_x_continuous(breaks = c(1,13,26, 111))+
    scale_y_continuous(label=math_format())+
    theme_bw()+
    theme(legend.position=c(0.8,0.2), legend.background = element_blank(),
          legend.key=element_rect(fill=NA),
          legend.text = element_text(size=13),
          axis.title = element_text(size=13), axis.text=element_text(size=13))+
    geom_line(data=allruns, aes(time, log10(Ps), group=sim.number, col="a"))+
    geom_line(data=avg.Ps, aes(time, log10(avg), col="b"))+
    geom_point(data=dataset.both, aes(day, log10(Ps), col="c"))+
    geom_line(data=soln.fp3, aes(time, log10(Ps), col="d"), lty="dashed")+
    geom_segment(aes(x=c(-1,11,24,67,109),y=averages.log$mean.log.Ps, 
                     xend=c(3,15,28,71,113),yend=averages.log$mean.log.Ps), color="red", lwd=1)+
    labs(x="Days after infection", y="Mtb Plasmid(spleen)", title="Plasmid/Free Sims(Starting CFU=267)")+
    scale_color_manual(name="", 
                       values=c("a"="gray", "b"="navy", "c"="black", "d"="darkgreen"),
                       guide=guide_legend(override.aes = 
                                            list(shape=c(NA,NA,16, NA), 
                                                 linetype=c("solid", "solid", "blank", "dashed"))),
                       labels=c("simulations(n=200)", "sim avg", "raw data", "ODE mod")) 
  
}

#Adaptivetau vs. ODE plots(spleen)
{
#spleen total
ggplot()+
  scale_x_continuous(breaks = c(1,13,26, 111))+
  scale_y_continuous(label=math_format())+
  theme_bw()+
  theme(legend.position=c(0.8,0.2), legend.background = element_blank(),
        legend.key=element_rect(fill=NA),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13), axis.text=element_text(size=13))+
  geom_line(data=allruns.ad, aes(time, log10(Ts), group=sim.number, col="a"))+
  geom_line(data=avg.Ts.ad, aes(time, log10(avg), col="b"))+
  geom_point(data=dataset.both, aes(day, log10(Ts), col="c"))+
  geom_line(data=soln.fp3, aes(time, log10(Ts), col="d"), lty="dashed")+
  geom_segment(aes(x=c(-1,11,24,67,109),y=averages.log$mean.log.Ts, 
                   xend=c(3,15,28,71,113),yend=averages.log$mean.log.Ts), color="red", lwd=1)+
  labs(x="Days after infection", y="Mtb Total(spleen)", title="Plasmid/Free Sims(Starting CFU=267)")+
  scale_color_manual(name="", 
                     values=c("a"="gray", "b"="navy", "c"="black", "d"="darkgreen"),
                     guide=guide_legend(override.aes = 
                                          list(shape=c(NA,NA,16, NA), 
                                               linetype=c("solid", "solid", "blank", "dashed"))),
                     labels=c("simulations(n=200)", "sim avg", "raw data", "ODE mod")) 

#spleen plasmid
ggplot()+
  scale_x_continuous(breaks = c(1,13,26, 111))+
  scale_y_continuous(label=math_format())+
  theme_bw()+
  theme(legend.position=c(0.8,0.2), legend.background = element_blank(),
        legend.key=element_rect(fill=NA),
        legend.text = element_text(size=13),
        axis.title = element_text(size=13), axis.text=element_text(size=13))+
  geom_line(data=allruns.ad, aes(time, log10(Ps), group=sim.number, col="a"))+
  geom_line(data=avg.Ps.ad, aes(time, log10(avg), col="b"))+
  geom_point(data=dataset.both, aes(day, log10(Ps), col="c"))+
  geom_line(data=soln.fp3, aes(time, log10(Ps), col="d"), lty="dashed")+
  geom_segment(aes(x=c(-1,11,24,67,109),y=averages.log$mean.log.Ps, 
                   xend=c(3,15,28,71,113),yend=averages.log$mean.log.Ps), color="red", lwd=1)+
  labs(x="Days after infection", y="Mtb Plasmid(spleen)", title="Plasmid/Free Sims(Starting CFU=267)")+
  scale_color_manual(name="", 
                     values=c("a"="gray", "b"="navy", "c"="black", "d"="darkgreen"),
                     guide=guide_legend(override.aes = 
                                          list(shape=c(NA,NA,16, NA), 
                                               linetype=c("solid", "solid", "blank", "dashed"))),
                     labels=c("simulations(n=200)", "sim avg", "raw data", "ODE mod")) 

}
