#######################################################################
####                  Reproducibility code for                     ####
# Modeling radiocarbon dynamics in soils: SoilR version 1.1           #
#######################################################################
#################################### Prepared by C.A. Sierra ##########
####################################################### April, 2014  #

library(SoilR)
library(FME)
library(forecast)

########## Example 1

#example("ThreepFeedbackModel14")

years=seq(1901,2009,by=0.5)
LitterInput=100
k1=1/2; k2=1/10; k3=1/50
a21=0.9*k1
a12=0.4*k2
a32=0.4*k2
a23=0.7*k3

Feedback=ThreepFeedbackModel14(
  t=years,
  ks=c(k1=k1, k2=k2, k3=k3),
  C0=c(100,500,1000),
  F0_Delta14C=c(0,0,0),
  In=LitterInput,
  a21=a21,
  a12=a12,
  a32=a32,
  a23=a23,
  inputFc=C14Atm_NH
)
F.R14m=getF14R(Feedback)
F.C14m=getF14C(Feedback)
F.C14t=getF14(Feedback)

Series=ThreepSeriesModel14(
  t=years,
  ks=c(k1=k1, k2=k2, k3=k3),
  C0=c(100,500,1000),
  F0_Delta14C=c(0,0,0),
  In=LitterInput,
  a21=a21,
  a32=a32,
  inputFc=C14Atm_NH
)
S.R14m=getF14R(Series)
S.C14m=getF14C(Series)
S.C14t=getF14(Series)

Parallel=ThreepParallelModel14(
  t=years,
  ks=c(k1=k1, k2=k2, k3=k3),
  C0=c(100,500,1000),
  F0_Delta14C=c(0,0,0),
  In=LitterInput,
  gam1=0.6,
  gam2=0.2,
  inputFc=C14Atm_NH
)
P.R14m=getF14R(Parallel)
P.C14m=getF14C(Parallel)
P.C14t=getF14(Parallel)

par(mfrow=c(3,2), mar=c(5,5,1,1),cex.lab=1.25)
plot(
  C14Atm_NH,
  type="l",
  xlab="Year",
  ylab=expression(paste(Delta^14,"C ","(\u2030)")),
  xlim=c(1940,2010)
) 
lines(years, P.C14t[,1], col=4)
lines(years, P.C14t[,2],col=4,lwd=2)
lines(years, P.C14t[,3],col=4,lwd=3)
legend(
  "topright",
  c("Atmosphere", "Pool 1", "Pool 2", "Pool 3"),
  lty=rep(1,4),
  col=c(1,4,4,4),
  lwd=c(1,1,2,3),
  bty="n"
)

plot(C14Atm_NH,type="l",xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlim=c(1940,2010)) 
lines(years,P.C14m,col=4)
lines(years,P.R14m,col=2)
legend("topright",c("Atmosphere","Bulk SOM", "Respired C"),
       lty=c(1,1,1), col=c(1,4,2),bty="n")

plot(C14Atm_NH,type="l",xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlim=c(1940,2010)) 
lines(years, S.C14t[,1], col=4)
lines(years, S.C14t[,2],col=4,lwd=2)
lines(years, S.C14t[,3],col=4,lwd=3)
legend("topright",c("Atmosphere", "Pool 1", "Pool 2", "Pool 3"),
       lty=rep(1,4),col=c(1,4,4,4),lwd=c(1,1,2,3),bty="n")

plot(C14Atm_NH,type="l",xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlim=c(1940,2010)) 
lines(years,S.C14m,col=4)
lines(years,S.R14m,col=2)
legend("topright",c("Atmosphere","Bulk SOM", "Respired C"),
       lty=c(1,1,1), col=c(1,4,2),bty="n")

plot(C14Atm_NH,type="l",xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlim=c(1940,2010)) 
lines(years, F.C14t[,1], col=4)
lines(years, F.C14t[,2],col=4,lwd=2)
lines(years, F.C14t[,3],col=4,lwd=3)
legend("topright",c("Atmosphere", "Pool 1", "Pool 2", "Pool 3"),
       lty=rep(1,4),col=c(1,4,4,4),lwd=c(1,1,2,3),bty="n")

plot(C14Atm_NH,type="l",xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),xlim=c(1940,2010)) 
lines(years,F.C14m,col=4)
lines(years,F.R14m,col=2)
legend("topright",c("Atmosphere","Bulk SOM", "Respired C"),
       lty=c(1,1,1), col=c(1,4,2),bty="n")


par(mfrow=c(1,1))

#Mean transit times

Af=new("ConstLinDecompOp",
        matrix(c(-k1,a12,0,
                 a21,-k2,a23,
                 0,a32,-k3),byrow=TRUE,3,3)
       )


TTf=getMeanTransitTime(Af, c(1,0,0))

As=new("ConstLinDecompOp",
       matrix(c(-k1,0,0,
                a21,-k2,0,
                0,a32,-k3),byrow=TRUE,3,3)
)

TTs=getMeanTransitTime(As, c(1,0,0))

Ap=new("ConstLinDecompOp",
       matrix(c(-k1,0,0,
                0,-k2,0,
                0,0,-k3),byrow=TRUE,3,3)
)


TTp=getMeanTransitTime(Ap, rep(1/3,3))

TT=data.frame(ModelStructure=c("Parallel","Series","Feedback"),TransitTime=round(c(TTp,TTs,TTf),1))
TT

########### Example 2

#We need to create the components of the model first.

#First, we define the points in time to run the model
time=C14Atm_NH$YEAR
t_start=min(time)
t_end=max(time)

#The input fluxes to the two-pools are created by this command
inputFluxes=BoundInFlux(
    function(t0){matrix(nrow=3,ncol=1,c(270,150,0))},
    t_start,
    t_end
)

#The initial amount of carbon is created by aggregating the organic and mineral pools from Sierra et al. (2012, Biogeosciences 9: 3013)
C0=c(390,220+390+1376,90+1800+560) 

# This function creates a model object using a set of parameters and returns the 
# Delta14C value of the respired carbon
Fc=BoundFc(C14Atm_NH,lag=0,format="Delta14C")
Mod1<-function(ks,pass=TRUE){ #teneille note - see where ks is used in At, this is later used for modFit
  At=new("BoundLinDecompOp",
         t_start,
         t_end,
         function(t0){
           matrix(nrow=3,ncol=3,byrow=TRUE,c(-ks[1],0,0,
                                             ks[4],-ks[2],0,
                                             ks[5],0,-ks[3]))
         }
  ) 
  mod=GeneralModel_14(time,At,C0,inputFluxes,initialValF=ConstFc(rep(0,3),"Delta14C"),Fc,pass=TRUE) 
  R14t=getF14R(mod)
  return(data.frame(time=time,R14t=R14t))
}

#The observed data needs to be orginazed in a dataframe of the form:
DataR14t=cbind(time=HarvardForest14CO2[,1],R14t=HarvardForest14CO2[,2],sd=sd(HarvardForest14CO2[,2]))


#Now we use the observed values and the model function to create a cost function
R14tCost <- function(pars){
  R14t <- Mod1(pars)
  return(modCost(model=R14t,obs=DataR14t,err="sd"))
}

#Fit the model to the observed data given some initial value for the parameters
#teneille note - how are initial parameters chosen?
Fit <- modFit(f=R14tCost,p=c(0.5,0.9,0.1,0.1,0.1))

# Now we run an MCMC using the variance and covariance results from the previous optimization
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled # The covariance matrix can be used for the jump, but wasn't used in this example.
MCMC <- modMCMC(f=R14tCost, p = Fit$par, niter = 2500, jump = NULL, var0 = var0, wvar0 = 0, #updatecov = 50,
                lower=c(0,0,0,0,0),upper=c(2,3,1,1,1))
summary(MCMC)
plot(MCMC)

sR=sensRange(func=Mod1, parInput=MCMC$par)

#Figure 4
par(mar=c(5,5,4,1))
plot(summary(sR),xlim=c(1950,2010),ylim=c(0,1000),xlab="Year",
     ylab=expression(paste(Delta^14,"C ","(\u2030)")),main="")
points(DataR14t,pch=20)
lines(C14Atm_NH,col=4)

#Figure 5
par(cex.axis=1.25)
pairs(MCMC,nsample=500)

##################################
# Example 3

years=seq(1966,2009,by=1/4) # A series of years by quarters
nz1=spline(Hua2013$NHZone1[,c(1,4)],xout=years) #Spline interpolation of the NH_Zone 1 dataset at a quaterly basis
nhz1=ts((nz1$y-1)*1000,start=1966,freq=4) #Transformation into a time-series object
m1=ets(nhz1) #Fits an exponential smoothing state space model to the time series
f1=forecast(m1,h=11*4) #Uses the fitted model to forecast 10 years into the future

par(mar=c(5,5,4,1))
plot(f1,main="",ylab=expression(paste(Delta^14,"C ","(\u2030)")), xlab="Years",
     xlim=c(2000,2020),ylim=c(-50,200))
abline(h=0,lty=3)






