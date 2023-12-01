library(mrdsalt)
#sim loglinear
par=c( 1.517424, 1.663175, -0.28, -0.5, 0.1311178)
W=10
N=1000
set.seed(1083821)
x=sim.loglinear(N,W,par=par)
#create data structure with 2 records for each observation
dd=make.design.data(x)
#fit model
results.1=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
            data=dd,method="bpi",meta.data=list(width=W),
            control=list(method="nlminb",debug=FALSE))
results.1$criterion
# %diff
100*(results.1$Nhat-N)/N
# unconditonal detection fcts
par(mfrow=c(3,2))
plot(results.1,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results.1,7:8,showlines=FALSE)

#fit model
results.2=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
            data=dd,method="loglinear",meta.data=list(width=W),
            control=list(method="nlminb",debug=FALSE))
results.2$criterion
# %diff
100*(results.2$Nhat-N)/N
# unconditonal detection fcts
par(mfrow=c(3,2))
plot(results.2,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results.2,7:8,showlines=FALSE)

#simbpi
x=sim.bpi(x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2),
                       distance=rep(runif(N/2,0,W),each=2)),
          par=par, pformula=~observer*distance,dformula=~-1+distance,PI=TRUE)
#fit model
results.3=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
            data=x,method="bpi",meta.data=list(width=W),
            control=list(method="nlminb",debug=FALSE))
results.3$criterion
# %diff
100*(results.3$Nhat-N)/N
# unconditional detection fcts
par(mfrow=c(3,2))
plot(results.3,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results.3,7:8,showlines=FALSE)

#' #fit model
results.4=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
            data=x,method="loglinear",meta.data=list(width=W),
            control=list(method="nlminb",debug=FALSE))
results.4$criterion
# %diff
100*(results.4$Nhat-N)/N
# unconditonal detection fcts
par(mfrow=c(3,2))
plot(results.4,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results.4,7:8,showlines=FALSE)



library(mrdsalt)
ptm <- proc.time()
reps=100
Nhat=matrix(NA,nrow=reps,ncol=4)
AIC=matrix(NA,nrow=reps,ncol=4)
convergence=matrix(NA,nrow=reps,ncol=4)
for(i in 1:reps)
{
cat("\n rep = ",i)
#sim loglinear
par=c( 1.5, 1.6, -0.5, -0.3, 0.1)
W=10
N=1000
x=sim.loglinear(N,W,par=par)
#create data structure with 2 records for each observation
dd=make.design.data(x)
#fit model
cat("\n     fit 1")
results.1=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
              data=dd,method="bpi",meta.data=list(width=W),
              control=list(method="nlminb",debug=FALSE))

#fit model
cat("\n     fit 2")
results.2=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
              data=dd,method="loglinear",meta.data=list(width=W),
              control=list(method="nlminb",debug=FALSE))

#simbpi
x=sim.bpi(x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2),
                       distance=rep(runif(N/2,0,W),each=2)),
          par=par, pformula=~observer*distance,dformula=~-1+distance,PI=TRUE)
#fit model
cat("\n     fit 3")
results.3=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
              data=x,method="bpi",meta.data=list(width=W),
              control=list(method="nlminb",debug=FALSE))


#' #fit model
cat("\n     fit 4")
results.4=ddf(mrmodel=list(pformula=~observer*distance,dformula=~-1+distance),
              data=x,method="loglinear",meta.data=list(width=W),
              control=list(method="nlminb",debug=FALSE))

Nhat[i,1]=results.1$Nhat
Nhat[i,2]=results.2$Nhat
Nhat[i,3]=results.3$Nhat
Nhat[i,4]=results.4$Nhat

AIC[i,1]=results.1$criterion
AIC[i,2]=results.2$criterion
AIC[i,3]=results.3$criterion
AIC[i,4]=results.4$criterion

convergence[i,1]=results.1$mod$convcode
convergence[i,2]=results.2$mod$convcode
convergence[i,3]=results.3$mod$convcode
convergence[i,4]=results.4$mod$convcode

}
proc.time() - ptm


library(splines)
N=1000
W=10
par=c(1,1,.4)
p0=c(1,1)
x=data.frame(observer=factor(rep(c(1,2),times=N)),object=rep(1:N,each=2), distance=rep(runif(N,0,W),each=2),cov=rep(rnorm(N,0,1),each=2))
dd=sim.general(x,halfnormal,W=10,par=par,p0=p0,scale.formula=~-1+observer+cov)
par(mfrow=c(2,1))
hist(dd$distance[dd$observer==1&dd$ch%in%c("11","10")],main="Observer 1",xlab="Distance")
hist(dd$distance[dd$observer==2&dd$ch%in%c("11","01")],main="Observer 2",xlab="Distance")
#Fit loglinear model to data
results=ddf(data=dd,meta.data=list(width=W),method="loglinear",mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance))
results$criterion
# %diff
100*(results$Nhat-N)/N
# unconditonal detection fcts
par(mfrow=c(3,2))
plot(results,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results,7:8,showlines=FALSE)
#Fit PI bpi model to data
results=ddf(data=dd,meta.data=list(width=W),method="bpi",control=list(indep=FALSE,PI=TRUE),mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance))
results$criterion
# %diff
100*(results$Nhat-N)/N
# unconditional detection fcts
par(mfrow=c(3,2))
plot(results,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results,7:8,showlines=FALSE)
#Fit full independence bpi model to data
results=ddf(data=dd,meta.data=list(width=W),control=list(indep=TRUE,PI=FALSE),method="bpi",mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~0))
results$criterion
# %diff
100*(results$Nhat-N)/N
# unconditional detection fcts
par(mfrow=c(3,2))
plot(results,1:6,showpoints=FALSE)
# conditional detection functions
par(mfrow=c(1,2))
plot(results,7:8,showlines=FALSE)

 N=250
 data(book.tee.data,package="mrds")
 x=book.tee.data$book.tee.dataframe
 results=ddf(data=x,mrmodel=list(pformula=~-1+observer+observer:distance,dformula=~-1+distance),
   meta.data=list(width=4),method="loglinear")
 results$criterion
 # %diff
 100*(results$Nhat-N)/N
 # unconditional detection fcts
 par(mfrow=c(3,2))
 plot(results,1:6,showpoints=FALSE)
 # conditional detection functions
 par(mfrow=c(1,2))
 plot(results,7:8,showlines=FALSE)

 results=ddf(data=x,mrmodel=list(pformula=~observer*distance,dformula=~-1+bs(distance)),
             meta.data=list(width=4),method="bpi")
 results$criterion
 results$Nhat
 # %diff
 100*(results$Nhat-N)/N
 # unconditional detection fcts
 par(mfrow=c(3,2))
 plot(results,1:6,showpoints=FALSE)
 # conditional detection functions
 par(mfrow=c(1,2))
 plot(results,7:8,showlines=FALSE)


