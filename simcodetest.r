#sim loglinear
par=c( 1.517424, 1.663175, -0.5, -0.3, 0.311178)
W=10
N=1000
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
# unconditonal detection fcts
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

results.1$criterion
results.2$criterion
results.3$criterion
results.4$criterion

