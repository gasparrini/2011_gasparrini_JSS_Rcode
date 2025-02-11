####################################################################
# Updated version of the code for the analysis in:
#
#   "Distributed lag linear and non-linear models in R: the package dlnm"
#   Gasparrini A
#   Journal of Statistical Software 2011
#   http://www.ag-myresearch.com/2011_gasparrini_jss.html
#
# Update: 11 January 2017
# * an updated version of this code, compatible with future versions of the
#   software, is available at:
#   https://github.com/gasparrini/2011_gasparrini_JSS_Rcode
####################################################################

####################################################################
# NB: THE CODE HAS BEEN ADAPTED TO THE NEW VERSION OF THE R PACKAGE dlnm
####################################################################

#options(prompt="R> ", continue="+  ", width=70, useFancyQuotes=FALSE)

library("dlnm")
vignette("dlnmOverview")

# CHECK VERSION OF THE PACKAGE
if(packageVersion("dlnm")<"2.2.0")
  stop("update dlnm package to version >= 2.2.0")

####################################################################
# NON-LINEAR AND DELAYED EFFECTS
####################################################################

# NB: THE FUNCTIONS mkbasis and mklagbasis HAVE BEEN REPLACED BY onebasis
# NB: CENTERING MOVED TO PREDICTION STAGE
onebasis(1:5, fun="bs", df=4, degree=2)
onebasis(1:5, fun="strata", breaks=c(2,4))

####################################################################
# SPECIFYING A DLNM
####################################################################

basis.o3 <- crossbasis(chicagoNMMAPS$o3, lag=10, 
  argvar=list(fun="thr",thr.value=40.3,side="h"),
  arglag=list(fun="strata",breaks=c(2,6)))

klag <- exp(((1+log(30))/4 * 1:3)-1)
basis.temp <- crossbasis(chicagoNMMAPS$temp, lag=30,
  argvar=list(fun="bs",degree=3,df=6), arglag=list(knots=klag))

summary(basis.temp)

library("splines")
model <- glm(death ~ basis.temp + basis.o3 + ns(time,7*14) + dow,
  family=quasipoisson(), chicagoNMMAPS)

####################################################################
# PREDICTING A DLNM
####################################################################

pred.o3 <- crosspred(basis.o3, model, at=c(0:65,40.3,50.3))
pred.temp <- crosspred(basis.temp, model, by=2, cen=25)

pred.o3$allRRfit["50.3"]
cbind(pred.o3$allRRlow,pred.o3$allRRhigh)["50.3",]

####################################################################
# REPRESENTING A DLNM
####################################################################

plot(pred.o3, var=50.3, type="p", pch=19, cex=1.5, ci="bars", col=2,
  ylab="RR",main="Lag-specific effects")
plot(pred.o3, "overall", ci="lines", ylim=c(0.95,1.25), lwd=2, col=4,
  xlab="Ozone", ylab="RR", main="Overall effect")

plot(pred.temp, xlab="Temperature", theta=240, phi=40, ltheta=-185,
  zlab="RR", main="3D graph")
plot(pred.temp, "contour", plot.title=title(xlab="Temperature",
  ylab="Lag", main="Contour graph"), key.title=title("RR"))

plot(pred.temp, var=-20, ci="n",ylim=c(0.95,1.22), lwd=1.5, col=2)
for(i in 1:2) lines(pred.temp, "slices", var=c(0,32)[i], col=i+2, lwd=1.5)
legend("topright",paste("Temperature =",c(-20,0,32)), col=2:4, lwd=1.5)

plot(pred.temp,var=c(-20,0,32), lag=c(0,5,20), ci.level=0.99, col=2,
  xlab="Temperature",ci.arg=list(density=20,col=grey(0.7)))

####################################################################
# MODELING STRATEGIES
####################################################################

basis.temp2 <- crossbasis(chicagoNMMAPS$temp, argvar=list(fun="poly",degree=6),
  arglag=list(knots=klag), lag=30)
model2 <- update(model, .~. - basis.temp + basis.temp2)
pred.temp2 <- crosspred(basis.temp2, model2, by=2, cen=25)

basis.temp3 <- crossbasis(chicagoNMMAPS$temp, argvar=list(fun="thr",
  thr.value=25,side="d"), arglag=list(knots=klag), lag=30)
model3 <- update(model, .~. - basis.temp + basis.temp3)
pred.temp3 <- crosspred(basis.temp3, model3, by=2)

plot(pred.temp, "overall", ylim=c(0.5,2.5), ci="n", lwd=1.5, col=2,
  main="Overall effect")
lines(pred.temp2, "overall", col=3, lty=2, lwd=2)
lines(pred.temp3, "overall", col=4, lty=4, lwd=2)
legend("top", c("natural spline","polynomial","double threshold"), col=2:4,
  lty=c(1:2,4), lwd=1.5, inset=0.1, cex=0.8)

plot(pred.temp, "slices", var=32, ylim=c(0.95,1.22), ci="n", lwd=1.5, col=2,
  main="Lag-specific effect")
lines(pred.temp2, "slices", var=32, col=3, lty=2, lwd=2)
lines(pred.temp3, "slices", var=32, col=4, lty=4, lwd=2)
legend("top", c("natural spline","polynomial","double threshold"), col=2:4,
  lty=c(1:2,4), inset=0.1, cex=0.8)

#
