rm(list=ls())

########################################################################
# Lecture: Kriging
########################################################################

# the meuse data set, see https://cran.r-project.org/web/packages/gstat/vignettes/gstat.pdf

###############################################################
# Exploratory data analysis
###############################################################
library(geoR)
library(sp)
library(fields)

data(meuse)

### Zinc concentrations ###########################

coordinates(meuse) = ~x+y  # convert to SpatialPointsDataFrame
bubble(meuse, "zinc",col="red",maxsize=1.5,pch=19,fill=T,main = "zinc concentrations (ppm)")


###############################################################
# Variogram fitting
###############################################################

meuse.geoR = as.geodata(meuse,data.col=4)  # NOTE: 4 corresponds to zinc

## No covariates #####################################

### 1) Obtain empirical variogram
lz.geoR = meuse.geoR
lz.geoR$data = log(meuse.geoR$data)
lz.v = variog(lz.geoR,max.dist=1500)
plot(lz.v,pch=19)


### 2) Fit variogram model
lz.fit.lin = variofit(lz.v,ini=c(.4,300),cov.model="linear",fix.nugget = F,nugget=.2,kap=2.5)
lz.fit.sph = variofit(lz.v,ini=c(.4,300),cov.model="spherical",fix.nugget = F,nugget=.2,kap=2.5)
lz.fit.mat = variofit(lz.v,ini=c(.4,300),cov.model="matern",fix.nugget = F,nugget=.2,kap=2.5)

### 3) Plot for visual inspection
plot(lz.v,pch=19,main="Semivariogram: log(zinc)",ylim=c(0,.8))
lines(lz.fit.lin,col="red",lwd=2)
lines(lz.fit.sph,col="green",lwd=2)
lines(lz.fit.mat,col="blue",lwd=2)


###############################################################
# Kriging
###############################################################

load("meuse_borders.RData")

## ordinary kriging

# pred_krig generates the 2D prediction grid, 
# using the borders from meuse.borders
# the by option is to have a slightly coarser grid

gr   =  pred_grid(meuse.borders,by=40)

# actual kriging
kc    = krige.control(type="ok",obj.mod = lz.fit.mat)
sk    = krige.conv(lz.geoR,krige=kc,loc=gr,borders = meuse.borders)


# points and borders
data(meuse.grid)
plot(meuse.grid$x,meuse.grid$y)
plot(meuse.borders$x,meuse.borders$y)

## plotting the kriging maps 

par(mfrow=c(1,2),mar=c(3,3,2,1))
image(sk, x.leg = c(178500, 181500), y.leg = c(329100, 329350), 
      ylim=c(329000,334000),ylab="",xlab="",
      main="Kriging: log(zinc)")
points(meuse.geoR,pch=,add=T)

## plotting kriging standard deviations
image(sk, val = sqrt(sk$krige.var), x.leg = c(178500, 181500), y.leg = c(329100, 329350), 
      col=terrain.colors(20), ylim=c(329000,334000),ylab="",xlab="",
      main="Kriging sd: log(zinc)")
points(meuse.geoR,pch=,add=T)


## lower 95% confidence interval
image(sk, val = sk$predict-1.96*sqrt(sk$krige.var), x.leg = c(178500, 181500), y.leg = c(329100, 329350), 
      col=terrain.colors(20), ylim=c(329000,334000),zlim=c(3.7,7.5),ylab="",xlab="",
      main="Kriging: lower 95%")
points(meuse.geoR, pch=, add=T)

## upper 95% confidence interval
image(sk, val = sk$predict+1.96*sqrt(sk$krige.var), x.leg = c(178500, 181500), y.leg = c(329100, 329350), 
      col=terrain.colors(20), ylim=c(329000,334000),zlim=c(3.7,7.5),ylab="",xlab="",
      main="Kriging: upper 95%")
points(meuse.geoR, pch=, add=T)


# no nugget?
kc    = krige.control(type="ok",obj.mod = lz.fit.mat,nugget=0)
sk.nonug    = krige.conv(lz.geoR,krige=kc,loc=gr,borders = meuse.borders)

## plotting the kriging maps 
image(sk, val = sqrt(sk$krige.var), x.leg = c(178500, 181500), y.leg = c(329100, 329350), 
      col=terrain.colors(20), ylim=c(329000,334000),ylab="",xlab="",
      main="Kriging sd: log(zinc)")
points(meuse.geoR,pch=,add=T)

## plotting kriging standard deviations
image(sk.nonug, val = sqrt(sk.nonug$krige.var), x.leg = c(178500, 181500), y.leg = c(329100, 329350), 
      col=terrain.colors(20), ylim=c(329000,334000),ylab="",xlab="",
      main="Kriging sd: log(zinc)")
points(meuse.geoR,pch=,add=T)
