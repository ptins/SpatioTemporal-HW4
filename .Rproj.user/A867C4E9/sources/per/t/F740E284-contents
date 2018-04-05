rm(list=ls())

library(geoR)
library(sp)
library(fields)

df = read.csv('precip.csv')
head(df)

coordinates(df) = ~LONGITUDE+LATITUDE
bubble(df, 'PRCP', col='blue', maxsize=1.5, pch=19, fill=T, main = 'Precipitation')

df.geoR = as.geodata(df,data.col=1)

lz.geoR = df.geoR
lz.geoR$data = log(df.geoR$data)
lz.v = variog(lz.geoR, max.dist=20)

par(mfcol=c(1,1))
plot(lz.v, pch=19)

## fit variogram models
lz.fit.lin = variofit(lz.v, ini=c(0.8,1), cov.model='linear', fix.nugget = F, nugget=0, kap=2.5)
lz.fit.sph = variofit(lz.v, ini=c(0.8,1), cov.model='spherical', fix.nugget = F, nugget=0, kap=2.5)
lz.fit.mat = variofit(lz.v, ini=c(0.8,1), cov.model='matern', fix.nugget = F, nugget=0, kap=2.5)

lz.fit.lin$value
lz.fit.sph$value
lz.fit.mat$value # we choose matern because minimal value of loss function

??variofit

## plot for visual inspection
plot(lz.v, pch=19, main='Semivariogram: log(PRCP)')
lines(lz.fit.lin, col='red', lwd=2)
lines(lz.fit.sph, col='green', lwd=2)
lines(lz.fit.mat, col='blue', lwd=2)

## create prediction grid 
x <- seq(min(lz.geoR$coords[,1]), max(lz.geoR$coords[,1]), length.out = 100)
y <- seq(min(lz.geoR$coords[,2]), max(lz.geoR$coords[,2]), length.out = 100)
gr <- expand.grid(x = x, y = y)

## actual kriging
kc = krige.control(type='ok', obj.mod = lz.fit.mat)
sk = krige.conv(lz.geoR, krige=kc, loc=gr)

## plotting the kriging maps 
par(mfrow=c(1,2))
image(sk,
      x.leg = c(min(df.geoR$coords[,1]), max(df.geoR$coords[,1])), y.leg = c(33, 34),
      ylim = c(32, 44),
      col=brewer.pal(10, "RdBu"),
      main="Kriging: log(PRCP)")
points(df.geoR, add=T)

## plotting kriging standard deviations
image(sk, val = sqrt(sk$krige.var), 
      x.leg = c(min(df.geoR$coords[,1]), max(df.geoR$coords[,1])), y.leg = c(33, 34),
      ylim = c(32, 44),
      col=rev(brewer.pal(10, "RdYlGn")),
      main="Kriging sd: log(PRCP)")
points(df.geoR,pch=,add=T)

par(mfrow=c(1,2))
## lower 95% confidence interval
image(sk, val = sk$predict-1.96*sqrt(sk$krige.var), 
      x.leg = c(min(df.geoR$coords[,1]), max(df.geoR$coords[,1])), y.leg = c(33, 34),
      ylim = c(32, 44),
      col=brewer.pal(10, "RdBu"),
      main="Kriging: lower 95%")
points(df.geoR,pch=,add=T)

## upper 95% confidence interval
image(sk, val = sk$predict+1.96*sqrt(sk$krige.var), 
      x.leg = c(min(df.geoR$coords[,1]), max(df.geoR$coords[,1])), y.leg = c(33, 34),
      ylim = c(32, 44),
      col=brewer.pal(10, "RdBu"),
      main = "Kriging: upper 95%")
points(df.geoR,pch=,add=T)

## interval size
par(mfrow=c(1,1))
image(sk, val = sk$predict+1.96*sqrt(sk$krige.var) - (sk$predict-1.96*sqrt(sk$krige.var)), 
      x.leg = c(min(df.geoR$coords[,1]), max(df.geoR$coords[,1])), y.leg = c(33, 34),
      ylim = c(32, 44),
      col=rev(brewer.pal(10, "PRGn")),
      main = "Kriging: 95% CI size")
points(df.geoR,pch=,add=T)

lz.geoR$data

