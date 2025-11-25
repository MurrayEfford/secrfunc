# polyCub experimental 2024-02-02

library(polyCub)

dfn <- function (xy, mu, lambda0 = 1, sigma) {
    xy <- matrix(xy, ncol=2)
    xy <- sweep(xy, MARGIN=2, STATS=mu, FUN='-')
    d2 <- apply(xy^2,1,sum)
    lambda0 * exp(-d2/2/sigma^2)
}

dfn2 <- function (xy, lambda0, sigma) {
    xy <- matrix(xy, ncol=2)
    d2 <- apply(xy^2,1,sum)
    lambda0 * exp(-d2/2/sigma^2)
}

intrfr <- function (R, sigma = 5)
{
    (1 - exp(-R^2/2/sigma^2))/2/pi
}

library(secr)
library(sf)

nc <- setNumThreads(7)    # adjust as required
setwd('d:/density communication/book')
polyexample0 <- read.traps(file = 'polygonexample0.txt', detector = 'polygon')
polyexample1 <- read.traps(file = 'polygonexample.txt', detector = 'polygon')
msk <- make.mask(polyexample0, buffer=150, spacing=10)
# points(-30,+60, pch = 16, col = 'blue')

sfpoly <- st_polygon(list(as.matrix(polyexample0)))
sfpoly1 <- st_polygon(list(as.matrix(polyexample1)))
plot(sfpoly)
system.time(h <- apply(msk,1,polyCub, polyregion = sfpoly, f = dfn, method = "SV", lambda0=1/(2*pi*100^2), sigma=100))
# user  system elapsed 
# 7.94    0.04   21.63 

system.time(h <- apply(msk,1,polyCub, polyregion = sfpoly, f = dfn, method = "SV", nGQ=4, lambda0=1/(2*pi*100^2), sigma=100))
# nGQ = 4
# user  system elapsed 
# 0.52    0.01    1.69 
# nGQ = 6
# user  system elapsed 
# 1.39    0.00    2.75 

# polygon with many vertices (137)
plot(sfpoly1)
system.time(h1 <- apply(msk,1,polyCub, polyregion = sfpoly1, f = dfn, method = "SV", nGQ=4, lambda0=1/(2*pi*100^2), sigma=100))
# nGQ = 4
# user  system elapsed 
# 4.29    0.00    8.61 


system.time(h <- apply(msk,1,polyCub, polyregion = sfpoly, f = dfn, method = "SV", nGQ=6, lambda0=1/(2*pi*100^2), sigma=100))

h <- apply(msk,1,polyCub, polyregion = sfpoly, method = "exact.Gauss", Sigma = diag(100^2))
# Error: 'polyCub.exact.Gauss' is currently unavailable.
# Contributions are welcome: <https://github.com/bastistician/polyCub/issues/2>

covariates(msk) <- data.frame(h=h)
plot(msk, cov='h')
plot(polyexample0, add=T, detpar=list(col=NA))

system.time(h <- apply(msk,1, function(cent) polyCub.iso(polyregion = sfpoly, intrfr = intrfr, sigma=100, center=cent)))
# user  system elapsed 
# 0.30    0.00    1.75 

system.time(h1 <- apply(msk,1, function(cent) polyCub.iso(polyregion = sfpoly1, intrfr = intrfr, sigma=100, center=cent)))
# user  system elapsed 
# 5.06    0.00   10.15 

covariates(msk) <- data.frame(h=h1)
plot(msk, cov='h')
plot(polyexample1, add=T, detpar=list(col=NA))
