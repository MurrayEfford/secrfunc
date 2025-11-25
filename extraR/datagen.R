library(secr)
trps <- make.poly()
trps[] <- trps*0.01
detlist <- list(lambda0 = 0.5, sigma = 0.2)
ch <- sim.capthist(trps, popn = list(D=1e5, buffer = 2), detectfn=14, 
                   detectpar = detlist, noccasions = 5,
                   seed = 123)

dettype <- secr:::secr_detectorcode(trps, MLonly = TRUE, noccasions = 5)
binomN <- secr:::secr_recodebinomN(dettype, 0, 0)
xy0 <- secr:::secr_getxy(dettype, ch)
dimw <- dim(ch)
w <- as.integer(as.array(ch))
xy <- as.matrix(xy0$xy)
start <- as.integer(xy0$start)

msk <- make.mask(trps, buffer = 1, spacing = 0.05)

gsbval <- matrix(unlist(detlist), nrow = 1)
cumk <- c(0,5)   # vertices of square detector area  

mask <- as.matrix(msk)
traps <- as.matrix(trps)
save(w, dimw, xy, start, mask, traps, gsbval, cumk, binomN, file = 'inst/extdata/testdata.RData')

plot(ch, tracks=T, border=0.1)
setNumThreads(18)
fit <- secr.fit(ch, mask = msk, detectfn = 14, start = detlist)