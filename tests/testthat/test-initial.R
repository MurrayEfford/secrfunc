## Started 2025-11-25
library(secrfunc)

## to avoid ASAN/UBSAN errors on CRAN, following advice of Kevin Ushey
## e.g. https://github.com/RcppCore/RcppParallel/issues/169
Sys.setenv(RCPP_PARALLEL_BACKEND = "tinythread")

###############################################################################
## 

# centred unit polygon 
traps <- matrix(c(0,0,1,1,0,0,1,1,0,0), ncol = 2) - 0.5
Tsk <- matrix(1, nrow = 1, ncol = 5)
markocc <- rep(1,5)
cumk <- c(0,5,0)   # zero-terminated
detectfn <- 14   # hazard half-normal

test_that("correctly computed hdot", {
    # central, edge and outer points
    xy01 <- matrix(c(0,-0.5,-5,0,0,-5), ncol = 2)
    h <- hdotpolycpp(xy01, traps, Tsk, markocc, cumk, detectfn, c(1,0.1), TRUE,2,1,2)
    expect_equal(h, c(5.0, 2.5, 0.0), tolerance = 1e-5, check.attributes = FALSE)
})

###############################################################################

# load pre-computed inputs from datagen.R
# w, dimw, xy, start, traps, mask, gsbval, cumk, binomN
datafilename <- system.file("extdata/testdata.RData", package = "secrfunc")
load(datafilename)  

detectfn <- 14   # hazard half-normal
dim <- 2
convex <- TRUE
grain <- 1
ncores <- 2
gkhk <- makegkPolygoncpp(detectfn, dim, convex, grain, ncores, gsbval, cumk, traps, mask)

test_that("correctly computed gkhk", {
    expect_equal(sum(gkhk$H), 0.125664, tolerance = 1e-5, check.attributes = FALSE)
    expect_equal(sum(gkhk$hk), 200, tolerance = 1e-5, check.attributes = FALSE)
})

nc      <- dimw[1]  # individuals
S       <- dimw[2]  # occasions
K       <- dimw[3]  # single polygon
minp    <- 1e-200
group   <- rep(0,nc)
pID     <- matrix(1, nrow = S, ncol = 1)
density <- matrix(1/nrow(mask), nrow(mask), 1)
PIA     <- as.integer(array (1, dim = c(nc, S, K)))
Tsk     <- matrix(1, nrow = K, ncol = S)
h       <- matrix(-1)
hindex  <- matrix(-1)
debug   <- FALSE

test_that("correctly computed prw", {
    prw <- polygonhistoriescpp(
        nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug)
    expect_equal(prw, c(0.00267240, 0.00407618, 0.02350034, 0.02881529, 0.00893066,
                        0.01947464, 0.00730938, 0.01595762, 0.02161582, 0.00832001,
                        0.02531177), tolerance = 1e-5, check.attributes = FALSE)
})

test_that("correctly computed fxi", {
    fxi <- polygonfxicpp(
        nc, detectfn, grain, ncores, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex)
    expect_equal(sum(fxi), 0.165984, tolerance = 1e-5, check.attributes = FALSE)
})

###############################################################################
