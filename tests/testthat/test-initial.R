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
