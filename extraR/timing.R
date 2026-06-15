library(microbenchmark)
microbenchmark(
    n2 = polygonhistoriescpp(
        nc, detectfn, grain, 2, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug),
    n4 = polygonhistoriescpp(
        nc, detectfn, grain, 4, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug),
    n6 = polygonhistoriescpp(
        nc, detectfn, grain, 6, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug),
    n8 = polygonhistoriescpp(
        nc, detectfn, grain, 8, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug),
    n10 = polygonhistoriescpp(
        nc, detectfn, grain, 10, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, debug)
)

    
microbenchmark(
    n2 = polygonhistories2cpp(
        nc, detectfn, grain, 2, safeLL, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, 
        mask_indices, mask_offsets, mask_id,
        debug),
    n4 = polygonhistories2cpp(
        nc, detectfn, grain, 4, safeLL, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, 
        mask_indices, mask_offsets, mask_id,
        debug),
    n6 = polygonhistories2cpp(
        nc, detectfn, grain, 6, safeLL, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, 
        mask_indices, mask_offsets, mask_id,
        debug),
    n8 = polygonhistories2cpp(
        nc, detectfn, grain, 8, safeLL, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, 
        mask_indices, mask_offsets, mask_id,
        debug),
    n10 = polygonhistories2cpp(
        nc, detectfn, grain, 10, safeLL, minp, binomN, w, xy, start, 
        group, gkhk$hk, gkhk$H, gsbval, pID, mask, density, PIA, Tsk, h, hindex, 
        mask_indices, mask_offsets, mask_id,
        debug)
)


