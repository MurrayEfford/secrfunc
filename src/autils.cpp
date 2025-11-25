#include "poly.h"

using namespace std;
using namespace Rcpp;

double minimumexp = -100;

//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k in 3D array
// stored in column-major order 

int i3 (int i, int j, int k, int ii, int jj) {
    return(ii * (jj * k + j) + i);
}
//--------------------------------------------------------------------------

// index to vector element corresponding to cell i,j,k,l in 4D array
// stored in column-major order 

int i4 (int i, int j, int k, int l, int ii, int jj, int kk) {
    return (ii *(jj*(kk*l + k) + j) + i);
}
//--------------------------------------------------------------------------

// customised dbinom 
double gbinom(int count, int size, double p)
{
    double x, q;
    int i;
    if ((count < 0) || (count > 0 && p <= 0)) {
        x = 0;
    }
    else if (count == 0) {
        q = 1 - p;
        x = q;
        for (i=1; i< size; i++) x *= q;
    }
    else {
        boost::math::binomial_distribution<> bin(size, p);
        x = boost::math::pdf(bin, count);
    }
    return (x);   
}
//--------------------------------------------------------------------------

// Calculate the length of intersection of a line segment and a circle
// Based on C code of Paul Bourke November 1992
// Line segment is defined from p1 to p2
// Circle is of radius r and centred at sc
// Two potential points of intersection given by
// p = p1 + mu1 (p2-p1)
// p = p1 + mu2 (p2-p1)
// Return 0 if line segment does not intersect circle

double SegCircle2 (
        double p1x, double p1y, 
        double p2x, double p2y, 
        double scx, double scy, 
        double r
) 
{
    double a,b,c;
    double bb4ac;
    double dpx;
    double dpy;
    
    double mu1;
    double mu2;
    int p1in;
    int p2in;
    
    double i1x;
    double i1y;
    double i2x;
    double i2y;
    
    int i1between;
    double d1,d2;
    double seg = 0;
    
    // case where both p1 and p2 inside circle 
    
    // Rprintf ("p1 %6.3f %6.3f\n", p1x, p1y);
    // Rprintf ("p2 %6.3f %6.3f\n", p2x, p2y);
    // Rprintf ("sc %6.3f %6.3f\n", scx, scy);
    // Rprintf ("r %6.3f \n", r);
    
    p1in = ((scx - p1x) * (scx - p1x) + 
        (scy - p1y) * (scy - p1y)) < (r * r);
    p2in = ((scx - p2x) * (scx - p2x) + 
        (scy - p2y) * (scy - p2y)) < (r * r);
    if (p1in && p2in) {        
        seg = std::sqrt ((p1x - p2x) * (p1x - p2x) + 
            (p1y - p2y) * (p1y - p2y));
        return (seg);
    }
    
    dpx = p2x - p1x;
    dpy = p2y - p1y;
    
    a = dpx * dpx + dpy * dpy;
    b = 2 * (dpx * (p1x - scx) + dpy * (p1y - scy));
    c = scx * scx + scy * scy;
    c += p1x * p1x + p1y * p1y;
    c -= 2 * (scx * p1x + scy * p1y);
    c -= r * r;
    bb4ac = b * b - 4 * a * c;
    
    // case of no intersection 
    if ((fabs(a) < 1e-10) || (bb4ac < 0)) {
        return (0);   
    }
    
    mu1 = (-b + std::sqrt(bb4ac)) / (2 * a);
    mu2 = (-b - std::sqrt(bb4ac)) / (2 * a);
    
    i1x = p1x + mu1 * (p2x - p1x);
    i1y = p1y + mu1 * (p2y - p1y);
    i2x = p1x + mu2 * (p2x - p1x);
    i2y = p1y + mu2 * (p2y - p1y);
    
    if (((mu1<0) && (mu2<0)) || ((mu1>1) && (mu2>1))) {
        // no intersection 
        seg = 0;
    }
    else {
        if (((mu1<0) && (mu2>1)) || ((mu1>1) && (mu2<0))) {
            // both inside 
            seg = std::sqrt ((p1x - p2x) * (p1x - p2x) + 
                (p1y - p2y) * (p1y - p2y));
        }
        else {
            if ((mu1>0) && (mu1<1) && (mu2>0) && (mu2<1)) {
                // two intersections 
                seg = std::sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
            }
            else {
                // one intersection 
                d1 = std::sqrt((i1x - p1x) * (i1x * p1x) + 
                    (i1y - p1y) * (i1y - p1y));
                d2 = std::sqrt((i1x - p2x) * (i1x * p2x) + 
                    (i1y - p2y) * (i1y - p2y));
                i1between = std::sqrt(a) < (d1 + d2 + 1e-10);
                if (p1in) {
                    if (i1between) {
                        i2x = p1x;
                        i2y = p1y;
                    }
                    else {
                        i1x = p1x;
                        i1y = p1y;
                    }
                }
                if (p2in) {
                    if (i1between) {
                        i2x = p2x;
                        i2y = p2y;
                    }
                    else {
                        i1x = p2x;
                        i1y = p2y;
                    }
                }
                seg = std::sqrt ((i1x - i2x) * (i1x - i2x) + 
                    (i1y - i2y) * (i1y - i2y));
            }
        }
    }
    return(seg);    
}

//----------------------------------------------------------------

// detect may take values -
// 0  multi-catch traps
// 1  binary proximity detectors
// 2  count  proximity detectors
// 3  exclusive polygon detector
// 4  exclusive transect detector
// 5  signal detector
// 6  polygon detector
// 7  transect detector
// 8  times  (undocumented)
// 9  cue    (undocumented) -- removed in secr 2.10.0
// 12 signalnoise


// probability of count for session s, detector k, animal i
// The argument 'g' is understood to be a cumulative hazard if binomN=0,
// a probability otherwise

double pski ( int binomN,
              int count,
              double Tski,
              double g,
              double pI) {
    
    double lambda;
    double result = 1.0;
    
    if (binomN == -1) {                              // binary proximity detectors : Bernoulli
        if (abs(Tski-1) > 1e-10) {                   // effort not unity; adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        if (count>0)                                 
            result = g*pI;  
        else 
            result = 1 - g*pI;
    }
    else if (binomN == 0) {                          // count detectors : Poisson 
        lambda = Tski * g * pI;
        if ((count < 0) || (count>0 && lambda<=0)) {         
            result = 0;
        }
        else if (count == 0) {
            result = exp(-lambda);            // routinely apply Tsk adjustment to cum. hazard 
        }
        else {
            boost::math::poisson_distribution<> pois(lambda);
            result = boost::math::pdf(pois,count);
        }
    }
    else if (binomN == 1) {                          // count detectors : Binomial, size from Tsk
        result = gbinom (count, round(Tski), g*pI); 
    }
    else if (binomN > 1) {                           // count detectors : Binomial, specified size 
        if (abs(Tski-1) > 1e-10) {                   // effort not unity, adjust g 
            g = 1 - pow(1 - g, Tski);
        }
        result = gbinom (count, binomN, g*pI);
    }
    else result = NAN; // Rcpp::stop("binomN < -1 not allowed");  // code multi -2 separately
    
    return (result);
}
//--------------------------------------------------------------------------

