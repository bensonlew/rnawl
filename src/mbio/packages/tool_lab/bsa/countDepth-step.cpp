#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector mdepths(NumericVector POS, NumericVector DEPTH, NumericVector pos1,NumericVector pos2) {
    unsigned int nout=POS.size(), i, nleft=0, nright=0,p;
    NumericVector out(nout);
    for( i=0; i < nout; i++ ) {
		double depth=0;
        while ((nright < nout-1) & (POS[nright + 1] <= pos2[i]))
            nright++;
        while (POS[nleft] < pos1[i])
			nleft++;
		for (p = nleft;p<=nright;p++ )
		{
			depth+=DEPTH[p];
		}
        out[i] = depth/(nright - nleft + 1 );
    }
    return out;
}
