#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rStraussProcess(int nsim, double muE, Function lambda) {
  
  //Simulate the number of points on space
  NumericVector N = Rcpp::rpois(1, muE);
  
  //Verify that there are points here:
  if (N[0] > 0){
    
    //Distribute the number of points uniformly in space
    NumericMatrix SimMatrix( N[0], 2 );
    
    //Simulate a homogeneous poisson process
    SimMatrix(_, 0) = Rcpp::runif(N[0], 0, 1);
    SimMatrix(_, 1) = Rcpp::runif(N[0], 0, 1);
    
    //Acceptance and rejection
    NumericVector pcriteria = Rcpp::runif(N[0],0,1);
    NumericMatrix pprocess (N[0] , 2);
    NumericVector lambdavec = lambda(SimMatrix(_, 0), SimMatrix(_, 1));
    int counter = 0;
    for (int i = 0; i < N[0]; i++){
      if (pcriteria[i] < lambdavec(i)/muE){
        pprocess(counter,_) = SimMatrix(i, _);
        counter += 1;
      }
    }
    if (counter > 0){
      return pprocess(Range(0,counter-1),_);  
    } else {
      NumericMatrix pprocess (0 , 2);
      return pprocess;
    }
  } else {
    NumericMatrix pprocess (0 , 2);
    return pprocess;
  }
}


/*** R
rStraussProcess(100, 600, function(x,y){300*(x^2 + y^2)})
*/
