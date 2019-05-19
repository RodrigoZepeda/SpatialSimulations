#include <RcppArmadillo.h>
#include <typeinfo>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rHawkesProcess(int nsim, double lambda, double rho) {
  
  //Simulate the number of centers of the clustering process
  int N = Rcpp::rpois(1, lambda)[0];
  Rcout << "N = " << N << std::endl;
  arma::mat Aux;
  arma::mat clusterC;
  arma::mat SimMatrix;
  int Nchild;
  
  //Verify that there are points here:
  if (N > 0){
    
    //Distribute the number of points uniformly in space
    SimMatrix = arma::mat(N, 2, arma::fill::none);

    //Simulate a homogeneous poisson process of centers
    SimMatrix.col(0) = as<arma::vec>(Rcpp::runif(N, 0, 1));
    SimMatrix.col(1) = as<arma::vec>(Rcpp::runif(N, 0, 1));

    //Generate the points procedurally generation by generation
    int tpoints = N; //Total number of points
    int gen = 1;       //STaring generation of children
    Rcout << "tpoints " << tpoints << std::endl;
    
    //Loop generation by generation iterating
    while (gen < tpoints){
      Rcout << "Loop " << gen << std::endl;
      
      //Select a center
      arma::subview_row<double> nextSim = SimMatrix.row(gen); 
      
      //Select number of children randomly
      Nchild = Rcpp::rpois(1, rho)[0];
      Rcout << "Nchild = " << Nchild << std::endl;
      
      //Repeat children matrix
      Aux = arma::repmat(nextSim, Nchild, 1);
      
      //Get coordinates (need to change this)
      clusterC = arma::mat(Nchild, 2, arma::fill::none);
      clusterC.col(0) = 0.02*as<arma::vec>(Rcpp::rnorm(Nchild));
      clusterC.col(1) = 0.02*as<arma::vec>(Rcpp::rnorm(Nchild));
      
      //Add to aux new cluster
      Aux = Aux + clusterC;
      
      //Update sim matrix
      SimMatrix = arma::join_cols(SimMatrix,Aux);
      
      //Update counters
      gen = gen + 1;
      tpoints = tpoints + Nchild;
    }
  }
    return SimMatrix;
    
}


/*** R
plot(rHawkesProcess(100, 30, 0.9))
*/
