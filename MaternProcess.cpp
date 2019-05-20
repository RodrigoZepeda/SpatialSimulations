#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List rMaternProcess(double kappa, double alpha, double r, arma::vec xlim, arma::vec ylim) {

  int NCluster;
  int cluster_tpoints;
  arma::mat SimMatrix;
  arma::vec candidate;
  arma::mat aux;
  
  //Get average number of points
  double mpoints = kappa*pow(1 + 2*r,2);

  //Simulate the number of clusters
  int N = Rcpp::rpois(1, mpoints)[0];

  //Verify that there are points here:
  if (N > 0){
    
    //Distribute the number of points uniformly in space
    SimMatrix = arma::mat(N, 2, arma::fill::none);
    
    //Simulate a homogeneous poisson process for the centers on window [-r, r+1]
    SimMatrix.col(0) = as<arma::vec>(Rcpp::runif(N, xlim(0)-r, xlim(1) + r));
    SimMatrix.col(1) = as<arma::vec>(Rcpp::runif(N, ylim(0)-r, ylim(1) + r));
    
    int tpoints = N;   //Total number of points

    //Loop through each cluster
    for (int i = 0; i < N; i++){

      //Simulate number of points per cluster
      NCluster = Rcpp::rpois(1, alpha)[0];
      cluster_tpoints = 0;
      
      //Create matrix of auxiliary points
      aux = arma::mat(NCluster, 2, arma::fill::none);
      
      //SImulate custer points via accept and reject
      while (cluster_tpoints < NCluster){
        
        //Simulate candidate point in ball of radius 1
        candidate    = arma::vec(2);
        candidate(0) = Rcpp::runif(2, xlim(0)-r, xlim(1) + r)[0];
        candidate(1) = Rcpp::runif(2, ylim(0)-r, ylim(1) + r)[1];
        if (arma::norm(candidate) < r){
          
          //Add point to matrix
          aux.row(cluster_tpoints) = candidate.t() + SimMatrix.row(i);
          
          //Update cluster and total points
          cluster_tpoints = cluster_tpoints + 1;
          tpoints = tpoints + 1;
        }
      }
      
      //Append to X
      SimMatrix = arma::join_cols(SimMatrix, aux);
      
    }
  }
  return Rcpp::List::create(Rcpp::Named("Simulations") = SimMatrix,
                            Rcpp::Named("N") = N);
  
}