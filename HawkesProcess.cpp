#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
List rHawkesProcess(double lambda, double rho, double alpha, arma::vec xlim, arma::vec ylim) {
  
  //Initialize the variables
  arma::mat Aux;
  arma::mat clusterC;
  arma::mat SimMatrix;
  int Nchild;
  
  //Simulate the number of centers of the clustering 
  int N = Rcpp::rpois(1, lambda)[0];
  
  //Verify that there are points to simulate:
  if (N > 0){
    
    //Distribute the number of points uniformly in space
    SimMatrix = arma::mat(N, 2, arma::fill::none);

    //Simulate a homogeneous poisson process of centers
    SimMatrix.col(0) = as<arma::vec>(Rcpp::runif(N, xlim(0), xlim(1)));
    SimMatrix.col(1) = as<arma::vec>(Rcpp::runif(N, ylim(0), ylim(1)));

    //Generate the points procedurally generation by generation
    int tpoints = N;   //Total number of points
    int gen = 1;       //Starting generation of children

    //Loop generation by generation iterating
    while (gen < tpoints){

      //Select a center
      arma::subview_row<double> nextSim = SimMatrix.row(gen); 
      
      //Select number of children randomly
      Nchild = Rcpp::rpois(1, rho)[0];
      
      if (Nchild > 0){
        
        //Create matrix of children
        Aux = arma::repmat(nextSim, Nchild, 1);
        
        //Get coordinates (need to change this)
        clusterC = arma::mat(Nchild, 2, arma::fill::none);
        clusterC.col(0) = alpha*as<arma::vec>(Rcpp::rnorm(Nchild));
        clusterC.col(1) = alpha*as<arma::vec>(Rcpp::rnorm(Nchild));
        
        //Add to aux new cluster
        Aux = Aux + clusterC;
        
        //Update sim matrix
        SimMatrix = arma::join_cols(SimMatrix,Aux);
        
        //Update total points
        tpoints = tpoints + Nchild;
      }
      
      //Update counters
      gen = gen + 1;
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("Simulations") = SimMatrix,
                            Rcpp::Named("N") = N);
    
}