#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector llrandpoiscpp(NumericVector y, NumericVector lp, NumericVector tau2, NumericMatrix gh) {

int nquad=gh.nrow();
int nobs=y.size();

NumericVector ll(nobs);
double thel;

for (int irow=0;irow<nobs;irow++) {
  thel=0.0;
  if (tau2[0]==0.0) ll(irow) = R::dpois(y(irow), std::exp(lp(irow)), true);
  else {
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+R::dpois(y(irow),std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0])),false)*gh(j,1);
    }
    ll(irow) = std::log(thel);
  }
}
return ll;
}

// [[Rcpp::export]]
NumericVector llrandnegbinomcpp(NumericVector y, NumericVector lp, NumericVector tau2,  NumericVector theta, NumericMatrix gh) {
  
  int nquad=gh.nrow();
  int nobs=y.size();
  
  NumericVector ll(nobs);
  double thel;
  
  for (int irow=0;irow<nobs;irow++) {
    thel=0.0;
    if (tau2[0]==0.0) ll(irow) = R::dnbinom_mu(y(irow), theta[0], std::exp(lp(irow)), true);
    else {
      thel=0.0;
      for (int j=0;j<nquad;j++) {
        thel=thel+R::dnbinom_mu(y(irow),theta[0], std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0])),false)*gh(j,1);
      }
      ll(irow) = std::log(thel);
    }
  }
  return ll;
}


  
// [[Rcpp::export]]
NumericVector llrandtruncpoiscpp(NumericVector y, NumericVector lp, NumericVector tau2, NumericMatrix gh) {
    
    int nquad=gh.nrow();
int nobs=y.size();

NumericVector ll(nobs);
double thel;

for (int irow=0;irow<nobs;irow++) {
  thel=0.0;
  /* ???? call actuar::dztpois
  */
  if (tau2[0]==0.0) ll(irow) = R::dpois(y(irow), std::exp(lp(irow)), true)-std::log(1-std::exp(-(std::exp(lp(irow)))));
  else {
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      double lambda = std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0]));
      thel=thel+R::dpois(y(irow), lambda, false)/(1-std::exp(-lambda))*gh(j,1);
    }
    ll(irow) = std::log(thel);
  }
}
return ll;
}

// [[Rcpp::export]]
NumericVector llrandgammacpp(NumericVector y, NumericVector lp, NumericVector tau2, NumericVector phi, NumericMatrix gh) {

int nquad=gh.nrow();
int nobs=y.size();

NumericVector ll(nobs);
double thel;

for (int irow=0;irow<nobs;irow++) {
  thel=0.0;
  if (tau2[0]==0.0) ll(irow) = R::dgamma(y(irow), 1.0/phi[0], 1.0/(phi[0]*std::exp(lp(irow))), true);
  else {
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+R::dgamma(y(irow), 1.0/phi[0], 1.0*(phi[0]*std::exp(lp(irow)+gh(j,0)*sqrt(tau2[0]))),false)*gh(j,1);
    }
    ll(irow) = std::log(thel);
  }
}
return ll;
}

// [[Rcpp::export]]
NumericVector llrandbinomcpp(NumericMatrix y, NumericVector lp, NumericVector tau2, NumericMatrix gh) {

int nquad=gh.nrow();
int nobs=y.nrow();

NumericVector ll(nobs);
double thel;

for (int irow=0;irow<nobs;irow++) {
  thel=0.0;
  if (tau2[0]==0.0) ll(irow) = R::dbinom(y(irow,0),y(irow,0)+y(irow,1), 1.0/(1.0+std::exp(-lp(irow))), true);
  else {
    thel=0.0;
    for (int j=0;j<nquad;j++) {
      thel=thel+R::dbinom(y(irow,0),y(irow,0)+y(irow,1),1.0/(1.0+std::exp(-lp(irow)-gh(j,0)*sqrt(tau2[0]))),false)*gh(j,1);
    }
    ll(irow) = std::log(thel);
  }
}
return ll;
}
