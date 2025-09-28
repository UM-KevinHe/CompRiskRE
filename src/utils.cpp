#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double CindexCR_cpp(const arma::vec& time,
                    const arma::vec& status,
                    const arma::vec& predicted,
                    const int Cause_int) {
  int n = time.n_elem;
  arma::vec Censoring = arma::conv_to<arma::vec>::from(status != 0);
  arma::vec Cause(n, arma::fill::ones);
  Cause.elem(arma::find(status == 2)).fill(2);
  arma::vec Prediction = -arma::log(predicted);
  double Time = time.max() + 1.0;
  
  arma::mat A(n,n,arma::fill::zeros);
  arma::mat B(n,n,arma::fill::zeros);
  arma::mat Q(n,n,arma::fill::zeros);
  arma::mat N_t(n,n,arma::fill::zeros);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (time[i] < time[j]) A(i,j) = 1.0;
      if ((time[i] >= time[j]) && (Cause[j] != Cause_int) && (Censoring[j] == 1)) B(i,j) = 1.0;
      if (Prediction[i] > Prediction[j]) Q(i,j) = 1.0;
    }
  }
  for (int i=0; i<n; i++) {
    if (time[i] <= Time && Cause[i] == Cause_int && Censoring[i] == 1) {
      N_t.row(i).ones();
    }
  }
  arma::mat Num_mat = (A + B) % Q % N_t;
  arma::mat Den_mat = (A + B) % N_t;
  double Num = arma::accu(Num_mat);
  double Den = arma::accu(Den_mat);
  return Num/Den;
}
