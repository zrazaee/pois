#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include<random>

using namespace Rcpp;
using arma::mat;
using arma::umat;
using arma::vec;
using arma::uvec;
using arma::urowvec;
using arma::ivec;
using arma::irowvec;
using arma::cube;
using arma::rowvec;
using arma::zeros;
using arma::ones;
using arma::span;
using arma::linspace;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// double sample_int(int n) {
//   vec temp = Rcpp::runif(1) * n;
//   return temp(0);
// }

// // [[Rcpp::export]]
// arma::vec seq_int(long int a, long int b){
//   long int d = std::abs(b-a)+1;
//
//   arma:
//
//   return arma::linspace(a, b, d);
// }

int sample_int(int N) {
  Rcpp::IntegerVector pool = Rcpp::seq(0, N-1);
  std::random_shuffle(pool.begin(), pool.end());
  return pool[0];
}


int my_sampler(arma::rowvec prob_vec) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::vector<double> weights = arma::conv_to< std::vector<double>>::from(prob_vec);

  std::discrete_distribution<> d(weights.begin(), weights.end());
  return( d(gen) );
}

//' @export
// [[Rcpp::export]]
arma::umat sample_pois(int n, arma::mat theta, arma::mat gam, int burn_in=1000, int spacing=100, bool verb=true) {
  int d = theta.n_rows;
  int r = theta.n_cols;

  umat data(n, d);
  if (verb) Rcout << "burn in = " << burn_in; // << '\n';
  // Rcout << "Generating ";
  int t = 0;
  int report_chunk = round(n/8);
  urowvec z = zeros<urowvec>(d);
  irowvec sig = ones<irowvec>(d); // initialize to 1
  // int counter_max = burn_in;
  for (int t=0; t < n; t++) {
    if (t % report_chunk == 0 && verb) Rcout << '.';
    // for (int iter=1; iter < counter_max; iter++) {
    for (int iter=0; iter < ((t == 0) ? burn_in : spacing); iter++) {
      // Rcout << '*';
      // Rcout << "\n" << z << sig;

      int i = sample_int(d);
      rowvec beta = exp(theta.row(i));
      double alpha = 0;
      for (int j=0; j < d; j++) {
        if (j != i) {
          alpha += gam(i,j) * sig(j);
        }
      }
      // alpha <- sum(gam(i,-i) * sig[-i])
      double f1 = exp(-2*alpha);
      double Const = 1 + f1*sum(beta);

      rowvec pr(r+1);
      pr(0) = 1/Const;
      pr(span(1,r)) = f1*beta/Const;

      z(i) = my_sampler(pr);
      //uvec temp = RcppArmadillo::sample(linspace<uvec>(0,r,r+1), 1, false, pr);
      sig(i) = (z(i) == 0) ? 1 : -1;
      //Rcout << z(i) << sig(i) << "\n";
    }
    // Rcout << '\n';
    data.row(t) = z;
    // counter_max = spacing;
    //data.row(t) = z.as_row();

    //if (t % 100 == 0 && t > 0 && verb) Rcout << '\n';

  }

  if (verb) Rcout << '\n';

  return data;
}



// // [[Rcpp::export]]
// arma::umat toxGenCpp(int n, arma::mat theta, arma::mat gam, int burn_in=1000, bool verb=true) {
//   int d = theta.n_rows;
//   int r = theta.n_cols;
//
//   umat data(n, d);
//   if (verb) Rcout << "burn in = " << burn_in << '\n';
//   // Rcout << "Generating ";
//     for (int t=0; t < n; t++) {
//
//         if (t % 10 == 0 && verb) Rcout << '.';
//         urowvec z = zeros<urowvec>(d);
//         irowvec sig = ones<irowvec>(d); // initialize to 1
//
//
//         for (int iter=1; iter < burn_in; iter++) {
//           // Rcout << "\n" << z << sig;
//
//           int i = sample_int(d);
//
//           // arma::rowvec temp0 = exp(theta.row(i));
//           // vec beta = temp0.as_col();
//           rowvec beta = exp(theta.row(i));
//
//           // beta = exp(theta.row(i));
//
//           double alpha = 0;
//           for (int j=0; j < d; j++) {
//             if (j != i) {
//               alpha += gam(i,j) * sig(j);
//             }
//           }
//           // alpha <- sum(gam(i,-i) * sig[-i])
//           double f1 = exp(-2*alpha);
//           double Const = 1 + f1*sum(beta);
//
//           rowvec pr(r+1);
//           pr(0) = 1/Const;
//           pr( span(1,r) ) = f1*beta/Const;
//           // for (int j=1; j <= r; j++) {
//           //  pr(j) = f1*beta(j)/Const;
//           // }
//
//           // Rcout << linspace<uvec>(0,r,r+1) << "\n";
//           // Rcout << pr << "--" << sum(pr) << "\n";
//
//           z(i) = my_sampler(pr);
//           //uvec temp = RcppArmadillo::sample(linspace<uvec>(0,r,r+1), 1, false, pr);
//           // z(i) = temp(0);
//           sig(i) = (z(i) == 0) ? 1 : -1;
//           //Rcout << z(i) << sig(i) << "\n";
//         }
//
//         data.row(t) = z;
//         //data.row(t) = z.as_row();
//         if (t % 100 == 0 && t > 0 && verb) Rcout << '\n';
//           // if (t %% 100 == 0) cat('\n')
//     }
//
//     // Rcout << '\n';
//
//     return data;
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
library(Matrix)
d <- 10 #number of toxicity
r <- 5  #toxicity levels
n <- 300 # number of patients

#load("realdata_powertest.RData")
#theta1 <- -4*matrix(runif(d*r),d)
theta1 <- cbind(matrix(rep(-1,d*2),nrow=d,2),matrix(rep(-Inf,d*(r-2)),nrow = d,r-2))
theta1 <- matrix(-Inf,nrow = d,r)
theta1[,2] <- -2
theta1[,3] <- -3
#theta1 <- do.call(rbind, lapply(split(theta1, row(theta1)), sort, decreasing = TRUE))

#theta2 <- -4*matrix(runif(d*r),d)
theta2 <- theta1
theta2[1:4,2:3] <- theta1[1:4,2:3]
theta2[5:10,1] = -2.5
theta2[1:4,4] = -3
theta2[1,5] = -4.5

#theta2 <- do.call(rbind, lapply(split(theta2, row(theta2)), sort, decreasing = TRUE))

#theta2 <-  theta1[,sample(ncol(theta1))]
#theta2 <- cbind(theta1[,1],theta1[,3],theta1[,2],theta1[,4],theta1[,5])

gam <- forceSymmetric(0.2*Matrix(rnorm(d^2),d))
diag(gam) <- 0


data1 <- toxGenCpp(200, theta1, as.matrix(gam), verb=F)
# toxGen(20, theta1, as.matrix(gam), 10)
# sample_int(10)
*/
