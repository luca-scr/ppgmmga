#include <RcppArmadillo.h>
#include <algorithm>

using namespace Rcpp;
//[[Rcpp::depends(RcppArmadillo)]]



// [[Rcpp::export]]
List LinTransf(arma::mat mean,
               arma::cube sigma,
               arma::mat B,
               arma::mat Z,
               int G,
               int d)
{

  arma::mat Bt = B.t();
  arma::mat m = Bt*mean;
  arma::mat sz = Bt*cov(Z)*B;
  arma::cube s(d,d,G); s.zeros();

  //arma::uword Brow =  B.n_rows;
  //arma::uword Bcol = B.n_cols;

  // if(Brow != d & Bcol != d)
  // {
  //
  //   B.save("B.csv", arma::csv_ascii);
  // }

  for(int i=0;i<G;i++)
  {
    s.slice(i) = Bt*sigma.slice(i)*B;
  }

  return List::create(Rcpp::Named("mean") = m,
                      Rcpp::Named("sigma") = s,
                      Rcpp::Named("sz") = sz);
}


// [[Rcpp::export]]
double EntropyGauss(arma::mat S, int d)
{
  return 0.5 * (d * (log(2 * arma::datum::pi) + 1) + log(arma::det(S)));
}


// [[Rcpp::export]]
arma::mat orth(arma::mat A, std::string method = "QR")
{
  double    eps = std::numeric_limits<double>::epsilon();
  arma::mat B;

  if(method == "SVD")
  {
    arma::mat U;
    arma::vec d;
    arma::mat V;
    double tol = 0;
    int    r = 0;
    arma::uword m = A.n_rows;
    arma::uword n = A.n_cols;
    arma::vec s(m); s.zeros();
    arma::svd_econ(U,d,V,A, "both", "dc");

    if(m >= 1)  s = d;

    double m1 = std::max(m,n);
    double m2 = arma::max(s);

    tol = m1 * m2 * eps;
    r = sum(s > tol) - 1;

    //Rcout << U.col(0) << "\n";
    //arma::mat out = U(arma::span(0,(m-1)), arma::span(0,r));
    //return U;
    //return out;
    
    B = U(arma::span(0,(m-1)), arma::span(0,r));

  } else if(method == "QR")
  {
    arma::mat Q, R;
    arma::qr_econ(Q,R,A);
    B = Q;
  }

  return B;
}

// [[Rcpp::export]]
NumericVector encode(NumericVector par, int p) 
{
  int n = p-1;
  NumericVector w(p, 1.0);

  if(p==2)
  {
    w(0) = sin(par(0));
    w(1) = cos(par(0));
  }
  else if(p==3)
  {
    w(0) = sin(par(1))*sin(par(0));
    w(1) = sin(par(1))*cos(par(0));
    w(2) = cos(par(1));
  } else 
  {
    int count = n-1;
    w(n) = cos(par(1));
    for(int i=1; i<n; i++)
    {
      w(0) *= sin(par(i));
      w(1) *= sin(par(i));
    }
    w(0) = w(0)*sin(par(0));
    w(1) = w(1)*cos(par(0));

    for(int i=2; i<p-1; i++)
    {
      for(int j=1; j <= n-i; j++)
      {
        w(i) *= sin(par(j));
      }
      w(i) = w(i)*cos(par(count));
      count--;
    }
  }

  return w;
}

// [[Rcpp::export]]
NumericMatrix encodebasis(NumericVector par, int d, int p)
{
  //par.attr("dim") = Dimension(d, p);
  NumericMatrix parMat(p-1, d, par.begin());
  NumericMatrix basis(p, d);
  for(int i=0;i<d;i++)
  {
    basis(_,i) = encode(parMat(_,i), p);
  }

  return basis;
}




