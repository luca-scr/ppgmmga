#include <RcppArmadillo.h>

using namespace Rcpp;

//[[Rcpp::depends(RcppArmadillo)]]

//////////////////////////////////////////////////
//            PRELIMINARY FUNCTIONS             //
//              internal functions              //
/////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////
//       KULBACK-LEIBNER DIVERGENCE FOR TWO GAUSSIAN DISTRIBUTIONS      //
//                                                                      //
//  ARGS:                                                               //
//                                                                      //
//  mean1 = mean first Normal distributions                             //
//  mean2 = mean second Normal distributions                            //
//  sigma1 = covariance matrix first Normal distributions               //
//  sigma2  = covariance matrix second Normal distributions             //
//////////////////////////////////////////////////////////////////////////


double KLMN(arma::vec mean1, arma::mat sigma1, arma::vec mean2, arma::mat sigma2){

  int d=sigma1.n_rows;
  arma::mat div; div.zeros();

  div = 0.5*(log(arma::det(sigma2)) - log(arma::det(sigma1)) + arma::trace( arma::inv(sigma2) *sigma1 ) - d + ((arma::trans((mean2-mean1)) * arma::inv(sigma2)) * (mean2-mean1)));

  return arma::as_scalar(div);
}

//////////////////////////////////////////////////////////////////////////
//       DENSITY MULTIVARIATE NORMAL DISTRIBUTIONS                      //
//            (for a single observation)                                //
//  ARGS:                                                               //
//                                                                      //
//  X = data (note that we need a rowvector)                            //
//  mean = mean   (note that we need a rowvector)                       //
//  sigma = covariance matrix                                           //
//  logd = logical value; if true density is returned in log scale      //
//////////////////////////////////////////////////////////////////////////

const double log2pi = std::log(2.0 * M_PI);
double dmvnrm(arma::rowvec x,
              arma::rowvec mean,
              arma::mat sigma,
              bool logd = false) 
{
  int    xdim = sigma.n_cols;
  
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  double constants = -(static_cast<double>(xdim)/2.0) * log2pi;

  arma::vec z = rooti * arma::trans( x - mean) ;
  double    out = constants - 0.5 * arma::sum(z%z) + rootisum;

  if (logd == false) 
  {
    out = exp(out);
  }
  
  return(out);
}

//////////////////////////////////////////////////////////////////////////
//       Logsumexp                                                      //
//////////////////////////////////////////////////////////////////////////

double logsumexp(arma::vec x) 
{
  double max = x.max();
  double lse = max + log(arma::sum(exp(x-max)));
  return lse;
}


//////////////////////////////////////////////////////////////////////////
//                         GMMs DENSITY                                 //
//                                                                      //
//  ARGS:                                                               //
//                                                                      //
//  data = data (note that we need a rowvector)                         //
//  G = Number of the GMM components                                    //
//  pro = mixing proportions                                            //
//  mean = mean   (note that we need a rowvector)                       //
//  sigma = covariance matrix                                           //
//  logarithm = if log density or normal scale                          //
//////////////////////////////////////////////////////////////////////////
double mixDensity(arma::rowvec data,
                  int G,
                  arma::vec pro,
                  arma::mat mean,
                  arma::cube sigma,
                  bool logarithm = false)
{
  
  arma::mat  m = mean.t();
  arma::vec  dens; dens.zeros(G);
  double out; out = 0;
  
  for(int i=0; i<G; i++)
  {
    dens(i) = log(pro(i)) + dmvnrm(data, m.row(i), sigma.slice(i), true);
  }
  
  out = logsumexp(dens);
    
    if(logarithm == false)
    {
      out = exp(out);
    }
    
    return out;
}

//////////////////////////////////////////////////////////////////////////
//                         GRADIENT OF GMM                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

arma::vec gradient(arma::vec data,
                   int G,
                   arma::vec pro,
                   arma::mat mean,
                   arma::cube sigma)
{

  arma::rowvec tradata = data.t();
  arma::mat traM = mean.t();

  unsigned int d = mean.n_rows;
  arma::vec output; output.zeros(d);

  for(int i=0; i<G; i++)
  {
    output = output +
             arma::inv_sympd(sigma.slice(i)) *
             (mean.col(i) - data) *
             pro(i) *
             dmvnrm(tradata, traM.row(i), sigma.slice(i));
  }

  return output;
}

//////////////////////////////////////////////////////////////////////////
//              INTERNAL FUNCTION TO COMPUTE SOTE                       //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

arma::mat F(arma::vec data,
            int G,
            arma::vec pro,
            arma::mat mean,
            arma::cube sigma)
{

  arma::rowvec tradata = data.t();
  arma::mat traM = mean.t();
  unsigned int d = mean.n_rows;
  arma::mat output; output.zeros(d,d);
  double MixtureDens = mixDensity(tradata,G, pro, mean,sigma);
  arma::vec grad = gradient(data, G, pro, mean, sigma);
  arma::mat I; I.eye( d, d );
  
  for(int i=0; i<G; i++)
  {
    output = output +
             pro(i) *
             arma::inv_sympd(sigma.slice(i)) *
             (((1/MixtureDens) *
             (data - mean.col(i)) *
             grad.t()) +
             ((data - mean.col(i)) *
             arma::trans(arma::inv_sympd(sigma.slice(i)) *
             (data - mean.col(i)))) - I) *
             dmvnrm(tradata, traM.row(i), sigma.slice(i));
  }

  return (1/MixtureDens)*output;
}

//////////////////////////////////////////////////////////////////////////
//                 UT APPROXIMATION ENTROPY GMMS                        //
//                                                                      //
//  ARGS:                                                               //
//                                                                      //
//  G = Number of the GMM components                                    //
//  pro = mixing proportions                                            //
//  mean = mean   (note that we need a rowvector)                       //
//  sigma = covariance matrix                                           //
//  d = dimension of the data                                           //
//////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]

double EntropyUT(int G,
                 arma::vec pro,
                 arma::mat mean,
                 arma::cube sigma,
                 double d)
{

  arma::mat u(d,d);  // for SVD function
  arma::mat v(d,d);  // for SVD function
  arma::vec D(d);    // for SVD function
  arma::vec m(d);
  double e = 0;
  arma::vec temp1(d);
  arma::vec temp2(d);
  arma::vec temp3(d);
  double en = 0;

  for(int i=0; i<G; i++)
  {
    arma::svd_econ(u,D,v,sigma.slice(i));
    m = mean.col(i);
    for(int j=0; j<d; j++)
    {
      temp3 = (sqrt(d * D(j)) * u.col(j));
      temp1 = m + temp3;
      temp2 = m - temp3;
      e =  e + (mixDensity(temp1.t(), G, pro, mean, sigma, true)) + (mixDensity(temp2.t(), G, pro, mean, sigma, true));
    }
    en = en + pro(i) * e;
    e = 0;
  }
  
  return -1.0/(2*d) * en;
}


//////////////////////////////////////////////////////////////////////////
//                 VAR APPROXIMATION ENTROPY GMMS                       //
//                                                                      //
//  ARGS:                                                               //
//                                                                      //
//  G = Number of the GMM components                                    //
//  pro = mixing proportions                                            //
//  mean = mean   (note that we need a rowvector)                       //
//  sigma = covariance matrix                                           //
//  d = dimension of the data                                           //
//////////////////////////////////////////////////////////////////////////

// [[Rcpp::export]]
double EntropyVAR(int G,
                 NumericVector pro,
                 NumericMatrix mean,
                 arma::cube sigma,
                 int d)
{

  double left = 0;
  double right = 0;
  double D=0;

  for (int i=0; i<G; i++)
  {
    left=0;
    for (int k=0; k<G; k++)
    {
      left = left + (pro[k] * exp(-KLMN(mean(_,i),sigma.slice(i),mean(_,k),sigma.slice(k))));
    }

    right = right + pro[i] * 0.5 * log(pow(2*M_PI*exp(1.0),d) * arma::det(sigma.slice(i)));
  	D = D + pro[i] * log(left);
  }

  return -(D - right);
}


//////////////////////////////////////////////////////////////////////////
//                  SOTE APPROXIMATION ENTROPY GMMS                     //
//                                                                      //
//  ARGS:                                                               //
//                                                                      //
//  data = data                                                         //
//  G = Number of the GMM components                                    //
//  pro = mixing proportions                                            //
//  mean = mean   (note that we need a rowvector)                       //
//  sigma = covariance matrix                                           //
//////////////////////////////////////////////////////////////////////////

//[[Rcpp::export]]
double EntropySOTE(arma::mat data,
                   int G,
                   arma::vec pro,
                   arma::mat mean,
                   arma::cube sigma)
{

  double output = 0;
  double H=0;
  arma::mat trandata = data.t();
  arma::mat tranmean  = mean.t();

  for(int i=0; i<G; i++)
  {
    H = H + pro(i) * mixDensity(tranmean.row(i), G, pro, mean, sigma, true);
    output = output +  (pro(i)/2) * (accu( F(mean.col(i), G, pro, mean, sigma) % sigma.slice(i)));
  }

  return (-H - output);
}


//[[Rcpp::export]]

List EntropyMCapprox(arma::mat data,
                     int G,
                     arma::vec pro,
                     arma::mat mean,
                     arma::cube sigma)
{
  int n = data.n_rows;
  arma::vec out; out.zeros(n);
  double    ent;
  double    se;

  for(int i=0; i<n; i++)
  {
    out(i) = mixDensity(data.row(i), G, pro, mean, sigma, true);
  }

  ent = -1.0*arma::mean(out);
  se  = sqrt(arma::var(out)/n);

  return List::create(Rcpp::Named("Entropy") = ent,
                      Rcpp::Named("se") = se);
}

