#ifndef COVMATRIX_FUNS_10_H
#define COVMATRIX_FUNS_10_H


// covariance functions
#include <RcppArmadillo.h>
#include <iostream>
#include <vector>

using namespace Rcpp;
using namespace arma;
//[[Rcpp::depends(RcppArmadillo)]]


//' Isotropic Matern covariance function, smoothness = 2.5
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs is a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range, nugget)
//' = \eqn{(\sigma^2,\alpha,\tau^2)}, and the covariance function is parameterized
//' as
//' \deqn{ M(x,y) = \sigma^2 (1 + || x - y ||/ \alpha + || x - y ||^2/3\alpha^2 ) exp( - || x - y ||/ \alpha )}
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern25_isotropic(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    int n = locs.n_rows;
    double nugget = covparms( 0 )*covparms( 2 );
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate covariances
    arma::mat covmat(n,n);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        // calculate distance
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        if( d == 0.0 ){
            covmat(i2,i1) = covparms(0);
        } else {
            // calculate covariance            
            covmat(i2,i1) = covparms(0)*(1 + d + d*d/3.0)*std::exp( -d );
        }
        // add nugget
        if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
        // fill in opposite entry
        else { covmat(i1,i2) = covmat(i2,i1); }
    }}

    return covmat;
}

//' @describeIn exponential_isotropic Derivatives of isotropic
//' matern covariance function with smoothness 2.5
// [[Rcpp::export]]
arma::cube d_matern25_isotropic(arma::vec covparms, arma::mat locs ){

    int dim = locs.n_cols;
    int n = locs.n_rows;
    //double nugget = covparms( 0 )*covparms( 2 );
    // create scaled locations
    mat locs_scaled(n,dim);
    for(int j=0; j<dim; j++){ 
        for(int i=0; i<n; i++){
            locs_scaled(i,j) = locs(i,j)/covparms(1);
        }
    }
    // calculate derivatives
    arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
    for(int i1=0; i1<n; i1++){ for(int i2=0; i2<=i1; i2++){
        double d = 0.0;
        for(int j=0; j<dim; j++){
            d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
        }
        d = pow( d, 0.5 );
        
        dcovmat(i1,i2,0) += (1 + d + d*d/3.0)*std::exp(-d);
        dcovmat(i1,i2,1) += 
            covparms(0)*std::exp(-d)*d*d/(3.0*covparms(1))*( d + 1.0 );
        if( i1 == i2 ){ // update diagonal entry
            dcovmat(i1,i2,0) += covparms(2);
            dcovmat(i1,i2,2) += covparms(0); 
        } else { // fill in opposite entry
            for(int j=0; j<covparms.n_elem; j++){
                dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
            }
        }
    }}

    return dcovmat;
}



// //' Matern covariance function, smoothess = 2.5, different range parameter for each dimension
// //'
// //' From a matrix of locations and covariance parameters of the form
// //' (variance, range_1, ..., range_d, nugget), return the square matrix of
// //' all pairwise covariances.
// //' @param locs A matrix with \code{n} rows and \code{d} columns.
// //' Each row of locs is a point in R^d.
// //' @param covparms A vector with covariance parameters
// //' in the form (variance, range_1, ..., range_d, nugget)
// //' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
// //' containing the covariance between observations at \code{locs[i,]} and
// //' \code{locs[j,]}.
// //' @section Parameterization:
// //' The covariance parameter vector is (variance, range_1, ..., range_d, nugget).
// //' The covariance function is parameterized as
// //' \deqn{ M(x,y) = \sigma^2 (1 + || D^{-1}(x - y) || + || D^{-1}(x - y) ||^2/3.0) exp( - || D^{-1}(x - y) || ) }
// //' where D is a diagonal matrix with (range_1, ..., range_d) on the diagonals.
// //' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
// //' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// // [[Rcpp::export]]
// arma::mat matern25_scaledim(arma::vec covparms, arma::mat locs ){
// 
//     int dim = locs.n_cols;
// 
//     if( covparms.n_elem - 2 != dim ){
//         stop("length of covparms does not match dim of locs");
//     }
//             
//     int n = locs.n_rows;
//     double nugget = covparms( 0 )*covparms( dim + 1 );
// 
//     // create scaled locations
//     mat locs_scaled(n,dim);
//     for(int j=0; j<dim; j++){ 
//         for(int i=0; i<n; i++){
//             locs_scaled(i,j) = locs(i,j)/covparms(1+j);
//         }
//     }
// 
//     // calculate covariances
//     arma::mat covmat(n,n);
//     for(int i1 = 0; i1 < n; i1++){ for(int i2 = 0; i2 <= i1; i2++){
//             
//         // calculate distance
//         double d = 0.0;
//         for(int j=0; j<dim; j++){
//             d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
//         }
//         d = pow( d, 0.5 );
//         
//         if( d == 0.0 ){
//             covmat(i2,i1) = covparms(0);
//         } else {
//             // calculate covariance            
//             covmat(i2,i1) = covparms(0)*(1 + d + d*d/3.0)*exp(-d);
//         }
//         // add nugget
//         if( i1 == i2 ){ covmat(i2,i2) += nugget; } 
//         // fill in opposite entry
//         else { covmat(i1,i2) = covmat(i2,i1); }
//     }}
//     return covmat;
// }

//' Matern covariance function, smoothess = 2.5, different range parameter for each dimension
//'
//' From a matrix of locations and covariance parameters of the form
//' (variance, range_1, ..., range_d, nugget), return the square matrix of
//' all pairwise covariances.
//' @param locs A matrix with \code{n} rows and \code{d} columns.
//' Each row of locs is a point in R^d.
//' @param covparms A vector with covariance parameters
//' in the form (variance, range_1, ..., range_d, nugget)
//' @return A matrix with \code{n} rows and \code{n} columns, with the i,j entry
//' containing the covariance between observations at \code{locs[i,]} and
//' \code{locs[j,]}.
//' @section Parameterization:
//' The covariance parameter vector is (variance, range_1, ..., range_d, nugget).
//' The covariance function is parameterized as
//' \deqn{ M(x_i,y_j) = \sigma^2 \prod_{l=1}^p (1 + \sqrt{5}|| \lambda_l^{-1}(x_{il} - y_{jl}) || + 5|| \lambda_l(x_{il} - y_{jl}) ||^{2}/3 ) exp( - \sqrt{5}|| \lambda_l^{-1}(x_{il} - y_{jl}) || ) }
//' where \eqn{ \boldsymbol{\lambda} } is a vector of range parameters \eqn{ (\lambda_1, ..., \lambda_p) }.
//' The nugget value \eqn{ \sigma^2 \tau^2 } is added to the diagonal of the covariance matrix.
//' NOTE: the nugget is \eqn{ \sigma^2 \tau^2 }, not \eqn{ \tau^2 }. 
// [[Rcpp::export]]
arma::mat matern25_scaledim(arma::vec covparms, arma::mat locs ){
 
 int dim = locs.n_cols;
 
 if( covparms.n_elem - 2 != dim ){
   stop("length of covparms does not match dim of locs");
 }
 
 int n = locs.n_rows;
 double nugget = covparms( 0 )*covparms( dim + 1 );
 
 // Initialize the final covariance matrix with ones
 arma::mat final_covmat = arma::ones(n, n);
 
 // create scaled locations
 mat locs_scaled(n,dim);
 for(int j=0; j<dim; j++){
   for(int i=0; i<n; i++){
     locs_scaled(i,j) = locs(i,j)/covparms(1+j);
   }
 }
 
 for(int j = 0; j < dim; j++){
   // create covariance matrix for each dimension and multiply elementwise
   arma::mat covmat(n,n);
   
   for(int i1 = 0; i1 < n; i1++){ for(int i2 = 0; i2 <= i1; i2++){
     // calculate distance
     double d = 0.0;
     d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
     d = pow( d, 0.5 );
     
     if( d == 0.0 ){
       covmat(i2,i1) = 1.0;
     } else {
       // calculate covariance
       const double sqrt_5 = sqrt(5.0);            
       // covmat(i2,i1) = covparms(0)*(1 + sqrt_3*d)*exp(-sqrt_3*d);
       covmat(i2,i1) = (1 + sqrt_5*d + 5.0*pow( d, 2.0)/3.0 )*exp(-sqrt_5*d);
     }
     // add nugget (do not add nugget here, add at end)
     if( i1 == i2 ){ covmat(i2,i2) += 0; } 
     // fill in opposite entry
     else { covmat(i1,i2) = covmat(i2,i1); }
   }}
   
   final_covmat %= covmat;
   
 }
 final_covmat *= covparms(0);
 final_covmat += nugget*arma::mat(n,n,arma::fill::eye);
 
 return (final_covmat);
}

// //' @describeIn matern25_scaledim Derivatives with respect to parameters
// // [[Rcpp::export]]
// arma::cube d_matern25_scaledim(arma::vec covparms, arma::mat locs ){
// 
//     int dim = locs.n_cols;
//     if( covparms.n_elem - 2 != dim ){
//         stop("length of covparms does not match dim of locs");
//     }
//             
//     int n = locs.n_rows;
//     double nugget = covparms( 0 )*covparms( dim + 1 );
// 
//     // create scaled locations
//     mat locs_scaled(n,dim);
//     for(int j=0; j<dim; j++){ 
//         for(int i=0; i<n; i++){
//             locs_scaled(i,j) = locs(i,j)/covparms(1+j);
//         }
//     }
// 
//     // calculate derivatives
//     arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
//     for(int i2=0; i2<n; i2++){ for(int i1=0; i1<=i2; i1++){
//         
//         double d = 0.0;
//         for(int j=0; j<dim; j++){
//             d += pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
//         }
//         d = pow( d, 0.5 );
//         
//         double cov;        
//         if( d == 0.0 ){
//             cov = covparms(0);
//             dcovmat(i1,i2,0) += 1.0;
//         } else {
//             cov = covparms(0)*(1 + d + d*d/3.0)*exp(-d);
//             // variance parameter
//             dcovmat(i1,i2,0) += cov/covparms(0);
//             // range parameters
//             for(int j=0; j<dim; j++){
//                 double dj2 = pow( locs_scaled(i1,j) - locs_scaled(i2,j), 2.0 );
//                 dcovmat(i1,i2,j+1) += 
//                     covparms(0)*exp(-d)*dj2/covparms(j+1)/3.0*( 1 + d );
//             }
//         }
//         if( i1 == i2 ){ // update diagonal entry
//             dcovmat(i1,i2,0) += covparms(dim+1);
//             dcovmat(i1,i2,dim+1) += covparms(0); 
//         } else { // fill in opposite entry
//             for(int j=0; j<covparms.n_elem; j++){
//                 dcovmat(i2,i1,j) = dcovmat(i1,i2,j);
//             }
//         }
//     }}
//     return dcovmat;
// }

//' @describeIn matern25_scaledim Derivatives with respect to parameters
 // [[Rcpp::export]]
 arma::cube d_matern25_scaledim(arma::vec covparms, arma::mat locs){
   
   int dim = locs.n_cols;
   if( covparms.n_elem - 2 != dim ){
     stop("length of covparms does not match dim of locs");
   }
   
   // if(covmat.empty()){
   //   arma::mat R = matern15_scaledim(covparms, locs);
   // }else{
   //   arma::mat R = covmat;
   // }
   
   int n = locs.n_rows;
   double nugget = covparms( 0 )*covparms( dim + 1 );
   
   arma::mat R = matern25_scaledim(covparms, locs);
   // need the covariance matrix without nugget or variance
   R = (R - nugget*arma::mat(n,n,arma::fill::eye))/covparms(0);
   
   // Create array of absolute difference matrices of locs
   arma::cube R0 = arma::cube(n,n,dim, fill::zeros);
   for(int k = 0; k < dim; k++){
     arma::vec col_k = locs.col(k);
     
     for(int i = 0; i < n; i++){
       for(int j = i; j < n; j++){
         double abs_diff = std::abs(col_k(i) - col_k(j));
         R0(i, j, k) = abs_diff;
         R0(j, i, k) = abs_diff;
       }
     }
   }
   
   // calculate derivatives
   arma::cube dcovmat = arma::cube(n,n,covparms.n_elem, fill::zeros);
   for(int k = 0; k < covparms.n_elem; k++){
     if(k == 0){
       dcovmat.slice(k) = R;
     } else if(k == covparms.n_elem - 1){
       // dcovmat.slice(k) = covparms(0)*arma::mat(n,n,arma::fill::eye);
       // For the gradient, omit the variance term to be consistent with RobustGaSP
       dcovmat.slice(k) = arma::mat(n,n,arma::fill::eye);
     } else{
       arma::mat R0_k;
       arma::mat R0_k_2;
       R0_k = R0.slice(k-1);
       R0_k_2 = arma::square(R0_k);
       
       const double sqrt_5 = sqrt(5.0);
       arma::mat matOnes = arma::mat(n,n,arma::fill::ones);
       arma::mat part1 = sqrt_5*R0_k + 10.0/3.0*R0_k_2/covparms(k);
       arma::mat part2 = matOnes + sqrt_5*R0_k/covparms(k) + 5.0/3.0*R0_k_2/pow( covparms(k), 2.0 );
       dcovmat.slice(k) = (part1/part2 - sqrt_5*R0_k)%(R*covparms(0))*(-1/pow(covparms(k), 2.0)); // need to divide by range parameter (chain rule)
     }
   }
   return dcovmat;
}





#endif

