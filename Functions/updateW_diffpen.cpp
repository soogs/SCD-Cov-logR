// updateW_diffpen.cpp //

// rcpp script for conditional update of W
// created following the update rules provided in the paper (Appendix)

// this is as according to the rcpp file
// updateW_faster_reference.cpp
// i let lasso and glasso be specified per component

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
double soft (double x, double lambda){	
	double x2;
	double result;
	double sgn;

	x2 = std::abs(x) - lambda;

	if (x2 < 0){ 
		x2 = 0;
	} 

	if (x < 0){
		sgn = -1;
	} else {
		sgn = 1;
	}

	result = sgn * x2;

	return result;
}

// [[Rcpp::export]]
double glasso_norm (
	arma::mat x,
	Rcpp::List blockindex,
	int R,
  arma::vec glasso){

	double ssq;
	arma::mat result; 
	int nblock;
	double l2norm;

	nblock = blockindex.size();
	
	int r;
	int i;

	l2norm = 0;

	for (r = 0; r < R; r++){
		for (i = 0; i < nblock; i++){
			Rcpp::NumericVector block;
			int block_min;
			int block_max;
			arma::mat inblock;
			double blocksize;

			block = blockindex(i);
			block_min = min(block);
			block_max = max(block);

			block_min = block_min - 1;
			block_max = block_max - 1;

			inblock = x(arma::span(block_min, block_max), r);

			ssq = sqrt(as_scalar(trans(inblock) * inblock));

			blocksize = block.size();

			l2norm = l2norm + glasso(r) * ssq * sqrt(blocksize);
		}
	}


	return l2norm;

}

// [[Rcpp::export]]
Rcpp::List checkit_reference (
	arma::mat X,
	arma::mat W,
	arma::mat Px,
  arma::mat Py,
  double py0,
	arma::vec lasso,
	arma::vec glasso,
  double ridge,
	Rcpp::List blockindex,
	int R,
  double alpha,
  arma::mat z,
  arma::mat qq){

  	arma::mat W0;
  	arma::mat W_new;
  	double nblock;

    arma::mat Xk;
  	arma::mat PXk;
  	arma::mat PXky;
  	arma::mat PXkPX;
  	arma::mat PXkPXk;
  	arma::mat Wk;
  	arma::mat PXkrk;
    arma::mat XkXk;
    arma::vec XTXdiag;
    arma::mat XTX;
    arma::vec sumPRESS;


    arma::mat log_res_k;
    arma::mat k_soft_log;

  	arma::mat Wk_new;
  	arma::mat k_soft_pca;
  	arma::mat k_soft;
  	arma::mat PX; 
  	arma::mat Xki;
  	arma::mat Xkirki;
	  int J_k;


    int W0h_index;
    double Wkh_diff;
    arma::mat z_py0;
    arma::mat PyXW0;
    Rcpp::IntegerVector block_index;



  	// define initial objects
  	W0 = W;
  	W_new = W;

   	nblock = blockindex.size();
   	PX = arma::kron(Px, X);

   	int r;
   	int k;

    XTX = trans(X) * X;

    XTXdiag = XTX.diag();

   	for (r = 0; r < R; r++){

      sumPRESS = X * (Px.col(r) - W_new.col(r));

   		for (k = 0; k < nblock; k++){
   			Rcpp::NumericVector block;
			  int block_min;
			  int block_max;

			  block = blockindex(k);
			  block_min = min(block);
			  block_max = max(block);

			  block_min = block_min - 1;
			  block_max = block_max - 1;

        Xk = X(arma::span(), arma::span(block_min, block_max));

			  // pre-define before the loop
   			// PXk = arma::kron(Px.col(r), Xk);

   			// PXky = trans(PXk) * arma::vectorise(X);

   			// PXkPX = trans(PXk) * PX;

   			XkXk = trans(Xk) * Xk;

   			J_k = block.size();
      
   			// group k: pca loss 
   			Wk = arma::vectorise(W0(arma::span(block_min, block_max), r));

   			W0 = arma::vectorise(W0);

        PXkrk = arma::vectorise(trans(sumPRESS) * Xk) + XkXk * Wk;

   			k_soft_pca = (2 * (1 - alpha)) * PXkrk;

        // group k: logistic regression loss

        arma::mat PyX;

        PyX = arma::kron(trans(Py), X);

        log_res_k = X * W_new * Py - Py(r,0) * (Xk * Wk);
        
        arma::mat inside_diag;

        inside_diag = arma::diagmat(z - log_res_k - py0);

        k_soft_log = (alpha) * trans(qq) * inside_diag * (Py(r,0) * Xk);

        k_soft_log = arma::vectorise(k_soft_log);

        // add the two losses and check if the group is 0 vector
        arma::mat k_to_soft;

        k_to_soft = k_soft_log + k_soft_pca;

        k_soft = k_to_soft;

        int k_soft_rows = k_soft.n_rows;

        // soft thresholding via for loop
        for (int i = 0; i < k_soft_rows; i++){

          double to_soft = k_to_soft(i,0);

          double softened = soft(to_soft, lasso(r));

          k_soft(i,0) = softened;
        }

        // calculate l2 norm 
        double groupnorm = sqrt(accu(square(k_soft)));

        if (groupnorm < (glasso(r) * sqrt(J_k))){
          // if the l2 norm is too small
          
          sumPRESS += Xk * Wk;

          Wk_new = arma::zeros( Wk.n_rows, Wk.n_cols );

          // update the entire W_new
          W_new(arma::span(block_min, block_max),r) = Wk_new;

          W0 = W_new;
        } else {
          
          // pre-define things prior to loop
          
          z_py0 = z - py0;
          PyXW0 = X * W_new * Py;

          block_index = Rcpp::seq(block_min, block_max);

          // coordinate descent      
          for (int h = 0; h < J_k; h++){

          arma::mat part1;
          arma::mat Xkh;
          arma::mat XkhXkh;

          Xkh = Xk.col(h);
          XkhXkh = trans(Xkh) * Xkh;
            
          part1 = - (alpha) * pow(Py(r,0),2) * trans(qq) * square(Xkh) - 
          (2 * (1-alpha)) * XkhXkh - 
          sqrt(J_k) * glasso(r) / sqrt(accu(square(Wk)));

          arma::mat part2;
          arma::mat part2_res_h;
          arma::mat part2_diag;


          part2_res_h = PyXW0 - Py(r,0) * Wk(h,0) * (Xkh);

          part2_diag = arma::diagmat(z_py0 - part2_res_h);

          part2 = alpha * trans(qq) * part2_diag * (Py(r,0) * Xkh);


          arma::mat PXkh;
          arma::mat part3_res_h;
          arma::mat part3;

          //PXkh = arma::kron(Px.col(r), Xkh);
         
          //part3_res_h = trans(PXkh) * X_PXW0 + XkhXkh * Wk(h,0);

          //t(sumPRESS) %*% Xk[,h] + Wk[h] * XTXdiag[(blockindex[[k]])[h]]

          part3_res_h = trans(sumPRESS) * Xkh + Wk(h,0) * XTXdiag(block_index(h));

          part3 = (2 * (1-alpha)) * part3_res_h;

          arma::mat bigger;
          arma::mat zeromat;

          zeromat = Wk;
          zeromat.fill(0);

          bool Wk_zero = approx_equal(Wk, zeromat, "absdiff", 1e-7);

          double soft_Wkh;
          double Wkh_new;
          double to_soft_Wkh = part2(0,0) + part3(0,0);

          if (Wk_zero){
            Wkh_new = 0;
          } else {
            soft_Wkh = soft(to_soft_Wkh, lasso(r));

            Wkh_new = soft_Wkh / (-1 * part1(0,0));
          }

          // updating the weights matrix
          Wk_new = Wk;

          Wk_new(h,0) = Wkh_new;
              
          W_new(arma::span(block_min, block_max),r) = Wk_new;


          // check convergence of W
          if (Wkh_new != 0){

            sumPRESS += (Wk(h,0) - Wkh_new) * Xkh;

            bool conv_ingroup = approx_equal(Wk, Wk_new, "absdiff", 1e-7);

          } else{

            sumPRESS += Xkh * Wk(h,0);

          }
           
          W0h_index = (r) * X.n_cols + block_index(h);

          Wkh_diff = Wkh_new - Wk(h,0);

          PyXW0 = PyXW0 + PyX.col(W0h_index) * Wkh_diff;
          
          Wk = arma::vectorise(W_new(arma::span(block_min, block_max), r));
          
        }

        W0 = W_new;

      }
    }
  }
  Rcpp::List result;

  result["Wnew"] = W_new;
  result["W0"] = W0;
  result["Wk"] = Wk;
  result["r"] = r;
  result["k"] = k;
  result["W0h_index"] = W0h_index;
  result["PyXW0"] = PyXW0;
  result["block_index"] = block_index;

  return result;
}
