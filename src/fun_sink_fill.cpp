#include <Rcpp.h>
using namespace Rcpp;
//' cpp wrapper function for filling of sinks
//' 
//' @param dem Digital elevation model
//' @param is_channel TRUE is a channel pixel
//' @param delta 3x3 matrix of minimum evelation drop to each adjacent pixel
//' @param max_iter maximum number of iterations
//'
//' @return matrix containing filled dem
//'
// [[Rcpp::export]]
NumericMatrix fun_sink_fill(NumericVector dem, LogicalVector is_channel,
			    NumericVector delta, int max_iter){

  
  int n_sink = 0;
  int n_finite = 0;
  int niter = 0;
  LogicalVector is_sink(dem.length(),true);

  double na_test_val = -10000; // test the for NAN against this - if NAN will return false
    
  // loop
  n_iter = 1;
  while( any(is_sink) && n_iter < max_iter ){
    n_sink=0;
    for(int i=0;i < dem.length(); i++){
      if( is_sink(i) ){
	if( dem(i) > na_test_val & not is_channel(i) ){
	  // work out index of neighbours
	  ngh(0) = i-nc-1;
	  ngh(1) = i-nc;
	  ngh(2) = i-nc+1;
	  ngh(3) = i-1;
	  ngh(4) = i+1;
	  ngh(5) = i+nc-1;
	  ngh(6) = i+nc;
	  ngh(7) = i+nc+1;
	  
	  // set the lowest neighbout value and number of finite neighbours
	  min_ngh = 1e32;
	  n_finite = 0;
	  // loop neighbours
	  for(int j=0; j < 8; j++){
	    if( ngh(j) < dem.length() &&
		ngh(j) > -1 &&
		dem(ngh(j)) > na_test_val &&
		dem(ngh(j)) < dem(i) ){
	      // increase number of finite values
	      n_finite = n_finite + 1;
	      // compute new lowest neighbour
	      n_lower(i) = n_lower(i) + 1;
	    }
	  }

	  // work out if changed

	  // set new value and neighbours to is_sink=true
	  
      if( (n_finite < 8) && (n_lower(i)==0) ){
	// potential edge drain set order to 1
	order(i) = 1;
      }
    }else{
      can_eval(i)=false;
    }
    
  for(int i=0;i < dem.nrow(); i++){
    for(int j=0; j < dem.ncol(); j++){
      // set to can't be filled
      can_eval(i,j) = false;
      // to fill it there must be 9 finite values in the block centred on i,j
      n_finite = 0;
      for(int ii=-1;ii<2;ii++){
	for(int jj=-1;jj<2;jj++){
	  if( i+ii > -1 && i+ii < dem.nrow() && 
	      j+jj > -1 && j+jj < dem.ncol() ){
	    if( dem(i+ii,j+jj) > na_test_val ){
	      n_finite = n_finite + 1;
	    }
	  }
	}
      }
      // see if there are 9 finite values
      if(n_finite==9){
	can_eval(i,j) = true;
      }
    }
  }
  
  // loop to fill
  n_sink = 1;
  while( n_sink > 0 && niter < max_iter){
    n_sink = 0;
    // remember cpp is 0 base
    for(int i=0;i < dem.nrow(); i++){
      for(int j=0; j < dem.ncol(); j++){
	// check it can be filled and isn't a channel
	if( can_eval(i,j)==true && is_channel(i,j)==false){
	  // presume it is a sink
	  bool is_sink = true;
	  double lowest_neighbour = 1e32;
	  
	  // check its status
	  for(int ii=-1;ii<2;ii++){
	    for(int jj=-1;jj<2;jj++){
	      if( i+ii > -1 && i+ii < dem.nrow() && 
		  j+jj > -1 && j+jj < dem.ncol() &&
		  !(ii == 0 && jj == 0) ){
		// value that dem(i,j) but excedd if the gradient in greater then min required
		double neighbour = dem(i+ii,j+jj) + delta(ii+1,jj+1);
		// see if this neighbour is low enough
		if( neighbour <= dem(i,j) ){
		  is_sink = false;
		}
		// if it is a sink then update the possible value
		if( is_sink==true && neighbour < lowest_neighbour ){
		  lowest_neighbour = neighbour;
		}
	      }
	    }
	  }
	  
	  // if it is still a sink
	  if(is_sink==true){
	    dem(i,j) = lowest_neighbour;
	    n_sink = n_sink + 1 ;
	  }
	}
	Rcpp::checkUserInterrupt();
      }
    }
    niter = niter + 1;
    Rcpp::checkUserInterrupt();
    Rcout << "The number of sinks handled in iteration " << niter << " is " << n_sink << std::endl;
    
  }
  
  return dem; 
}
