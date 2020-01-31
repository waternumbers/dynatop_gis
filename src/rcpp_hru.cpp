#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// #include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

List fN(int k, NumericVector rst_prop){
  //Rcout << "in fN ffss" << "\n";
  
  int nc = rst_prop(1);
  int nr = rst_prop(0);
  double dx = rst_prop(2);
  double dxd = dx*sqrt(2);
  double cl = dx /(1+sqrt(2));
  
  int row = (k/nc); // row - should be integer division
  int col = k - row*nc; //column

  IntegerVector new_col = {col-1,col,col+1,col-1,col+1,col-1,col,col+1};
  IntegerVector new_row = {row-1,row-1,row-1,row,row,row+1,row+1,row+1};
  IntegerVector new_idx = new_row*nc + new_col;
  NumericVector new_dx = {dxd,dx,dxd,dx,dx,dxd,dx,dxd};
  NumericVector new_cl = {cl,cl,cl,cl,cl,cl,cl,cl}; // this needs changing

  LogicalVector is_valid = (new_col>-1) & (new_col < nc) & (new_row>-1) & (new_col<nr);

  //Rcout << new_col.length() << "\n";
  //Rcout << new_row.length() << "\n";
  //Rcout << new_idx.length() << "\n";
  //Rcout << new_dx.length() << "\n";
  //Rcout << new_cl.length() << "\n";
  //Rcout << is_valid.length() << "\n";
  
  return List::create(Named("idx") = new_idx[is_valid],
		      Named("dx") = new_dx[is_valid],
		      Named("cl") = new_cl[is_valid]);
}


//' cpp wrapper function for computation of redistribution matrices
//' 
//' @param dem Digital elevation model
//' @param grad average downslope gradient
//' @param land_area Area of land surface in pixel
//' @param channel_area Area of channel surface in pixel
//' @param channel_id UID of channel present in the pixel
//' @param hillslope_id hillslope class of each pixel
//' @param offset - difference between cell index of adjacent cells and current cell index - clockwise from top left
//' @param dx distance between cell centres - from top left in clockwise direction
//' @param cl contour length - from top left in a clockwise direction. The 9th value is used for cells split beteen land and channel
//' @param max_index maximum value of the hillslope and channel id's
//' @return list of hillslope and channel properties
//'
// [[Rcpp::export]]
void rcpp_hru(NumericVector dem,
	      IntegerVector hillslope_id,
	      IntegerVector channel_id,
	      NumericVector land_area,
	      NumericVector channel_area,
	      NumericVector grad,
	      NumericVector atb,
	      arma::sp_mat W,
	      NumericVector area,
	      NumericVector s_bar,
	      NumericVector atb_bar,
	      NumericVector rst_prop){

  //Rcout << "got here" << "\n";

  for( int i=0; i<dem.length(); i++){

    
    //Rcpp::Rcout << dem(i) << "\n";
    if( !(NumericVector::is_na(dem(i))) ){
      Rcout << i << "\n";
      // then a hillslope or channel cell
      bool is_channel = !(IntegerVector::is_na( channel_id(i) ));
      bool has_land = land_area(i) > 0 ;
      
      //Rcpp::Rcout << is_channel << "\n";
      //Rcpp::Rcout << channel_id(i) << "\n";
      //Rcpp::Rcout << has_land << "\n";
      //Rcpp::Rcout << land_area(i) << "\n";
      if( is_channel ){
	//Rcout << "in channel" << "\n";
	area( channel_id(i)-1 ) += channel_area(i);
	if( has_land ){
	  area( hillslope_id(i)-1 ) += land_area(i);
	  s_bar( hillslope_id(i)-1 ) += grad(i)*land_area(i);
	  atb_bar( hillslope_id(i)-1 ) += atb(i)*land_area(i);
	  //Rcout << "assigning to W" << "\n";
	  W(channel_id(i)-1,hillslope_id(i)-1) += land_area(i);
	  
	}
      }else{
	//Rcout << "in the not channel part" << "\n";
	
	area( hillslope_id(i)-1 ) += land_area(i);
	s_bar( hillslope_id(i)-1 ) += grad(i)*land_area(i);
	atb_bar( hillslope_id(i)-1 ) += atb(i)*land_area(i);
	
	List j = fN(i,rst_prop);
	
	//Rcout << "Got neighbours" << "\n";
	if( i == 1946138 ){
	  IntegerVector idx = j["idx"];
	  Rcout << idx << "\n";
	}

	NumericVector cl = j["cl"];
	NumericVector dx = j["dx"];
	IntegerVector idx = j["idx"];
	NumericVector ddn = dem[idx];
	NumericVector ddi (idx.length(),dem(i));

	//Rcout << idx << "\n";
	//Rcout << dx << "\n";
	//Rcout << cl << "\n";

	//Rcout << ddn << "\n";
	//Rcout << ddi << "\n";
	
	NumericVector grd = (ddi - ddn) / dx ;
	//Rcout << grd << "\n";
		
	LogicalVector is_lower = (grd>0) ;
	is_lower[is_na(is_lower)] = 0; // set na values to false
	
	//LogicalVector is_shit = is_na(grd);
	//Rcout << is_shit << "\n";				      
	//Rcout << is_lower << "\n";
	if( any(is_lower).is_true() ){

	  NumericVector trim_cl = cl[is_lower];
	  double sum_cl = Rcpp::sum( trim_cl );
	  NumericVector frc = grd*cl/sum_cl;
	  
	  for(int k=0; k< idx.length(); k++){
	    if( is_lower(k) ){
	      if( land_area(idx(k)) > 0 ){
		W(hillslope_id(idx(k))-1,hillslope_id(i)-1) += frc(k)*land_area(i);
	      }else{
		W(channel_id(idx(k))-1,hillslope_id(i)-1) += frc(k)*land_area(i);
	      }
	    }
	  }
	}
      }
    }
    Rcpp::checkUserInterrupt();
  }
  
  //return area;
}

    
	  
	    
				     
//       // then a hillslope element
//       area( hillslope_id(i) ) += land_area(i); // add area
//       av_grad( hillslope_id(i) ) += land_area(i)*grad(i); //weighted sum of average gradient
//       av_atanb( hillslope_id(i) ) += land_area(i)*atanb(i); // weighted sum of a/tanb
//       // redistribute downstream
//       if( !(IntegerVector::is_na(channel_id(i))) ){
// 	// then part cell with channel - all flow to channel
// 	W( channel_id(i), hillslope_id(i) ) += land_area(i);
//       }else{
// 	// work out downslope gradients
// 	IntegerVector ngh = offset + i;
// 	LogicalVector in_range = (ngh<dem.length()) & (ngh>-1);
// 	NumericVector gl(8,0.0); // gradient
// 	for(int j=0;j<ngh.length();j++){	
// 	  if( in_range(j) ){
// 	    if( !(NumericVector::is_na(dem(ngh(j)))) &&
// 		(land_area(i) > 0) &&
// 		dem(ngh(j)) < dem(i) ){
// 	      gl(j) = cl(j)*( (dem(i) - dem(ngh(j))) / dx(j) ); //tan(beta)*L
// 	    }
// 	  }
// 	}
// 	double sum_gl = sum(gl);
// 	// distribute downstream
// 	for(int j=0;j<ngh.length();j++){	
// 	  if( in_range(j) ){
// 	    // drain to hillslope if next downstream has hillslope_id
// 	    if( !(IntegerVector::is_na(hillslope_id(ngh(j)))) &&
// 		dem(ngh(j)) < dem(i) ){
// 	      W( hillslope_id(ngh(j)), hillslope_id(i) ) += gl(j)*land_area(i)/sum_gl;
// 	    }else{
// 	      // if not got a hillslope_id then drain to channel if available
// 	      if( !(IntegerVector::is_na(channel_id(ngh(j)))) &&
// 		  dem(ngh(j)) < dem(i) ){
// 		W( channel_id(ngh(j)), hillslope_id(i) ) += gl(j)*land_area(i)/sum_gl;
// 	      }
// 	    }
// 	  }
// 	}
//       }
//     }
//     // handle if has a channel_id
//     if( !(IntegerVector::is_na(channel_id(i))) ) {
//       // then a hillslope element
//       area( channel_id(i) ) += channel_area(i); // add area
//     }
//   }
  
//   return List::create(Named("area") = area,
// 		      Named("av_grad") = av_grad,
// 		      Named("av_atanb") = av_atanb
// 		      Named("W") = W);
// }
