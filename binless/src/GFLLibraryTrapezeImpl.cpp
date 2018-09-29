#include <Rcpp.h>
#include <vector>
#include <numeric>

#include "GFLLibraryTrapezeImpl.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm

std::vector<std::vector<int> > trapezoidal_grid_chain(int nrows, int maxdiag) {
  if (maxdiag > (nrows+1)/2) Rcpp::stop("maxdiag is too large for this implementation!\n");
  if (maxdiag <= 1) Rcpp::stop("maxdiag must be at least 2!\n");
  int nfull = nrows-maxdiag+1; //number of rows with maxdiag elements
  //int ntotal = nfull*maxdiag + (maxdiag-1)*maxdiag/2; //total number of elements in trapeze
  std::vector<std::vector<int> > chains;
  // full rows
  for (int rowno=1; rowno <= nfull; ++rowno) {
    //fill with increasing sequence
    std::vector<int> current(maxdiag);
    std::iota(current.begin(), current.end(), (rowno-1)*maxdiag); //sequence starts at 0
    chains.push_back(current);
  }
  //last rows with decreasing size
  unsigned csz = maxdiag-1;
  std::vector<int> current;
  for (int i=nfull*maxdiag; csz > 1; ++i) {
    current.push_back(i);
    if (current.size() == csz) {
      chains.push_back(current);
      current.clear();
      --csz;
    }
  }
  //first columns of increasing size
  //with Ui+1 = Ui + maxdiag-1 with U1 from 2 to maxdiag
  for (int U1=2; U1<=maxdiag; ++U1) {
    int Ui=U1;
    std::vector<int> current(1,Ui-1);
    for (int i=1; i<U1; ++i) {
      int Uip1 = Ui + maxdiag - 1;
      current.push_back(Uip1-1);
      Ui=Uip1;
    }
    chains.push_back(current);
  }
  //middle columns of constant size maxdiag
  //with Ui+1 = Ui + maxdiag-1
  //with colno from maxdiag+1 to nrows-(maxdiag-2)
  for (int colno=maxdiag+1; colno<=nrows-(maxdiag-2); ++colno) {
    int Ui = chains.back()[1] + 1;
    std::vector<int> current(1,Ui);
    for (int i=1; i<maxdiag; ++i) {
      int Uip1 = Ui + maxdiag - 1;
      current.push_back(Uip1);
      Ui=Uip1;
    }
    chains.push_back(current);
  }
  //final columns of constant size maxdiag
  //with colno from 1 to maxdiag-2 (counting from the end)
  for (int colno=maxdiag-2; colno>=1; --colno) {
    //first maxdiag-1 elements are taken from previous column
    std::vector<int> previous(chains.back());
    std::vector<int> current;
    for (unsigned i=1; i<previous.size(); ++i) current.push_back(previous[i]+1);
    //last element is computed
    current.push_back(current.back()+colno);
    chains.push_back(current);
  }
  return(chains);
}

void GFLLibraryTrapezeImpl::store_trails(int nrows, int maxdiag) {
    const std::vector<std::vector<int> > chains = trapezoidal_grid_chain(nrows, maxdiag);
    trails_.clear();
    breakpoints_.clear();
    for (std::vector<std::vector<int> >::const_iterator it = chains.begin() ;
         it != chains.end(); ++it) {
        if (trails_.size()>0) breakpoints_.push_back(trails_.size());
        trails_.insert(trails_.end(), it->begin(), it->end());
    }
    if (trails_.size()>0) breakpoints_.push_back(trails_.size());
    ntrails_ = breakpoints_.size();
    tsz_ = trails_.size();
}

void GFLLibraryTrapezeImpl::reset() {
    //setup initial values for a cold start
    set_ninner(0);
    beta_impl_ = std::vector<double>(N_,0);
    z_ = std::vector<double>(tsz_,0);
    u_ = std::vector<double>(tsz_,0);
}

std::vector<double> GFLLibraryTrapezeImpl::extract_trapeze(const std::vector<double>& vec) const {
  std::vector<double> ret;
  //Rcpp::Rcout << "max_idx_base:" << (((unsigned)nrows_-1)*((unsigned)nrows_)/2+((unsigned)nrows_-1)+(unsigned)maxdiag_) << "\n";
  ret.reserve(N_);
  for(unsigned i=0; i<(unsigned)nrows_; ++i) {
    unsigned idx_base = i*(i+1)/2 + i*(nrows_-i);
    auto beg_itr = vec.begin();
    auto end_itr = vec.begin();
    std::advance( beg_itr, idx_base );
    std::advance( end_itr, idx_base+std::min<unsigned>((unsigned)maxdiag_,nrows_-i));
    ret.insert(ret.end(), beg_itr, end_itr);
  }
  return(ret);
}

std::vector<double> GFLLibraryTrapezeImpl::fill_triangle(const std::vector<double>& y, const std::vector<double>& w) const {
  //compute average value
  double sum_wy=0;
  double sum_w=0;
  for (unsigned i=0; i<(unsigned)nrows_; ++i) {
    unsigned idx_base = i*(i+1)/2 + i*(nrows_-i);
    for (unsigned d=(unsigned)maxdiag_; i+d<(unsigned)nrows_; ++d) {
      sum_wy += w[idx_base+d]*y[idx_base+d];
      sum_w += w[idx_base+d];
    }
  }
  double avg = sum_wy/sum_w;
  //create return vector, filled with average value
  std::vector<double> ret((unsigned)nrows_*((unsigned)nrows_+1)/2, avg);
  //fill in lasso solution on trapeze
  unsigned j=0;
  for (unsigned i=0; i<(unsigned)nrows_; ++i) {
    unsigned idx_base = i*(i+1)/2 + i*(nrows_-i);
    for (unsigned d=0; d<std::min<unsigned>((unsigned)maxdiag_,(unsigned)nrows_-i); ++d) {
      ret[idx_base+d] = beta_impl_[j++];
    }
  }
  Rcpp::Rcout << "average = " << avg << " beta[1,3] = " << beta_impl_[2] << " " << ret[2] << "\n";
  return(ret);
}

void GFLLibraryTrapezeImpl::optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
    //perform optimization on the C side
    const std::vector<double> y_trap = extract_trapeze(y);
    const std::vector<double> w_trap = extract_trapeze(w);
    double* py = const_cast<double*>(&y_trap[0]);
    double* pw = const_cast<double*>(&w_trap[0]);
    double alpha = get_alpha();
    int counter = graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, &trails_[0], &breakpoints_[0],
                                               lambda2, &alpha, get_inflate(), get_ninner_max(), converge,
                                               &beta_impl_[0], &z_[0], &u_[0]);
    set_alpha(alpha);
    set_beta(fill_triangle(y,w));
    set_ninner(get_ninner() + counter);
}

GFLLibraryTrapezeImpl::GFLState_t GFLLibraryTrapezeImpl::get_state() const {
    return Rcpp::List::create(_["z"]=z_, _["u"]=u_,  _["alpha"]=get_alpha(),  _["beta"]=beta_impl_,
                              _["counter"]=get_ninner(), _["maxdiag"]=maxdiag_);
}

void GFLLibraryTrapezeImpl::set_state(const GFLState_t& state) {
    if (state.containsElementNamed("u") && Rcpp::as<std::vector<double> >(state["u"]).size() == tsz_
          && state.containsElementNamed("maxdiag") && Rcpp::as<int>(state["maxdiag"])==maxdiag_ ) {
        z_ = Rcpp::as<std::vector<double> >(state["z"]);
        u_ = Rcpp::as<std::vector<double> >(state["u"]);
        beta_impl_ = Rcpp::as<std::vector<double> >(state["beta"]);
        set_alpha(Rcpp::as<double>(state["alpha"]));
        set_ninner(Rcpp::as<unsigned>(state["counter"]));
    }
}
