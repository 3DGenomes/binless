#include <Rcpp.h>
#include <vector>
#include <numeric>

#include "GFLLibrary_trapezoidal.hpp"
#include "gfl_graph_fl.h" //graph_fused_lasso_weight_warm

std::vector<std::vector<int> > trapezoidal_grid_chain(int nrows, int maxdiag) {
  if (maxdiag > (nrows+1)/2) Rcpp::stop("maxdiag is too large for this implementation!\n");
  if (maxdiag <= 1) Rcpp::stop("maxdiag must be at least 2!\n");
  int nfull = nrows-maxdiag+1; //number of rows with maxdiag elements
  int ntotal = nfull*maxdiag + (maxdiag-1)*maxdiag/2; //total number of elements in trapeze
  std::vector<std::vector<int> > chains;
  // full rows
  for (int rowno=1; rowno <= nfull; ++rowno) {
    //fill with increasing sequence
    std::vector<int> current(maxdiag);
    std::iota(current.begin(), current.end(), (rowno-1)*maxdiag); //sequence starts at 0
    chains.push_back(current);
  }
  //last rows with decreasing size
  int csz = maxdiag-1;
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
    for (int i=1; i<previous.size(); ++i) current.push_back(previous[i]+1);
    //last element is computed
    current.push_back(current.back()+colno);
    chains.push_back(current);
  }
  return(chains);
}

void GFLLibrary_trapezoidal::store_trails(int nrows, int maxdiag) {
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

void GFLLibrary_trapezoidal::reset() {
    //setup initial values for a cold start
    counter_ = 0;
    beta_trap_ = std::vector<double>(N_,0);
    z_ = std::vector<double>(tsz_,0);
    u_ = std::vector<double>(tsz_,0);
}

std::vector<double> GFLLibrary_trapezoidal::extract_trapeze(const std::vector<double>& vec) const {
  std::vector<double> ret;
  ret.reserve(N_);
  for (unsigned i=0; i<nrows_; ++i) {
    unsigned idx_base = i*(i+1)/2 + i*(nrows_-i);
    ret.insert(ret.end(), vec.begin()+idx_base, vec.begin()+idx_base+std::min(maxdiag_,nrows_-i));
  }
  return(ret);
}

std::vector<double> GFLLibrary_trapezoidal::fill_triangle(const std::vector<double>& y, const std::vector<double>& w) const {
  //compute average value
  double sum_wy=0;
  double sum_w=0;
  for (unsigned i=0; i<nrows_; ++i) {
    unsigned idx_base = i*(i+1)/2 + i*(nrows_-i);
    for (unsigned d=maxdiag_; i+d<nrows_; ++d) {
      sum_wy += w[idx_base+d]*y[idx_base+d];
      sum_w += w[idx_base+d];
    }
  }
  double avg = sum_wy/sum_w;
  //create return vector, filled with average value
  std::vector<double> ret(nrows_*(nrows_+1)/2, avg);
  //fill in lasso solution on trapeze
  unsigned j=0;
  for (unsigned i=0; i<nrows_; ++i) {
    unsigned idx_base = i*(i+1)/2 + i*(nrows_-i);
    for (unsigned d=0; d<std::min(maxdiag_,nrows_-i); ++d) {
      ret[idx_base+d] = beta_trap_[j++];
    }
  }
  return(ret);
}

void GFLLibrary_trapezoidal::optimize(const std::vector<double>& y, const std::vector<double>& w, double lambda2, double converge) {
    //perform optimization on the C side
    const std::vector<double> y_trap = extract_trapeze(y);
    const std::vector<double> w_trap = extract_trapeze(w);
    double* py = const_cast<double*>(&y_trap[0]);
    double* pw = const_cast<double*>(&w_trap[0]);
    int counter = graph_fused_lasso_weight_warm (N_, py, pw, ntrails_, &trails_[0], &breakpoints_[0],
                                               lambda2, &alpha_, inflate_, ninner_, converge,
                                               &beta_trap_[0], &z_[0], &u_[0]);
    beta_tri_ = fill_triangle(y,w);
    //Rcpp::Rcout << "GFLLibrary_trapezoidal: " << counter << " steps\n";
    counter_ += counter;
}

GFLLibrary_trapezoidal::GFLState_t GFLLibrary_trapezoidal::get_state() const {
    return Rcpp::List::create(_["z"]=z_, _["u"]=u_,  _["alpha"]=alpha_,  _["beta_trap"]=beta_trap_,  _["counter"]=counter_);
}

void GFLLibrary_trapezoidal::set_state(const GFLState_t& state) {
    if (state.containsElementNamed("u") && Rcpp::as<std::vector<double> >(state["u"]).size() == tsz_) {
        z_ = Rcpp::as<std::vector<double> >(state["z"]);
        u_ = Rcpp::as<std::vector<double> >(state["u"]);
        beta_trap_ = Rcpp::as<std::vector<double> >(state["beta_trap"]);
        alpha_ = Rcpp::as<double>(state["alpha"]);
        counter_ = Rcpp::as<unsigned>(state["counter"]);
    }
}

void test_trap(unsigned nrows, unsigned maxdiag) {
  GFLLibrary_trapezoidal gfl(nrows, maxdiag);
  //std::vector<std::vector<int> > ret = trapezoidal_grid_chain(nrows, maxdiag);
  std::vector<double> beta_tri;
  int k = 0;
  Rcpp::Rcout << "\n";
  for (unsigned i=0; i<nrows; ++i) {
    for (unsigned j=i; j<nrows; ++j) {
      beta_tri.push_back(k++);
      Rcpp::Rcout << k << " ";
    }
    Rcpp::Rcout << "\n";
  }
  std::vector<double> w(beta_tri.size(),1.);
  gfl.beta_trap_ = gfl.extract_trapeze(beta_tri);
  std::vector<double> ret = gfl.fill_triangle(beta_tri,w);
  Rcpp::Rcout << "ret has length " << ret.size() << "\n";
  for (unsigned i=0; i<ret.size(); ++i) {
    Rcpp::Rcout << ret[i]+1 << " ";
    Rcpp::Rcout << "\n";
  }
}

