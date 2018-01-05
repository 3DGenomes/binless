
template<typename Derived>
FastDataCore<Derived>::FastDataCore(const DataFrame& obs, unsigned nbins) :
name_(Rcpp::as<std::vector<unsigned> >(obs["name"])),
bin1_(Rcpp::as<std::vector<unsigned> >(obs["bin1"])),
bin2_(Rcpp::as<std::vector<unsigned> >(obs["bin2"])),
observed_(Rcpp::as<std::vector<unsigned> >(obs["observed"])),
distance_(Rcpp::as<std::vector<double> >(obs["distance"])),
nbins_(nbins), ncells_(nbins*(nbins+1)/2), ndatasets_(*std::max_element(name_.begin(), name_.end())), N_(ndatasets_*ncells_),
log_biases_(std::vector<double>(nbins_,0)),
beta_(std::vector<double>(N_,0)),
betahat_(std::vector<double>(N_,0)),
weights_(std::vector<double>(N_,0)),
exposures_(std::vector<double>(ndatasets_,0)) {
    if (name_.size() != N_) Rcpp::stop("Input size should be N(N+1)/2 dense matrix of observed counts for each dataset");
    if (*std::min_element(name_.begin(), name_.end()) != 1) Rcpp::stop("Name column values should start at 1");
    if (*std::min_element(bin1_.begin(), bin1_.end()) != 1 || *std::min_element(bin2_.begin(), bin2_.end()) != 1)
        Rcpp::stop("bins should range from 1 to the number of bins");
    if (*std::max_element(bin1_.begin(), bin1_.end()) != nbins_ || *std::max_element(bin2_.begin(), bin2_.end()) != nbins_)
        Rcpp::stop("bins should range from 1 to the number of bins");
}
