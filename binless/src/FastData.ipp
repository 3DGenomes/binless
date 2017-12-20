
template<typename Derived>
FastDataCore<Derived>::FastDataCore(const DataFrame& obs, unsigned nbins) :
name_(Rcpp::as<std::vector<unsigned> >(obs["name"])),
bin1_(Rcpp::as<std::vector<unsigned> >(obs["bin1"])),
bin2_(Rcpp::as<std::vector<unsigned> >(obs["bin2"])),
observed_(Rcpp::as<std::vector<unsigned> >(obs["observed"])),
distance_(Rcpp::as<std::vector<double> >(obs["distance"])),
nbins_(nbins), ncells_(nbins*(nbins+1)/2), ndatasets_(*std::max_element(name_.begin(), name_.end())), N_(ndatasets_*ncells_),
log_biases_(std::vector<double>(nbins_,0)),
log_decay_(std::vector<double>(nbins_,0)),
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

template<typename Derived>
std::vector<double> FastDataCore<Derived>::get_log_expected() const {
    std::vector<double> log_expected;
    log_expected.reserve(get_N());
    std::vector<double> dbeta = static_cast<const Derived*>(this)->get_log_signal();
    for (unsigned i=0; i<get_N(); ++i) {
        unsigned bin1 = bin1_[i]-1; //offset by 1 for vector indexing
        unsigned bin2 = bin2_[i]-1;
        double bi = log_biases_[bin1];
        double bj = log_biases_[bin2];
        double decay = log_decay_[bin2-bin1];
        double beta = dbeta[i];
        unsigned name = name_[i]-1;
        double exposure = exposures_[name];
        log_expected.push_back(bi + bj + decay + beta + exposure);
    }
    return log_expected;
}
