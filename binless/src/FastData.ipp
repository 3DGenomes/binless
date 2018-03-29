
template<typename Derived>
FastDataCore<Derived>::FastDataCore(const DataFrame& obs, unsigned nbins) :
name_(Rcpp::as<std::vector<unsigned> >(obs["name"])),
bin1_(Rcpp::as<std::vector<unsigned> >(obs["bin1"])),
pos1_(Rcpp::as<std::vector<unsigned> >(obs["pos1"])),
bin2_(Rcpp::as<std::vector<unsigned> >(obs["bin2"])),
pos2_(Rcpp::as<std::vector<unsigned> >(obs["pos2"])),
name_levels_(Rcpp::as<Rcpp::IntegerVector>(obs["name"]).attr("levels")),
bin1_levels_(Rcpp::as<Rcpp::IntegerVector>(obs["bin1"]).attr("levels")),
bin2_levels_(Rcpp::as<Rcpp::IntegerVector>(obs["bin2"]).attr("levels")),
observed_(Rcpp::as<std::vector<unsigned> >(obs["observed"])),
nobs_(Rcpp::as<std::vector<unsigned> >(obs["nobs"])),
distance_(Rcpp::as<std::vector<double> >(obs["distance"])),
nbins_(nbins), ncells_(nbins*(nbins+1)/2), ndatasets_(*std::max_element(name_.begin(), name_.end())), N_(ndatasets_*ncells_),
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
