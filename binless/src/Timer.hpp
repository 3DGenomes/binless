#ifndef TIMER_HPP
#define TIMER_HPP

#include <Rcpp.h>
#include <vector>
#include <string>
#include <ctime>
#include <unordered_map>

#include "util.hpp"

class ClockTimer {
public:
    ClockTimer(const std::string& timer_name) : timer_name_(timer_name), timer_on_(false) {}
    
    void stop_timer() {
        if (timer_on_) {
            times_[current_sample_] += double(std::clock() - start_)/CLOCKS_PER_SEC;
            timer_on_ = false;
        }
    }
    
    void print_times() {
        stop_timer();
        for (auto it = times_.cbegin(); it != times_.cend(); ++it) {
            Rcpp::Rcout << timer_name_ << " timer : " << it->first << " : " << it->second << " s\n";
        }
    }
    
    void start_timer(const std::string& sample_name) {
        stop_timer();
        current_sample_ = sample_name;
        timer_on_ = true;
        start_ = std::clock();
    }
    
private:
    const std::string timer_name_;
    bool timer_on_;
    std::unordered_map<std::string, double> times_;
    std::unordered_map<std::string, unsigned> order_;
    std::string current_sample_;
    std::clock_t start_;
};

class NoTimer {
public:
    NoTimer(const std::string&) {}
    void stop_timer() {}
    void print_times() {}
    void start_timer(const std::string&) {}
};



typedef NoTimer Timer;
//typedef ClockTimer Timer;


#endif

