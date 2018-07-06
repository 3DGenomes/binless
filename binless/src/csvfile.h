#ifndef USE_BOOST
#define USE_BOOST 0
#endif
#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#if USE_BOOST
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/gzip.hpp>

namespace bio = boost::iostreams;

#endif

class csvfile;

inline static csvfile& endrow(csvfile& file);
inline static csvfile& flush(csvfile& file);

class csvfile
{
#if USE_BOOST
  bio::filtering_ostream fs_;
#else
  std::ofstream fs_;
#endif
  bool is_first_;
  const std::string separator_;
  const std::string escape_seq_;
  const std::string special_chars_;
public:
  csvfile(const std::string filename, const std::string separator = ",")
    : fs_()
  , is_first_(true)
  , separator_(separator)
  , escape_seq_("\"")
  , special_chars_("\"")
  {
#if USE_BOOST
    fs_.push(bio::gzip_compressor());
    fs_.push(bio::file_descriptor_sink(filename));
#else
    fs_.open(filename);
#endif
    fs_.exceptions(std::ios::failbit | std::ios::badbit);
    
    
  }
  
  ~csvfile()
  {
    flush();
#if !USE_BOOST
    fs_.close();
#endif
  }
  
  void flush()
  {
    fs_.flush();
  }
  
  void endrow()
  {
    fs_ << std::endl;
    is_first_ = true;
  }
  
  csvfile& operator << ( csvfile& (* val)(csvfile&))
  {
    return val(*this);
  }
  
  csvfile& operator << (const char * val)
  {
    return write(escape(val));
  }
  
  csvfile& operator << (const std::string & val)
  {
    return write(escape(val));
  }
  
  template<typename T>
  csvfile& operator << (const T& val)
  {
    return write(val);
  }
  
private:
  template<typename T>
  csvfile& write (const T& val)
  {
    if (!is_first_)
    {
      fs_ << separator_;
    }
    else
    {
      is_first_ = false;
    }
    fs_ << std::fixed << std::setprecision(4) << val;
    return *this;
  }
  
  std::string escape(const std::string & val)
  {
    std::ostringstream result;
    result << '"';
    std::string::size_type to, from = 0u, len = val.length();
    while (from < len &&
           std::string::npos != (to = val.find_first_of(special_chars_, from)))
    {
      result << val.substr(from, to - from) << escape_seq_ << val[to];
      from = to + 1;
    }
    result << val.substr(from) << '"';
    return result.str();
  }
};


inline static csvfile& endrow(csvfile& file)
{
  file.endrow();
  return file;
}

inline static csvfile& flush(csvfile& file)
{
  file.flush();
  return file;
}