// sequence.h

#ifndef __SEQUENCE_H_INCLUDED__ 
#define __SEQUENCE_H_INCLUDED__

// [[Rcpp::depends(BH)]]
//#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include <string>
#include <iostream>
#include <vector>
#include <memory>
#include <sstream>

// Needs apparently C++11
//#include <cstdint>

#include <boost/shared_ptr.hpp>
#include <boost/cstdint.hpp>

#define SEQUENCE_BITS_PER_BASE 0x3
#define SEQUENCE_BIT_MASK 0x7

class Sequence { 

  private:
    Sequence(); // no default constructor

  protected:
    static const uint64_t  Triplets[];

    static uint64_t parse(const std::string&);
    
    unsigned int countGC() const;
    
    uint64_t      value_;
    size_t        n_;

  public:

    static const char Bases[];
    static const uint64_t BASE_MISREAD;
    static std::vector<uint64_t> REAL_BASES;

    Sequence(const uint64_t &value, const size_t &n);
    Sequence(const std::string&);

    virtual ~Sequence();

    virtual unsigned int at(const size_t index) const;
    virtual bool operator==(const Sequence &other) const;

    Sequence append(const Sequence&) const;
    Sequence append(const uint64_t) const;
    Sequence substitute(const size_t index, const uint64_t base) const;
    Sequence insert(const size_t index, const uint64_t base) const;
    Sequence remove(const size_t index) const;
    Sequence truncate(const size_t new_n) const;

    size_t        length() const;
    uint64_t      value() const;
    bool          isGCContentRight() const;
    bool          containsTriples() const;

    bool          isSelfComplementary() const;
    bool          operator!=(const Sequence &other) const;
    bool          operator<(const Sequence&) const;
    bool          operator() (const Sequence& x, const Sequence& y) const;

    std::string   asString() const;

    friend std::ostream&  operator<<(std::ostream&, Sequence const&);

};

class SequenceComparator {
  public:
    bool operator()(const boost::shared_ptr<Sequence>& s1, const boost::shared_ptr<Sequence>& s2) const
    {
      return *s1 < *s2;
    };

    bool operator()(const Sequence& s1, const Sequence& s2) const
    {
      return s1 < s2;
    }
};

/*
static bool comparePtrToSequence(const boost::shared_ptr<Sequence> &a, const boost::shared_ptr<Sequence> &b) {
    return (*a < *b);
}
*/

#endif 
