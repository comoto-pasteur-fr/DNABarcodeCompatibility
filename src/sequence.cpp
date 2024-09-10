// in sequence.cpp

#include <boost/foreach.hpp>
#include <boost/cstdint.hpp>
#include <boost/assign/std/vector.hpp>
#include <boost/assign/list_of.hpp>

#include <cstdlib>
#include <string>
#include <memory>
#include <limits>

#include "sequence.h"

using namespace boost::assign; 

//#include "Rcpp.h"
 
const char Sequence::Bases[]    = { 'A', 'G', 'N', 'X', 'X', 'N', 'C', 'T' };
const uint64_t Sequence::Triplets[] = { parse("AAA"), parse("TTT"), parse("GGG"), parse("CCC") };
const uint64_t Sequence::BASE_MISREAD = parse("N");

std::vector<uint64_t> Sequence::REAL_BASES = list_of(parse("A"))(parse("G"))(parse("C"))(parse("T"));

Sequence::Sequence(const uint64_t &value, const size_t &n) : value_(value), n_(n) {
}

Sequence::Sequence(const std::string &sequence)  {
  n_      = sequence.length();
  value_  = parse(sequence);
}

Sequence::~Sequence(){};

uint64_t Sequence::parse(const std::string &sequence)  {
  uint64_t  value = 0;

  for (int seq_position = (sequence.length() - 1); seq_position >= 0; seq_position--) {
    char base = sequence[seq_position];

    // Find position
    int base_int = -1;
    for (size_t i = 0; i < 8; i++) {
      if (Bases[i] == base) {
        base_int = i;
      }
    }

    if (base_int == -1) {
      //Rcpp::Rcerr << "Wrong base submitted" << std::endl;
      return(0);
    }
    //std::cerr << "Parsed base " << base << " to " << base_int << std::endl;
    value = (value << SEQUENCE_BITS_PER_BASE) | base_int;
  }
  return (value);
}

Sequence Sequence::append(const Sequence &appended_seq) const {
  uint64_t val_append = appended_seq.value();
  size_t   n_append   = appended_seq.length();

  uint64_t val_new  = value_ | (val_append << (SEQUENCE_BITS_PER_BASE * n_));
  size_t   n_new    = n_ + n_append;

  return Sequence(val_new, n_new);
}

Sequence Sequence::append(const uint64_t base) const {
  uint64_t val_new  = value_ | (base << (SEQUENCE_BITS_PER_BASE * n_));

  return Sequence(val_new, n_ + 1);
}

Sequence Sequence::substitute(const size_t index, const uint64_t base) const {
  if (index >= n_)
    return Sequence(value_, n_); // return a copy of this sequence

  size_t index_in_val = SEQUENCE_BITS_PER_BASE * index;

  uint64_t mask = ~(SEQUENCE_BIT_MASK << index_in_val);

  return Sequence((value_ & mask) | (base << index_in_val), n_);
}

Sequence Sequence::insert(const size_t index, const uint64_t base) const {
  if (index > n_) // Important, we are saying "> n" (not ">= n") because we may want to insert a base at the last position
    return Sequence(value_, n_);

  size_t index_in_val = SEQUENCE_BITS_PER_BASE * index;

  int mask = (~0) << index_in_val;
  return Sequence(((value_ & mask) << SEQUENCE_BITS_PER_BASE) | (value_ & (~mask)) | (base << index_in_val),n_ + 1);
}

Sequence Sequence::remove(const size_t index) const {
  if (index >= n_)
    return Sequence(value_, n_);

  size_t index_in_val = SEQUENCE_BITS_PER_BASE * index;

  uint64_t mask = (~0) << index_in_val;

  return Sequence((value_ & ~(mask)) | ((value_ & (mask << SEQUENCE_BITS_PER_BASE)) >> SEQUENCE_BITS_PER_BASE),n_ - 1);
}

Sequence Sequence::truncate(const size_t new_n) const {
  if (n_ > new_n) {
    // Need to chip the bits away, that are not needed anymore
    uint64_t new_value = value_;

    uint64_t mask = ~((~0) << (SEQUENCE_BITS_PER_BASE * new_n));
    return(Sequence(new_value & mask,new_n));
  } else {
    return(Sequence(value_,n_));
  }
}

size_t Sequence::length() const {
  return n_;
}

uint64_t Sequence::value() const{
  return value_;
}

unsigned int Sequence::at(const size_t index) const {
  return ((value_ >> (index * SEQUENCE_BITS_PER_BASE)) & SEQUENCE_BIT_MASK);
}

unsigned int Sequence::countGC() const {
  unsigned int count = 0;

  for(size_t index = 0; index < n_; index++) {
    unsigned int base = ((value_ >> (SEQUENCE_BITS_PER_BASE * index)) & SEQUENCE_BIT_MASK);

    if ((base == 0x1) || (base == 0x6))
      count++;
  }

  //std::cerr << "GC count: " << count << std::endl;

  return count;
}

bool Sequence::isGCContentRight() const {
  double gcPercentage = ((double) countGC() / (double) n_);
  return (((0.4 - gcPercentage) < std::numeric_limits<double>::epsilon()) && ((gcPercentage-0.6) < std::numeric_limits<double>::epsilon()));
}

bool Sequence::containsTriples() const {

  if (n_ <= 2) {
    return false;
  }

  // a mask of 6 bits
  unsigned int mask = 0x1FF;

  size_t moves = n_ - 3;

  for (size_t i = 0; i <= moves; i++) {
    uint64_t tmp = (value_ >> (i * SEQUENCE_BITS_PER_BASE) & mask);

    BOOST_FOREACH(uint64_t triplet, Triplets ) {
      if (tmp  == triplet) 
        return true;
    }
  }

  return false;
}

bool Sequence::isSelfComplementary() const {
  uint64_t  mask = ~0 << (SEQUENCE_BITS_PER_BASE * n_);

  // Invert this sequence
  uint64_t  old_value = ~value_ & ~mask;
  uint64_t  new_value = 0;

  // reverse the bases
  for (size_t i = 0; i < n_; i++) {
    // Get base at this position
    unsigned int base = (old_value >> (SEQUENCE_BITS_PER_BASE * i)) & SEQUENCE_BIT_MASK;
    // Put it into new position
    new_value = (new_value << SEQUENCE_BITS_PER_BASE) | base;
  }

  return (value_ == new_value);
}


std::ostream& operator<<(std::ostream& os, Sequence const& seq) {
  os << seq.asString();

  return os;
}

bool Sequence::operator==(const Sequence &other) const {
  return (value_ == other.value() && n_ == other.length());
}

bool Sequence::operator!=(const Sequence &other) const {
  return !(*this == other);
}


bool Sequence::operator<(const Sequence& other) const {
  //return (n_ < other.length() || htobe64(value_) < htobe64(other.value()));
  return (n_ < other.length() || (n_ == other.length() && (value_) < (other.value())));
}

bool Sequence::operator()(const Sequence& x, const Sequence& y) const {
  return x < y;
}

std::string Sequence::asString() const {
  std::stringstream out;

  for (size_t i = 0; i < n_; i++) {
    int base_index = (value_ >> (SEQUENCE_BITS_PER_BASE * i)) & SEQUENCE_BIT_MASK;
    out << Sequence::Bases[base_index];
  }
  
  return out.str();
}

