#pragma once

#include "relic.h"
#include "relic_conf.h"
#include <cassert>
#include <string>
#include <vector>
#include <tuple>
#include <iostream>

class relic_int {
 private:
  bn_t bn;

 public:
  relic_int() {
      bn_new(bn);
      bn_zero(bn);
  }

  relic_int(const relic_int& b) {
    bn_new(bn);
    bn_copy(bn, b.bn);
  }
  
  relic_int(const bn_t& b) {
    bn_new(bn);
    bn_copy(bn, b);
  }

  relic_int& operator=(const relic_int& b) {
      bn_zero(bn);
    bn_copy(bn, b.bn);
    return (*this);
  }
  
  relic_int& operator=(int64_t b) {
    uint64_t a = std::abs(b);
    bool sign = (b < 0);
    bn_new(bn);
    bn_set_dig(bn, a);
    if (sign) bn_neg(bn, bn);
    return (*this);
  }

  // Initialize from a decimal string
  relic_int& operator=(const std::string& str) {
    bn_new(bn);
    char* cstr = new char[str.length() + 1];
    strcpy(cstr, str.c_str());
    bn_read_str(bn, cstr, str.length(), 10);
    delete[] cstr;
    return *this;
  }

  ~relic_int() { bn_clean(bn); }

  relic_int(int64_t n) {
    uint64_t a = std::abs(n);
    bool sign = (n < 0);
    bn_new(bn);
    bn_set_dig(bn, a);
    if (sign) bn_neg(bn, bn);
  }

  // Initialize from a decimal string
  relic_int(const std::string& str) {
    bn_new(bn);
    char* cstr = new char[str.length() + 1];
    strcpy(cstr, str.c_str());      
    bn_read_str(bn, cstr, str.length(), 10);
    delete[] cstr;
  }

  // Initialize from bytes
  relic_int(const uint8_t* bytes, int len) {
    bn_new(bn);
    bn_read_bin(bn, bytes, len);
  }

  // Initialize from std::vector of bytes
  relic_int(const std::vector<uint8_t> bytes) {
    bn_new(bn);
    bn_read_bin(bn, &bytes[0], bytes.size());
  }

  
  int bitSize() const { return bn_bits(bn); }
  bool get_bit(int i) const { return bn_get_bit(bn, i);}
  void get_val(bn_t& v) const { bn_copy(v,bn); }

  int to_int() const {
    dig_t a;
    bn_get_dig(&a, bn);
    return (int)a;
  }
  
  uint64_t getUint64() const {
    uint8_t res;
    bn_write_bin(&res, 1, bn);
    return (uint64_t)res;
  }

  void SetDecimal(const std::string& str) {
    const char* psz = str.c_str();
    bn_read_str(bn, psz, str.size(), 10);
  }

  std::string ToString(int nBase = 10) const {
    char c_str[1000]; // perhaps excessive?
    bn_write_str(c_str, sizeof(c_str), bn, nBase);
    std::string str(c_str);
    return str;
  }

  std::vector<uint8_t> Serialize() const {
    std::string s = ToString(10);
    std::vector<uint8_t> vec(s.begin(), s.end());
    return vec;
  }
  
  void setHex(const std::string& str) { bn_read_str(bn, str.c_str(), str.size(), 16); }
  std::string getHex() const { return ToString(16); }
  std::string getDec() const { return ToString(10); }

  // exponentiation
  relic_int pow(const relic_int& e) const {
    relic_int res;
    bn_new(res.bn);
    bn_mxp(res.bn, bn, e.bn, e.bn);
    return res;
  }

  // multiplication: (this * b)
  relic_int mult(const relic_int& b) const {
    relic_int ret;
    bn_mul(ret.bn, bn, b.bn);
    return ret;
  }

  // Rabin primality test
  bool isPrime() const { return bn_is_prime(bn); }

  relic_int& operator+=(const relic_int& b) {
    bn_add(bn, bn, b.bn);
    return *this;
  }

  relic_int& operator-=(const relic_int& b) {
    bn_sub(bn, bn, b.bn);
    return *this;
  }

  relic_int& operator*=(const relic_int& b) {
    bn_mul(bn, bn, b.bn);
    return *this;
  }

  relic_int& operator/=(const relic_int& b) {
    *this = *this / b;
    return *this;
  }

  relic_int& operator<<=(uint32_t shift) {
    bn_lsh(bn, bn, shift);
    return *this;
  }

  relic_int& operator>>=(uint32_t shift) {
    bn_rsh(bn, bn, shift);
    return *this;
  }

  relic_int& operator++() {
    // prefix operator
    relic_int b(1);
    bn_add(bn, bn, b.bn);
    return *this;
  }

  const relic_int operator++(int) {
    // postfix operator
    const relic_int ret = *this;
    ++(*this);
    return ret;
  }

  relic_int& operator--() {
    // prefix operator
    relic_int b(1);
    bn_sub(bn, bn, b.bn);
    return *this;
  }

  const relic_int operator--(int) {
    // postfix operator
    const relic_int ret = *this;
    --(*this);
    return ret;
  }

  friend inline const relic_int operator-(const relic_int& a);
  friend inline const relic_int operator+(const relic_int& a, const relic_int& b);
  friend inline const relic_int operator-(const relic_int& a, const relic_int& b);
  friend inline const relic_int operator/(const relic_int& a, const relic_int& b);
  friend inline const relic_int operator*(const relic_int& a, const relic_int& b);
  friend inline const relic_int operator<<(const relic_int& a, uint32_t shift);
  friend inline const relic_int operator-(const relic_int& a);
  friend inline bool operator==(const relic_int& a, const relic_int& b);
  friend inline bool operator!=(const relic_int& a, const relic_int& b);
  friend inline bool operator<=(const relic_int& a, const relic_int& b);
  friend inline bool operator>=(const relic_int& a, const relic_int& b);
  friend inline bool operator<(const relic_int& a, const relic_int& b);
  friend inline bool operator>(const relic_int& a, const relic_int& b);
  friend inline relic_int operator%(const relic_int& a, const relic_int& b);
  friend inline relic_int gcd(const relic_int& a, const relic_int& b);
  friend inline relic_int sqr(const relic_int& a);
  friend inline std::tuple<relic_int, relic_int, relic_int> extended_gcd(const relic_int& a, const relic_int& b);
    

};

inline const relic_int operator+(const relic_int& a, const relic_int& b) {
  relic_int r;
  bn_add(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator-(const relic_int& a, const relic_int& b) {
  relic_int r;
  bn_sub(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator-(const relic_int& a) {
  relic_int r;
  bn_neg(r.bn, a.bn);
  return r;
}

inline const relic_int operator*(const relic_int& a, const relic_int& b) {
  relic_int r;
  bn_mul(r.bn, a.bn, b.bn);
  //  std::cout << "For * a size = " << a.bitSize() << " * " << b.bitSize() << " res = " << r.bitSize() << "\n";
  return r;
}

inline const relic_int operator/(const relic_int& a, const relic_int& b) {
  relic_int r;
  bn_div(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator<<(const relic_int& a, uint32_t shift) {
  relic_int r;
  bn_lsh(r.bn, a.bn, shift);
  return r;
}

inline const relic_int operator>>(const relic_int& a, uint32_t shift) {
  relic_int r = a;
  r >>= shift;
  return r;
}

inline relic_int operator%(const relic_int& a, const relic_int& b) {
  relic_int c;
  bn_mod(c.bn, a.bn, b.bn);
  return c;
}
inline bool operator==(const relic_int& a, const relic_int& b) { return (bn_cmp(a.bn, b.bn) == CMP_EQ); }
inline bool operator!=(const relic_int& a, const relic_int& b) { return (bn_cmp(a.bn, b.bn) != CMP_EQ); }
inline bool operator<=(const relic_int& a, const relic_int& b) {
  return ((bn_cmp(a.bn, b.bn) == CMP_LT) || (bn_cmp(a.bn, b.bn) == CMP_EQ));
}
inline bool operator>=(const relic_int& a, const relic_int& b) {
  return ((bn_cmp(a.bn, b.bn) == CMP_GT) || (bn_cmp(a.bn, b.bn) == CMP_EQ));
}
inline bool operator<(const relic_int& a, const relic_int& b) { return (bn_cmp(a.bn, b.bn) == CMP_LT); }
inline bool operator>(const relic_int& a, const relic_int& b) { return (bn_cmp(a.bn, b.bn) == CMP_GT); }
inline std::ostream& operator<<(std::ostream& strm, const relic_int& b) { return strm << b.ToString(10); }

// Functions

inline bool isPrime(const relic_int& a) { return a.isPrime(); }
inline relic_int pow(const relic_int& a, const relic_int& b) { return a.pow(b); }

inline relic_int sqr(const relic_int& a) {
  relic_int r(a);
  bn_sqr(r.bn, a.bn);
  return r;
}

inline std::pair<relic_int, relic_int> divmod(const relic_int& a, const relic_int& b) {
  //std::cout << "divmod in = " << a << "," << b << "\n";
  return std::make_pair(a / b, a % b);
}

inline relic_int gcd(const relic_int& a, const relic_int& b) {
  relic_int c;
  bn_gcd(c.bn, a.bn, b.bn);
  return c;
}
inline std::tuple<relic_int, relic_int, relic_int> extended_gcd(const relic_int& a, const relic_int& b) {
  // Return r, s, t such that gcd(a, b) = r = a * s + b * t
  relic_int r0 = a;
  relic_int r1 = b;
  relic_int s0 = 1;
  relic_int s1 = 0;
  relic_int t0 = 0;
  relic_int t1 = 1;

  if (r0 > r1) {
    r0 = b;
    r1 = a;
    s0 = 0;
    s1 = 1;
    t0 = 1;
    t1 = 0;
  }
  relic_int q,r;
  relic_int ss1,tt1;

  // 1/2 the time we get 2 loops and 1/2 about 140
  while (r1 > 0) {
    q = r0 / r1;
      /* This should be as fast as below but isn't!
      std::tie(r0,r1) = std::make_tuple(r1, r0 % r1);
      std::tie(t0,t1) = std::make_tuple(t1, t0 -q*t1);
      std::tie(s0,s1) = std::make_tuple(s1, s0 -q*s1);
*/
    r = r0 % r1;
    tt1 = t0 - q * t1;
    ss1 = s0 - q * s1;
    r0 = r1;
    s0 = s1;
    t0 = t1;
    
    r1 = r;
    s1 = ss1;
    t1 = tt1;

  }
    /* bn_gcd_ext may be broken
    std::cout << "r0 " << r0 << " s0 " << s0 << " t0 = " << t0 << "\n";
    bn_gcd_ext(r0.bn, s0.bn, t0.bn, a.bn, b.bn);
    std::cout << "r0 " << r0 << " s0 " << s0 << " t0 = " << t0 << "\n";
     */
  return std::make_tuple(r0, s0, t0);
}

inline std::pair<relic_int, relic_int> solve_mod(const relic_int& a, const relic_int& b, const relic_int& m) {
  //std::cout << "mod solve in = " << a << "," << b << "," << m << "\n";
  relic_int g, d, e;
  std::tie(g, d, e) = extended_gcd(a, m);
  // std::cout << "mod solve g,d,e = " << g <<"," << d << "," << e << "\n";
  relic_int q, r;
  std::tie(q, r) = divmod(b, g);
  //std::cout << "q,r = "<< q << "," << r << "\n";
  if (r != 0) { throw std::runtime_error("no solution..."); }
  assert(b == q * g);
  return std::make_pair((q * d) % m, m / g);
}
