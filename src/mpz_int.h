#pragma once

#include <gmp.h>
#include <cassert>
#include <string>
#include <iostream>
#include <vector>

class relic_int {
 private:
  mpz_t bn;

 public:
  relic_int() { mpz_init(bn); }

  relic_int(const relic_int& b) {
    mpz_init(bn);
    mpz_set(bn, b.bn);
  }

  relic_int& operator=(const relic_int& b) {
    mpz_set(bn, b.bn);
    return (*this);
  }

  ~relic_int() { mpz_clear(bn); }

  relic_int(int64_t n) {
    uint64_t a = std::abs(n);
    bool sign = (n < 0);
    mpz_init(bn);
    mpz_set_ui(bn, a);
    if (sign) mpz_neg(bn, bn);
  }

  // Initialize from a hex String
  relic_int(const std::string& str) {
    mpz_init(bn);
    char* cstr = new char[str.length() + 1];
    strcpy(cstr, str.c_str());      
    mpz_set_str(bn, cstr, 10);
    delete[] cstr;
  }

  // Initialize from bytes
  relic_int(const uint8_t* bytes, int len) {
    mpz_init(bn);
    mpz_set_str(bn, bytes, len);
  }

  // Initialize from std::vector of bytes
  relic_int(const std::vector<uint8_t> bytes) {
    mpz_init(bn);
    mpz_set_str(bn, &bytes[0], bytes.size());
  }

  
  bool get_bit(int i) const {
    char* c_str = mpz_get_str(nullptr, 2, bn);
    return (c_str[i] == '1');
  }
  
  int bitSize() const { return mpz_sizeinbase(bn,2); }

  int to_int() const {
    int n = mpz_get_ui(bn);
    mpz_t zero;
    mpz_init(zero);
    mpz_set_ui(zero, 0);
    if (mpz_cmp(bn, zero) > 0) n *= -1;
    return n;
  }
  
  std::string ToString(int nBase = 10) const {
    char* c_str = mpz_get_str(nullptr, nBase, bn);
    std::string str(c_str);
    return str;
  }

  std::vector<uint8_t> Serialize() const {
    std::string s = ToString(10);
    std::vector<uint8_t> vec(s.begin(), s.end());
    return vec;
  }
  
  //void setHex(const std::string& str) { bn_read_str(bn, str.c_str(), str.size(), 16); }
  std::string getHex() const { return ToString(16); }
  std::string getDec() const { return ToString(10); }

  // exponentiation
  relic_int pow(const relic_int& e) const {
    relic_int res;
    uint64_t ei = mpz_get_ui(e.bn);
    mpz_pow_ui(res.bn, bn, ei);
    return res;
  }

  // multiplication: (this * b)
  relic_int mult(const relic_int& b) const {
    relic_int ret;
    mpz_mul(ret.bn, bn, b.bn);
    return ret;
  }

  // the inverse of this element mod m. ????
  relic_int inverse(const relic_int& m) const {
    relic_int ret;
    mpz_invert(ret.bn, bn, m.bn);
    return ret;
  }

  // Rabin primality test
  bool isPrime(const int checks=15) const { return mpz_probab_prime_p(bn, checks);}

  relic_int& operator+=(const relic_int& b) {
    mpz_add(bn, bn, b.bn);
    return *this;
  }

  relic_int& operator-=(const relic_int& b) {
    mpz_sub(bn, bn, b.bn);
    return *this;
  }

  relic_int& operator*=(const relic_int& b) {
    mpz_mul(bn, bn, b.bn);
    return *this;
  }

  relic_int& operator/=(const relic_int& b) {
    *this = *this / b;
    return *this;
  }

  relic_int& operator<<=(uint32_t shift) {
    mpz_mul_2exp(bn, bn, shift);
    return *this;
  }

  relic_int& operator>>=(uint32_t shift) {
    mpz_div_2exp(bn, bn, shift);
    return *this;
  }

  relic_int& operator++() {
    // prefix operator
    relic_int b(1);
    mpz_add(bn, bn, b.bn);
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
    mpz_sub(bn, bn, b.bn);
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
};

inline const relic_int operator+(const relic_int& a, const relic_int& b) {
  relic_int r;
  mpz_add(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator-(const relic_int& a, const relic_int& b) {
  relic_int r;
  mpz_sub(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator-(const relic_int& a) {
  relic_int r;
  mpz_neg(r.bn, a.bn);
  return r;
}

inline const relic_int operator*(const relic_int& a, const relic_int& b) {
  relic_int r;
  mpz_mul(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator/(const relic_int& a, const relic_int& b) {
  relic_int r;
  mpz_div(r.bn, a.bn, b.bn);
  return r;
}

inline const relic_int operator<<(const relic_int& a, uint32_t shift) {
  relic_int r;
  mpz_mul_2exp(r.bn, a.bn, shift);
  return r;
}

inline const relic_int operator>>(const relic_int& a, uint32_t shift) {
  relic_int r(a);
  r >>= shift;
  return r;
}

inline relic_int operator%(const relic_int& a, const relic_int& b) {
  relic_int c;
  mpz_mmod(c.bn, a.bn, b.bn);
  return c;
}
inline bool operator==(const relic_int& a, const relic_int& b) { return (mpz_cmp(a.bn, b.bn) == 0); }
inline bool operator!=(const relic_int& a, const relic_int& b) { return (mpz_cmp(a.bn, b.bn) != 0); }
inline bool operator<=(const relic_int& a, const relic_int& b) { return (mpz_cmp(a.bn, b.bn) <= 0); }
inline bool operator>=(const relic_int& a, const relic_int& b) {return  (mpz_cmp(a.bn, b.bn) >= 0); }
inline bool operator<(const relic_int& a, const relic_int& b) { return (mpz_cmp(a.bn, b.bn) <0); }
inline bool operator>(const relic_int& a, const relic_int& b) { return (mpz_cmp(a.bn, b.bn) >0); }
inline std::ostream& operator<<(std::ostream& strm, const relic_int& b) { return strm << b.ToString(10); }

// Functions

inline bool isPrime(const relic_int& a) { return a.isPrime(); }
inline relic_int inverse(const relic_int& a, const relic_int& b) { return a.inverse(b); }
inline relic_int pow(const relic_int& a, const relic_int& b) { return a.pow(b); }

inline std::pair<relic_int, relic_int> divmod(const relic_int& a, const relic_int& b) {
  //std::cout << "divmod in = " << a << "," << b << "\n";
  return std::make_pair(a / b, a % b);
}

inline relic_int gcd(const relic_int& a, const relic_int& b) {
  relic_int c;
  mpz_gcd(c.bn, a.bn, b.bn);
  return c;
}
inline std::tuple<relic_int, relic_int, relic_int> extended_gcd(const relic_int& a, const relic_int& b) {
  // Return r, s, t such that gcd(a, b) = r = a * s + b * t
  auto r0 = a;
  auto r1 = b;
  relic_int s0 = 1;
  relic_int s1 = 0;
  relic_int t0 = 0;
  relic_int t1 = 1;

  if (r0 > r1) {
    auto tt0 = s0;
    auto tt1 = s1;
    auto ss0 = t0;
    auto ss1 = t1;
    auto rr0 = r1;
    auto rr1 = r0;
    r0 = rr0;
    r1 = rr1;
    s0 = ss0;
    s1 = ss1;
    t0 = tt0;
    t1 = tt1;
  }
  while (r1 > 0) {
    relic_int q;
    relic_int r;
    //std::cout <<"divmod r0,r1 " << r0 << "," << r1 << "\n";

    std::tie(q, r) = divmod(r0, r1);
    //std::cout << "extended_gcd q,r = " << q << "," << r << "\n";
 
    auto tt1 = t0 - q * t1;
    auto tt0 = t1;
    auto ss1 = s0 - q * s1;
    auto ss0 = s1;
    auto rr1 = r;
    auto rr0 = r1;
    r0 = rr0;
    r1 = rr1;
    s0 = ss0;
    s1 = ss1;
    t0 = tt0;
    t1 = tt1;
  }
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
