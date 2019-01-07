#pragma once
#include "relic_int.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

class Classgroup {
 private:
  relic_int a;
  relic_int b;
  relic_int c;
  relic_int d;  // for discriminant;

 public:
  Classgroup() = default;
  Classgroup(const std::tuple<relic_int, relic_int, relic_int> tup) {
    a = std::get<0>(tup);
    b = std::get<1>(tup);
    c = std::get<2>(tup);
  }
  Classgroup(const relic_int& aa, const relic_int& bb, const relic_int& cc) {
    a = aa;
    b = bb;
    c = cc;
    d = discriminant();
  }
  Classgroup(const relic_int& aa, const relic_int& bb, const relic_int& cc, const relic_int& dd) {
    a = aa;
    b = bb;
    c = cc;
    d = dd;
  }
  Classgroup(const Classgroup& x) {
    a = x.a;
    b = x.b;
    c = x.c;
    d = x.d;
  }

  Classgroup operator*=(const Classgroup b) {
    this->multiply(b);
    return *this;
  }

  int a_bits() const { return a.bitSize(); }
  int b_bits() const { return b.bitSize(); }
  int c_bits() const { return c.bitSize(); }
  int d_bits() const { return d.bitSize(); }

  void get_a(bn_t& aa) const { a.get_val(aa); }
  void get_b(bn_t& bb) const { b.get_val(bb); }
  void get_c(bn_t& cc) const { c.get_val(cc); }
  void get_d(bn_t& dd) const { d.get_val(dd); }
  
  relic_int discriminant() {
    if (d == 0) d = b * b - 4 * a * c;
    return d;
  }

  // http://www.luschny.de/math/seq/binaryqf/IPYNB-BinaryQF.html
  // reduce_real()
  // Avoid a,b, and c from getting too large
  Classgroup reduced() {
    normalized();
    //std::cout << "After 1st normalized : a,b,c = " << *this << "\n";
      
    while (a > c or (a == c and b < 0)) {
      auto s = (c + b) / (c + c);  //
      auto aa = c;
      auto bb = -b + 2 * s * c;
      c = c * s * s - b * s + a;
      a = aa;
      b = bb;
    }
    //std::cout << "Before 2nd normalized : a,b,c = " << *this << "\n";

    normalized();
    //  std::cout << "After 2nd normalized : a,b,c = " << *this << "\n";

    return *this;
  }

  std::vector<uint8_t> Serialize() const {
    std::vector<uint8_t> bytes = a.Serialize();
    std::vector<uint8_t> bs = b.Serialize();
    std::vector<uint8_t> cs = c.Serialize();
    bytes.insert(bytes.end(), bs.begin(), bs.end());
    bytes.insert(bytes.end(), cs.begin(), cs.end());
    return bytes;
  }
    
  
  void normalized() {
    if ((-a < b) && (b <= a)) return;
    auto r = (a - b) / (2 * a);  // Make sure this rounds down (floor)
    //std::cout << "Normalized r = " << r << "\n";
    auto bb = b + 2 * r * a;
    //std::cout << "Normalized a in = " << a << "\n";
    //std::cout << "Normalized b in = " << b << "\n";
    //std::cout << "Normalized 2*r*a = " << 2*r*a << "\n";
    //std::cout << "Normalized (b)bb = " << bb << "\n";
    c = a * r * r + b * r + c;
    //std::cout << "Normalized c = " << c << "\n";
    b = bb;
  }

  Classgroup identity() {
    auto d = discriminant();
    c = (b * b - d) / (4 * a);
    reduced();
    return *this;
  }

  Classgroup pow(int n) {
    if (n == 0) return identity();
    Classgroup cg(*this);
    bool items_set=false;
    Classgroup items_prod;
    int largestbit = std::log2(n);
    int mask = 1;
    for (int i = 0; i < largestbit; i++) {
      if (n & mask) {
        if (not items_set) items_prod = cg;
        else items_prod *= cg;
        items_set = true;
      }
      mask <<= 1;
      cg.square();
      // std::cout << "Sq # " << i << " a = " << cg << "\n";
    }
    if (items_set) cg *= items_prod;
    //std::cout << "cg*items " << cg << "\n";
    return cg.reduced();
  }

  Classgroup pow(const relic_int& n) {
    if (n == 0) return identity();
    Classgroup cg(*this);
    bool items_set=false;
    Classgroup items_prod;
    int largestbit = n.bitSize();
    int mask = 1;
    for (int i = 0; i < largestbit; i++) {
      if (n.get_bit(i)) {
        if (not items_set) items_prod = cg;
        else items_prod *= cg;
        items_set = true;
      }
      mask <<= 1;
      cg.square();
      // std::cout << "Sq # " << i << " a = " << cg << "\n";
    }
    if (items_set) cg *= items_prod;
    //std::cout << "cg*items " << cg << "\n";
    return cg.reduced();
  }

  Classgroup multiply(const Classgroup& other) {
    Classgroup bb = other;
    Classgroup re = bb.reduced();

    auto g = (re.b + b) >> 1;  //
    auto h = (re.b - b) >> 1;
      
    //std::cout << *this << "\n";
    //std::cout << "re = " << re << "\n";

    relic_int w = gcd(a, gcd(re.a, g));
      
    //std::cout << "g,h,w = " << g << "," << h << "," << w << "\n";

    auto j = w;
    auto r = 0;
    auto s = a / w;
    auto t = re.a / w;
    auto u = g / w;
      
    //std::cout << "j,r,s,t,u = " << j << "," << r << "," << s << "," << t << "," << u << "\n";


    // solve these equations for k, l, m
    // k * t - l * s = h
    // k * u - m * s = c2
    // l * u - m * t = c1

    // solve
    //(tu)k - (hu + sc) = 0 mod st
    // k = (- hu - sc) * (tu)^-1
      
    //std::cout << " mod inputs " << t*u << "," << h*u+s*c << "," << s*t << "\n";
      
    auto [k_temp, constant_factor] = solve_mod(t * u, h * u + s * c, s * t);
    auto [n, constant_factor_2] = solve_mod(t * constant_factor, h - t * k_temp, s);

    auto k = k_temp + constant_factor * n;
    auto l = (t * k - h) / s;
    auto m = (t * u * k - h * u - s * c) / (s * t);

     a = s * t - r * u;
    b = (j * u + m * r) - (k * t + l * s);
    c = k * l - j * m;
    return reduced();
  }

  Classgroup square() {
    reduced(); // possibly unncecessary

    auto g = b;
    auto h = 0;
      
    relic_int w = gcd(a, g);
    //std::cout << "g,h,w = " << g << "," << h << "," << w << "\n";

    auto j = w;
    auto r = 0;
    auto s = a / w;
    auto t = s;
    auto u = g / w;
      
    //std::cout << "j,r,s,t,u = " << j << "," << r << "," << s << "," << t << "," << u << "\n";

    // solve these equations for k, l, m
    // k * t - l * s = h
    // k * u - m * s = c2
    // l * u - m * t = c1

    // solve
    //(tu)k - (hu + sc) = 0 mod st
    // k = (- hu - sc) * (tu)^-1
      
    //std::cout << " solve_mod inputs " << t*u << "," << h*u+s*c << "," << s*t << "\n";
      
    auto [k_temp, constant_factor] = solve_mod(t * u, h * u + s * c, s * t);
    auto [n, constant_factor_2] = solve_mod(t * constant_factor, h - t * k_temp, s);

    //std::cout << "k_temp = " << k_temp << "\n";
    //std::cout << "constant_factor = " << constant_factor << "\n";

    //std::cout << "n = " << n << "\n";
    //std::cout << "constant_factor_2 = " << constant_factor_2 << "\n";
    
    auto k = k_temp + constant_factor * n;
    auto m = (t * u * k - h * u - s * c) / (s * t);
    auto l = (t * m + c) / u;

    //std::cout << "k = " << k << "\n";
    //std::cout << "kt = " << k*t << "\n";

    //std::cout << "****** m = " << m << "\n";
    //std::cout << "tm = " << t*m << "\n";
    //std::cout << "tmc = " << t*m + c << "\n";
    //std::cout << "u = " << u << "\n";
    //std::cout << "l = " << l << "\n";
    //std::cout << "s = " << s << "\n";

    //std::cout << "ls = " << l*s << "\n";

    
    //std::cout << "jumr = " << j*u + m*r << "\n";
    //std::cout << "ktls = " << k*t + l*s << "\n";
    

    
    a = s * t - r * u;
    b = (j * u + m * r) - (k * t + l * s);
    c = k * l - j * m;
    /*
    std::cout << "Before Reduced : a,b,c = ";
    std::cout << "a = " << a << "\n";
    std::cout << "b = " << b << "\n";
    std::cout << "c = " << c << "\n";
    */
    return reduced();
  }

  Classgroup simpler_square() {
    Classgroup re = reduced();
    auto g = b;
    auto w = gcd(a, g);
    auto j = w;
    auto s = a / w;
    auto u = g / w;
    auto s2 = sqr(s);
    auto [k_temp, constant_factor] = solve_mod(s * u, s * c, s2);
    auto [n, constant_factor_2] = solve_mod(s * constant_factor, - s * k_temp, s);
    auto k = k_temp + constant_factor * n;
    auto m = (s * u * k - s * c) / (s2);
    auto l = (s * m + c) / u;
    a = s2;
    b = (j * u) - (k * s + l * s);
    c = k * l - j * m;
    return reduced();
  }

  // Initialize from String
  //  Classgroup(const std::string& str) {    SetHexBool(str); }

  // exponentiation
  // Classgroup pow(const int e) const { return this->pow(Classgroup(e)); }

  // inverse
  Classgroup inverse() {
    b = -b;
    return *this;
  }
  /*
  Classgroup& operator+=(const Classgroup& b) { return *this; }
  Classgroup& operator/=(const Classgroup& b) { return *this; }
  Classgroup& operator%=(const Classgroup& b) { return *this; }
  Classgroup& operator<<=(const Classgroup& b) { return *this; }
  Classgroup& operator>>=(const Classgroup& b) { return *this; }
  Classgroup& operator++() { return *this; }
  Classgroup& operator--() { return *this; }
  */
  friend inline bool operator==(const Classgroup& a, const Classgroup& b);
  friend inline bool operator!=(const Classgroup& a, const Classgroup& b);
  friend inline bool operator<=(const Classgroup& a, const Classgroup& b);
  friend inline bool operator>=(const Classgroup& a, const Classgroup& b);
  friend inline bool operator<(const Classgroup& a, const Classgroup& b);
  friend inline bool operator>(const Classgroup& a, const Classgroup& b);
  friend inline Classgroup inverse(const Classgroup& m);
  friend inline std::ostream& operator<<(std::ostream& strm, const Classgroup& b);
};

inline std::ostream& operator<<(std::ostream& strm, const Classgroup& b) {
  std::stringstream s;
  //  s << "{" << b.a.ToString(16) << ",\n" << b.b.ToString(16) << ",\n" << b.c.ToString(16) << "}\n";
  s << "{" << b.a.ToString(10) << "," << b.b.ToString(10) << "," << b.c.ToString(10) << "}";
  return strm << s.str();
}

inline Classgroup from_ab_discriminant(const relic_int& aa, const relic_int& bb, const relic_int& discriminant) {
  auto c = (bb * bb - discriminant) / (4 * aa);  // check /?
  Classgroup p(aa, bb, c, discriminant);
  return p.reduced();
}

inline Classgroup identity_for_discriminant(const relic_int& d) {
  relic_int one(1);
  return from_ab_discriminant(one, one, d);
}

inline Classgroup normalized(const Classgroup& a) {
  Classgroup b = a;
  b.normalized();
  return b;
}

inline Classgroup square(const Classgroup& a) {
  Classgroup b = a;
  b.square();
  return b;

}

inline Classgroup inverse(const Classgroup& m) {
  Classgroup ret(m);
  ret.b = -m.b;
  return ret;
}

inline Classgroup pow(const Classgroup& x, int m) {
  Classgroup ret(x);
  return ret.pow(m);
}

inline Classgroup pow(const Classgroup& x, relic_int m) {
  Classgroup ret(x);
  return ret.pow(m);
}

inline const Classgroup operator+(const Classgroup& a, const Classgroup& b) {
  Classgroup r;
  return r;
}

inline const Classgroup operator-(const Classgroup& a, const Classgroup& b) {
  Classgroup r;
  return r;
}

inline const Classgroup operator-(const Classgroup& a) {
  Classgroup r;
  return r;
}

inline const Classgroup operator*(const Classgroup& a, const Classgroup& b) {
  Classgroup r(a);
  r.multiply(b);
  return r;
}

inline const Classgroup operator/(const Classgroup& a, const Classgroup& b) {
  Classgroup r(a);
  // r.divide(b);
  return r;
}
inline bool operator!=(const Classgroup& a, const Classgroup& b) {
  return (a.a != b.a) || (a.b != b.b) || (a.c != b.c); // skip d???
}

inline bool operator==(const Classgroup& a, const Classgroup& b) { return !(a != b); }

/* Do these make sense?
inline const Classgroup operator<<(const Classgroup& a, uint32_t shift) {
  Classgroup r;
  return r;
}

inline const Classgroup operator>>(const Classgroup& a, uint32_t shift) {
  Classgroup r = a;
  return r;
}
*/
/*
inline bool operator<=(const Classgroup& a, const Classgroup& b) { return (a.a <= b.a); }
inline bool operator>=(const Classgroup& a, const Classgroup& b) { return (a.a >= b.a); }
inline bool operator<(const Classgroup& a, const Classgroup& b) { return (a.a < b.a); }
inline bool operator>(const Classgroup& a, const Classgroup& b) { return (a.a > b.a); }
*/
