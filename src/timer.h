#pragma once
#include <chrono>

class timer {
  std::chrono::high_resolution_clock::time_point t1;
  std::chrono::high_resolution_clock::time_point t2;

 public:
  void start() {    t1 = std::chrono::high_resolution_clock::now(); }
  void stop() {     t2 = std::chrono::high_resolution_clock::now(); }

  double duration_in_milliseconds(int div_factor=1) const {
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count();
    return static_cast<double>((duration/1000.0) / div_factor );
  }
};
