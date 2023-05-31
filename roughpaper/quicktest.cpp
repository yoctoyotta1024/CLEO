// /opt/homebrew/bin/g++-13 quicktest.cpp && ./a.out 

#include <iostream> 
#include <limits>

int main()
{ 
  const auto sz = std::numeric_limits<size_t>::max();
  const auto u = std::numeric_limits<unsigned int>::max();
  const auto ul = std::numeric_limits<unsigned long>::max();
  const auto ull = std::numeric_limits<unsigned long long>::max();
  
  std::cout << "szt max: " << sz << ", bytes:" << sizeof(sz) <<"\n";
  std::cout << "uuu max: " << u << ", bytes:" << sizeof(u) <<"\n";
  std::cout << "luu max: " << ul << ", bytes:" << sizeof(ul) <<"\n";
  std::cout << "llu max: " << ull << ", bytes:" << sizeof(ull) <<"\n";

  return 0;
}