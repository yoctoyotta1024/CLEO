// /opt/homebrew/bin/g++-13 quicktest.cpp && ./a.out 

#include <iostream> 
#include <cmath>

int main()
{ 
  const double delt = 1.0;
  const double subdelt = 1.0;
   
  const unsigned int nsubs = std::ceil(delt/subdelt);
  std::cout << "subs: " << nsubs << "\n";

  return 0;
}