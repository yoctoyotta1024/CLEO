// /opt/homebrew/bin/g++-13 quicktest.cpp && ./a.out 

#include <iostream> 
#include <functional>

struct IncFunctional
{
    double w;
    std::function<double(const unsigned int)> get_wvel;

    IncFunctional(const unsigned int SDnspace)
      : w(3.0), get_wvel()
    {
      if (SDnspace > 0)
      {
        get_wvel = [&](const unsigned int ii) {return w;}; 
      }
      else
      {
       get_wvel =  [](const unsigned int ii) {return 0.0;};
      }
    }
};

int main()
{
  IncFunctional testit(0);
  double a = testit.get_wvel(3); 
  std::cout << a;

  IncFunctional testit2(1);
  a = testit2.get_wvel(3); 
  std::cout << a;

  return 0;
}