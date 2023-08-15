// /opt/homebrew/bin/g++-13 quicktest.cpp --std=c++20 && ./a.out 
// g++ quicktest.cpp --std=c++20 && ./a.out 

#include <iostream> 
#include <limits>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

struct Superdrop
{
  unsigned long long eps;
  double radius;
};

Superdrop &assign_superdroplet_old(Superdrop &dropA, Superdrop &dropB,
                                  const unsigned int whichdrop)
/* compare dropA.eps with dropB.eps and return either
drop1 or drop2 such that drop1.eps is always > drop2.eps */
{
  if (dropA.eps > dropB.eps)
  {
    if (whichdrop == 1) // "drop1"
    {
      return dropA;
    }
    else // "drop2"
    {
      return dropB;
    }
  }

  else
  {
    if (whichdrop == 1) // "drop1"
    {
      return dropB;
    }
    else // "drop2"
    {
      return dropA;
    }
  }
}

std::pair<Superdrop&, Superdrop&>
assign_superdroplet(const Superdrop &dropA, const Superdrop &dropB)
/* compare dropA.eps with dropB.eps and return (non-const)
references to dropA and dropB in a pair {drop1, drop2}
such that drop1.eps is always > drop2.eps */
{
  auto comp = [](const Superdrop &dropA, const Superdrop &dropB)
  {
    return dropA.eps < dropB.eps; // returns true if epsA < epsB
  };

  auto [drop2, drop1] = std::minmax(dropA, dropB, comp); // drop2.eps =< drop1.eps

  return {const_cast<Superdrop&>(drop1), const_cast<Superdrop&>(drop2)};
}

int main()
{ 
  Superdrop dropA{30000, 1.0}; 
  Superdrop dropB{10000, 1.0}; 

  Superdrop &drop1o = assign_superdroplet_old(dropA, dropB, 1);
  Superdrop &drop2o = assign_superdroplet_old(dropA, dropB, 2);

  std::cout << dropA.eps << ", " << dropB.eps << "\n";
  std::cout << drop1o.eps << ", " << drop2o.eps << "\n";

  auto [drop1, drop2] = assign_superdroplet(dropA, dropB);

  std::cout << dropA.eps << ", " << dropB.eps << "\n";
  std::cout << drop1.eps << ", " << drop2.eps << "\n";

  drop1.eps = 4;
  drop2.eps = 2;

  std::cout << dropA.eps << ", " << dropB.eps << "\n";
  std::cout << drop1.eps << ", " << drop2.eps << "\n";

  std::cout << "\n---------------------------------\n";

  const double d1(std::pow(1.4015e-05*2.0, 3.0));
  const double d2(std::pow(6.2649e-05*2.0, 3.0));

  const double dratio(d1 * d2 / (d1 + d2));
  std::cout << 998.0 * std::numbers::pi / 12.0 * dratio * (-7.0817e-01) * (-7.0817e-01) << "\n"; //cke
  std::cout << 7.28e-2 * std::numbers::pi * (1.4015e-05*2.0) * (1.4015e-05*2.0) << "\n"; //surft of min radius


  std::cout << "\n----------------------------------\n";
 
  unsigned long long eps(1);
  std::cout << "eps-1 = " << eps - 1 << "\n";
  const auto lim(std::numeric_limits<unsigned long long>::max());
  unsigned long long val(std::min(eps-1, lim-1));
  std::cout << "val = " << eps - 1 << "\n";
  std::cout << "new eps = " << val + 1 << "\n";

  std::cout << "\n-----------favsdca-----------------------\n";
  unsigned long long eps00(1);
  std::cout << "eps-1 = " << eps00/2 << ", " <<  eps00  - eps00/2 << "\n";

  return 0;
}