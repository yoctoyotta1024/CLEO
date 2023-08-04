// /opt/homebrew/bin/g++-13 quicktest.cpp --std=c++20 && ./a.out 

#include <iostream> 
#include <limits>
#include <memory>
#include <vector>
#include <algorithm>

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
    return dropA.eps < dropB.eps; //returns true if epsA < epsB
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

  return 0;
}