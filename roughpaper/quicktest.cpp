// /opt/homebrew/bin/g++-13 quicktest.cpp && ./a.out 

#include <iostream> 
#include <functional>
#include <concepts>


class DefaultDeviceType
{
};

// template <typename P, class DT = Kokkos::DefaultExecutionSpace, typename... Args>
template <typename P, class DT, typename... Args>
concept SdmProcess = requires(P p, const int currenttimestep)
/* concept SdmProcess is all types that meet requirements
(constraints) of 2 timstepping functions called "on_step"
and "next_step" and have a "run_step" function */
{
  {
    p.next_step(currenttimestep)
    } -> std::convertible_to<int>;
  {
    p.on_step(currenttimestep)
    } -> std::convertible_to<bool>;
  {
    p.template run_step<DT>(currenttimestep)
  };
};

struct ConstTStep
{
  int interval;
  
  int next_step(const int t) const
  {
    return ((t / interval) + 1) * interval;
  }

  bool on_step(const int t) const
  {
    return t % interval == 0;
  }

  template <class DT = DefaultDeviceType>
  void run_step(const int t) {}
};

int main()
{ 
   
  SdmProcess<> auto a = ConstTStep{1};

  return 0;
}