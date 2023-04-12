// Author: Clara Bayley
// File: "sdmprocess.hpp"
/* SDM Process Concept as well as structures for
various different SDM processes that may occur
in the superdroplet model, eg. condensation or
collision-coalescence (see ConstTstepProcess struct) */

#ifndef SDMPROCESS
#define SDMPROCESS

#include <concepts>
#include <span>
#include <random>
#include <limits>

#include "./superdrop.hpp"
#include "./thermostate.hpp"
#include "./collisionsmethod.hpp"
#include "./condensationmethod.hpp"

template <typename F>
concept StepFunc = requires(F f, int currenttimestep,
                            std::span<SuperdropWithGbxindex> span4SDsinGBx,
                            ThermoState state,
                            std::mt19937 gen)
/* concept StepFunc is all (function-like) types
(ie. types that can be called with some arguments)
that have the same signature as the "run_step"
function (see below in SdmProcess) */
{
  {f(currenttimestep, span4SDsinGBx, state, gen)};
};

template <typename P, typename... Args>
concept SdmProcess = requires(P p, int currenttimestep,
                              std::span<SuperdropWithGbxindex> span4SDsinGBx,
                              ThermoState state,
                              std::mt19937 gen)
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
  {p.run_step(currenttimestep, span4SDsinGBx, state, gen)};
};

template <SdmProcess A, SdmProcess B>
struct CombinedSdmProcess
{
  A a;
  B b;

  CombinedSdmProcess(A a, B b) : a(a), b(b) {}

  int next_step(const int currenttimestep) const
  /* for combination of 2 SDM proceses, the next step
  is the smaller out of the two possible next steps */
  {
    int s1 = a.next_step(currenttimestep);
    int s2 = b.next_step(currenttimestep);
    if (s1 < s2)
    {
      return s1;
    }
    else
    {
      return s2;
    }
  }

  bool on_step(const int currenttimestep) const
  /* for combination of 2 SDM proceses, a tstep
  is on_step true when either process is on_step true */
  {
    return a.on_step(currenttimestep) || b.on_step(currenttimestep);
  }

  void run_step(const int currenttimestep,
                std::span<SuperdropWithGbxindex> span4SDsinGBx,
                ThermoState &state,
                std::mt19937 &gen) const
  /* for combination of 2 SDM proceses, each process is
  run if it's respective on_step returns true */
  {
    if (a.on_step(currenttimestep))
    {
      a.run_step(currenttimestep, span4SDsinGBx, state, gen);
    }

    if (b.on_step(currenttimestep))
    {
      b.run_step(currenttimestep, span4SDsinGBx, state, gen);
    }
  }
};

auto operator>>(SdmProcess auto a, SdmProcess auto b)
/* define ">>" operator that combines
two Superdroplet Model Processes */
{
  return CombinedSdmProcess(a, b);
}

struct NullProcess
/* NullProcess does nothing at all
(is defined for a Monoid Structure) */
{
  int next_step(int currenttimestep) const
  {
    return std::numeric_limits<int>::max();
  }

  bool on_step(int currenttimestep) const
  {
    return false;
  }

  void run_step(const int currenttimestep,
                std::span<SuperdropWithGbxindex> span4SDsinGBx,
                ThermoState &state,
                std::mt19937 &gen) const {}
};

template <StepFunc F>
struct ConstTstepProcess
/* this structure is a type that satisfies the concept of an SDM
Process and has a constant tstep 'interval'. It can be used to create
SDM processes from a constant timestep and a method whose type
satisfies the StepFunc concept */
{
  int interval;
  F run_step;

  int next_step(int t) const
  {
    return ((t / interval) + 1) * interval;
  }

  bool on_step(int t) const
  {
    return t % interval == 0;
  }
};

#endif // SDMPROCESS