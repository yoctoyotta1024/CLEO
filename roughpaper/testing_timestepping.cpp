// Author: Clara Bayley
// File: testing_timestepping.cpp
/* This file shows model of how timestepping works.
/opt/homebrew/bin/g++-12 testing_timestepping.cpp && ./a.out
*/

#include <iostream>
#include <algorithm>
#include <vector>
#include <limits>

int next_step(int t, int step)
{
  return ((t / step) + 1) * step;
}

void start_step(const int t_mdl)
{
  std::cout << t_mdl << " -> observe\n";
}

void run_driverstep(const int t_mdl)
{
  std::cout << t_mdl << " -> run driver \n";
}

int proceed_tonextstep(int t_mdl, const int couplstep)
{
  std::cout << t_mdl+couplstep << " -> couple \n";
  
  return t_mdl += couplstep;
}

int nextt_coupl_or_motion(const int t_sdm, const int couplstep,
                               const int motionstep)
/* given current time, t_sdm, work out which event (motion or coupling)
is next to occur and return the time of the sooner event */
{
  const int next_motion = ((t_sdm / motionstep) + 1) * motionstep; // t of next xchange
  const int next_coupl = ((t_sdm / couplstep) + 1) * couplstep;             // t of next output

  return std::min(next_motion, next_coupl);
}

void run_sdmstep(const int t_mdl, const int couplstep,
                 const int motionstep, const int sdmstep)
{
  int t_sdm = t_mdl;
  while (t_sdm < t_mdl + couplstep)
  {
    const int nextt = nextt_coupl_or_motion(t_sdm, couplstep, motionstep);
    
    if (t_sdm % motionstep == 0)
    {
      std::cout << t_sdm <<  " --> motion step\n";
    }

    for (int subt = t_sdm; subt < nextt;
         subt = next_step(subt, sdmstep))
    {
      if (subt % sdmstep == 0)
      {
        std::cout << subt << " ---> process step\n";
      }
    }


    t_sdm = nextt;
  }
}


int main()
{
  const int t_end = 10;
  const int couplstep = 8; // outstep
  const int motionstep = 5;
  const int sdmstep = 6;
  
  int t_mdl = 0; // time that is incremented by length 'outstep' between coupling communication
  while (t_mdl <= t_end)
  {
    start_step(t_mdl);
    run_sdmstep(t_mdl, couplstep, motionstep, sdmstep);
    run_driverstep(t_mdl);
    t_mdl = proceed_tonextstep(t_mdl, couplstep);
  }


  return 0;
}
