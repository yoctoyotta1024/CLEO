/*
 * ----- CLEO -----
 * File: initialconditions.hpp
 * Project: runcleo
 * Created Date: Tuesday 17th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 29th October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 * concept and struct obeying concept for
 * use in initial conditions of superdroplets
 * and gridboxes in CLEO SDM
 */

#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

#include <concepts>
#include <vector>
#include <utility>

template <typename IC>
concept InitialConditions = requires(IC ic, unsigned int t,
                                     const viewh_constgbx h_gbxs)
/* concept InitialConditions is all types that have
initsupers and initgbxs structs which can call functions listed */
{
  {
    ic.initsupers.get_nspacedims()
  } -> std::convertible_to<unsigned int>;
  {
    ic.initsupers.get_totnsupers()
  } -> std::convertible_to<size_t>;
  {
    ic.initsupers.get_size()
  } -> std::convertible_to<size_t>;
  {
    ic.initsupers.sdgbxindex()
  } -> std::convertible_to<std::vector<unsigned int>>;
  {
    ic.initsupers.coord3()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initsupers.coord1()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initsupers.coord2()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initsupers.radius()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initsupers.msol()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initsupers.xi()
  } -> std::convertible_to<std::vector<unsigned long long>>;

  {
    ic.initgbxs.get_ngbxs()
  } -> std::convertible_to<size_t>;
  {
    ic.initgbxs.get_size()
  } -> std::convertible_to<size_t>;
  {
    ic.initgbxs.volume()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initgbxs.press()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initgbxs.temp()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initgbxs.qvap()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initgbxs.qcond()
  } -> std::convertible_to<std::vector<double>>;
  {
    ic.initgbxs.wvel()
  } -> std::convertible_to<std::vector<std::pair<double, double>>>;
  {
    ic.initgbxs.uvel()
  } -> std::convertible_to<std::vector<std::pair<double, double>>>;
  {
    ic.initgbxs.vvel()
  } -> std::convertible_to<std::vector<std::pair<double, double>>>;
};

template <typename SuperdropInitConds, typename GbxInitConds>
struct InitConds
/* struct for functions to generate
intial conditions for CLEO */
{
  const SuperdropInitConds initsupers; // initial conditions for creating superdroplets
  const GbxInitConds initgbxs;         // initial conditions for creating gridboxes

  InitConds(const SuperdropInitConds initsupers,
            const GbxInitConds initgbxs)
      : initsupers(initsupers),
        initgbxs(initgbxs) {}
};

#endif // INITIALCONDITIONS_HPP