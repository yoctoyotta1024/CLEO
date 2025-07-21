/*
 * Copyright (c) 2024 MPI-M, Clara Bayley
 *
 *
 * ----- CLEO -----
 * File: state.hpp
 * Project: superdrops
 * Created Date: Friday 13th October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors: Tobias Kölling (TK)
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * File Description:
 * Header file for functions and structures related to the 'State' of a certain volume, meaning
 * e.g. the thermodynamics required for SDM microphysical processes.
 */

#ifndef LIBS_SUPERDROPS_STATE_HPP_
#define LIBS_SUPERDROPS_STATE_HPP_

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>

/**
 * @brief Struct representing the state of a certain volume at a given time.
 *
 * Variables which define the state of a certain volume (e.g., inside a gridbox) at a given time.
 * State variables are, for example, thermodynamic variables required for SDM such as temperature,
 * pressure, etc.
 */
struct State {
 private:
  double volume;

 public:
  double press; /**< Pressure defined at the center of volume. */
  double temp;  /**< Temperature defined at the center of volume. */
  double qvap;  /**< Vapor mass mixing ratio at the center of volume. */
  double qcond; /**< Liquid mass mixing ratio at the center of volume. */
  Kokkos::pair<double, double> wvel;
  /**< Vertical velocity (coord3 direction) defined on {lower, upper} coord3 faces of volume. */
  Kokkos::pair<double, double> uvel;
  /**< Eastwards velocity (coord1 direction) defined on {lower, upper} coord1 faces of volume. */
  Kokkos::pair<double, double> vvel;
  /**< Northwards velocity (coord2 direction) defined on {lower, upper} coord2 faces of volume. */

  /**
   * @brief Default constructor requirement for use of State in Kokkos View
   */
  State() = default;

  /**
   * @brief Default destructor requirement for use of State in Kokkos View
   */
  ~State() = default;

  /**
   * @brief Constructor with parameters.
   *
   * @param volume The volume of the state.
   * @param press The pressure defined at the center of the volume.
   * @param temp The temperature defined at the center of the volume.
   * @param qvap The vapor mass mixing ratio defined at the center of the volume.
   * @param qcond The liquid (condensate) mass mixing ratio defined at the center of the volume.
   * @param wvel The vertical velocity (coord3 direction) defined on {lower, upper}
   * coord3 faces of the volume.
   * @param uvel The horizontal eastwards velocity (coord1 direction) defined on {lower, upper}
   * coord1 faces of the volume.
   * @param vvel The horizontal northwards velocity (coord2 direction) defined on {lower, upper}
   * coord2 faces of the volume.
   */
  KOKKOS_INLINE_FUNCTION
  State(const double volume, const double press, const double temp, const double qvap,
        const double qcond, const Kokkos::pair<double, double> wvel,
        const Kokkos::pair<double, double> uvel, const Kokkos::pair<double, double> vvel)
      : volume(volume),
        press(press),
        temp(temp),
        qvap(qvap),
        qcond(qcond),
        wvel(wvel),
        uvel(uvel),
        vvel(vvel) {}

  /**
   * @brief Get the volume of the state.
   *
   * @return The volume of the state.
   */
  KOKKOS_INLINE_FUNCTION
  auto get_volume() const { return volume; }

  /**
   * @brief Get the vertical velocity (in the coord3 direction) defined at the center of volume.
   *
   * @return The vertical velocity component 'w' defined at the center of volume.
   */
  KOKKOS_INLINE_FUNCTION
  double wvelcentre() const { return (wvel.first + wvel.second) / 2; }

  /**
   * @brief Get the horizontal eastwards velocity (in the coord1 direction) defined at
   * the center of volume.
   *
   * @return The horizontal eastwards velocity component 'u' defined at the center of volume.
   */
  KOKKOS_INLINE_FUNCTION
  double uvelcentre() const { return (uvel.first + uvel.second) / 2; }

  /**
   * @brief Get the horizontal northwards velocity (in the coord2 direction) defined at
   * the center of volume.
   *
   * @return The horizontal northwards velocity component 'v' defined at the center of volume.
   */
  KOKKOS_INLINE_FUNCTION
  double vvelcentre() const { return (vvel.first + vvel.second) / 2; }
};

#endif  // LIBS_SUPERDROPS_STATE_HPP_
