// Author: Clara Bayley
// File: "detectors_ptr.hpp"
/* Header file for concept that defines
how shared pointer to detectors instance 
is created in a gridbox */

#ifndef DETECTORS_PTR_HPP 
#define DETECTORS_PTR_HPP 

#include <memory>
#include <concepts>

#include "./detectors.hpp"
#include "./logbooks.hpp"

template <typename F>
concept CreateDetectorsPtr = requires(F f, const unsigned int ii)
/* concept CreateDetectors is all
(function-like) types (ie. types that can be
called with some arguments) which given an
unsigned int return a shared pointer to a
Detectors instance */
{
  {
    f(ii)
  } -> std::same_as<std::shared_ptr<Detectors>>;
};

struct NullDetectorsPtr
/* operator() returns a smart pointer to
default instantiated Detectors */
{
  std::shared_ptr<Detectors> operator()(const unsigned int gbxindex) const
  {
    return std::make_shared<Detectors>();
  }
};

struct PrecipDetectorsPtr
/* operator() returns a smart pointer to a
detectors instance that may modify data in
vectors pointed to by logbooks */
{
private:
  const DetectorLogbooks &logbooks;
  const Maps4GridBoxes &gbxmaps;

  std::shared_ptr<Detectors>
  install_precipitation_detectors(const std::shared_ptr<Detectors> detectors,
                                  const unsigned int gbxindex) const;
  /* if upper z boundary of gbx is <= precip_zlim install
  a detector to detect accumulated precipitation */

  std::shared_ptr<Detectors>
  install_detectors(std::shared_ptr<Detectors> detectors,
                    const unsigned int gbxindex) const;
  /* operator installs certain types of detector in
  detectors struct given its pointer */

public:
  PrecipDetectorsPtr(const DetectorLogbooks &logbooks,
                     const Maps4GridBoxes &gbxmaps)
      : logbooks(logbooks), gbxmaps(gbxmaps) {}

  std::shared_ptr<Detectors> operator()(const unsigned int gbxindex) const
  /* operator creates a smart pointer to a detectors struct
  and installs certain types of detector in it according
  to install_detectors function */
  {
    auto detectors = std::make_shared<Detectors>();

    return install_detectors(detectors, gbxindex);
  }
};

#endif // DETECTORS_PTR_HPP 