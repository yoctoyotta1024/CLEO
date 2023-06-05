// Author: Clara Bayley
// File: "logbooks_intostore.hpp"
/* structures obeying the ObserveLbks
concept for various ways of observing
logbooks which ends up writing data
into a (zarr) store on disk */

#ifndef LOGBOOKS_INTOSTORE_HPP
#define LOGBOOKS_INTOSTORE_HPP

#include <vector>
#include <memory>

#include "sdmgridboxes/logbooks.hpp"
#include "zarrstorage/logbooksstorage.hpp"
#include "zarrstorage/singlevarstorage.hpp"

struct ObservePrecip
/* satisfies ObserveLbks concept and
writes precipation data to zarr storage */
{
private:
  LogbooksStorage<double> &zarr;

public:  
  ObservePrecip(LogbooksStorage<double> &zarr) : zarr(zarr)
  {
    zarr.is_name("surfpp");
    zarr.is_dim1(NOTSETVALUE(), "logbooktags");
  }

  void prepare(const DetectorLogbooks &logbooks) const
  {
    zarr.prepare(logbooks.surfpp -> get_size());
  }

  void observe_surfprecip(const std::shared_ptr<Logbook<double>> &logbook) const
  {
    std::vector<double> record = logbook -> get_and_reset_record(0.0);
    zarr.value_to_storage(record);
    ++zarr.nobs;
  }

  void operator()(const DetectorLogbooks &logbooks) const
  {
    observe_surfprecip(logbooks.surfpp); 
  }
};

#endif // LOGBOOKS_INTOSTORE_HPP