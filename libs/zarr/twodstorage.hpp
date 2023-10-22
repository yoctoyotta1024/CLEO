/*
 * ----- CLEO -----
 * File: twodstorage.hpp
 * Project: zarr
 * Created Date: Sunday 22nd October 2023
 * Author: Clara Bayley (CB)
 * Additional Contributors:
 * -----
 * Last Modified: Sunday 22nd October 2023
 * Modified By: CB
 * -----
 * License: BSD 3-Clause "New" or "Revised" License
 * https://opensource.org/licenses/BSD-3-Clause
 * -----
 * Copyright (c) 2023 MPI-M, Clara Bayley
 * -----
 * File Description:
 */


#ifndef TWODSTORAGE_HPP 
#define TWODSTORAGE_HPP 

#include "./storehelpers.hpp"
#include "./singlevarstorage.hpp"

template <typename T>
struct TwoDStorage : SingleVarStorage<T>
/* 2D storage with dimensions [time, dim1] where
ntime is number of observation events (no. time outputs)
and ndim1 is the number of elements in 1st dimension
of 2-D data i.e. no. elements observed for each time.
For example, ndim1 could equal the number of gridboxes
an observer observes during 1 observation. Data for values
of time and dim1 could be output using a CoordinateStorage */
{
private:
  const std::string dim1name; // name of 1st dimension (e.g. "gbxindex")
  unsigned int ndim1;         // number elements in 1st dimensin (e.g. number of gridboxes that are observed)

  void writechunk()
  /* write data in buffer to a chunk in store alongside metadata jsons */
  {
    const std::string chunknum = std::to_string(this->chunkcount) + ".0";
    std::tie(this->chunkcount, this->bufferfill) = storehelpers::
        writebuffer2chunk(this->store, this->buffer,
                          this->name, chunknum,
                          this->chunkcount);

    writejsons();
  }

  void writejsons()
  /* write strictly required metadata to decode chunks (MUST).
  Assert also check 2D data dimensions is as expected */
  {
    assert((this->ndata == nobs * ndim1) &&
           "1D data length matches 2D array size");
    assert((this->get_chunksize() % ndim1 == 0.0) &&
           "chunks are integer multiple of 1st dimension of 2-D data");

    const auto n1str = std::to_string(ndim1);
    const auto nobstr = std::to_string(nobs);
    const auto nchstr = std::to_string(this->get_chunksize() / ndim1);

    const auto shape("[" + nobstr + ", " + n1str + "]");
    const auto chunks("[" + nchstr + ", " + n1str + "]");
    const std::string dims = "[\"time\", \"" + dim1name + "\"]";
    this->zarrayjsons(shape, chunks, dims);
  }

protected:
  void set_ndim1(const unsigned int i_ndim1)
  {
    if (ndim1 != storehelpers::NOTSETVALUE)
    {
      const std::string err("ndim1 already set; it cannot be changed");
      throw std::invalid_argument(err);
    }
    ndim1 = i_ndim1;
  }

public:
  unsigned int nobs; // number of output times that have been observed

  TwoDStorage(FSStore &store, const unsigned int maxchunk,
              const std::string name, const std::string dtype,
              const std::string units, const double scale_factor,
              const std::string i_dim1name, const unsigned int i_ndim1)
      : SingleVarStorage<T>(store, storehelpers::good2Dchunk(maxchunk, i_ndim1),
                            name, dtype, units, scale_factor),
        dim1name(i_dim1name), ndim1(i_ndim1), nobs(0) {}

  ~TwoDStorage()
  /* upon destruction write any data leftover in buffer
  to a chunk and write array's metadata to a .json file */
  {
    if (this->bufferfill != 0)
    {
      writechunk();
    }
  }

  void is_dim1(const size_t goodndim1,
               const std::string &goodname) const
  {
    if ((size_t)ndim1 != goodndim1)
    {
      const std::string errmsg("ndim1 is" +
                               std::to_string(ndim1) +
                               ", but should be " +
                               std::to_string(goodndim1));
      throw std::invalid_argument(errmsg);
    }

    if (dim1name != goodname)
    {
      const std::string errmsg("name of dim1 is " + dim1name +
                               ", but should be " + goodname);
      throw std::invalid_argument(errmsg);
    }
  }
};

#endif // TWODSTORAGE_HPP