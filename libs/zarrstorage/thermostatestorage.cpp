// Author: Clara Bayley
// File: thermostatestorage.cpp
/* structs/classes used to create
ObserveThermoState that writes data
from thermostate into orthogonal
multidimensional (2D) array(s) */

#include "thermostatestorage.hpp"

unsigned int ThermoIntoStore::
    copy2buffers(const ThermoState &state, unsigned int j)
/* copy press, temp, qvap and qcond data in the state to buffers at index j */
{
  storagehelper::val2buffer(state.press, pressbuffer, j);
  storagehelper::val2buffer(state.temp, tempbuffer, j);
  storagehelper::val2buffer(state.qvap, qvapbuffer, j);
  storagehelper::val2buffer(state.qcond, qcondbuffer, j);

  return ++j;
}

unsigned int ThermoIntoStore::
    writechunks(FSStore &store, unsigned int chunkcount)
/* write buffer vector into attr's store at chunkcount
and then replace contents of buffer with numeric limit */
{
  const std::string chunknum = std::to_string(chunkcount) + ".0";

  storagehelper::writebuffer2chunk(store, pressbuffer, "press",
                                   chunknum, chunkcount);
  storagehelper::writebuffer2chunk(store, tempbuffer, "temp",
                                   chunknum, chunkcount);
  storagehelper::writebuffer2chunk(store, qvapbuffer, "qvap",
                                   chunknum, chunkcount);
  storagehelper::writebuffer2chunk(store, qcondbuffer, "qcond",
                                   chunknum, chunkcount);

  return ++chunkcount;
}

void ThermoIntoStore::writejsons(FSStore &store,
                                 const std::string &metadata) const
/* write same .zarray metadata to a json file for each
thermostate array in store alongside distinct .zattrs json files */
{
  const std::string dims("[\"time\", \"gbxindex\"]");

  const auto press_a(storagehelper::arrayattrs(dims, "hPa", dlc::P0 / 100));
  storagehelper::writezarrjsons(store, "press", metadata, press_a);

  const auto temp_a(storagehelper::arrayattrs(dims, "K", dlc::TEMP0));
  storagehelper::writezarrjsons(store, "temp", metadata, temp_a);

  const auto qvap_a(storagehelper::arrayattrs(dims));
  storagehelper::writezarrjsons(store, "qvap", metadata, qvap_a);

  const auto qcond_a(storagehelper::arrayattrs(dims));
  storagehelper::writezarrjsons(store, "qcond", metadata, qcond_a);
}