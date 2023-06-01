// Author: Clara Bayley
// File: "logbooks.hpp"
/* Header file for functions and
structures related to logbooks
which are used by detectors to store
detection data (e.g. of SDM processes)*/

#ifndef LOGBOOKS_HPP
#define LOGBOOKS_HPP

#include <limits>
#include <memory>
#include <vector>

#include "../superdrop_solver/superdrop.hpp"

template <typename T>
struct Logbook
/* generic logbook that stores data of type 'T'
at positions in vector 'record' such that a
stored T value at position 'idx' is associated
with the tag in 'tags' at the same
positon (idx) */
{
private:
  std::vector<unsigned int> tags; // please note: a tag should never be modified once in tags
  std::vector<T> record;

public:
  size_t new_entry(const unsigned int tag)
  /* creates new entry in logbook and returns index
  'idx' of that new entry whereby value at
  record.at(idx) is value associated with
  tag at tags.at(idx). */
  {
    tags.push_back(tag);
    record.push_back(0.0);

    return record.size() - 1;
  }

  void update_entry(const size_t idx, const T val)
  /* copies value 'val' to record at position 'idx' */
  {
    record.at(idx) = val;
  }

  unsigned int get_tag(const size_t idx) const
  /* returns gbxindex associated with
  value at position 'idx' in record */
  {
    return gbxindexes.at(idx);
  }

  T get_from_record(const size_t idx) const
  /* returns value in record at
  position 'idx' in record */
  {
    return record.at(idx);
  }
};

template <typename T>
struct EntryInLogbook
/* generic struct for handling an
entry at position 'idx' in a
logbook given by the shared pointer */
{
private:
  std::shared_ptr<Logbook<T>> logbook;
  size_t idx;

  void new_entry_in_logbook(
      const std::shared_ptr<Logbook<T>> ilogbook,
      const unsigned int tag)
  /* create a new entry in logbook 'ilogbook'
  with tag 'tag' and change 'idx' to
  be position of entry in that logbook */
  {
    logbook = ilogbook;
    idx = logbook->new_entry(tag);
  }

public:
  EntryInLogbook()
      : logbook(),
        idx(std::numeric_limits<size_t>::max()) {}
  /* create nullptr for logbook and
  max value for position 'idx' */

  EntryInLogbook(const std::shared_ptr<Logbook<T>> logbook,
                 const unsigned int gbxindex)
      : logbook(ilogbook),
        idx(logbook->new_entry(gbxindex)) {}
  /* create an entry in logbook 'ilogbook'
  with tag 'gbxindex' and make 'idx' to
  be position of entry in that logbook */

  auto get_logbook() const { return logbook; }

  void operator()(const T val) const
  /* copy value 'val' into logbook
  entry at position 'idx' */
  {
    logbook->update_entry(idx, val);
  }
};

struct DetectionLogbooks
/* DetectionLogbooks makes and stores
shared pointers to various logbook
instances (for detectors to use) */
{
  std::shared_ptr<Logbook<double>> accpp; // logbook for accumulated precipitation

  DetectionLogbooks()
      : accpp(std::make_shared<Logbook<double>>()) {}
};

#endif // LOGBOOKS_HPP