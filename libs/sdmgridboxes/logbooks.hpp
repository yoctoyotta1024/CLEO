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

  size_t get_size() const
  /* get size of record (= size of tags)
  (ie. number of entries in logbook) */
  {
    return tags.size();
  }

  unsigned int get_tag(const size_t idx) const
  /* returns tag associated with
  value at position 'idx' in record */
  {
    return tags.at(idx);
  }

  T get_entry(const size_t idx) const
  /* returns value in record at
  position 'idx' in record */
  {
    return record.at(idx);
  }

  void set_entry(const size_t idx, const T val)
  /* copies value 'val' to record at position 'idx' */
  {
    record.at(idx) = val;
  }

  void increment_entry(const size_t idx, const T val)
  /* increment value in record at
  position 'idx' by value 'val' */
  {
    record.at(idx) += val;
  }

  std::vector<T> get_and_reset_record(const T reset_val)
  /* replaces all the values in record with 'reset_val'
  and returns a copy of record before */
  {
    std::vector<T> record_copy(record);
    
    record.assign(record.size(), reset_val);

    return record_copy;
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
      const unsigned int itag)
  /* create a new entry in logbook 'ilogbook'
  with tag 'tag' and change 'idx' to
  be position of entry in that logbook */
  {
    logbook = ilogbook;
    idx = logbook->new_entry(itag);
  }

public:
  KOKKOS_INLINE_FUNCTION ~EntryInLogbook() = default; // Default destructor is Kokkos requirement for a (dual)View

  KOKKOS_INLINE_FUNCTION
  EntryInLogbook()
      : logbook(),
        idx(std::numeric_limits<size_t>::max()) {} // Default constructor is also Kokkos requirement for a (dual)View
  /* create nullptr for logbook and
  max value for position 'idx' */

  KOKKOS_INLINE_FUNCTION
  EntryInLogbook(const std::shared_ptr<Logbook<T>> ilogbook,
                 const unsigned int itag)
      : logbook(ilogbook),
        idx(logbook->new_entry(itag)) {}
  /* create an entry in logbook 'ilogbook'
  with tag 'itag' and make 'idx' to
  be position of entry in that logbook */

  auto get_logbook() const { return logbook; }

  void set_to(const T val) const
  /* copy value 'val' into logbook
  entry at position 'idx' */
  {
    logbook->set_entry(idx, val);
  }

  void increment_by(const T val) const
  /* increment value at position 'idx'
  in logbook entry by 'val' */
  {
    logbook->increment_entry(idx, val);
  }
};

struct DetectorLogbooks
/* DetectorLogbooks makes and stores
shared pointers to various logbook instances
(e.g. for a Detectors instance to use) */
{
  std::shared_ptr<Logbook<double>> surfpp; // logbook for surface precipitation

  DetectorLogbooks() : surfpp(std::make_shared<Logbook<double>>()) {}
};

#endif // LOGBOOKS_HPP