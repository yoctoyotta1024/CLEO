// /opt/homebrew/bin/g++-13 quicktest.cpp --std=c++20 && ./a.out 

#include <iostream> 
#include <limits>
#include <memory>
#include <vector>

void dtype_sizes()
{ 
  const auto sz = std::numeric_limits<size_t>::max();
  const auto u = std::numeric_limits<unsigned int>::max();
  const auto ul = std::numeric_limits<unsigned long>::max();
  const auto ull = std::numeric_limits<unsigned long long>::max();
  
  const auto dbl = std::numeric_limits<double>::max();
  
  std::cout << "szt max: " << sz << ", bytes:" << sizeof(sz) <<"\n";
  std::cout << "uuu max: " << u << ", bytes:" << sizeof(u) <<"\n";
  std::cout << "luu max: " << ul << ", bytes:" << sizeof(ul) <<"\n";
  std::cout << "llu max: " << ull << ", bytes:" << sizeof(ull) <<"\n";
  
  std::cout << "dbl max: " << ull << ", bytes:" << sizeof(dbl) <<"\n";
}

template <typename T>
struct LogbookWithGBxIndexes
{
private:
  std::vector<unsigned int> gbxindexes;
  std::vector<T> record;

public:
  size_t new_entry(const unsigned int gbxindex)
  {
    gbxindexes.push_back(gbxindex);
    record.push_back(0.0);

    return record.size() - 1;
  }

  void update_entry(const size_t idx, const T val)
  /* copies value 'val' to record at position 'idx' */
  {
    record.at(idx) = val;
  }

  size_t get_size()
  {
    return record.size(); 
  }

  unsigned int get_gbxindex(const size_t idx)
  /* returns gbxindex associated with
  value at position 'idx' in record */
  {
    return gbxindexes.at(idx);
  }

  T get_from_record(const size_t idx)
  /* returns value in record at
  position 'idx' in record */
  {
    return record.at(idx);
  }
};

template <typename T>
struct ManageEntryInLogbook
{
private:
  std::shared_ptr<LogbookWithGBxIndexes<T>> logbook;
  size_t recordidx = std::numeric_limits<size_t>::max();

public:
  void create_entry_in_logbook(
      const std::shared_ptr<LogbookWithGBxIndexes<T>> ilogbook,
      const unsigned int gbxindex)
  {
    logbook = ilogbook;
    recordidx = logbook->new_entry(gbxindex);
  }

  auto get_logbook() { return logbook; }

  void operator()(const T val)
  {
    logbook -> update_entry(recordidx, val);
  }
};

struct AccumPrecipDetector
{
  ManageEntryInLogbook<double> manage_entry;

  void operator()(const double drop)
  {
    if (manage_entry.get_logbook())
    {
      double mass_precipitation(drop);
      manage_entry(mass_precipitation);
    }
  }
};

struct DetectionLogbooks
{
  std::shared_ptr<LogbookWithGBxIndexes<double>> accprecip;

  DetectionLogbooks()
      : accprecip(std::make_shared<LogbookWithGBxIndexes<double>>()) {}
};

class Detectors
{
private:
  const DetectionLogbooks &logbooks;

public:
  AccumPrecipDetector accprecip_dtr;
  
  Detectors(const DetectionLogbooks &logbooks)
      : logbooks(logbooks) {}

  void install_precip_detector(const unsigned int gbxindex)
  {
    accprecip_dtr.manage_entry.create_entry_in_logbook(logbooks.accprecip,
                                                       gbxindex);

    std::cout << "New entry in accumprecip logbook at for gbx" << gbxindex
              << '\n';
  }
};

std::unique_ptr<Detectors> create_detectors(const DetectionLogbooks &logbooks,
                                            const unsigned int gbxindex,
                                            const double zbound)
{
  auto detectors = std::make_unique<Detectors>(logbooks);

  if (zbound < 1.0)
  {
    detectors -> install_precip_detector(gbxindex);
  }

  return detectors;
}

struct GridBox
{
  unsigned int gbxindex;
  std::unique_ptr<Detectors> detectors;

  GridBox(const unsigned int gbxindex,
          const DetectionLogbooks &logbooks,
          const double zbound)
      : gbxindex(gbxindex),
        detectors(create_detectors(logbooks, gbxindex, zbound))
  {}
};

void move_superdrops(GridBox &gbx)
{
  const double drop(2.4);
  gbx.detectors -> accprecip_dtr(drop);
}

int main()
{ 
  std::cout <<"--- let's make detectors --- \n";

  const DetectionLogbooks logbooks;

  const double zbound(0.9);

  for (unsigned int gbxindex=0; gbxindex < 3; ++gbxindex)
  {
    GridBox gbx(gbxindex, logbooks, zbound);
    move_superdrops(gbx);
  }
  
  for (size_t idx=0; idx < logbooks.accprecip -> get_size(); ++idx)
  {
    std::cout << "At idx " << idx
              << " gbx=" << logbooks.accprecip -> get_gbxindex(idx)
              << ", Record states "
              << logbooks.accprecip -> get_from_record(idx)
              << "\n";
  }

  return 0;
}