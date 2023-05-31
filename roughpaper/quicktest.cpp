// /opt/homebrew/bin/g++-13 quicktest.cpp && ./a.out 

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

  void operator()(const size_t idx, const T val)
  /* returns value in record at
  position 'idx' in record */
  {
    record.at(idx) = val;
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

struct DetectionLogbooks
{
  LogbookWithGBxIndexes<double> accumprecip;    
};

class Detectors
{
private:
  std::shared_ptr<DetectionLogbooks> logbooks;

  size_t accumprecip_index = std::numeric_limits<size_t>::max();

public:
  Detectors(const std::shared_ptr<DetectionLogbooks> logbooks)
      : logbooks(logbooks) {}

  void new_precip_detector(const unsigned int gbxindex)
  {
    accumprecip_index = logbooks->accumprecip.new_entry(gbxindex);
    std::cout << "New entry in accumprecip logbook at "
              << accumprecip_index << " for gbx"
              << logbooks->accumprecip.get_gbxindex(accumprecip_index)
              << '\n';
  }

  void precip(const double testval)
  {
    if (accumprecip_index != std::numeric_limits<size_t>::max())
    {
      double precip_measured(testval);
      logbooks->accumprecip(accumprecip_index, precip_measured);
    }
  }

  void printtest(const size_t gbxindex, const double testval)
  {
    std::cout << "At gbx" << gbxindex << ", " << testval
              << " was detected. Record states "
              << logbooks->accumprecip.get_from_record(accumprecip_index)
              << "\n";
  }
};

struct GridBox
{
  unsigned int gbxindex;
  Detectors detectors;

  GridBox(const unsigned int gbxindex,
          const std::shared_ptr<DetectionLogbooks> logbooks,
          const double zbound)
      : gbxindex(gbxindex), detectors(logbooks)
  {
    if (zbound < 1.0)
    {
      detectors.new_precip_detector(gbxindex);
    }
  }
};

void move_superdrops(GridBox &gbx)
{

  const double testval(5.4);
  gbx.detectors.precip(testval);

  gbx.detectors.printtest(gbx.gbxindex, testval);
}

int main()
{ 
  std::cout <<"--- let's make detectors --- \n";

  const auto logbooks(std::make_shared<DetectionLogbooks>());

  const double zbound(0.5);

  for (unsigned int gbxindex=0; gbxindex < 3; ++gbxindex)
  {
    GridBox gbx(gbxindex, logbooks, zbound);
    move_superdrops(gbx);
  }
  
  return 0;
}