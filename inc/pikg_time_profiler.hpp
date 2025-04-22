
#ifndef HPP_PIKG_TIME_PROFILER
#define HPP_PIKG_TIME_PROFILER

#include<iostream>
#include<iomanip>

#include<chrono>

#include<string>
#include<map>
#include<atomic>

#include<cassert>
#define PROFILE_PIKG

namespace PIKG{
class TimeProfiler{
#ifdef PROFILE_PIKG
  class TimerData{
    public:
      std::string name;
      std::chrono::high_resolution_clock::time_point s,e;
      double ns,max,min;
      uint32_t level;
      uint32_t ncall;
      bool isMeasuring;

      TimerData() : ns(0.0), level(0), ncall(0), isMeasuring(false) {}
      TimerData(const std::string& _name,const uint32_t _level) : name(_name), ns(0.0), max(0.0), min(1e32), level(_level), ncall(0), isMeasuring(false) {}
      TimerData(const TimerData& t) : name(t.name), s(t.s), e(t.e), ns(t.ns), level(t.level), ncall(t.ncall), isMeasuring(t.isMeasuring) {}
      ~TimerData() {}
  };
  class Index{
    int i;
    public:
    Index() : i (-1) {}
    Index(const int i) : i(i) {}
    template<class T>
    operator T() {return (T)i;}
    Index operator=(const int _i){ i = _i; return *this; }
  };
  std::map<std::string,Index> index;
  TimerData* data;
  std::atomic<int> item_count;
  uint32_t current_level;
#endif
public:
  TimeProfiler(){
#ifdef PROFILE_PIKG
    data = new TimerData[512];
    item_count = 0;
    current_level = 0;
#endif
  }
  ~TimeProfiler(){
#ifdef PROFILE_PIKG
    dump(); delete[] data;
#endif
  }
  void start(const std::string& name){
#ifdef PROFILE_PIKG
    int i = index[name]; 
    if(i == -1){
      assert(item_count < 512);
      i = item_count++;
      index[name] = i; // std::map::insert is not thread safe. do not call first start in multi threaded environment
      data[i] = TimerData(name,current_level);
    }
    if(data[i].isMeasuring){
      std::cerr << "error: TimeProfiler::start() is called again while measuring " << name << std::endl;
      exit(0);
    }
    data[i].s = std::chrono::high_resolution_clock::now();
    data[i].isMeasuring = true;
#endif
  }
  void end(const std::string& name){
#ifdef PROFILE_PIKG
    const int i = index[name]; 
    if(i == -1){
      std::cerr << "error: TimeProfiler: " << name << " is not initialized" << std::endl;
      exit(0);
    }
    if(!data[i].isMeasuring){
      std::cerr << "error: TimeProfiler::end() is called not while measuring " << name << std::endl;
    }
    data[i].e = std::chrono::high_resolution_clock::now();
    data[i].isMeasuring = false;
    double elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(data[i].e-data[i].s).count();
    data[i].ns += elapsed;
    data[i].max = std::max(data[i].max,elapsed);
    data[i].min = std::min(data[i].min,elapsed);
    data[i].ncall++;
#endif
  }
  void increaseLevel(){ 
#ifdef PROFILE_PIKG
    current_level++; 
#endif
  }
  void decreaseLevel(){ 
#ifdef PROFILE_PIKG
    assert(current_level > 0);
    current_level--; 
#endif
  }
  void clearAll(){
    for(int i=0;i<item_count;i++){
      data[i].ns = 0.0;
      data[i].max = 0.0;
      data[i].min = 1e32;
      data[i].ncall = 0;
    }
  }
  void dump(std::ostream& os = std::cout){
#ifdef PROFILE_PIKG
    std::cout << "#item time(sec) ncall max min" << std::endl;
    for(int i=0;i<item_count;i++){
      for(size_t j=0;j<data[i].level;j++) os << "    ";
      os << std::setw(40) << std::left << data[i].name << "    " << std::setw(10) << std::fixed << std::setprecision(8) << data[i].ns * 0.000000001 << "    " << data[i].ncall << "   " << data[i].max*0.000000001 << "   " << data[i].min*0.000000001 << std::endl;
    }
#endif
  }
};

} // namespace PIKG
#endif
