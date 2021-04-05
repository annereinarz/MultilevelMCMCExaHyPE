#ifndef __InitialData_CLASS_HEADER__
#define __InitialData_CLASS_HEADER__

#if defined(USE_ASAGI)
namespace easi {
  class YAMLParser;
  class Component;
};
class AsagiReader;
#endif

class InitialData {
 public:
  InitialData();
  InitialData(int scenario, char* filename);  
  ~InitialData();
  void getInitialData(const double* const x,double* Q);
  
 private:
  int scenario;
#if defined(USE_ASAGI)
  easi::YAMLParser* parser;
  easi::Component* model;
  AsagiReader* asagiReader;
#endif

  void readAsagiData(const double* const x,double* Q);
  void readAsagiData_nobath(const double* const x,double* Q);
  void readAsagiData_onlybath(const double* const x,double* Q);
};

#endif // __InitialData_CLASS_HEADER__
