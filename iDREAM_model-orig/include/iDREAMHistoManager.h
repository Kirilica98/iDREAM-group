#ifndef iDREAMHistoManager_H
#define iDREAMHistoManager_H 1

#include <string>

class iDREAMHistoManager
{
  public:
 //   iDREAMHistoManager(double, double, int);
    iDREAMHistoManager(unsigned int);
    ~iDREAMHistoManager();
    
    void fill(double);
    void save(std::string, std::string);
    void statistics(double& mean, double& rms);
        
  private:
    inline double bin(int);
        double min, max, h;
    int nbins;
    double* hist;
    unsigned int k;
};

#endif

