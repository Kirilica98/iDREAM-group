#include "iDREAMHistoManager.h"
#include "iDREAMRunAction.hh"
#include <fstream>
#include <cmath>
using namespace std;

/*iDREAMHistoManager::iDREAMHistoManager(double mi, double ma, int n)
: min(mi), max(ma), nbins(n)
{
  h = (max-min)/nbins;
  hist = new int[nbins];
  for (int i=0; i<nbins; i++) hist[i] = 0;
}*/
iDREAMHistoManager::iDREAMHistoManager(unsigned int n) //0 - 4 294 967 295
: nbins(n)
{
  //h = (max-min)/nbins;
  hist = new double[nbins];
  for (int i=0; i<nbins; i++) hist[i] = 0;
  k=0;
}

iDREAMHistoManager::~iDREAMHistoManager()
{
  delete [] hist;
}

void iDREAMHistoManager::fill(double x)
{
  //if ((min<=x)&&(x<=max)) {
  //  int i = int((x - min)/h);
    hist[k] = x;
    k++;
  //}
}

void iDREAMHistoManager::save(string fname, string banner)
{
  ofstream f(fname.data());
  //f << banner << "\n";
  for (int i=0; i<nbins; i++)
   // f << bin(i) << ", " << hist[i] << "\n";
    f << hist[i] << "\n";
  f.close();
}

void iDREAMHistoManager::statistics(double& mean, double& rms)
{
  double N = 0.;
	double Sx = 0.;
	double Sx2 = 0.;
	for (int i=0; i<nbins; i++) {
	  N   += hist[i];
	  Sx  += hist[i]*bin(i);
    Sx2 += hist[i]*bin(i)*bin(i);
  }
	
	mean = Sx/N;
	rms = sqrt(Sx2/N - mean*mean);
}

double iDREAMHistoManager::bin(int i)
{
  return min + (i + .5)*h;
}

