#include <string>
#include <vector>
#include <cmath>
#include "gnuplot-iostream.h"

#pragma once //Replacement for IFNDEF

class FiniteFunction{

public:
  FiniteFunction(); //Empty constructor
  FiniteFunction(double range_min, double range_max, std::string outfile); //Variable constructor
  ~FiniteFunction(); //Destructor
  double rangeMin(); //Low end of the range the function is defined within
  double rangeMax(); //High end of the range the function is defined within
  double integral(int Ndiv = 1000); 
  std::vector< std::pair<double,double> > scanFunction(int Nscan = 1000); //Scan over function to plot it (slight hack needed to plot function in gnuplot)
  void setRangeMin(double RMin);
  void setRangeMax(double RMax);
  void setOutfile(std::string outfile);
  void plotFunction(); //Plot the function using scanFunction
  
  //Plot the supplied data points (either provided data or points sampled from function) as a histogram using NBins
  void plotData(std::vector<double> &points, int NBins, bool isdata=true); //NB! use isdata flag to pick between data and sampled distributions
  std::vector<double> metropolis(std::vector<double> sample);

  virtual void printInfo(); //Dump parameter info about the current function (Overridable)
  virtual double callFunction(double x); //Call the function with value x (Overridable)

//Protected members can be accessed by child classes but not users
protected:
  double m_RMin;
  double m_RMax;
  double m_Integral;
  int m_IntDiv = 0; //Number of division for performing integral
  std::string m_FunctionName;
  std::string m_OutData; //Output filename for data
  std::string m_OutPng; //Output filename for plot
  std::vector< std::pair<double,double> > m_data; //input data points to plot
  std::vector< std::pair<double,double> > m_samples; //Holder for randomly sampled data 
  std::vector< std::pair<double,double> > m_function_scan; //holder for data from scanFunction (slight hack needed to plot function in gnuplot)
  bool m_plotfunction = false; //Flag to determine whether to plot function
  bool m_plotdatapoints = false; //Flag to determine whether to plot input data
  bool m_plotsamplepoints = false; //Flag to determine whether to plot sampled data 
  double integrate(int Ndiv);
  std::vector< std::pair<double, double> > makeHist(std::vector<double> &points, int Nbins); //Helper function to turn data points into histogram with Nbins
  void checkPath(std::string outstring); //Helper function to ensure data and png paths are correct
  void generatePlot(Gnuplot &gp); 
  
private:
  double invxsquared(double x); //The default functional form
};

class NormalDistribution : public FiniteFunction{
  public:
  NormalDistribution();
  NormalDistribution(double range_min, double range_max, std::string outfile, double mean, double std_dev); //Variable constructor

  double callFunction(double x);
  void printInfo();
  protected:
  
  double nd_mean;
  double nd_std_dev;

  private:
  double normal_pmf(double x,double mean, double std_dev);
};


class CauchyLorentzDistribution : public FiniteFunction{
  public:
  CauchyLorentzDistribution();
  CauchyLorentzDistribution(double range_min, double range_max, std::string outfile, double mean, double gamma); //Variable constructor

  double callFunction(double x);
  void printInfo();
  protected:
  
  double cl_mean;
  double cl_gamma;

  private:
  double cauchy_pmf(double x,double mean, double gamma);
};

class NegCrystalBallDistribution : public FiniteFunction{
  public:
  NegCrystalBallDistribution();
  NegCrystalBallDistribution(double range_min, double range_max, std::string outfile, double mean, double std_dev, double n, double alpha); 

  double callFunction(double x);
  void printInfo();
  protected:
  
  double cb_mean;
  double cb_std_dev;
  double cb_n;
  double cb_alpha;

  private:
  double crystalball_pmf(double x, double mean, double std_dev, double n, double alpha);
};


