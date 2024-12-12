#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include "FiniteFunctions.h"
#include <filesystem> //To check extensions in a nice way
#include <random>
#include "gnuplot-iostream.h" //Needed to produce plots (not part of the course) 

using std::filesystem::path;

//Empty constructor
FiniteFunction::FiniteFunction(){
  m_RMin = -5.0;
  m_RMax = 5.0;
  this->checkPath("DefaultFunction");
  m_Integral = NULL;
}

//initialised constructor
FiniteFunction::FiniteFunction(double range_min, double range_max, std::string outfile){
  m_RMin = range_min;
  m_RMax = range_max;
  m_Integral = NULL;
  this->checkPath(outfile); //Use provided string to name output files
}

//Plots are called in the destructor
//SUPACPP note: They syntax of the plotting code is not part of the course
FiniteFunction::~FiniteFunction(){
  Gnuplot gp; //Set up gnuplot object
  this->generatePlot(gp); //Generate the plot and save it to a png using "outfile" for naming 
}


NormalDistribution::NormalDistribution() : FiniteFunction() {
  nd_mean = 1.0;
  nd_std_dev = 1.4;
  std::cout << "Derived default constructor\n";
}

NormalDistribution::NormalDistribution(double range_min, double range_max, std::string outfile, double mean, double std_dev)  
    : FiniteFunction(range_min, range_max, outfile), nd_mean(mean), nd_std_dev(std_dev) {
    std::cout << "Derived parameterized constructor\n";
}


CauchyLorentzDistribution::CauchyLorentzDistribution() : FiniteFunction() {
  cl_mean = 1.0;
  cl_gamma = 1.4;
  std::cout << "Derived default constructor\n";
}

CauchyLorentzDistribution::CauchyLorentzDistribution(double range_min, double range_max, std::string outfile, double mean, double gamma) 
    : FiniteFunction(range_min, range_max, outfile), cl_mean(mean), cl_gamma(gamma) {
    if (gamma <= 0) {
        throw std::invalid_argument("Gamma must be positive.");
    }
    std::cout << "Derived parameterized constructor\n";
}

NegCrystalBallDistribution::NegCrystalBallDistribution() : FiniteFunction() {
  cb_mean = -1.5;
  cb_std_dev = 1.5;
  cb_n = 2.0;
  cb_alpha = 0.85;
  std::cout << "Derived default constructor\n";
}

NegCrystalBallDistribution::NegCrystalBallDistribution(double range_min, double range_max, std::string outfile, double mean, double std_dev, double n, double alpha) 
    : FiniteFunction(range_min, range_max, outfile), cb_mean(mean), cb_std_dev(std_dev),cb_n(n),cb_alpha(alpha)  {
    if ((n <= 1)||(alpha<=0)) {
        throw std::invalid_argument("fix arguments (n>1, alpha>0) :)");
    }
    std::cout << "Derived parameterized constructor\n";
}

/*
###################
//Setters
###################
*/ 
void FiniteFunction::setRangeMin(double RMin) {m_RMin = RMin;};
void FiniteFunction::setRangeMax(double RMax) {m_RMax = RMax;};
void FiniteFunction::setOutfile(std::string Outfile) {this->checkPath(Outfile);};

/*
###################
//Getters
###################
*/ 
double FiniteFunction::rangeMin() {return m_RMin;};
double FiniteFunction::rangeMax() {return m_RMax;};
/*
###################
//Function eval
###################
*/ 
double FiniteFunction::invxsquared(double x) {return 1/(1+x*x);};
double FiniteFunction::callFunction(double x) {return this->invxsquared(x);}; //(overridable)

double NormalDistribution::normal_pmf(double x,double mean, double std_dev) {
    double exponent = -pow(x - mean, 2) / (2 * pow(std_dev, 2));
    double pmf = (1 / (std_dev * (sqrt(2 * M_PI)))) * exp(exponent);
    return pmf;
}

double NormalDistribution::callFunction(double x) {return this->NormalDistribution::normal_pmf(x,nd_mean,nd_std_dev);}; 

double CauchyLorentzDistribution::cauchy_pmf(double x, double mean, double gamma) {
    // Probability density function for Cauchy-Lorentz distribution
    double numerator = gamma;
    double denominator = M_PI * (gamma * gamma + pow(x - mean, 2));
    return numerator / denominator;
}
double CauchyLorentzDistribution::callFunction(double x) {
    return this->cauchy_pmf(x, cl_mean, cl_gamma);
}

double NegCrystalBallDistribution::crystalball_pmf(double x, double mean, double std_dev, double n, double alpha) {
    // Probability density function for crystal ball
    double exponent = -pow(fabs(alpha), 2) / 2;
    double A_coeff = pow(n/fabs(alpha),n)*exp(exponent);
    double B_coeff = n/fabs(alpha) - fabs(alpha);
    double D_coeff = sqrt(M_PI/2)*(1+std::erf(fabs(alpha)/sqrt(2)));
    double C_coeff = (n/alpha)*(1/(n-1))*exp(exponent);
    double N_coeff = 1.0/(std_dev*(C_coeff+D_coeff));
    double condition_val = (x-mean)/std_dev;

    if (condition_val>-alpha) {
      double condition_greater_exponent = -pow(x-mean,2)/(2*pow(std_dev,2));
      return N_coeff*exp(condition_greater_exponent);}
    else if (condition_val<=-alpha) {
      return N_coeff*A_coeff*pow((B_coeff-(x-mean)/std_dev),-n);}
    else throw std::invalid_argument("Alpha is nonsensical");
}

double NegCrystalBallDistribution::callFunction(double x) {
    return this->crystalball_pmf(x, cb_mean, cb_std_dev, cb_n, cb_alpha);
}

/*
###################
Integration by hand (output needed to normalise function when plotting)
###################
*/ 
double FiniteFunction::integrate(int Ndiv){ //private
  double count;
  double range = m_RMax-m_RMin;
  double step = range/(double)Ndiv;
  double x = m_RMin;
  for (unsigned int i=0; i<Ndiv; i++){
    double integral_mini;
    integral_mini=callFunction(x);
    x += range/Ndiv;
    count+=integral_mini*step;
  }
  //get bin contents over the for loop AND TIMES BY STEP SIZE :-)
  return count;
}
double FiniteFunction::integral(int Ndiv) { //public
  if (Ndiv <= 0){
    std::cout << "Invalid number of divisions for integral, setting Ndiv to 1000" <<std::endl;
    Ndiv = 1000;
  }
  if (m_Integral == NULL || Ndiv != m_IntDiv){
    m_IntDiv = Ndiv;
    m_Integral = this->integrate(Ndiv);
    return m_Integral;
  }
  else return m_Integral; //Don't bother re-calculating integral if Ndiv is the same as the last call
}




std::vector<double> FiniteFunction::metropolis(std::vector<double> sample){
    int n_random = 20000;
    int seed = 6667;
    std::mt19937 mtEngine(seed);  // Fixed seed
    std::uniform_real_distribution<double> dicePDF(m_RMin, m_RMax);
    std::uniform_real_distribution<double> distribution_TVar(0.0, 1.0);

    double x_i = dicePDF(mtEngine);  // Initial random value

    for (int i = 0; i < n_random; i++) {
        double std_dev_i = 1.5; /*values within [0.7,8] seem to work properly -
         outside this limit impede the strength of the sampling technique*/
        std::normal_distribution<double> norm(x_i, std_dev_i);
        double y = norm(mtEngine);

        double fx = this->callFunction(x_i);
        double fy = this->callFunction(y);

        double A = std::min(fy / fx, 1.0);
        double t = distribution_TVar(mtEngine);

        if (t < A) {
            x_i = y;  
        }
        sample.push_back(x_i);//fill sample
    }
    return sample;
}


/*
###################
//Helper functions 
###################
*/
// Generate paths from user defined stem
void FiniteFunction::checkPath(std::string outfile){
 path fp = outfile;
 m_FunctionName = fp.stem(); 
 m_OutData = m_FunctionName+".data";
 m_OutPng = m_FunctionName+".png";
}

//Print (overridable)
void FiniteFunction::printInfo(){
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

void NormalDistribution::printInfo() {
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "Mean: " << nd_mean << std::endl;
  std::cout << "standard deviation: " << nd_std_dev << std::endl; 
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}
void CauchyLorentzDistribution::printInfo() {
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "Mean: " << cl_mean << std::endl;
  std::cout << "Gamma: " << cl_gamma << std::endl; 
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}
void NegCrystalBallDistribution::printInfo() {
  std::cout << "rangeMin: " << m_RMin << std::endl;
  std::cout << "rangeMax: " << m_RMax << std::endl;
  std::cout << "Mean: " << cb_mean << std::endl;
  std::cout << "Standard Deviation: " << cb_std_dev << std::endl; 
  std::cout << "n: " << cb_n << std::endl;
  std::cout << "Alpha: " << cb_alpha << std::endl; 
  std::cout << "integral: " << m_Integral << ", calculated using " << m_IntDiv << " divisions" << std::endl;
  std::cout << "function: " << m_FunctionName << std::endl;
}

/*
###################
//Plotting
###################
*/

//Hack because gnuplot-io can't read in custom functions, just scan over function and connect points with a line... 
void FiniteFunction::plotFunction(){
  m_function_scan = this->scanFunction(10000);
  m_plotfunction = true;
}

//Transform data points into a format gnuplot can use (histogram) and set flag to enable drawing of data to output plot
//set isdata to true (default) to plot data points in black, set to false to plot sample points in blue
void FiniteFunction::plotData(std::vector<double> &points, int Nbins, bool isdata){
  if (isdata){
    m_data = this->makeHist(points,Nbins);
    m_plotdatapoints = true;
  }
  else{
    m_samples = this->makeHist(points,Nbins);
    m_plotsamplepoints = true;
  }
}


/*
  #######################################################################################################
  ## SUPACPP Note:
  ## The three helper functions below are needed to get the correct format for plotting with gnuplot
  ## In theory you shouldn't have to touch them
  ## However it might be helpful to read through them and understand what they are doing
  #######################################################################################################
 */

//Scan over range of function using range/Nscan steps (just a hack so we can plot the function)
std::vector< std::pair<double,double> > FiniteFunction::scanFunction(int Nscan){
  std::vector< std::pair<double,double> > function_scan;
  double step = (m_RMax - m_RMin)/(double)Nscan;
  double x = m_RMin;
  //We use the integral to normalise the function points
  if (m_Integral == NULL) {
    std::cout << "Integral not set, doing it now" << std::endl;
    this->integral(Nscan);
    std::cout << "integral: " << m_Integral << ", calculated using " << Nscan << " divisions" << std::endl;
  }
  //For each scan point push back the x and y values 
  for (int i = 0; i < Nscan; i++){
    function_scan.push_back( std::make_pair(x,this->callFunction(x)/m_Integral));
    x += step;
  }
  return function_scan;
}

//Function to make histogram out of sampled x-values - use for input data and sampling
std::vector< std::pair<double,double> > FiniteFunction::makeHist(std::vector<double> &points, int Nbins){

  std::vector< std::pair<double,double> > histdata; //Plottable output shape: (midpoint,frequency)
  std::vector<int> bins(Nbins,0); //vector of Nbins ints with default value 0 
  int norm = 0;
  for (double point : points){
    //Get bin index (starting from 0) the point falls into using point value, range, and Nbins
    int bindex = static_cast<int>(floor((point-m_RMin)/((m_RMax-m_RMin)/(double)Nbins)));
    if (bindex<0 || bindex>Nbins){
      continue;
    }
    bins[bindex]++; //weight of 1 for each data point
    norm++; //Total number of data points
  }
  double binwidth = (m_RMax-m_RMin)/(double)Nbins;
  for (int i=0; i<Nbins; i++){
    double midpoint = m_RMin + i*binwidth + binwidth/2; //Just put markers at the midpoint rather than drawing bars
    double normdata = bins[i]/((double)norm*binwidth); //Normalise with N = 1/(Ndata*binwidth)
    histdata.push_back(std::make_pair(midpoint,normdata));
  }
  return histdata;
}

//Function which handles generating the gnuplot output, called in destructor
//If an m_plot... flag is set, the we must have filled the related data vector
//SUPACPP note: They syntax of the plotting code is not part of the course
void FiniteFunction::generatePlot(Gnuplot &gp){

  if (m_plotfunction==true && m_plotdatapoints==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotdatapoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_data);
  }
  else if (m_plotfunction==true && m_plotsamplepoints==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title '"<<m_FunctionName<<"', '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_function_scan);
    gp.send1d(m_samples);
  }
  else if (m_plotfunction==true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "set style line 1 lt 1 lw 2 pi 1 ps 0\n";
    gp << "plot '-' with linespoints ls 1 title 'function'\n";
    gp.send1d(m_function_scan);
  }

  else if (m_plotdatapoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 1 lc rgb 'black' pt 7 title 'data'\n";
    gp.send1d(m_data);
  }

  else if (m_plotsamplepoints == true){
    gp << "set terminal pngcairo\n";
    gp << "set output 'Outputs/png/"<<m_FunctionName<<".png'\n"; 
    gp << "set xrange ["<<m_RMin<<":"<<m_RMax<<"]\n";
    gp << "plot '-' with points ps 2 lc rgb 'blue' title 'sampled data'\n";
    gp.send1d(m_samples);
  }
}
