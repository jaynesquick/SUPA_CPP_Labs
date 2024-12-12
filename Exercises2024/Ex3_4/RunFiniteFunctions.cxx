#include "FiniteFunctions.h"
#include "InputFuncs.h"
#include<string>
#include<vector>

int main(){
    FiniteFunction standard = FiniteFunction();

    standard.printInfo();
    // standard.plotData();
    std::string fileloc = "func001";
    FiniteFunction gauss = FiniteFunction(-10.0,10.0, fileloc);
    gauss.printInfo();

    std::string filename = "/workspaces/SUPA_CPP_Labs/Exercises2024/Ex3_4/Outputs/data/MysteryData21210.txt";
    std::vector<double> place;
    place = readfile(filename);
    
    // double integral = gauss.integral(25);
    gauss.plotData(place, 25, true);
    gauss.plotFunction();
    NormalDistribution norm = NormalDistribution(-10.0,10.0,"normal",-0.65,2.0);
    norm.plotFunction();
    norm.plotData(place,25,true);
    norm.printInfo();
    std::vector<double> samples_norm;
    samples_norm = norm.metropolis(samples_norm);
    norm.plotData(samples_norm,100,false);



    CauchyLorentzDistribution cldist = CauchyLorentzDistribution(-10.0,10.0,"cauchylorentz",-0.3,2.0);
    cldist.printInfo();
    cldist.plotFunction();
    cldist.plotData(place,25,true);
    std::vector<double> samples_cauch;
    samples_cauch = cldist.metropolis(samples_cauch);
    cldist.plotData(samples_cauch,100,false);


    NegCrystalBallDistribution crystal = NegCrystalBallDistribution(-10.0,10.0,"crystalball",-1.2,1.7,8,0.9);
    crystal.printInfo();
    crystal.plotFunction();
    crystal.plotData(place,25,true);
    std::vector<double> samples_crystal;
    samples_crystal = crystal.metropolis(samples_crystal);
    crystal.plotData(samples_crystal,100,false);

    return 0;
}