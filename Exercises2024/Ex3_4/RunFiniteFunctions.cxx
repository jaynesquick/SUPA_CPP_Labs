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
    NormalDistribution norm = NormalDistribution(-10.0,10.0,"normal",5,1.0);
    norm.plotFunction();
    norm.plotData(place,25,true);
    norm.printInfo();

    CauchyLorentzDistribution cldist = CauchyLorentzDistribution(-10.0,10.0,"cauchylorentz",3.6,0.9);
    cldist.printInfo();
    cldist.plotFunction();
    cldist.plotData(place,25,true);


    NegCrystalBallDistribution crusty = NegCrystalBallDistribution(-10.0,10.0,"crystalball",-1.2,1.25,2.4,0.75);
    crusty.printInfo();
    crusty.plotFunction();
    crusty.plotData(place,25,true);

    return 0;
}