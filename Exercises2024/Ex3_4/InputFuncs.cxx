#include"InputFuncs.h"
#include<iostream>
#include <fstream>
#include <limits>
#include<vector>
#include<sstream>
#include<cmath>
#include <numeric>


std::vector<double> readfile(std::string filename) {
    std::ifstream MyReadFile(filename); 
    if (MyReadFile.fail()) {
        std::cerr << "File not found" << std::endl;
        exit(1);
    }
    std::string myText;
    std::getline(MyReadFile, myText);

    std::vector<double> list;

    while (std::getline(MyReadFile, myText)) {
        double val = 0;
        val = std::stof(myText);
        if (val) {
            list.push_back(val);
        }
    }

    MyReadFile.close();
    return list;
}