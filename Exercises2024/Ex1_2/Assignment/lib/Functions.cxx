#include<iostream>
#include <fstream>
#include <limits>
#include<vector>
#include<sstream>
#include<cmath>
#include <numeric>
#include "../share/vector_utils.h"
#include "../share/Functions.h"

float mod(float x,float y){
    float modulo= std::sqrt(x*x+y*y);
    return modulo;
}

std::vector<float> access_column(std::vector<std::vector<float>> data, int column_index){
    std::vector<float> column;
    for (const auto& row : data) {
        if (column_index < row.size()) { // Ensure the column exists in this row
            column.push_back(row[column_index]);
        }
    }
    return column;
}

std::vector<std::vector<float>> readfile(std::string filename, int n) {
    std::ifstream MyReadFile(filename); 
    if (MyReadFile.fail()) {
        std::cerr << "File not found" << std::endl;
        exit(1);
    }

    std::string myText;
    std::getline(MyReadFile, myText);

    std::vector<std::vector<float>> list_vecs;

    int line_count = 0;
    while (line_count < n && std::getline(MyReadFile, myText)) {
        std::vector<float> vec;
        std::istringstream ss(myText);
        std::string token;

        while (std::getline(ss, token, ',')) {
            try {
                float val = std::stof(token); // Convert string to float
                vec.push_back(val);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid val in line: '" << token << "'" << std::endl;
                exit(1);
            }
        }

        // Add the parsed vector to the list
        if (!vec.empty()) {
            list_vecs.push_back(vec);
        }

        line_count++;
    }

    MyReadFile.close();
    return list_vecs;
}

void printing(std::string text, float val) {
    std::cout<<text<<val<<std::endl;
}

void printing(std::vector<float> results) {
    for (float res : results){
        std::cout<<res<<std::endl;
    }
}

void printing(std::vector<std::vector<float>> coords) {
    for (std::vector<float> coord : coords){
        std::cout<<coord<<std::endl;
        //overloaded for datatype of interest - see vector_utils.h for inserrtion operator overloading procedure.
    }
}

void printing(int n, std::vector<std::vector<float>> coords){
    if (n > coords.size()){
        std::cout<<"Number of coordinates specified exceeds length of dataset.\n Setting n=25."<<std::endl;
        n=25;
    }
    std::vector<float> expos;
    std::vector<float> moduli;

    for (int i=0; i<n; i++){
        std::cout<<"\n\n\n"<<std::endl;
        std::vector<float> vec = coords[i];
        std::cout<< "Line "<< i+1<< "`s coordinates are: " <<vec <<std::endl;
        float modulus = mod(vec[0],vec[1]);
        std::cout<<"Magnitude of "<< vec << " is: " <<modulus<<"."<<std::endl;
        moduli.push_back(modulus);
        float expo = exponentiate (vec[0],vec[1]);
        std::cout<<"x^y for "<< vec << " is: " <<expo<<std::endl;
        expos.push_back(expo);
    }
    std::string file_vecs = "../output/vectors.txt";
    writing(coords,file_vecs);
    std::string file_moduli = "../output/moduli.txt";
    writing(moduli,file_moduli);
    std::string file_expos = "../output/exponents.txt";
    writing(expos,file_expos);

    //overloaded printing function to allow custom number of lines of dataset to be printed.
}


float exponentiate(float x,float y){
    if (x == 0.0) return 0.0;
    if (y == 0.0) return 1.0;
    int Y = (int) round(y); // convert to int to allow recursive multiplication
    if (Y<0) return (1/x) * exponentiate (x, Y+1);
    else return x * exponentiate (x , Y-1);
}

//choose column to sum
float summation(const std::vector<std::vector<float>>& coords, int col) {
    float total = 0.0f;
    for (const std::vector<float>& coord : coords) {
        if (col < coord.size()) {  // Ensure the column index is valid
            total += coord[col];
        } else {
            std::cerr << "Invalid column index: " << col << std::endl;
        }
    }
    return total;
}

//sum different columns
float summation_product(const std::vector<float>& vec1, const std::vector<float>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vector shape mismatch.");
    }

    float total = 0.0f;
    for (size_t i = 0; i < vec1.size(); ++i) {
        total += vec1[i] * vec2[i];  
    }
    return total;
}

//sum same columns
float summation_square(const std::vector<float>& vec) {
    float total = 0.0f;
    for (size_t i = 0; i < vec.size(); ++i) {
        total += vec[i] * vec[i];  
    }
    return total;
}

float p_calc (const float sumx, const float sumy, const float sumxsquared, const float sumxy, const int n){
    float p = (n*sumxy-sumx*sumy)/(n*sumxsquared-sumx*sumx);
    std::cout<<"\n\n\n"<<std::endl;
    std::cout<<"Gradient is:"<<p<<std::endl;
    return p;
}

float q_calc (const float sumx, const float sumy, const float sumxsquared, const float sumxy, const int n){
    float q = (sumxsquared*sumy-sumxy*sumx)/(n*sumxsquared-sumx*sumx);
    std::cout<<"Intercept is:"<<q<<std::endl;
    return q;
}

std::vector<float> y_fit(std::vector<float> x, float p, float q){
    std::vector<float> yi;
    for (float xi: x){
        yi.push_back(p*xi+q);
    }
    return yi;
}


float chisquare(std::vector<float> y_ex, std::vector<float> y_ob, std::vector<float> y_err){
    float chi2;
    for (int i=0; i<y_ex.size();i++){
        chi2 += std::pow(y_ob[i]-y_ex[i],2)/std::pow(y_err[i],2);
    }
    return chi2;
}



int query(std::string options) {
    int n;
    std::cout << options << std::endl;
    while (true) {
        std::cin >> n;

        // Check if the input fails or if the input is less than 1
        if (std::cin.fail() || n < 1) {
            std::cout << "Error: INVALID INT - Try again" << std::endl;
            std::cin.clear();

            // Ignore the invalid input left in the buffer
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        } else {
            return n;
        }
    }
}

void writing(std::vector<float> input, std::string filename){
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "File could not be opened. Please try again" << std::endl;
        return;
    }
    outFile << input << "\n";
    outFile.close();
}

void writing(std::vector<std::vector<float>> input, std::string filename){
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "File could not be opened. Please try again" << std::endl;
        return;
    }
    outFile << input << "\n";
    outFile.close();
}

void writing(float input, std::string filename){
    std::ofstream outFile(filename);
    if (!outFile) {
        std::cerr << "File could not be opened. Please try again" << std::endl;
        return;
    }
    outFile << input << "\n";
    outFile.close();
}