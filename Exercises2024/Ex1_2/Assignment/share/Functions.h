#include<iostream>
#include <fstream>
#include<vector>
#include<sstream>
#include<cmath>


float mod(float x,float y);
std::vector<float> access_column(std::vector<std::vector<float>> data, int column_index);
std::vector<std::vector<float>> readfile(std::string filename, int n);
void printing(std::string text, float val);
void printing(std::vector<float> results);
void printing(std::vector<std::vector<float>> coords);
void printing(int n, std::vector<std::vector<float>> coords);
float summation(const std::vector<std::vector<float>>& coords, int col);
float summation_product(const std::vector<float>& vec1, const std::vector<float>& vec2);
float exponentiate(float x,float y);
float p_calc (const float sumx, const float sumy, const float sumxsquared, const float sumxy, const int n);
float q_calc (const float sumx, const float sumy, const float sumxsquared, const float sumxy, const int n);
std::vector<float> y_fit(std::vector<float> x, float p, float q);
float chisquare(std::vector<float> y_ex, std::vector<float> y_ob, std::vector<float> y_err);
int query(std::string options);
void writing(std::vector<float> input, std::string filename);
void writing(std::vector<std::vector<float>> input, std::string filename);
void writing(float input, std::string filename);
