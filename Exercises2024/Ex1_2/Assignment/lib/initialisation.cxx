#include<iostream>
#include<vector>
#include<cmath>
#include<string>

using namespace std;

float modulo (float x, float y){
    float mod = sqrt(pow(x,2)+pow(y,2));
    return mod;
}

int main() {
    float y;
    float x;
    cout<<"Enter x:"<<endl;
    cin>>x;
    cout<<"Enter y:"<<endl;
    cin>>y;
    float eee;
    eee = modulo(x,y);
    cout<<"Magnitude of vector is:"<<eee<<"!!!"<<endl;

    return 0;
}