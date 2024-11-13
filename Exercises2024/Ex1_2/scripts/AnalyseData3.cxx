#include <iostream>
#include <fstream>
#include<vector>
#include<sstream>
#include<cmath>
using namespace std;


float modulo (float x, float y){
    float mod = sqrt(pow(x,2)+pow(y,2));
    return mod;
}



int main () {
    string myText;

    // Read from the text file
    ifstream MyReadFile("../input2D_float.txt");
    int n;
    cout<<"Enter n:"<<endl;
    cin>>n;
    vector<vector<float>> list_vecs;
    vector<float> mags;





    // cout<<"2"<<endl;
    for(int i=0; i<=n;i++){
        
        
        getline (MyReadFile, myText);
        if(myText=="x,y") continue;

        vector<float> vec;
        istringstream ss(myText);
        string token;
        float val;
        while(std::getline(ss, token, ',')) {


            val = stof(token);
            vec.push_back(val);
        }
        cout <<"x: " <<vec[0]<<",   y: "<<vec[1]<<endl;
        float mag;
        mag = modulo(vec[0],vec[1]);
        cout << "Magnitude is: " << mag << endl;
        list_vecs.push_back(vec);
        mags.push_back(mag);
        

        
    }
    // Close the file
    MyReadFile.close();
    return 0;
}