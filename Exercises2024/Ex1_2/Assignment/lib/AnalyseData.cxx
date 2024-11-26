#include <iostream>
#include <fstream>
#include<vector>
#include<sstream>
#include<cmath>
#include"../share/Functions.h"

int main () {
    std::string options = "Enter an option from the list below:\n 1 - Print dataset with modulus and exponents of entries\n 2 - Print data and perform a linear regression fit between x and y coordinates\n 3 - All of the above, with the addition of a reduced chi-square calculation for the fit\n 4 - Perform nothing and exit";
    int option = query(options);

    if (option==1){
        std::cout<<"Option 1 has been selected"<<std::endl;
        std::string inputs = "../data/input2D_float.txt";
        std::string errors = "../data/error2D_float.txt";
        std::string num = "Enter number of lines to be examined (i.e used in I/O):";
        int n = query(num);
        // std::cout<<n<<std::endl;
        //read in file and store vec of vec
        std::vector<std::vector<float>> list_vec = readfile(inputs, n);
        // std::cout<<"test"<<std::endl;
        std::vector<float> xvals = access_column(list_vec,0);
        std::vector<float> yvals = access_column(list_vec,1);

        //print everything
        printing(n, list_vec);
        }    
    else if(option ==2){
        std::cout<<"Option 2 has been selected"<<std::endl;
        std::string inputs = "../data/input2D_float.txt";
        std::string errors = "../data/error2D_float.txt";
        std::string num = "Enter number of lines to be examined (i.e used in I/O):";
        int n = query(num);
        // std::cout<<n<<std::endl;
        //read in file and store vec of vec
        std::vector<std::vector<float>> list_vec = readfile(inputs, n);
        std::vector<float> xvals = access_column(list_vec,0);
        std::vector<float> yvals = access_column(list_vec,1);

        //print everything
        printing(n, list_vec);
        //calc summations 
        float xsum = summation(list_vec,0);
        float ysum = summation(list_vec,1);
        float xysum = summation_product(xvals,yvals);
        float xxsum = summation_product(xvals,xvals);
        //calc p and q
        float p = p_calc (xsum, ysum, xxsum, xysum, n);
        float q = q_calc (xsum, ysum, xxsum, xysum, n);


        //fit for each val in the x_vec and check if n same for input and errors
        std::vector<float> y_ex = y_fit(xvals,p, q);
        std::vector<float> y_obs = yvals;
        
        std::string file_yex = "../output/y_ex.txt";
        writing(y_ex,file_yex);
        std::cout<<"\n\n\n"<<std::endl;
        std::cout<<"Printing expected y values from the regression model."<<std::endl;
        printing(y_ex);
    }
    // code block
    else if(option ==3){
        std::cout<<"Option 3 has been selected"<<std::endl;
        std::string inputs = "../data/input2D_float.txt";
        std::string errors = "../data/error2D_float.txt";
        std::string num = "Enter number of lines to be examined (i.e used in I/O):";
        int n = query(num);
        // std::cout<<n<<std::endl;
        //read in file and store vec of vec
        std::vector<std::vector<float>> list_vec = readfile(inputs, n);
        std::vector<float> xvals = access_column(list_vec,0);
        std::vector<float> yvals = access_column(list_vec,1);

        //print everything
        printing(n, list_vec);
        //calc summations 
        float xsum = summation(list_vec,0);
        float ysum = summation(list_vec,1);
        float xysum = summation_product(xvals,yvals);
        float xxsum = summation_product(xvals,xvals);
        //calc p and q
        float p = p_calc (xsum, ysum, xxsum, xysum, n);
        float q = q_calc (xsum, ysum, xxsum, xysum, n);

        //fit for each val in the x_vec and check if n same for input and errors
        std::vector<float> y_ex = y_fit(xvals,p, q);
        std::vector<float> y_obs = yvals;
        
        std::string file_yex = "../output/y_ex.txt";
        writing(y_ex,file_yex);
        std::cout<<"\n\n\n"<<std::endl;
        std::cout<<"Printing expected y values from the regression model."<<std::endl;
        printing(y_ex);

        std::vector<std::vector<float>> list_errors = readfile(errors, n);
        std::vector<float> yerr = access_column(list_errors, 1);
        // std::cout<<yerr.size()<<std::endl;
        float chi2 = chisquare(y_ex, y_obs, yerr);
        std::string textchi2 = "The chisquare value for this fit is: ";
        std::cout<<"\n\n\n"<<std::endl;

        printing(textchi2, chi2);
        float ndeg = static_cast<float>(y_ex.size()-2);
        float chi2red = chi2/ndeg;
        std::string textchi2red = "The reduced chisquare value for this fit is: ";
        printing(textchi2red, chi2red);

        std::string file_chi2 = "../output/reduced-chisquare.txt";
        writing(chi2red,file_chi2);
        }
    else if(option ==4){
        std::cout<<"Exiting code"<<std::endl;
        return 4;
    }

    else{
        std::cout<<"Invalid option - selection option 3"<<std::endl;
        std::string inputs = "../data/input2D_float.txt";
        std::string errors = "../data/error2D_float.txt";
        std::string num = "Enter number of lines to be examined (i.e used in I/O):";
        int n = query(num);
        // std::cout<<n<<std::endl;
        //read in file and store vec of vec
        std::vector<std::vector<float>> list_vec = readfile(inputs, n);
        std::vector<float> xvals = access_column(list_vec,0);
        std::vector<float> yvals = access_column(list_vec,1);

        //print everything
        printing(n, list_vec);
        //calc summations 
        float xsum = summation(list_vec,0);
        float ysum = summation(list_vec,1);
        float xysum = summation_product(xvals,yvals);
        float xxsum = summation_product(xvals,xvals);
        //calc p and q
        float p = p_calc (xsum, ysum, xxsum, xysum, n);
        float q = q_calc (xsum, ysum, xxsum, xysum, n);

        //fit for each val in the x_vec and check if n same for input and errors
        std::vector<float> y_ex = y_fit(xvals,p, q);
        std::vector<float> y_obs = yvals;
        
        std::string file_yex = "../output/y_ex.txt";
        writing(y_ex,file_yex);
        std::cout<<"\n\n\n"<<std::endl;
        std::cout<<"Printing expected y values from the regression model."<<std::endl;
        printing(y_ex);

        std::vector<std::vector<float>> list_errors = readfile(errors, n);
        std::vector<float> yerr = access_column(list_errors, 1);
        // std::cout<<yerr.size()<<std::endl;
        float chi2 = chisquare(y_ex, y_obs, yerr);
        std::string textchi2 = "The chisquare value for this fit is: ";
        printing(textchi2, chi2);
        float ndeg = static_cast<float>(y_ex.size()-2);
        float chi2red = chi2/ndeg;
        std::string textchi2red = "The reduced chisquare value for this fit is: ";
        printing(textchi2red, chi2red);

        std::string file_chi2 = "../output/reduced-chisquare.txt";
        writing(chi2red,file_chi2);
    }

    return 0;
}






