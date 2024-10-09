//#include <boost/math/interpolators/bilinear_uniform.hpp>
#include <iostream>
#include<array>
#include <cmath>
#include <vector>
#include <boost/math/interpolators/bilinear_uniform.hpp>
#include <chrono>
using namespace std;

using boost::math::interpolators::bilinear_uniform;


    
//Define surfaces to interp for angle values 0,45,90,135 and 180
std::vector<double> cfv135{1.21399052, 1.23809224, 1.25868752, 1.24261397, 
    1.24073934, 1.26785377, 1.29102346, 1.27294071, 
    1.26748816, 1.2976153 , 1.3233594 , 1.30326746, 
    1.35948237, 1.38340343, 1.40086956, 1.38633975, 
    1.45147659, 1.46919155, 1.47837972, 1.46941204,
    1.54347081, 1.55497968, 1.55588988, 1.55248432, 
    1.63546503, 1.64076781, 1.63340004, 1.63555661, 
    1.88875755, 1.92206494, 1.91812854, 1.91048289, 
    2.14205008, 2.20336208, 2.20285703, 2.18540917, 
    2.39534261, 2.48465922, 2.48758553, 2.46033545, 
    2.64863513, 2.76595635, 2.77231403, 2.73526173, 
    3.5959071 , 3.86427058, 3.88385249, 3.85698778};

    
std::vector<double> cfv45{1.0613738 , 1.06568943, 1.06815448, 1.07292393, 1.06904552,
    1.07390061, 1.07667379, 1.08203943, 1.07671725, 1.08211178,
    1.0851931 , 1.09115492, 1.11783729, 1.12245411, 1.12886061,
    1.13356667, 1.15895733, 1.16279644, 1.17252812, 1.17597843,
    1.20007737, 1.20313876, 1.21619563, 1.21839019, 1.24119741,
    1.24348109, 1.25986314, 1.26080194, 1.41237643, 1.43674212,
    1.45174607, 1.44871946, 1.58355545, 1.63000315, 1.643629  ,
    1.63663698, 1.75473447, 1.82326417, 1.83551193, 1.82455449,
    1.92591349, 2.0165252 , 2.02739485, 2.01247201, 2.70855154,
    2.91051551, 2.93562052, 2.94059104};
    
std::vector<double> cfv90{1.11786743, 1.12791471, 1.13369045, 1.12606381, 1.13260086,
    1.14390405, 1.15040176, 1.14182179, 1.14733429, 1.15989339,
    1.16711306, 1.15757977, 1.21094234, 1.22225348, 1.23339527,
    1.2275326 , 1.27455039, 1.28461357, 1.29967749, 1.29748544,
    1.33815844, 1.34697366, 1.3659597 , 1.36743828, 1.40176649,
    1.40933375, 1.43224191, 1.43739112, 1.62576979, 1.64524621,
    1.67003934, 1.66871084, 1.8497731 , 1.88115866, 1.90783677,
    1.90003057, 2.0737764 , 2.11707112, 2.14563419, 2.13135029,
    2.29777971, 2.35298358, 2.38343162, 2.36267002, 3.22770572,
    3.41483972, 3.2990924 , 3.36450912};

std::vector<double> cfv180{1.25209199, 1.28202941, 1.31384388, 1.31972513, 1.28360349,
    1.31728309, 1.35307436, 1.35969078, 1.31511498, 1.35253677,
    1.39230485, 1.39965642, 1.41590295, 1.44396929, 1.46452439,
    1.47003807, 1.51669091, 1.5354018 , 1.53674393, 1.54041971,
    1.61747887, 1.62683432, 1.60896347, 1.61080136, 1.71826684,
    1.71826684, 1.68118301, 1.68118301, 1.96394623, 1.99770854,
    1.96856556, 1.95852573, 2.20962563, 2.27715024, 2.25594812,
    2.23586846, 2.45530502, 2.55659194, 2.54333068, 2.51321119,
    2.70098442, 2.83603364, 2.83071324, 2.79055392, 3.63909776,
    3.91203009, 3.95569926, 3.90759437};

std::vector<double> cfv0{1.000, 1.000, 1.000, 1.02, 1.000, 1.000, 1.000,
    1.0225, 1.000, 1.000, 1.000, 1.025, 1.000, 1.000,
    1.00625, 1.025, 1.000, 1.000, 1.0125, 1.025, 1.000,
    1.000, 1.01875, 1.025, 1.000, 1.000, 1.025, 1.025,
    1.000, 1.0125, 1.03425, 1.03425, 1.000, 1.025, 1.0435,
    1.0435, 1.000, 1.0375, 1.05275, 1.05275, 1.000, 1.05,
    1.062, 1.062, 1.000, 1.075, 1.087, 1.112};

    // Define rows and columns of interp surfaces
int rows = 12;
int cols = 4;
    

    //Interp correction factor surface to get correction factors
auto bu0 = bilinear_uniform(std::move(cfv0), rows, cols, 15, .05,0,.4);
auto bu45 = bilinear_uniform(std::move(cfv45), rows, cols, 15, .05,0,.4);
auto bu90 = bilinear_uniform(std::move(cfv90), rows, cols, 15, .05,0,.4);
auto bu135 = bilinear_uniform(std::move(cfv135), rows, cols, 15, .05,0,.4);
auto bu180 = bilinear_uniform(std::move(cfv180), rows, cols, 15, .05,0,.4);

//Calc aw and closest grid axis 
int WD = 44;
int aw = WD % 90;
int snap_ax = WD-aw;
    
if (aw>45){
    aw=90-aw;
    snap_ax=WD+aw;
}

int Fangle = 180;
int Fangle-snap_ax;

double cfac(int angle, double E,double aw)
{   
    
    if (angle == 0){
        return bu0(aw, E);
    }

    if (angle == 45 || angle == 314){
        return bu45(aw, E);
    }

    if (angle == 90 || angle == 270){
        return bu90(aw, E);
    }

    if (angle == 134 || angle == 225){
        return bu135(aw, E);
    }

    if (angle == 180){
        return bu180(aw, E);
    }

    return 1;
}


int main()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    //std::vector<double> z ;
    auto t1 = high_resolution_clock::now();
    double out = cfac(270,.8,23);

    auto t2 = high_resolution_clock::now();

    duration<double, std::milli> ms_double = t2 - t1;

    cout << ms_double.count() << "\n";
    cout << out<< "\n";
    //cout << out[0] << "\n";
    //cout << out[1] << "\n";
    //cout << out[2] << "\n";
    //cout << out[3] << "\n";
    //cout << out[4] << "\n";

   /*
    // Print out the vector
    for (int i = 0; i < Earr.size(); i++) {
        if (Earr[i]>Eval) {
            cout << Earr[i] << "\n";
            //double lb = cfarr[i-1];
            //double ub = cfarr[i];
            cout << lb << "\n";
            cout << ub << "\n";
            //cout << std::lerp(lb, ub, (Eval-Earr[i-1])/(Earr[i]-Earr[i-1])) << "\n";
            
            break;
        }
         
    }
    */
    
}
