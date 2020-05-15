#include <iostream>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <random>
#include <vector>

#define c2 1.136847234338557L

using namespace std;

/**** Parameter Initialization ****/
// //// Here we create global Variables. The values are assigned in the function "readFile()" //// //
int seed = 0;
int samples = 0, Iterations = 0, binSize = 0; // Statistics
double dt; // Step size
double zeta = 0, ks = 0, g = 0; // g = 2*zeta*kb*T // Dynamics

/**** [FUNCTION] Read "data.txt" file [FUNCTION] ****/
void readFile(){

    /**** Reads File with Variables in use ****/
    ifstream infile("data.txt"); // Read "data.txt" file
    vector<double> arguments;
    double line = 0;

    while(infile >> line){
        arguments.push_back(line);
    }

    /**** Assign Values in "arguments" to Parameters ****/
    // //// Seed //// //
    seed = arguments[0];
    // //// Prop to Simul Time - Time, Samples, Step Size //// //
    samples = arguments[1];
    Iterations = arguments[2];
    dt = arguments[3];
    // //// Prop to File Size - How many iterations do I record //// //
    binSize = arguments[4];
    // //// Parameters of Interactions //// //
    ks = arguments[5];
    zeta = arguments[6];
    g = arguments[7];

}

int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    /**** SRK2 Vars ****/
    double pos = 0;
    double posTemp = 0;

    const int arraySizeSimul = Iterations/binSize; // gives the # of entries on the files that plot the
                                        // e2e distance in respect to time
    /**** For Loop Variables ****/
    int spl = 0, t = 0, tt = 0;    
    const double timeJump = dt*binSize;
    double tS = timeJump;

    /**** Statistics ****/
    vector<double> savePos(arraySizeSimul);

    /**** Random stuff ****/
    double al = 0;
    double std = sqrt(g/dt);
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0, std);

    /**** Current Values ****/
    double ksOzetaTdt = (ks/zeta)*dt;
    double dtOzeta = dt/zeta;

    /**** Keep track of used parameters ****/
    cout << "    *** Global Parameters *** \n \nSamples: " << samples << "\n" << "Number of Iterations: " <<
    Iterations << "\n" << "Time Step: " << dt << "\n\n" << "Bin Size: " << binSize <<
    "\n\n" << "ks = " << ks << "\n" << "zeta = " << zeta << "\n" << "g = " << g << "\n" <<  endl;

    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        tt = 0;
        pos = 0;

        /**** Calculate stuff for some time ****/
        for(t = 0; t < Iterations; t++){

            /**** Set Random Motion ****/
            // do{ 
            //     al = normal(generator)*c2;
            // }
            // while(al*al > 4*std*c2*c2);

            // //// DISTRIBUTION NOT TRUNCATED - FOR TESTING ONLY //// //
            al = normal(generator);
            //        //////////// //////////// ////////////        //

            /****** Langevin ******/
            posTemp = (1-ksOzetaTdt)*pos + dtOzeta*al;
            pos += -0.5*ksOzetaTdt*(pos + posTemp) + dtOzeta*al;

            /**** Statistics ****/
            if((t + 1) % (binSize) == 0){
                // //// Saves the position value at iteration X //// //
                savePos[tt] += pos*pos; // keeps adding to then perform an avg
                tt++; // counter to run the array above
            }
        } 
    }

    
    /**** Output file with pos of bead ****/
    // /// !! Calculo o Max no R !! /// //
    ofstream sp;
        sp.open("out1.dat");
    for(auto&& i : savePos){ 
        sp << i/samples << "\t" << tS << endl;
        tS += timeJump; // time passed
    }
    // tS = 0; // set to 0 again
    sp.close();
    //     /////// /////// ///////     //

    return 0;

}


