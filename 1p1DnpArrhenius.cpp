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
int samples = 0, Iterations = 0; // Statistics
double dt = 0; // Step size
double zeta = 0, g = 0; // g = 2*zeta*kb*T // Dynamics

// //// Potential Parameters //// //
// Const //
const int sigma = 1;
const int epsilon = 1;
const long double rMin = 1.12246204830937298*sigma;
const long double rCut = 2.5*rMin;
const double rc2 = rCut*rCut;
const double s2 = sigma*sigma;
const double s6 = s2*s2*s2;
const double s12 = s6*s6;
// On data.txt file//
double peak = 0, a = 0, b = 0;


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
    // //// Parameters of Interactions //// //
    zeta = arguments[4];
    g = arguments[5];
    // //// Parameters of the Potential //// //
    peak = arguments[6];
    a = arguments[7];
    b = arguments[8];

}


int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    const double totalTime = Iterations*dt;

    double elapsedTime = 0;

    /**** Random stuff ****/
    double al = 0;
    double std = sqrt(g/dt);
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0, std);

    /**** Positions ****/
    double rx = rMin, rxTemp = 0, dist = 0, d6 = 0, d12 = 0;
    double force = 0, forceTemp = 0;

    /**** For Loop Variables ****/
    int spl = 0, breakCount = 0, p = 0;
    double t = 0;

    /**** Initialization of Repeated Calculations ****/
    double distMpeak = 0;
    const double b2 = b*b;
    const double dtOzeta = dt/zeta;

    /**** Statistics ****/
    double avgElapsTime = 0;

    cout << "\t*** Global Parameters *** \n \nSamples: " << samples << "\n"  << "Number of Iterations: " << Iterations << "\n" 
    << "Time Step: " << dt << "\n\n" << "\t*** Potential Parameters *** \n\nsigma = " << sigma << "\n" << "epsilon = " << epsilon << "\n" 
    << "a = " << a << "\n" << "b = " << b << "\n" << "Peak of Potential Wall at " << peak << "\n\n\t*** Other Parameters *** \n\n" << "zeta = " << zeta << "\n\n" << 
    "\t*** Parameter in Study *** \n\ng = " << g << endl;

    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        /**** Initialize Variables ****/
        t = 0;
        rx = rMin;

        /**** Calculate stuff for some time ****/
        while(t < totalTime){

            // //////// ///////////// //////// //
            /**** Force Evaluation - Begin ****/
            dist = rx; // Distance Calculation

            if(dist < peak){

                distMpeak = dist - peak;
                d6 = dist*dist*dist*dist*dist*dist;
                d12 = d6*d6;

                force = (48*epsilon/dist)*(s12/d12 - 0.5*s6/d6) + ((2*a*distMpeak)/b2)*exp(-(distMpeak*distMpeak)/b2);
            }
            else{ // if "dist" is greater than "peak" the particle escaped the potential well

                elapsedTime += t - dt; // save the time it took to happen and add it to previous samples
                breakCount++; // count the number of times I added so I can avg the result
                break; // break while loop

            }
            /**** Force Evaluation - End ****/
            // /////// //////////// /////// //

            /**** Set Random Motion ****/
            // do{ 
            //     al = normal(generator)*c2;
            // }
            // while(al*al > 4*std*c2*c2);

            // //// DISTRIBUTION NOT TRUNCATED - FOR TESTING ONLY //// //
            al = normal(generator);
            //        //////////// //////////// ////////////        //

            /**** Langevin - Temporary ****/
            rxTemp = rx + dtOzeta*(al + force);

            // //////// ///////////// //////// //
            /**** Force Evaluation - Begin ****/
            dist = rxTemp; // Distance Calculation

            if(dist < peak){

                distMpeak = dist - peak;
                d6 = dist*dist*dist*dist*dist*dist;
                d12 = d6*d6;

                forceTemp = (48*epsilon/dist)*(s12/d12 - 0.5*s6/d6) + ((2*a*distMpeak)/b2)*exp(-(distMpeak*distMpeak)/b2);
            }
            /**** Force Evaluation - End ****/
            // /////// //////////// /////// //

            /**** Langevin - Final ****/
            rx += dtOzeta*(al + 0.5*(force + forceTemp));

            t +=  dt;
        }

    }

    avgElapsTime = elapsedTime/breakCount;

    ofstream Arr;
        Arr.open("Arrhenius.dat"); 
        Arr << avgElapsTime << "\t" << g << endl;
    Arr.close();

    return 0;
}






