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
int videoFrames = 0;
double dt = 0; // Step size
double zeta = 0, ks = 0, g = 0; // g = 2*zeta*kb*T // Dynamics

/**** [FUNCTION] Read "data.txt" file [FUNCTION] ****/
void readFile(){

    /**** Reads File with Variables in use ****/
    ifstream infile("data4.txt"); // Read "data.txt" file
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
    // videoFrames = arguments[5];
    // //// Parameters of Interactions //// //
    ks = arguments[5];
    zeta = arguments[6];
    g = arguments[7];

}


int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    /**** SRK2 Vars ****/
    // !!!! FALTA METER AQUI AS CONDIÇÕES INICIAIS !!!! //
    vector<double> posX(2);
    vector<double> posY(2);
    vector<double> posXTemp(2);
    vector<double> posYTemp(2);
    vector<double> fx(4);
    vector<double> fy(4);
    double dx = 0;
    double dy = 0;

    const int arraySizeSimul = Iterations/binSize; // gives the # of entries on the files that plot the
                                        // e2e distance in respect to time
    /**** For Loop Variables ****/
    int spl = 0, t = 0, tt = 0, p = 0;
    const double timeJump = dt*binSize;
    double tS = timeJump;

    /**** Statistics ****/
    vector<double> savePos(arraySizeSimul);
    // vector<double> savePos(500);

    /**** Random stuff ****/
    vector<double> alx(2);
    vector<double> aly(2);
    double std = sqrt(g/dt);
    std::mt19937 generator(seed);
    std::normal_distribution<double> normal(0, std);
 
    /**** Calc constant shit before loops ****/
    double ksOzetaTdt = (ks/zeta)*dt;
    double dtOzeta = dt/zeta;
    const double truncLimit = 4*std*std*c2*c2;

    /**** Keep track of used parameters ****/
    cout << "\n\t*** Global Parameters *** \n \nSamples: " << samples << "\n" << "Number of Iterations: " <<
    Iterations << "\n" << "Time Step: " << dt << "\n\n" << "Bin Size: " << binSize << "\n\n\t*** Other Parameters *** \n\n" <<
    "ks = " << ks << "\n" << "zeta = " << zeta << "\n" << "g = " << g << "\n" << endl;

    /**** Create Video File ****/ 
    // int recordVideoHere = samples/2; // Gives the Recorded Sample
    // int videoTimeJump = Iterations/videoFrames; // jump 'videoTimeJump' iterations
    // std::string atom = "C"; // make pretty
    // std::string videoFileName = "videoFile.xyz";
    // ofstream video;
    //     video.open(videoFileName);
    //     video.precision(8);
    //     video.setf(ios::fixed | ios::showpoint);

    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        tt = 0;

        /**** Initial Conditions ****/
        posX[0] = 0;
        posY[0] = 0;
        posX[1] = 0;
        posY[1] = 0;

        /**** Calculate stuff for some time ****/
        for(t = 0; t < Iterations; t++){

            // //////// ///////////// //////// //
            /**** Force Evaluation - Begin ****/
            // //////// ///////////// //////// //

            /**** Distance Between Monomers ****/
            dx = posX[1] - posX[0];
            dy = posY[1] - posY[0];

            /**** Force Evaluation ****/
            fx[0] = ks*dx;
            fx[1] = -fx[0];
            fy[0] = ks*dy;
            fy[1] = -fy[0];
            
            // /////// //////////// /////// //
            /**** Force Evaluation - End ****/
            // /////// //////////// /////// //

            /**** Statistics ****/
            if((t + 1) % (binSize) == 0){
            // if(t < 500){
                // //// Saves the position value at iteration X //// //
                savePos[tt] += dx*dx + dy*dy; // keeps adding to then perform an avg
                tt++; // counter to run the array above
            }

            /**** Langevin - Temporary ****/
            for(p = 0; p <= 1; p++){

                /**** Set Random Motion ****/
                do{ 
                    alx[p] = normal(generator)*c2;
                }
                while(alx[p]*alx[p] > truncLimit);
                do{
                    aly[p] = normal(generator)*c2;
                }
                while(aly[p]*aly[p] > truncLimit);

                // //// DISTRIBUTION NOT TRUNCATED - FOR TESTING ONLY //// //
                // alx[p] = normal(generator);
                // aly[p] = normal(generator);
                //        //////////// //////////// ////////////        //

                /**** Position ****/
                posXTemp[p] = posX[p] + dtOzeta*(alx[p] + fx[p]);
                posYTemp[p] = posY[p] + dtOzeta*(aly[p] + fy[p]);

            }

            // //////// ///////////// //////// //
            /**** Force Evaluation - Begin ****/
            // //////// ///////////// //////// //

            /**** Distance Between Monomers ****/
            dx = posXTemp[1] - posXTemp[0];
            dy = posYTemp[1] - posYTemp[0];

            /**** Force Evaluation ****/
            fx[2] = ks*dx;
            fx[3] = -fx[0];
            fy[2] = ks*dy;
            fy[3] = -fy[0];

            // /////// //////////// /////// //
            /**** Force Evaluation - End ****/
            // /////// //////////// /////// //

            /**** Create Video of Particles Moving - Parameters ****/
            // if(spl == recordVideoHere && t % videoTimeJump == 0){
            //     video << 2 << "\n" << "t = " << t << endl;
            // }

            /**** Langevin - Final ****/
            for(p = 0; p <= 1; p++){

                // /**** Create Video of Particles Moving - Position ****/
                // if(spl == recordVideoHere && t % videoTimeJump == 0){
                //     video << atom << "\t" << posX[p] << "\t"<< posY[p] << "\t" << 0.0 << endl;
                // }
                /**** Positions - Langevin ****/
                posX[p] += dtOzeta*(alx[p] + 0.5*(fx[p] + fx[p + 2]));
                posY[p] += dtOzeta*(aly[p] + 0.5*(fy[p] + fy[p + 2]));

            }
        }// time loop
    }// samples loop

    /**** Output file with pos of bead ****/
    // /// !! Calculo o Max no R !! /// //
    ofstream sp;
        sp.open("spring.dat");
    for (auto&& i : savePos){ 
        sp << i/samples << "\t" << tS << endl;
        tS += timeJump; // time passed
        // tS++;
    } 
    sp.close();

    // tS = timeJump; // set to "timeJump" again
    //     /////// /////// ///////     //

    return 0;

}