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
int seed = 0, samples = 0, Iterations = 0, binSize = 0, startRecordingAt = 0; // Statistics
int beadNumber = 0; // bead number
double dt = 0; // Step size
double zeta = 0, gama = 0, ks = 0, g = 0; // g = 2*zeta*kb*T // Dynamics

/**** Int to give position in loops ****/
const int inTemp = 0;
const int inFinal = 1;

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
    // //// Statistics //// //
    binSize = arguments[4];
    startRecordingAt = arguments[5];
    // //// Prop to Simul Time - Number of Beads //// //
    beadNumber = arguments[6];
    // //// Parameters of Interactions //// //
    gama = arguments[7];
    ks = arguments[8];
    zeta = arguments[9];
    g = arguments[10];

}

/**** [FUNCTION] Evaluate Force [FUNCTION] ****/
inline void EvalForce(int loopPos, int beadnumber, std::vector<double> &rx, std::vector<double> &ry,
                std::vector<double> &fx, std::vector<double> &fy){ 
    // // //////// loopPos -> 0 it's temp, 1 it's not //////// //

    double sig = 0;

    /**** Loop Variables ****/
    int n = 0;

    /**** Distance and Force Variables ****/
    double dxL = 0, dyL = 0, dxR = 0, dyR = 0, dL = 0, dL2 = 0, dR = 0, dR2 = 0, f = 0;

    if(loopPos == inTemp){

        for(n = 0; n < beadnumber; n++){

            // ///////////////////////////////////////    /////////////////////////////////////// //
            // We sum the number of beads here so we can deal with the second part of the array //
            // that corresponds to the temporary positions                                      //
            // ///////////////////////////////////////    /////////////////////////////////////// //

            /**** Distance between 2 Temporary Points ****/ 
            // ///// First Bead ///// //
            if(n == 0){
                dxR = rx[n + 1 + beadnumber] - rx[n + beadnumber]; // Right bead
                dyR = ry[n + 1 + beadnumber] - ry[n + beadnumber]; // Right bead
                dxL = 0; // No Left Bead
                dyL = 0; // No Left Bead
            }
            // ///// Last Bead ///// //
            else if(n == beadnumber - 1){
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = 0; // No Right Bead
                dyR = 0; // No Right Bead
            }
            // ///// Rest of the Chain ///// //
            else{
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = rx[n + 1 + beadnumber] - rx[n + beadnumber]; // Right bead
                dyR = ry[n + 1 + beadnumber] - ry[n + beadnumber]; // Right bead
            }

            /**** Force Calculation ****/
            fx[n + beadnumber] = ks*(dxL + dxR);
            fy[n + beadnumber] = ks*(dyL + dyR);

        } //end for
    } // end if for temporary position

    else if(loopPos == inFinal){

        for(n = 0; n < beadnumber; n++){

            /**** Distance between 2 Points ****/
            // ///// First Bead ///// //
            if(n == 0){
                dxR = rx[n+1] - rx[n]; // Right bead
                dyR = ry[n+1] - ry[n]; // Right bead
                dxL = 0; // No Left Bead
                dyL = 0; // No Left Bead
            }
            // ///// Last Bead ///// //
            else if(n == beadnumber - 1){
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = 0; // No Right Bead
                dyR = 0; // No Right Bead
            }
            // ///// Rest of the Chain ///// //
            else{
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = rx[n+1] - rx[n]; // Right bead
                dyR = ry[n+1] - ry[n]; // Right bead
            }

            /**** Force Calculation ****/
            fx[n] = ks*(dxL + dxR);
            fy[n] = ks*(dyL + dyR);

        } // end for

    } // end else if for final positions

} // end function


int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    /**** SRK2 Vars ****/
    vector<double> rx(2*beadNumber); // from the begining to (startBeadNumber - 1) is position
    vector<double> ry(2*beadNumber); // from startBeadNumber to the end is temporary position
    vector<double> fx(2*beadNumber); // from the begining to (startBeadNumber - 1) is force
    vector<double> fy(2*beadNumber); // from startBeadNumber to the end is temporary force
    vector<double> alx(beadNumber); // random portion for x
    vector<double> aly(beadNumber); // random portion for y

    const int arraySizeSimul = (Iterations - startRecordingAt)/binSize; // gives the # of entries on the files that plot the
                                          // e2e distance in respect to time

    /**** For Loop Variables ****/
    int pos = 0, t = 0, p = 0, spl = 0, InitCheck;
    const double timeJump = dt*binSize;
    double tS = (startRecordingAt + binSize)*dt; // we ignore the first iterations in order to reach a state near equilibrium (what matters)

    /**** Others ****/
    int tt = 0;

    /**** Random stuff ****/
    double std = sqrt(g/dt); // usually it's "g/dt" but since dt will change I divide after
    std::mt19937 generator (seed);
    std::normal_distribution<double> normal(0, std);

    /**** File Names ****/
    char filename_007[100];

    /**** Statistics ****/
    vector<double> sig(arraySizeSimul);
    double sigDouble = 0;
    int ifCheck = 0;

    /**** Calc constant shit before loops ****/
    const double dtOzeta = dt/zeta;
    const double dtTgamma = dt*gama;

    /**** Keep track of used parameters ****/
    cout << "    *** Global Parameters *** \n \nSamples: " << samples << "\n" << "Number of Iterations: " <<
    Iterations << "\n" << "Time Step: " << dt << "\n\n" << " Iterations" << "\n\n" << "Recorded last " << 
    Iterations - startRecordingAt << " Iterations" << "\n" << "Bin Size: " << binSize << "\n\n" << "# of Beads: " 
    << beadNumber << "\n\n" << "gamma = " << gama <<"\n" << "ks = " << ks << "\n" << "zeta = " << zeta << "\n" <<
    "g = " << g << "\n\n    *** Variable Parameters *** \n" <<  endl;


    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        /**** Assign Initial Positions from the Previous Sample ****/
        std::fill(rx.begin(), rx.end(), 0);
        std::fill(ry.begin(), ry.end(), 0);
        
        tt = 0;

        /**** Calculate Stuff for some Time ****/
        for(t = 0; t < Iterations; t++){

            /**** Evaluation of Force ***/
            EvalForce(inFinal, beadNumber, rx, ry, fx, fy);
            
            /**** Calculate Temporary Position ****/
            for(p = 0; p < beadNumber; p++){

                /**** Set Random Motion ****/
                // do{ 
                //     alx[p] = normal(generator)*c2;
                // }
                // while(alx[p]*alx[p] > 4*std*c2*c2);
                // do{
                //     aly[p] = normal(generator)*c2;
                // }
                // while(aly[p]*aly[p] > 4*std*c2*c2);
                
                // //// DISTRIBUTION NOT TRUNCATED - FOR TESTING ONLY // ////
                alx[p] = normal(generator);
                aly[p] = normal(generator);
                //        //////////// //////////// ////////////          //

                /**** Temporary Position ****/
                rx[p + beadNumber] = rx[p] + dtTgamma*ry[p] + dtOzeta*(alx[p] + fx[p]);
                ry[p + beadNumber] = ry[p] + dtOzeta*(aly[p] + fy[p]);

            }

            /**** Evaluation of Force - Temp ***/
            EvalForce(inTemp, beadNumber, rx, ry, fx, fy);

            /**** Calculate Position ****/
            for(p = 0; p < beadNumber; p++){

                rx[p] += dtTgamma*ry[p] + dtOzeta*(alx[p] + 0.5*(fx[p] + fx[p + beadNumber])); // f           +       f_temporary
                ry[p] += dtOzeta*(aly[p] + 0.5*(fy[p] + fy[p + beadNumber])); //            0 to (beadNumber - 1)    beadNumber to end

                /**** Statistics ****/
                if((t > startRecordingAt) && ((t + 1) % binSize == 0)){ // record at every "bin" timesteps after we leave the begining of simul
                    sig[tt] -= fx[p]*ry[p];
                    ifCheck = 1;
                }
                else{ifCheck = 0;}
            }

            if(ifCheck == 1){ // checks if we recorded something in this loop or not - if yes, then +1 so we can go foward in the array
                tt++;
            }

        } // Closes iteration loop

    } // Closes samples loop

    /**** Output file with pos of bead ****/
    // /// !! Calculo o Max no R !! /// //
    ofstream sp;
        sp.open("gcShear.dat");
    for (auto&& i : sig){ 
        sp << i/samples << "\t" << tS << endl;
        tS += timeJump; // time passed
        // tS++;
    } 
    sp.close();

    return 0;
}