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

// ///////////////////////////////////////////////////////////////////////////////// //
// //////////////////////////// /**** GLOBAL VARS ****/ //////////////////////////// //
// ///////////////////////////////////////////////////////////////////////////////// //

/**** Parameter Initialization ****/
// //// Here we create global Variables. The values are assigned in the function "readFile()" //// //
int seed = 0, samples = 0, Iterations = 0, binSize = 0; // Statistics
int beadNumber = 0, chainNumber = 0; // Each chain has 'beadNumber' beads and there are 'chainNumber' chains
double density = 0;
double dt = 0; // Step size
double zeta = 0, ks = 0, g = 0; // g = 2*zeta*kb*T // Dynamics
// //// Potential Parameters //// //
double sigma = 0, epsilon = 0, peak = 0, a = 0, b = 0;


// ///////////////////////////////////////////////////////////////////////////////// //
// ///////////////////////////// /**** FUNCTIONS ****/ ///////////////////////////// //
// ///////////////////////////////////////////////////////////////////////////////// //

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
    // //// Chains //// //
    beadNumber = arguments[5]; // # of beads per chain
    chainNumber = arguments[6]; // # of chains
    // //// Box And Initial Conditions  //// //
    density = arguments[7];
    // //// Parameters of Interactions //// //
    ks = arguments[8];
    zeta = arguments[9];
    g = arguments[10];
    // //// Parameters of the Potential //// //
    sigma = arguments[11];
    epsilon = arguments[12];
    peak = arguments[13];
    a = arguments[14];
    b = arguments[15];

}

inline void nic(double& x, double& y, const double& BoxLength, const double& HalfLength){

    if(x > HalfLength){
        x -= BoxLength;
    }
    else if(x < -HalfLength){
        x += BoxLength;
    }
    if(y > HalfLength){
        y -= BoxLength;
    }
    else if(y < -HalfLength){
        y += BoxLength;
    }

}

inline void pbc(double& x, double& y, const double& BoxLength){

    if(x < 0){
        x += BoxLength;
    }
    else if(x > BoxLength){
        x -= BoxLength;
    }
    if(y < 0){
        y += BoxLength;
    }
    else if(y > BoxLength){
        y -= BoxLength;
    }

}

void setInitialPositions(const double& BoxLength, const double& HalfLength,
                std::vector<double> &rx, std::vector<double> &ry){

    /**** Distances ****/
    double dx = 0, dy = 0;
    double dr = 0;
    double dr2 = 0;
    double minAprovedDist = 4; // Used in setting the initial positions

    /**** Random Stuff ****/
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uniBox(1, BoxLength - 1);
    std::uniform_real_distribution<double> uniSmall(-0.2, 0.2);

    /**** Loop Variables ****/
    int chain = 0, bead = 0, i = 0, j = 0, n = 0, ErrorCount = 0;

    /**** Points ****/
    double pointX = 0;
    double pointY = 0;

    /**** Others ****/
    int check = 1;
    int bPc = 0;

    check = 1;

    while(chain < chainNumber){

        /**** Assign a Point ****/
        pointX = uniBox(generator);
        pointY = uniBox(generator);

        if(chain == 0){
            for(bead = 0; bead < beadNumber; bead++){

                bPc = bead + chain*beadNumber;

                rx[bPc] = pointX + uniSmall(generator);
                ry[bPc] = pointY + uniSmall(generator);

            }
            chain++;
        }
            
        else{

            check = 1;

            for(j = 0; j < chain*beadNumber; j++){

                /**** Calculate Distance in 1D to Point ****/
                dx = pointX - rx[j];
                dy = pointY - ry[j];

                /**** Nearest Image Convention ****/
                nic(dx, dy, BoxLength, HalfLength);

                /**** Calculate Distance in 2D to Point ****/
                dr2 = dx*dx + dy*dy;

                /**** Check if it's viable ****/
                if(dr2 < minAprovedDist){
                    check *= 0; // it's not
                }
                else{
                    check *= 1; // it is
                }

            }

            if(check != 0){
            
                for(bead = 0; bead < beadNumber; bead++){

                    bPc = bead + chain*beadNumber;

                    rx[bPc] = pointX + uniSmall(generator);
                    ry[bPc] = pointY + uniSmall(generator);

                }

                chain++;

            }
        }
         
    }
}


/**** [FUNCTION] Evaluate Force [FUNCTION] ****/
inline void evalForce(double BoxLength, double HalfLength, int totalBeads,
                std::vector<double> &X, std::vector<double> &Y,
                std::vector<double> &springFx, std::vector<double> &springFy,
                std::vector<double> &repelFx, std::vector<double> &repelFy,
                std::vector<double> &bondFx, std::vector<double> &bondFy,
                double rc2, double rMin2, double peak2, double s6, double s12, double b2){

    /**** Loop Variables ****/
    int i = 0, j = 0;

    /**** Intra Force Variables ****/
    double dxL = 0, dyL = 0, dxR = 0, dyR = 0;

    /**** Inter Force Variables ****/
    double dx = 0, dy = 0, d = 0, d2 = 0, d6 = 0, d12 = 0, f = 0;
    double dMpeak = 0;

    /**** Set Force to 0 because we sum over it ****/
    std::fill(repelFx.begin(), repelFx.end(), 0);
    std::fill(bondFx.begin(), bondFx.end(), 0);
    std::fill(repelFy.begin(), repelFy.end(), 0);
    std::fill(bondFy.begin(), bondFy.end(), 0);
        
    for(i = 0; i < totalBeads-1; i++){

        // ///// First Bead ///// //
        if(i%beadNumber == 0){

            for(j = i + 1; j < totalBeads; j++){

                dx = X[j] - X[i];
                dy = Y[j] - Y[i];

                /**** Nearest Image Convention ****/
                nic(dx, dy, BoxLength, HalfLength);

                d2 = dx*dx + dy*dy;

                if(j == i+1){
                    
                    // Save position to use next
                    dxR = dx;
                    dyR = dy; 

                    springFx[i] = ks*dxR;
                    springFy[i] = ks*dyR;

                }

                if(i/beadNumber != j/beadNumber){

                    if(j%beadNumber == 0 || j%beadNumber == beadNumber-1){

                        /**** Check if Beads are close enough to Interact ****/
                        if(d2 < rc2){
                            
                            /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                            d = sqrt(d2);
                            d6 = d2*d2*d2;
                            d12 = d6*d6;
                            dMpeak = d - peak;

                            /**** Calculate Force ****/
                            f = (48*epsilon/d)*(s12/d12 - 0.5*s6/d6) + ((2*a*(dMpeak))/(b2))*exp(-(dMpeak*dMpeak)/(b2));

                            bondFx[j] += dx*f/d;
                            bondFx[i] -= dx*f/d;
                            bondFy[j] += dy*f/d;
                            bondFy[i] -= dy*f/d;

                        } 
                    } // end Attractive 'if'

                    else{

                        /**** Check if Beads are close enough to Interact ****/
                        if(d2 < rMin2){
                            
                            /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                            d6 = d2*d2*d2;
                            d12 = d6*d6;

                            /**** Calculate Force ****/
                            f = (48*epsilon)*(s12/d12 - 0.5*s6/d6);

                            repelFx[j] += dx*f/d2;
                            repelFx[i] -= dx*f/d2;
                            repelFy[j] += dy*f/d2;
                            repelFy[i] -= dy*f/d2;

                        }
                    } // end Repulsive 'if'
                } // end 'if' that checks if we are not on the same chain 
            } // end 'for' loop for every bead on the right
        } // end 'if' that checks if its the FIRST bead

        // ///// Last Bead ///// //
        else if(i % beadNumber == beadNumber - 1){

            dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
            dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right

            springFx[i] = ks*dxL;
            springFy[i] = ks*dyL;

            for(j = i + 1; j < totalBeads; j++){

                dx = X[j] - X[i];
                dy = Y[j] - Y[i];

                /**** Nearest Image Convention ****/
                nic(dx, dy, BoxLength, HalfLength);

                d2 = dx*dx + dy*dy;

                if(i/beadNumber != j/beadNumber){

                    if(j%beadNumber == 0 || j%beadNumber == beadNumber-1){

                        /**** Check if Beads are close enough to Interact ****/
                        if(d2 < rc2){
                            
                            /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                            d = sqrt(d2);
                            d6 = d2*d2*d2;
                            d12 = d6*d6;
                            dMpeak = d - peak;

                            /**** Calculate Force ****/
                            f = (48*epsilon/d)*(s12/d12 - 0.5*s6/d6) + ((2*a*(dMpeak))/(b2))*exp(-(dMpeak*dMpeak)/(b2));

                            bondFx[j] += dx*f/d;
                            bondFx[i] -= dx*f/d;
                            bondFy[j] += dy*f/d;
                            bondFy[i] -= dy*f/d;

                        } 
                    } // end Attractive 'if'

                    else{

                        /**** Check if Beads are close enough to Interact ****/
                        if(d2 < rMin2){
                            
                            /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                            d6 = d2*d2*d2;
                            d12 = d6*d6;

                            /**** Calculate Force ****/
                            f = (48*epsilon)*(s12/d12 - 0.5*s6/d6);

                            repelFx[j] += dx*f/d2;
                            repelFx[i] -= dx*f/d2;
                            repelFy[j] += dy*f/d2;
                            repelFy[i] -= dy*f/d2;

                        }
                    } // end Repulsive 'else'
                } // end 'if' that checks if we are not on the same chain 
            } // end 'for' loop for every bead on the right
        } // end 'if' that checks if its the LAST bead

        else{

            for(j = i + 1; j < totalBeads; j++){

                dx = X[j] - X[i];
                dy = Y[j] - Y[i];

                /**** Nearest Image Convention ****/
                nic(dx, dy, BoxLength, HalfLength);

                d2 = dx*dx + dy*dy;

                if(j == i+1){
                    
                    // The distance to the left bead is equal to the previous distance to the right
                    dxL = -dxR; 
                    dyL = -dyR;
                    
                    // Save the value for the distance on the right so it can be used for the left later ^^
                    dxR = dx;
                    dyR = dy; 

                    springFx[i] = ks*(dxL + dxR);
                    springFy[i] = ks*(dyL + dyR);

                }

                if(i/beadNumber != j/beadNumber){

                    /**** Check if Beads are close enough to Interact ****/
                    if(d2 < rMin2){
                        
                        /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                        d6 = d2*d2*d2;
                        d12 = d6*d6;

                        /**** Calculate Force ****/
                        f = (48*epsilon)*(s12/d12 - 0.5*s6/d6);

                        repelFx[j] += dx*f/d2;
                        repelFx[i] -= dx*f/d2;
                        repelFy[j] += dy*f/d2;
                        repelFy[i] -= dy*f/d2;

                    }
                } // end 'if' that checks if we are not on the same chain 
            } // end 'for' loop for every bead on the right
        } // end 'else' that runs comands for the middle chain
    } // end 'for' that runs over every bead-1

    /****************************************************************************/
    /**** the last bead to be checked has had all it's forces calculated in
    * respect to the others leaving only the harmonic potential to be calculated.
    * We already have the distance to the one on the left
    * from the last bead that was calculated inside the loop ****/
    dxL = -dxR;
    dyL = -dyR; 

    springFx[totalBeads-1] = ks*dxL;
    springFy[totalBeads-1] = ks*dyL;

    /****************************************************************************/
}

/**** [FUNCTION] Check Bonds [FUNCTION] ****/
inline void Mi6(int t, double BoxLength, double HalfLength, int totalBeads, int numOfPossibleBonds,
        std::vector<double> &X, std::vector<double> &Y, double peak2,
        vector<int>& bondGPSprev, vector<int>& formedVector, vector<int>& brokenVector, vector<int>& numBondsVector){

    /**** Loop Variables ****/
    int i = 0, j = 0, n = 0; 

    /**** Distance ****/
    double dx = 0, dy = 0, d2 = 0;

    /**** Statistics ****/
    int formedBonds = 0, brokenBonds = 0, numberOfBonds = 0;       

    /**** This doesn't need to exit the function ****/
    vector<int> bondGPSnow(numOfPossibleBonds);

    /* ERROR!!! ERROR!!! ERROR!!! */
    for(i = 0; i < totalBeads-1; i++){

        if(i%beadNumber == 0 || i%beadNumber == beadNumber-1){

            for(j = i+1; j < totalBeads; j++){
            
                if(i/beadNumber != j/beadNumber){

                    if(j%beadNumber == 0 || j%beadNumber == beadNumber-1){

                        dx = X[j] - X[i];
                        dy = Y[j] - Y[i];

                        /**** Nearest Image Convention ****/
                        nic(dx, dy, BoxLength, HalfLength);

                        d2 = dx*dx + dy*dy;

                        /**** Check if there is a bond or not ****/
                        if(d2 < peak2){
                            bondGPSnow[n] = 1;
                            numberOfBonds++;
                        }
                        else{
                            bondGPSnow[n] = 0;
                        }

                        /**** Compare the current and previous times and see if a bond has formed or broken ****/
                        if(bondGPSprev[n] < bondGPSnow[n]){
                            formedBonds++;
                        }
                        else if(bondGPSprev[n] > bondGPSnow[n]){
                            brokenBonds++;
                        }

                        bondGPSprev[n] = bondGPSnow[n];

                        n++;

                    } // end of 'if' to see if 'j' interacts or not 
                } // end of 'if' to see if 'j' and 'i' are in different chains 
            } // end 'j' loop
        } // end of 'if' to see if 'i' interacts or not 
    } // end 'i' loop

    formedVector[t] += formedBonds; // keeps adding to then perform an avg
    brokenVector[t] += brokenBonds;
    numBondsVector[t] += numberOfBonds;

}

/**** [FUNCTION] Make Video [FUNCTION] ****/
// inline void Cannes(std::ofstream& video, int t, int attractingBeads, int totalBeads, int numOfPossibleBonds, std::vector<int> James, std::vector<double> X, std::vector<double> Y){

//     std::vector<int> bond(attractingBeads);

//     /**** Loop Variables ****/
//     int n = 0, r = 0, var1 = 0, var2 = 0, var3 = 0, atom = 1;

//     for(var1 = 0; var1 < attractingBeads-1; var1++){

//         if(var1 % 2 == 0){

//             for(var2 = var1+2; var2 < attractingBeads; var2++){
            
//                 if(James[n] == 1){
//                     bond[var1]++;
//                     bond[var2]++;
//                     n++;
//                 }
//                 else{n++;}
                
//             }

//         }

//         else{

//             for(var2 = var1+1; var2 < attractingBeads; var2++){
            
//                 if(James[n] == 1){
//                     bond[var1]++;
//                     bond[var2]++;
//                     n++;
//                 }
//                 else{n++;}
                
//             }

//         }
//     }

//     video << totalBeads << "\n" << "t = " << t << endl;

//     for(r = 0; r < totalBeads; r++){

//         if(r%beadNumber == 0 || r%beadNumber == beadNumber-1){
//             atom = bond[var3]+2;
//             var3++;
//         }
//         else{
//             atom = 1;
//         }

//         video << atom << "\t" << X[r] << "\t"<< Y[r] << "\t" << 0.0 << endl;
//     }
// }

// //////////////////////////////////////////////////////////////////////////////// //
// /////////////////////////////// /**** MAIN ****/ /////////////////////////////// //
// //////////////////////////////////////////////////////////////////////////////// //

int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    const int totalBeads = chainNumber*beadNumber;
    const int attractingBeads = chainNumber*2;

    /**** Potential Constants ****/
    const long double rMin = 1.12246204830937298*sigma;
    const long double rMin2 = rMin*rMin; 
    const long double rCut = 2.5*rMin;
    const double rc2 = rCut*rCut;
    const double peak2 = peak*peak;
    const double s2 = sigma*sigma;
    const double s6 = s2*s2*s2;
    const double s12 = s6*s6;
    const double b2 = b*b;

    const double BoxLength = sqrt(attractingBeads/density);
    const double HalfLength = BoxLength*0.5;
    /* this goes in main and then we send to the functions so that
    so that I don't have to calculate this over and over again */

    /**** Position ****/
    vector<double> rx(totalBeads); 
    vector<double> ry(totalBeads); 
    vector<double> rx_temp(totalBeads); 
    vector<double> ry_temp(totalBeads); 

    /**** Force ****/
    vector<double> fSpring_x(totalBeads);
    vector<double> fSpring_y(totalBeads);
    vector<double> fSpring_x_temp(totalBeads);
    vector<double> fSpring_y_temp(totalBeads);
    vector<double> fBond_x(totalBeads);
    vector<double> fBond_y(totalBeads);
    vector<double> fBond_x_temp(totalBeads);
    vector<double> fBond_y_temp(totalBeads);
    vector<double> fRepel_x(totalBeads);
    vector<double> fRepel_y(totalBeads);
    vector<double> fRepel_x_temp(totalBeads);
    vector<double> fRepel_y_temp(totalBeads);

    /**** Thermal Noise ****/
    vector<double> alx(totalBeads); 
    vector<double> aly(totalBeads);

    const int arraySizeSimul = Iterations/binSize;

    /**** For Loop Variables ****/
    int pos = 0, t = 0, p = 0, j = 0, i = 0, alberto = 0, spl = 0, arrayPos = 0;
    const double timeJump = dt*binSize;
    double tS = 0; // I'm starting at 0 time... might change

    /**** # of Possible Bonds (Non Repeated) ****/
    for(i = 1; i < attractingBeads; i++){
        alberto += i;
    }
    const int numOfPossibleBonds = alberto; 

    /**** Video ****/
    // int recSample = 1;

    /**** Random stuff ****/
    double std = sqrt(g/dt); // usually it's "g/dt" but since dt will change I divide after
    std::mt19937 generator (seed);
    std::normal_distribution<double> normal(0, std);

    /**** Statistics ****/
    vector<int> bondGPSprev(numOfPossibleBonds);
    vector<int> formedVector(arraySizeSimul);
    vector<int> brokenVector(arraySizeSimul);
    vector<int> numBondsVector(arraySizeSimul);

    /**** Calc constant shit before loops ****/
    const double dtOzeta = dt/zeta;
    const double samplesThalfBeads = samples*attractingBeads*0.5;
    const double truncLimit = 4*std*std*c2*c2;

    /**** Keep track of used parameters ****/
    cout << "\n\t*** Global Parameters *** \n \nSamples: " << samples << "\n" << "Number of Iterations: " <<
    Iterations << "\n" << "Time Step: " << dt << "\n\n" << "Bin Size: " << binSize << "\n\n" << "Number of Chains: " << chainNumber << "\n" 
    << "Number of Beads per Chain: " << beadNumber << "\n" << "Box Side Length: " << BoxLength << "\n\n\t*** Potential Parameters *** \n\nsigma = " << sigma << "\n" <<
    "epsilon = " << epsilon << "\n" << "a = " << a << "\n" << "b = " << b << "\n" << "Peak of Potential Wall at " << peak <<
    "\n\n\t*** Other Parameters *** \n\n" << "zeta = " << zeta << "\n" << "ks = " << ks << "\n" << "g = " << g << "\n\n" << endl;

    // ofstream video;
    //     video.open("repel1.imd");
    //     video.precision(8);
    //     video.setf(ios::fixed | ios::showpoint);

    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        /**** Assign Initial Positions ****/
        setInitialPositions(BoxLength, HalfLength, rx, ry);

        /**** Set this to 0 so it doesn't compare the begining
         * of the sample with the end of the previous one ****/
        std::fill(bondGPSprev.begin(), bondGPSprev.end(), 0);
        
        arrayPos = 0;

        /**** Calculate Stuff for some Time ****/
        for(t = 0; t < Iterations; t++){

            /**** Evaluation of Force ***/
            evalForce(BoxLength, HalfLength, totalBeads, rx, ry, 
                    fSpring_x, fSpring_y, fRepel_x, fRepel_y, fBond_x, fBond_y,
                    rc2, rMin2, peak2, s6, s12, b2);

            // ////////////////// /**** Statistics ****/ ////////////////// //
            /*  Saves the number of formed and broken bonds as well as the number of bonds. 

                Need to divide the first 2 by the # of samples and elapsed time (binSize). 
                The #Bonds needs to be divided by the # of samples, elapsed time and # of beads

                This will be done in the end    */

            if(t % binSize == 0){
                Mi6(arrayPos, BoxLength, HalfLength, totalBeads, numOfPossibleBonds, rx, ry, 
                    peak2, bondGPSprev, formedVector, brokenVector, numBondsVector);
                arrayPos++;
            }
            
            /**** Calculate Temporary Position ****/
            for(p = 0; p < totalBeads; p++){

                /**** Set Random Motion ****/
                do{ 
                    alx[p] = normal(generator)*c2;
                }
                while(alx[p]*alx[p] > truncLimit);
                do{
                    aly[p] = normal(generator)*c2;
                }
                while(aly[p]*aly[p] > truncLimit);
                
                // //// DISTRIBUTION NOT TRUNCATED - FOR TESTING ONLY // ////
                // alx[p] = normal(generator);
                // aly[p] = normal(generator);
                //        //////////// //////////// ////////////          //

                /**** Temporary Position ****/
                rx_temp[p] = rx[p] + dtOzeta*(alx[p] + fSpring_x[p] + fRepel_x[p] + fBond_x[p]);
                ry_temp[p] = ry[p] + dtOzeta*(aly[p] + fSpring_y[p] + fRepel_y[p] + fBond_y[p]);

                pbc(rx_temp[p], ry_temp[p], BoxLength);

            }

            /**** Evaluation of Force - Temp ***/
            evalForce(BoxLength, HalfLength, totalBeads, rx_temp, ry_temp, 
                    fSpring_x_temp, fSpring_y_temp, fRepel_x_temp, fRepel_y_temp, fBond_x_temp, fBond_y_temp,
                    rc2, rMin2, peak2, s6, s12, b2);

            /**** Create Video of Particles Moving - Parameters ****/
            // if(spl == recSample && t % binSize == 0){
            //     Cannes(video, t, attractingBeads, totalBeads, numOfPossibleBonds, bondGPSprev, rx, ry);
            // }

            /**** Calculate Position ****/
            for(p = 0; p < totalBeads; p++){

                rx[p] += dtOzeta*(alx[p] + 0.5*(fSpring_x[p] + fRepel_x[p] + fBond_x[p] + fSpring_x_temp[p] + fRepel_x_temp[p] + fBond_x_temp[p]));
                ry[p] += dtOzeta*(aly[p] + 0.5*(fSpring_y[p] + fRepel_y[p] + fBond_y[p] + fSpring_y_temp[p] + fRepel_y_temp[p] + fBond_y_temp[p]));

                pbc(rx[p], ry[p], BoxLength);

            }

        } // Closes iteration loop

    } // Closes samples loop

    // video.close();

    ofstream file1;
        file1.open("out1.dat");

    for(j = 0; j < arraySizeSimul; j++){ 
        file1 << ((double)formedVector[j])/((double)samples) << "\t" << ((double)brokenVector[j])/((double)samples) <<
            "\t" << numBondsVector[j]/samplesThalfBeads << "\t" << tS << endl;
        tS += timeJump; // time passed
    } 

    file1.close();
    
    return 0;
}