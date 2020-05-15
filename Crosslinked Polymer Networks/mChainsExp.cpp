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
int chainSize = 0, chainNumber = 0; // Each chain has 'chainSize' beads and there are 'chainNumber' chains
double dt = 0; // Step size
double density = 0;
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
    chainSize = arguments[5]; // # of beads per chain
    chainNumber = arguments[6]; // # of chains
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

/**** [FUNCTION] Evaluate Force [FUNCTION] ****/
inline void evalForce(const int& t, const int& numOfPossibleBonds, const double& BoxLength, const double& HalfLength, 
                const int& middleOfChain, const int& totalBeads,
                std::vector<double> &rx, std::vector<double> &ry,
                std::vector<double> &intraFx, std::vector<double> &intraFy,
                std::vector<double> &interFx, std::vector<double> &interFy,
                const double& rc2, const double& peak2, const double& s6, const double& s12, const double& b2, 
                vector<int>& bondGPSprev, int& formedBonds, int& brokenBonds, int& numberOfBonds){

    int bPc = 0;

    /**** Loop Variables ****/
    int chain = 0, bead = 0, next = 0, nb = 0;

    /**** Intra Force Variables ****/
    double dxL = 0, dyL = 0, dxR = 0, dyR = 0, dL = 0, dL2 = 0, dR = 0, dR2 = 0;

    /**** Inter Force Variables ****/
    double dx = 0, dy = 0, d = 0, d2 = 0, d6 = 0, d12 = 0, f = 0;
    double dMpeak = 0;

    /**** This doesn't need to exit the function ****/
    vector<int> BondGPSnow(numOfPossibleBonds);

    /**** Set Force to 0 because we sum interFx and interFy ****/
    std::fill(interFx.begin(), interFx.end(), 0);
    std::fill(interFy.begin(), interFy.end(), 0);

        
    for(chain = 0; chain < chainNumber; chain ++){
        for(bead = 0; bead < chainSize; bead++){

            bPc = bead + chain*chainSize;

            // //////////////////////////////////////// //
            // ///// /**** INTRA Bead Force ****/ ///// //
            // //////////////////////////////////////// //

            /**** Distance between 2 Points ****/
            // ///// First Bead ///// //
            if(bead == 0){
                dxR = rx[bPc + 1] - rx[bPc]; // Right bead
                dyR = ry[bPc + 1] - ry[bPc]; // Right bead
                dxL = 0; // No Left Bead
                dyL = 0; // No Left Bead
            }
            // ///// Last Bead ///// //
            else if(bead == chainSize - 1){
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = 0; // No Right Bead
                dyR = 0; // No Right Bead
            }
            // ///// Rest of the Chain ///// //
            else{
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = rx[bPc + 1] - rx[bPc]; // Right bead
                dyR = ry[bPc + 1] - ry[bPc]; // Right bead
            }

            /**** Nearest Image Convention ****/
            nic(dxR, dyR, BoxLength, HalfLength);
            nic(dxL, dyL, BoxLength, HalfLength);

            /**** Force Calculation ****/
            intraFx[bPc] = ks*(dxL + dxR);
            intraFy[bPc] = ks*(dyL + dyR);

            // //////////////////////////////////////// //
            // ///// /**** INTER Bead Force ****/ ///// //
            // //////////////////////////////////////// //

            if(bead == middleOfChain){ // middle bead has LJ interactions

                for(next = bPc + chainSize; next < totalBeads; next += chainSize){ // only works for 1 bead in the middle
        
                    /**** Distance between 2 points ****/
                    dx = rx[bPc] - rx[next];
                    dy = ry[bPc] - ry[next];

                    /**** Nearest Image Convention ****/
                    nic(dx, dy, BoxLength, HalfLength);

                    d2 = dx*dx + dy*dy;

                    /**** Check if Beads are close enough to Interact ****/
                    if(d2 < rc2){
                        
                        /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                        d = sqrt(d2);
                        d6 = d2*d2*d2;
                        d12 = d6*d6;
                        dMpeak = d - peak;

                        /**** Calculate Force and Potential ****/
                        f = (48*epsilon/d)*(s12/d12 - 0.5*s6/d6) + ((2*a*(dMpeak))/(b2))*exp(-(dMpeak*dMpeak)/(b2));

                        interFx[bPc] += dx*f/d;
                        interFx[next] -= dx*f/d;
                        interFy[bPc] += dy*f/d;
                        interFy[next] -= dy*f/d;

                        /**** MI6 Statistics ****/
                        if(d2 < peak2){
                            BondGPSnow[nb] = 1;
                            numberOfBonds++;
                        }
                        else{
                            BondGPSnow[nb] = 0;
                        }

                    } // end distance if

                    else{
                        BondGPSnow[nb] = 0;
                    }

                    /**** Compare the current and previous times and see if a bond has formed or broken ****/
                    if(t != 0){ // no point in comparing at time 0

                        if(bondGPSprev[nb] < BondGPSnow[nb]){
                            formedBonds++;
                        }
                        else if(bondGPSprev[nb] > BondGPSnow[nb]){
                            brokenBonds++;
                        }
                    }

                    bondGPSprev[nb] = BondGPSnow[nb];
                    nb++;

                } // end for loop used to check next interacting beads 
            } // end if to see if the beads are of the interacting kind
        } // end for Bead
    } // end for Chain

}

/**** [FUNCTION] Evaluate Force [FUNCTION] ****/
inline void evalForce_temp(const double& BoxLength, const double& HalfLength, const int& middleOfChain, const int& totalBeads,
                std::vector<double> &rx_temp, std::vector<double> &ry_temp,
                std::vector<double> &intraFx_temp, std::vector<double> &intraFy_temp,
                std::vector<double> &interFx_temp, std::vector<double> &interFy_temp,
                const double& rc2, const double& peak2, const double& s6, const double& s12, const double& b2){

    int bPc = 0;

    /**** Loop Variables ****/
    int chain = 0, bead = 0, next = 0;

    /**** Intra Force Variables ****/
    double dxL = 0, dyL = 0, dxR = 0, dyR = 0, dL = 0, dL2 = 0, dR = 0, dR2 = 0;

    /**** Inter Force Variables ****/
    double dx = 0, dy = 0, d = 0, d2 = 0, d6 = 0, d12 = 0, f = 0;
    double dMpeak = 0;

    /**** Set Force to 0 because we sum interFx and interFy ****/
    std::fill(interFx_temp.begin(), interFx_temp.end(), 0);
    std::fill(interFy_temp.begin(), interFy_temp.end(), 0);

        
    for(chain = 0; chain < chainNumber; chain ++){
        for(bead = 0; bead < chainSize; bead++){

            bPc = bead + chain*chainSize;

            // //////////////////////////////////////// //
            // ///// /**** INTRA Bead Force ****/ ///// //
            // //////////////////////////////////////// //

            /**** Distance between 2 Points ****/
            // ///// First Bead ///// //
            if(bead == 0){
                dxR = rx_temp[bPc + 1] - rx_temp[bPc]; // Right bead
                dyR = ry_temp[bPc + 1] - ry_temp[bPc]; // Right bead
                dxL = 0; // No Left Bead
                dyL = 0; // No Left Bead
            }
            // ///// Last Bead ///// //
            else if(bead == chainSize - 1){
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = 0; // No Right Bead
                dyR = 0; // No Right Bead
            }
            // ///// Rest of the Chain ///// //
            else{
                dxL = -dxR; // The distance to the left bead is equal to the previous distance to the right
                dyL = -dyR; // The distance to the left bead is equal to the previous distance to the right
                dxR = rx_temp[bPc + 1] - rx_temp[bPc]; // Right bead
                dyR = ry_temp[bPc + 1] - ry_temp[bPc]; // Right bead
            }

            /**** Nearest Image Convention ****/
            nic(dxR, dyR, BoxLength, HalfLength);
            nic(dxL, dyL, BoxLength, HalfLength);

            /**** Force Calculation ****/
            intraFx_temp[bPc] = ks*(dxL + dxR);
            intraFy_temp[bPc] = ks*(dyL + dyR);

            // //////////////////////////////////////// //
            // ///// /**** INTER Bead Force ****/ ///// //
            // //////////////////////////////////////// //

            if(bead == middleOfChain){ // middle bead has LJ interactions

                for(next = bPc + chainSize; next < totalBeads; next += chainSize){ // only works for 1 bead in the middle
        
                    /**** Distance between 2 points ****/
                    dx = rx_temp[bPc] - rx_temp[next];
                    dy = ry_temp[bPc] - ry_temp[next];

                    /**** Nearest Image Convention ****/
                    nic(dx, dy, BoxLength, HalfLength);

                    d2 = dx*dx + dy*dy;

                    /**** Check if Beads are close enough to Interact ****/
                    if(d2 < rc2){
                        
                        /**** Calculate Distances but only if it's in range, otherwise we don't bother ****/
                        d = sqrt(d2);
                        d6 = d2*d2*d2;
                        d12 = d6*d6;
                        dMpeak = d - peak;

                        /**** Calculate Force and Potential ****/
                        f = (48*epsilon/d)*(s12/d12 - 0.5*s6/d6) + ((2*a*(dMpeak))/(b2))*exp(-(dMpeak*dMpeak)/(b2));

                        interFx_temp[bPc] += dx*f/d;
                        interFx_temp[next] -= dx*f/d;
                        interFy_temp[bPc] += dy*f/d;
                        interFy_temp[next] -= dy*f/d;

                    } // end distance if
                } // end for loop used to check next interacting beads 
            } // end if to see if the beads are of the interacting kind
        } // end for Bead
    } // end for Chain

}

void setInitialPositions(const double& BoxLength, const double& HalfLength, const int& middleOfChain,
                std::vector<double> &rx, std::vector<double> &ry){

    /**** Distances ****/
    double dx = 0, dy = 0;
    double dr = 0;
    double dr2 = 0;
    double minAprovedDist = 2; // Used in setting the initial positions

    /**** Random Stuff ****/
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uniBox(1, BoxLength - 1);
    std::uniform_real_distribution<double> uniSmall(-1, 1);

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
            for(bead = 0; bead < chainSize; bead++){

                bPc = bead + chain*chainSize;

                rx[bPc] = pointX + uniSmall(generator);
                ry[bPc] = pointY + uniSmall(generator);

            }
            chain++;
        }
            
        else{

            check = 1;

            for(j = 0; j <= chain*chainSize; j++){

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
            
                for(bead = 0; bead < chainSize; bead++){

                    bPc = bead + chain*chainSize;

                    rx[bPc] = pointX + uniSmall(generator);
                    ry[bPc] = pointY + uniSmall(generator);

                }

                chain++;

            }
        }
         
    }
}

/**** [FUNCTION] Make Video [FUNCTION] ****/
inline void Cannes(std::ofstream& video, int t, int attractingBeads, int totalBeads, int numOfPossibleBonds, std::vector<int> James, std::vector<double> X, std::vector<double> Y){

    std::vector<int> bond(attractingBeads);

    /**** Loop Variables ****/
    int n = 0, r = 0, var1 = 0, var2 = 0, var3 = 0, atom = 1;

    for(var1 = 0; var1 < attractingBeads-1; var1++){
        for(var2 = var1+1; var2 < attractingBeads; var2++){
            
            if(James[n] == 1){
                bond[var1]++;
                bond[var2]++;
                n++;
            }
            else{n++;}
            
        }
    }

    video << totalBeads << "\n" << "t = " << t << endl;

    for(r = 0; r < totalBeads; r++){

        if(r%chainSize == (chainSize-1)/2){
            atom = bond[var3]+2;
            var3++;
        }
        else{
            atom = 1;
        }

        video << atom << "\t" << X[r] << "\t"<< Y[r] << "\t" << 0.0 << endl;
    }
}


// //////////////////////////////////////////////////////////////////////////////// //
// /////////////////////////////// /**** MAIN ****/ /////////////////////////////// //
// //////////////////////////////////////////////////////////////////////////////// //


int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    const int totalBeads = chainNumber*chainSize;
    const int interactBeads = chainNumber;
    const int middleOfChain = chainSize/2;

    /**** Potential Constants ****/
    const long double rMin = 1.12246204830937298*sigma;
    const long double rCut = 2.5*rMin;
    const double rc2 = rCut*rCut;
    const double peak2 = peak*peak;
    const double s2 = sigma*sigma;
    const double s6 = s2*s2*s2;
    const double s12 = s6*s6;
    const double b2 = b*b;

    const double BoxLength = sqrt(interactBeads/density);
    cout << BoxLength << endl;
    const double HalfLength = BoxLength*0.5;
    /* this goes in main and then we send to the functions so that
    so that I don't have to calculate this over and over again */

    /**** Position ****/
    vector<double> rx(totalBeads); 
    vector<double> ry(totalBeads); 
    vector<double> rx_temp(totalBeads); 
    vector<double> ry_temp(totalBeads); 

    /**** Force ****/
    vector<double> intraFx(totalBeads);
    vector<double> intraFy(totalBeads);
    vector<double> intraFx_temp(totalBeads);
    vector<double> intraFy_temp(totalBeads);
    vector<double> interFx(totalBeads);
    vector<double> interFy(totalBeads);
    vector<double> interFx_temp(totalBeads);
    vector<double> interFy_temp(totalBeads);

    /**** Thermal Noise ****/
    vector<double> alx(totalBeads); 
    vector<double> aly(totalBeads);

    const int arraySizeSimul = Iterations/binSize; 

    const int numOfPossibleBonds = interactBeads*(interactBeads + 1)/2; // number of possible bonds (non repeated)

    /**** For Loop Variables ****/
    int pos = 0, t = 0, p = 0, spl = 0, InitCheck = 0;
    const double timeJump = dt*binSize;
    double tS = timeJump;

    /**** Others ****/
    int tt = 0, j = 0;

    /**** Video ****/
    const int recSample = 0; // Gives the Recorded Sample

    /**** Random stuff ****/
    double std = sqrt(g/dt); // usually it's "g/dt" but since dt will change I divide after
    std::mt19937 generator (seed);
    std::normal_distribution<double> normal(0, std);

    /**** Statistics ****/
    vector<int> bondGPSprev(numOfPossibleBonds);
    int formedBonds = 0;
    int brokenBonds = 0;
    int numberOfBonds = 0;
    vector<int> formedVector(arraySizeSimul);
    vector<int> brokenVector(arraySizeSimul);
    vector<int> numBondsVector(arraySizeSimul);

    /**** Calc constant shit before loops ****/
    const double dtOzeta = dt/zeta;
    const double samplesTbinSizeTparticleNumber = samples*binSize*interactBeads*0.5;
    const double truncLimit = 4*std*std*c2*c2;

    /**** Keep track of used parameters ****/
    cout << "\n\t*** Global Parameters *** \n \nSamples: " << samples << "\n" << "Number of Iterations: " <<
    Iterations << "\n" << "Time Step: " << dt << "\n\n" << "Bin Size: " << binSize << "\n\n" << "Number of Chains: " << chainNumber << "\n" 
    << "Number of Beads per Chain: " << chainSize << "\n\n\t*** Potential Parameters *** \n\nsigma = " << sigma << "\n" <<
    "epsilon = " << epsilon << "\n" << "a = " << a << "\n" << "b = " << b << "\n" << "Peak of Potential Wall at " << peak <<
    "\n\n\t*** Other Parameters *** \n\n" << "zeta = " << zeta << "\n" << "ks = " << ks << "\n" << "g = " << g << "\n\n" << endl;

    ofstream video;
        video.open("mChains.imd");
        video.precision(8);
        video.setf(ios::fixed | ios::showpoint);

    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        /**** Assign Initial Positions ****/
        setInitialPositions(BoxLength, HalfLength, middleOfChain, rx, ry);
        
        /**** Check Initial Positions... cause idk what I'm doing ****/ 
        // if(spl==0){

        //     ofstream file_2;
        //         file_2.open("out3.dat");

        //     for(j = 0; j < totalBeads; j++){ 
        //         file_2 << rx[j] << "\t" << ry[j] << endl;
        //     } 

        //     file_2.close();
        // }
        
        tt = 0;

        /**** Calculate Stuff for some Time ****/
        for(t = 0; t < Iterations; t++){

            /**** Evaluation of Force ***/
            evalForce(t, numOfPossibleBonds, BoxLength, HalfLength, middleOfChain, totalBeads, rx, ry, intraFx, intraFy,
                     interFx, interFy, rc2, peak2, s6, s12, b2, bondGPSprev, formedBonds, brokenBonds, numberOfBonds);

            // ////////////////// /**** Statistics ****/ ////////////////// //
            /*  Saves the number of formed and broken bonds as well as the number of bonds. 

                Need to divide the first 2 by the # of samples and elapsed time (binSize). 
                The #Bonds needs to be divided by the # of samples, elapsed time and # of beads

                This will be done in the end    */

            if((t + 1) % (binSize) == 0){

                formedVector[tt] += formedBonds; // keeps adding to then perform an avg
                brokenVector[tt] += brokenBonds;
                numBondsVector[tt] += numberOfBonds;

                formedBonds = 0;
                brokenBonds = 0;
                numberOfBonds = 0;
                tt++; // counter to run the array above
            }
            // ////////////////// //////////////////// ////////////////// //

            
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
                rx_temp[p] = rx[p] + dtOzeta*(alx[p] + intraFx[p] + interFx[p]);
                ry_temp[p] = ry[p] + dtOzeta*(aly[p] + intraFy[p] + interFy[p]);

                pbc(rx_temp[p], ry_temp[p], BoxLength);

            }

            /**** Evaluation of Force - Temp ***/
            evalForce_temp(BoxLength, HalfLength, middleOfChain, totalBeads, rx_temp, ry_temp, intraFx_temp, intraFy_temp,
                     interFx_temp, interFy_temp, rc2, peak2, s6, s12, b2);

            /**** Create Video of Particles Moving - Parameters ****/
            if(spl == recSample && t % binSize == 0){
                Cannes(video, t, interactBeads, totalBeads, numOfPossibleBonds, bondGPSprev, rx, ry);
            }

            /**** Calculate Position ****/
            for(p = 0; p < totalBeads; p++){

                rx[p] += dtOzeta*(alx[p] + 0.5*(intraFx[p] + interFx[p] + intraFx_temp[p] + interFx_temp[p]));
                ry[p] += dtOzeta*(aly[p] + 0.5*(intraFy[p] + interFy[p] + intraFy_temp[p] + interFy_temp[p]));

                pbc(rx[p], ry[p], BoxLength);

            }

        } // Closes iteration loop

    } // Closes samples loop

    // video.close();

    ofstream file;
        file.open("out1.dat");

    for(j = 0; j < arraySizeSimul; j++){ 
        // file << formedVector[j]/samplesTbinSize << "\t" << brokenVector[j]/samplesTbinSize << "\t" <<
        file << numBondsVector[j]/samplesTbinSizeTparticleNumber << "\t" << tS << endl;
        tS += timeJump; // time passed
    } 

    file.close();

    return 0;
}
