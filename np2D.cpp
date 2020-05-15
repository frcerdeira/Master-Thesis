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
int samples = 0, Iterations = 0, binSize = 0;
int frames = 0;
int particleNumber = 0;
double density = 0;
double dt = 0; // Step size
double zeta = 0, g = 0; // g = 2*zeta*kb*T // Dynamics
// //// Potential Parameters //// //
double sigma = 0, epsilon = 0, peak = 0, a = 0, b = 0;


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
    particleNumber = arguments[3];
    density = arguments[4];
    dt = arguments[5];
    // //// Statistics //// //
    binSize = arguments[6];
    // //// Parameters of Interactions //// //
    zeta = arguments[7];
    g = arguments[8];
    // //// Parameters of the Potential //// //
    sigma = arguments[9];
    epsilon = arguments[10];
    peak = arguments[11];
    a = arguments[12];
    b = arguments[13];
    frames = arguments[14];

}

/**** [FUNCTION] Calculates Forces [FUNCTION] ****/
inline void evalForce_temp(const double& BoxLength, const double& HalfLength, 
                            vector<double>& rx_temp, vector<double>& ry_temp, vector<double>& fx_temp, vector<double>& fy_temp, 
                            const double& rc2, const double& peak2, const double& s6, const double& s12, const double& b2){

    /**** For Loop Variables ****/
    int p1 = 0, p2 = 0, nb = 0;

    /**** Distance and Force Variables ****/
    double dx = 0, dy = 0, d = 0, d2 = 0, d6 = 0, d12 = 0, f = 0;

    double dMpeak = 0;

    /**** Initialize Variables Again Because In Previous Codes Everything Went To Shit If I Didn't Do So ****/
    dx = 0;
    dy = 0;
    d = 0;
    d2 = 0;
    d6 = 0;
    d12 = 0;
    f = 0;

    /**** Set Force to 0 because we sum fx and fy ****/
    std::fill(fx_temp.begin(), fx_temp.end(), 0);
    std::fill(fy_temp.begin(), fy_temp.end(), 0);

    nb = 0;

    /**** Force and Bond Calculation ****/
    for(p1 = 0; p1 < particleNumber-1; p1++){
        for(p2 = p1 + 1; p2 < particleNumber; p2++){
        
            /**** Distance between 2 points ****/
            dx = rx_temp[p1] - rx_temp[p2];
            dy = ry_temp[p1] - ry_temp[p2];

            /**** PBC ****/
            if(dx > HalfLength){
                dx -= BoxLength;
            }
            else if(dx < -HalfLength){
                dx += BoxLength;
            }
            if(dy > HalfLength){
                dy -= BoxLength;
            }
            else if(dy < -HalfLength){
                dy += BoxLength;
            }

            d2 = dx*dx + dy*dy;
            d = sqrt(d2);
            d6 = d2*d2*d2;
            d12 = d6*d6;

            dMpeak = d - peak;

            /**** Check if Partcles Interact ****/
            if(d2 < rc2){

                /**** Calculate Force and Potential ****/
                f = (48*epsilon/d)*(s12/d12 - 0.5*s6/d6) + ((2*a*(dMpeak))/(b*b))*exp(-(dMpeak*dMpeak)/(b*b));

                fx_temp[p1] += dx*f/d;
                fx_temp[p2] -= dx*f/d;
                fy_temp[p1] += dy*f/d;
                fy_temp[p2] -= dy*f/d;

            }
        }
    } // Closes Force loop

}

/**** [FUNCTION] Calculates Forces and Bond Statistics [FUNCTION] ****/
inline void evalForce(const int& t, const int& numOfPossibleBonds, const double& BoxLength, const double& HalfLength, 
                            vector<double>& rx, vector<double>& ry, vector<double>& fx, vector<double>& fy, 
                            const double& rc2, const double& peak2, const double& s6, const double& s12, const double& b2, 
                            vector<int>& BondGPSprev, int& formedBonds, int& brokenBonds, int& numberOfBonds){

    /**** For Loop Variables ****/
    int p1 = 0, p2 = 0, nb = 0;

    /**** Distance and Force Variables ****/
    double dx = 0, dy = 0, d = 0, d2 = 0, d6 = 0, d12 = 0, f = 0;

    double dMpeak = 0;

    /**** Initialize Variables Again Because In Previous Codes Everything Went To Shit If I Didn't Do So ****/
    dx = 0;
    dy = 0;
    d = 0;
    d2 = 0;
    d6 = 0;
    d12 = 0;
    f = 0;

    /**** This doesn't need to exit the function ****/
    vector<int> BondGPSnow(numOfPossibleBonds);

    /**** Set Force to 0 because we sum fx and fy ****/
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(BondGPSnow.begin(), BondGPSnow.end(), 0); // Just to make sure XD
    nb = 0;

    /**** Force and Bond Calculation ****/
    for(p1 = 0; p1 < particleNumber-1; p1++){
        for(p2 = p1 + 1; p2 < particleNumber; p2++){
        
            /**** Distance between 2 points ****/
            dx = rx[p1] - rx[p2];
            dy = ry[p1] - ry[p2];

            /**** PBC ****/
            if(dx > HalfLength){
                dx -= BoxLength;
            }
            else if(dx < -HalfLength){
                dx += BoxLength;
            }
            if(dy > HalfLength){
                dy -= BoxLength;
            }
            else if(dy < -HalfLength){
                dy += BoxLength;
            }

            d2 = dx*dx + dy*dy;
            d = sqrt(d2);
            d6 = d2*d2*d2;
            d12 = d6*d6;

            dMpeak = d - peak;

            /**** Check if Partcles Interact ****/
            if(d2 < rc2){

                /**** Calculate Force and Potential ****/
                f = (48*epsilon/d)*(s12/d12 - 0.5*s6/d6) + ((2*a*(dMpeak))/b2)*exp(-(dMpeak*dMpeak)/b2);

                fx[p1] += dx*f/d;
                fx[p2] -= dx*f/d;
                fy[p1] += dy*f/d;
                fy[p2] -= dy*f/d;

                /**** MI6 Statistics ****/
                if(d2 < peak2){
                    BondGPSnow[nb] = 1;
                    numberOfBonds++;
                }
                else{
                    BondGPSnow[nb] = 0;
                }

            }
            else{
                BondGPSnow[nb] = 0;
            }

            /**** Compare the current and previous times and see if a bond has formed or broken ****/
            if(t != 0){ // no point in comparing at time 0

                if(BondGPSprev[nb] < BondGPSnow[nb]){
                    formedBonds++;
                }
                else if(BondGPSprev[nb] > BondGPSnow[nb]){
                    brokenBonds++;
                }
            }

            BondGPSprev[nb] = BondGPSnow[nb];
            nb++;

        }
    } // Closes Force loop

}

/**** [FUNCTION] Randomly Assign Positions [FUNCTION] ****/
inline void setInitialPositions(const int& seed, const double& BoxLength, const double& HalfLength, 
            vector<double>& rx, vector<double>& ry){

    /**** Distances ****/
    double dx = 0, dy = 0;
    double dr = 0;
    double dr2 = 0;
    double minAprovedDist = 0.7; // Used in setting the initial positions

    /**** Random Stuff ****/
    std::mt19937 generator(seed);
    std::uniform_real_distribution<double> uni(1, BoxLength - 1);

    /**** Loop Variables ****/
    int i = 0, j = 0, ErrorCount = 0;

    /**** Points ****/
    double pointX = 0;
    double pointY = 0;

    /**** Others ****/
    int check = 1;

    check = 1;

    while(i < particleNumber){

        /**** Assign a Point ****/
        pointX = uni(generator);
        pointY = uni(generator);
        
        if(i == 0){
            rx[i] = pointX;
            ry[i] = pointY;
            i++;
        }

        else{
            
            check = 1;

            for(j = 1; j <= i; j++){

                /**** Calculate Distance in 1D to Point ****/
                dx = pointX - rx[j];
                dy = pointY - ry[j];

                /**** PBC ****/
                if(dx > HalfLength){
                    dx -= BoxLength;
                }
                else if(dx < -HalfLength){
                    dx += BoxLength;
                }
                if(dy > HalfLength){
                    dy -= BoxLength;
                }
                else if(dy < -HalfLength){
                    dy += BoxLength;
                }

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
                rx[i] = pointX;
                ry[i] = pointY;
                i++;
            }
        }
    }
}


int main(int argc, char const *argv[]){

    /**** Read File with Parameters ****/
    readFile();

    /**** Potential Constants ****/
    const long double rMin = 1.12246204830937298*sigma;
    const long double rCut = 2.5*rMin;
    const double rc2 = rCut*rCut;
    const double peak2 = peak*peak;
    const double s2 = sigma*sigma;
    const double s6 = s2*s2*s2;
    const double s12 = s6*s6;
    const double b2 = b*b;

    const double BoxLength = sqrt(particleNumber/density);
    cout << BoxLength << endl;
    const double HalfLength = BoxLength*0.5;
    /* this goes in main and then we send to the functions so that
    so that I don't have to calculate this over and over again */


    /**** SRK2 Vars ****/
    vector<double> rx(particleNumber); // position x coordinate
    vector<double> ry(particleNumber); // position y coordinate
    vector<double> fx(particleNumber); // force x coordinate
    vector<double> fy(particleNumber); // force y coordinate
    vector<double> rx_temp(particleNumber); // temporary position x coordinate
    vector<double> ry_temp(particleNumber); // temporary position y coordinate
    vector<double> fx_temp(particleNumber); // temporary force x coordinate
    vector<double> fy_temp(particleNumber); // temporary force y coordinate
    vector<double> alx(particleNumber); // random portion x coordinate
    vector<double> aly(particleNumber); // random portion y coordinate

    const int arraySizeSimul = Iterations/binSize; /* gives the # of entries on the files that plot the
                                           e2e distance in respect to time */

    const int numOfPossibleBonds = particleNumber*(particleNumber + 1)/2; // number of possible bonds (non repeated)

    /**** For Loop Variables ****/
    int t = 0, p = 0, spl = 0, tt = 0, j = 0;
    const double timeJump = dt*binSize;
    double tS = timeJump; // we ignore the first iterations in order to reach a state near equilibrium (what matters)

    /**** Video ****/
    string atom = "C";
    const int recSample = 1; // Gives the Recorded Sample
    const int VideoTimeJump = Iterations/frames; 

    /**** Statistics ****/
    vector<int> BondGPSprev(numOfPossibleBonds);
    int formedBonds = 0;
    int brokenBonds = 0;
    int numberOfBonds = 0;
    vector<int> formedVector(arraySizeSimul);
    vector<int> brokenVector(arraySizeSimul);
    vector<int> numBondsVector(arraySizeSimul);


    /**** Random stuff ****/
    double std = sqrt(g/dt);
    std::mt19937 generator (seed);
    std::normal_distribution<double> normal(0, std);


    /**** Keep track of used parameters ****/
    cout << "\n\t*** Global Parameters *** \n \nSamples: " << samples << "\n" << "Number of Iterations: " <<
    Iterations << "\n" << "Time Step: " << dt << "\n\n" << "Bin Size: " << binSize << "\n\n" << "Number of Particles Used: " 
    << particleNumber << "\n\n\t*** Potential Parameters *** \n\nsigma = " << sigma << "\n" << "epsilon = " << epsilon << "\n" 
    << "a = " << a << "\n" << "b = " << b << "\n" << "Peak of Potential Wall at " << peak << "\n\n\t*** Other Parameters *** \n\n" <<
    "zeta = " << zeta << "\n" << "g = " << g << "\n" << endl;

    /**** Calc constant shit before loops ****/
    const double dtOzeta = dt/zeta;
    const double samplesTbinSize = samples*binSize;
    const double samplesTbinSizeTparticleNumber = samples*binSize*particleNumber*0.5;
    const double truncLimit = 4*std*std*c2*c2;

    /**** Open Video File ****/
    ofstream video;
        video.open("np2D.xyz");
        video.precision(8);
        video.setf(ios::fixed | ios::showpoint);

    /**** Repeat stuff ****/
    for(spl = 0; spl < samples; spl++){

        /**** Assign Initial Positions ****/
        setInitialPositions(seed, BoxLength, HalfLength, rx, ry);

        formedBonds = 0;
        brokenBonds = 0;
        numberOfBonds = 0;
        tt = 0;

        /**** Calculate Stuff for some Time ****/
        for(t = 0; t < Iterations; t++){

            /**** Evaluation of Force ***/
            evalForce(t, numOfPossibleBonds, BoxLength, HalfLength, rx, ry, fx, fy, rc2, 
                                    peak2, s6, s12, b2, BondGPSprev, formedBonds, 
                                    brokenBonds, numberOfBonds);

                                        /**** Statistics ****/
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
            
            /**** Calculate Temporary Position ****/
            for(p = 0; p < particleNumber; p++){

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
                rx_temp[p] = rx[p] + dtOzeta*(alx[p] + fx[p]);
                ry_temp[p] = ry[p] + dtOzeta*(aly[p] + fy[p]);

                /**** PBC ****/
                if(rx_temp[p] < 0){
                    rx_temp[p] += BoxLength;
                }
                else if(rx_temp[p] > BoxLength){
                    rx_temp[p] -= BoxLength;
                }
                if(ry_temp[p] < 0){
                    ry_temp[p] += BoxLength;
                }
                else if(ry_temp[p] > BoxLength){
                    ry_temp[p] -= BoxLength;
                }

            }

            /**** Evaluation of Force - Temp ***/
            evalForce_temp(BoxLength, HalfLength, rx_temp, ry_temp, fx_temp, fy_temp, rc2,
                            peak2, s6, s12, b2);

            /**** Create Video of Particles Moving - Parameters ****/
            if(spl == recSample && t % VideoTimeJump == 0){
                video << particleNumber << "\n" << "t = " << t << endl;
            }

            /**** Calculate Position ****/
            for(p = 0; p < particleNumber; p++){

                /**** Create Video of Particles Moving - Position ****/
                if(spl == recSample && t % VideoTimeJump == 0){
                    video << atom << "\t" << rx[p] << "\t"<< ry[p] << "\t" << 0.0 << endl;
                }

                rx[p] += dtOzeta*(alx[p] + 0.5*(fx[p] + fx_temp[p])); 
                ry[p] += dtOzeta*(aly[p] + 0.5*(fy[p] + fy_temp[p]));

                /**** PBC ****/
                if(rx[p] < 0){
                    rx[p] += BoxLength;
                }
                else if(rx[p] > BoxLength){
                    rx[p] -= BoxLength;
                }
                if(ry[p] < 0){
                    ry[p] += BoxLength;
                }
                else if(ry[p] > BoxLength){
                    ry[p] -= BoxLength;
                }

            }

        } // Closes iteration loop

    } // Closes samples loop

    video.close();

    ofstream file;
        file.open("out1.dat");

    for(j = 0; j < arraySizeSimul; j++){ 

        file << formedVector[j]/samplesTbinSize << "\t" << brokenVector[j]/samplesTbinSize << "\t" << numBondsVector[j]/samplesTbinSizeTparticleNumber << "\t" << tS << endl;
        tS += timeJump; // time passed

    } 

    file.close();


    return 0;
}