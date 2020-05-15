#include <iostream>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <random>
#include <tuple>
#include <vector>

#define samples 20
#define frames 5000
#define nParticles 60
#define simul 25000
#define dt 0.001 // tem que ser da ordem 10^-5 though.... rip
#define sigma 1
#define gamma 1
#define epsilon 1
#define m 1
#define c2 1.136847234338557L

using namespace std;

/**** Constants ****/
double rMin = 1.12246204830937298*sigma;
double rCut = 2.5*rMin;

double BoxLength = nParticles*0.6;
double HalfLength = BoxLength*0.5;

/**** Bin ****/
int binSize = 25;

double rc2 = rCut*rCut;
double rc6 = rc2*rc2*rc2;
double rc12 = rc6*rc6;
double s2 = sigma*sigma;
double s6 = sigma*sigma*sigma*sigma*sigma*sigma;
double s12 = s6*s6;
double Ucut = 4*epsilon*(s12/rc12 - s6/rc6);

/**** SRK2 Definitive Vars ****/
vector<double> rx(nParticles);
vector<double> ry(nParticles);
vector<double> fx(nParticles);
vector<double> fy(nParticles);

/**** SRK2 Temp Vars ****/
vector<double> rxTemp(nParticles);
vector<double> ryTemp(nParticles);
vector<double> fxTemp(nParticles);
vector<double> fyTemp(nParticles);

/**** Bond Variables ****/
int bpossible = nParticles*(nParticles + 1)/2;
vector<int> BondGPSprev(bpossible);
vector<double> brokenB(simul);
vector<double> formedB(simul);
vector<double> numBonds(simul);


/**** [FUNCTION] Calculate Temp Force [FUNCTION] ****/
void EvalForceTemp(int timeStamp){

    /**** For Loop Variables ****/
    int f1 = 0, f2 = 0, forInit = 0;

    /**** Distance and Force Variables ****/
    double dx = 0, dy = 0, d2 = 0, d6 = 0, d12 = 0, f = 0, u = 0;

    /**** Initialize Variables Again Because In Previous Codes Everything Went To Shit If I Didn't Do So ****/
    dx = 0;
    dy = 0;
    d2 = 0;
    d6 = 0;
    d12 = 0;
    f = 0;
    u = 0;

    /**** Set Force to 0 because we sum fx and fy ****/
    std::fill(fxTemp.begin(), fxTemp.end(), 0);
    std::fill(fyTemp.begin(), fyTemp.end(), 0);

    /**** Force Calculation ****/
    for(f1 = 0; f1 < nParticles-1; f1++){
        for(f2 = f1 + 1; f2 < nParticles; f2++){
        
            /**** Distance between 2 points ****/
            dx = rxTemp[f1] - rxTemp[f2];
            dy = ryTemp[f1] - ryTemp[f2];

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
            d6 = d2*d2*d2;
            d12 = d6*d6;

            /**** Check if Partcles Interact ****/
            if(d2 < rc2){

                /**** Calculate Force and Potential ****/
                f = (48*epsilon)*(s12/d12 - 0.5*s6/d6);

                fxTemp[f1] += dx*f/d2;
                fxTemp[f2] -= dx*f/d2;
                fyTemp[f1] += dy*f/d2;
                fyTemp[f2] -= dy*f/d2;
            }
        }
    } // Closes Force Loop

}

std::vector<int> EvalForceBond(vector<int> &BondGPSprev, int timeStamp){

    /**** For Loop Variables ****/
    int f1 = 0, f2 = 0, forInit = 0, nb = 0;

    /**** Distance and Force Variables ****/
    double dx = 0, dy = 0, d2 = 0, d6 = 0, d12 = 0, f = 0, u = 0;

    /**** This doesn't need to exit the function ****/
    vector<int> BondGPSnow(bpossible);

    /**** Set Force to 0 because we sum fx and fy ****/
    std::fill(fx.begin(), fx.end(), 0);
    std::fill(fy.begin(), fy.end(), 0);
    std::fill(BondGPSnow.begin(), BondGPSnow.end(), 0); // Just to make sure XD
    nb = 0;

    /**** Force and Bond Calculation ****/
    for(f1 = 0; f1 < nParticles-1; f1++){
        for(f2 = f1 + 1; f2 < nParticles; f2++){
        
            /**** Distance between 2 points ****/
            dx = rx[f1] - rx[f2];
            dy = ry[f1] - ry[f2];

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
            d6 = d2*d2*d2;
            d12 = d6*d6;

            /**** Check if Partcles Interact ****/
            if(d2 < rc2){

                /**** Calculate Force and Potential ****/
                f = (48*epsilon)*(s12/d12 - 0.5*s6/d6);

                fx[f1] += dx*f/d2;
                fx[f2] -= dx*f/d2;
                fy[f1] += dy*f/d2;
                fy[f2] -= dy*f/d2;

                /**** MI6 Statistics ****/
                BondGPSnow[nb] = 1; // if dx2 < rc2 then we have a bond between the particles x and y (facCount) at current time
                numBonds[timeStamp]++;

            }

            else{
                BondGPSnow[nb] = 0; // if dx2 > rc2 then we don't have a bond between the particles x and y (facCount) at current time
            }

            /**** Compare the current and previous times and see if a bond has formed or broken ****/
            if(timeStamp != 0){ // no point in comparing at time 0

                if(BondGPSprev[nb] < BondGPSnow[nb]){
                    formedB[timeStamp]++;
                }
                else if(BondGPSprev[nb] > BondGPSnow[nb]){
                    brokenB[timeStamp]++;
                }
            }

            BondGPSprev[nb] = BondGPSnow[nb];
            nb++;

        }
    } // Closes Force loop

    return BondGPSnow;

}

/**** [FUNCTION] Randomly Assign Positions [FUNCTION] ****/
void initP(){

    /**** Distances ****/
    double dx = 0, dy = 0;
    double dr = 0;
    double dr2 = 0;
    double min2 = 0.8; // rMin*rMin;

    /**** Random Stuff ****/
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    std::uniform_real_distribution<double> uni(1, BoxLength - 1);

    /**** Loop Variables ****/
    int i = 0, j = 0, ErrorCount = 0;

    /**** Points ****/
    double pointX = 0;
    double pointY = 0;

    /**** Others ****/
    int check = 1;

    ErrorCount = 0;
    check = 1;

    while(i < nParticles){

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
                if(dr2 < min2){
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

        ErrorCount++;

        if(ErrorCount == nParticles*5000){
            cout << "ERROR - INFINITE LOOP (probs XD)" << endl;
        }
    }

}

/**** [FUNCTION] Bins The Data [FUNCTION] ****/
void bin(int it){

    /**** Statistics ****/
    double bondTemp, formedTemp, brokenTemp, bondAvg, formedAvg, brokenAvg = 0;

    /**** Time ****/
    double tSimul = 0, avgTime = 0;

    /**** Others ****/
    int i, count = 0;
    char filename_007[100];

    sprintf(filename_007, "/home/cerdeira/Documents/Data/Code Output/SRK2/LJ2D/BOND_%d.txt", it);
    ofstream james;
        james.open(filename_007);

    /**** Initialize Variables.... for some reason idk, doesn't work otherwise... ****/
    tSimul = 0;
    avgTime = 0;
    bondTemp = 0;
    formedTemp = 0;
    brokenTemp = 0;
    count = 0;

    for(i = 0; i < simul; i++){
        
        if(i % binSize != 0 || i == 0){
            bondTemp += numBonds[i];
            formedTemp += formedB[i];
            brokenTemp += brokenB[i];
            count++;
        }
        else{

            bondAvg = bondTemp/(samples*count*nParticles);
            formedAvg = formedTemp/(samples*count*nParticles);
            brokenAvg = brokenTemp/(samples*count*nParticles);

            avgTime = tSimul - (count/2)*dt;
            james << bondAvg << "\t" << formedAvg << "\t" << brokenAvg << "\t" << avgTime << "\t" << count << endl;

            bondTemp = numBonds[i];
            formedTemp = formedB[i];
            brokenTemp = brokenB[i];
            count = 1;
 
        }

        tSimul += dt; 

    }

    james.close();

}

int main(int argc, char const *argv[]){

    /**** For loops variables ****/
    int sample, errorFor, InitCheck, t, forInit, fp1, fp2, p = 0;

    /**** Progress Bar Variables ****/
    double partNum = 1;
    double partDen = 6;
    int samplesD = samples/partDen;
    int part = samplesD;

    /**** Others ****/
    int i, count, ifExpression = 0;
    int c = 1;

    /**** Video ****/
    string atom = "C";
    int VideoRep = samples/2; // Gives the Recorded Sample
    int VideoTimeJump = simul/frames; 

    /**** Time ****/
    double TSimul = 0;
    double dtTrak = 0;

    /**** Bonds - Error Check ****/
    double bondError = 0;
    double formedError = 0;
    double brokenError = 0;
    
    /**** Random stuff ****/
    double alx[nParticles] = {};
    double aly[nParticles] = {};
    double gStart = 0.4;
    double g = gStart; // g = 2*gamma*kb*T
    auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);

    /**** File Names ****/
    char filename_M[100];
    char filename_Vid[100];

    /**** Open file to keep track of used parameters ****/
    ofstream par;
        par.open("/home/cerdeira/Documents/Data/Code Output/SRK2/LJ2D/parametersLJ2D.txt");

        par << "    *** Global Parameters *** \n \nSamples: " << samples << "\n" << "Particles: " << nParticles << "\n" 
            << "Box Length: " << BoxLength << "\n" << "Number of Iterations: " << simul << "\n" << "Time Step: " << dt << "\n\n" 
            << "sigma = " << sigma << "\n" << "epsilon = " << epsilon << "\n" << "gamma = " << gamma 
            << "\n\n    *** Variable Parameters *** \n" <<  endl;

    while (g < 1.7){

        double std = sqrt(g/dt);

        std::uniform_real_distribution<double> uni(1, BoxLength - 1);
        std::normal_distribution<double> normal(0, std);

        /**** Initialize variables ****/
        dtTrak = 0;
        part = samplesD;
        partNum = 1;

        for(i = 0; i < simul; i++){ 
            numBonds[i] = 0;
            formedB[i] = 0;
            brokenB[i] = 0;
        }

        /**** Open Video File ****/
        sprintf(filename_Vid, "/home/cerdeira/Documents/Data/Code Output/SRK2/LJ2D/LJ2D_VIDEO_%d.xyz", c);
        ofstream video;
            video.open(filename_Vid);
            video.precision(8);
            video.setf(ios::fixed | ios::showpoint);

        /**** Repeat stuff ****/
        for(sample = 0; sample < samples; sample++){

            /**** Assign Initial Positions ****/
            initP();

            /**** Check if Initialized Shit ain't Shit ****/
            if(sample == 0 && g == gStart){

                ofstream Position;
                    Position.open("/home/cerdeira/Documents/Data/Code Output/SRK2/LJ2D/INIT_Check.txt");

                for(InitCheck = 0; InitCheck < nParticles; InitCheck++){
                    Position << rx[InitCheck] << "\t" << ry[InitCheck] << endl;
                }

                Position.close();
                cout << "File LJ1D_INIT_Check.txt Created" << "\n" << endl;
            }
            
            /**** Calculate Stuff for some Time ****/
            for(t = 0; t < simul; t++){

                /**** Calculate Temporary Position ****/
                for(p = 0; p < nParticles; p++){

                    /**** Set Random Motion ****/
                    do{ 
                        alx[p] = normal(generator)*c2;
                    }
                    while(alx[p]*alx[p] > 4*std*c2*c2);
                    do{
                        aly[p] = normal(generator)*c2;
                    }
                    while(aly[p]*aly[p] > 4*std*c2*c2);

                    /**** Temporary Position ****/
                    rxTemp[p] = rx[p] + (dt/gamma)*(alx[p] + fx[p]);
                    ryTemp[p] = ry[p] + (dt/gamma)*(aly[p] + fy[p]);

                    /**** PBC ****/
                    if(rxTemp[p] < 0){
                        rxTemp[p] += BoxLength;
                    }
                    else if(rxTemp[p] > BoxLength){
                        rxTemp[p] -= BoxLength;
                    }
                    if(ryTemp[p] < 0){
                        ryTemp[p] += BoxLength;
                    }
                    else if(ryTemp[p] > BoxLength){
                        ryTemp[p] -= BoxLength;
                    }
                }

                /**** Evaluation of Force - Temp ***/
                EvalForceTemp(t);

                /**** Create Video of Particles Moving - Parameters ****/
                if(sample == VideoRep && t % VideoTimeJump == 0){
                    video << nParticles << "\n" << "t = " << t << endl;
                }

                /**** Calculate Position ****/
                for(p = 0; p < nParticles; p++){

                    /**** Create Video of Particles Moving - Position ****/
                    if(sample == VideoRep && t % VideoTimeJump == 0){
                        video << atom << "\t" << rx[p] << "\t"<< ry[p] << "\t" << 0.0 << endl;
                    }

                    rx[p] += (dt/gamma)*(alx[p] + 0.5*(fx[p] + fxTemp[p]));
                    ry[p] += (dt/gamma)*(aly[p] + 0.5*(fy[p] + fyTemp[p]));

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

                /**** Evaluation of Force and Bond Statistics ****/
                // std::vector<int> isto antes do "BondGPSprev =" ou não ???
                BondGPSprev = EvalForceBond(BondGPSprev, t);

            } // Closes simul loop

            /**** Progress Bar (sort of...) ****/
            if(sample == part){
                cout << (partNum/partDen)*100 << "% done for g = " << g << endl;
                partNum++;
                part += samplesD;
            }
            
        } // Closes repeat loop
        
        /**** Bins the Values of the Bonds Staistics ****/
        bin(c);

        sprintf(filename_M, "/home/cerdeira/Documents/Data/Code Output/SRK2/LJ2D/Error_%d.txt", c);
        ofstream error;
            error.open(filename_M);

        for(errorFor = 0; errorFor < 1500; errorFor++){

            bondError = numBonds[errorFor]/samples;
            formedError = formedB[errorFor]/samples;
            brokenError = brokenB[errorFor]/samples;

            error << bondError << "\t" << formedError << "\t" << brokenError << "\t" << errorFor << endl;

        }

        error.close();
        video.close();

        cout << "Files Created" << "\t" << "g = " << g << endl;
        par << "File #" << c << "\t" << "g = " << g << endl;
        cout << "\n" << endl;

        g += 0.6;
        c++;

    }

    par.close();
    
    return 0;

}
