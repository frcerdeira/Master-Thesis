#include <iostream>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <fstream>
#include <bits/stdc++.h>
#include <random>
#include <chrono>

#define M_PIl 3.141592653589793238462643383279502884L /* pi */
#define st 0.198748043098799L
#define str 0.291125094772793L
#define c1 1.853361676622186L
#define c2 1.136847234338557L
#define c3 1.01360419764221L

#define samples 500000
#define dt 0.001
#define sigma 1
#define gamma 1
#define g 2.0 // sigma e não a var
#define avg 0 // media da dist
#define epsilon 1
#define m 1
#define dt 0.001

using namespace std; 

int main(int argc, char const *argv[])
{
    
    long double al = 0, alTr = 0;
    int i = 0;
    int cut = g; // factor by which we cut the dist --> can be 2*std, 4*std, etc where cut is the number

    std::mt19937 generator (12345); // constant number instead of seed to get repeated results
    std::normal_distribution<long double> distribution(avg, sqrt(g/dt));

    double root = sqrt(str);

    // double *lista;
    // lista = new double [samples];

    ofstream normal;
	    normal.open("normalRNG.txt");
    ofstream tru;
	    tru.open("truRNG.txt");
    // ofstream myway;
	//     myway.open("mywayRNG.txt");

    for(i = 0; i < samples; i++){

        /**** Generation of Numbers ****/
        al = distribution(generator); 
        alTr = distribution(generator)*c2;

        normal << al << endl;

        if(alTr*alTr < 4*(g/dt)*c2*c2){
            tru << alTr << endl;
        }
        // if(al*al < g*g){
        //     myway << al/root << endl;
        // }


    }

    normal.close();
    tru.close();
    // myway.close();
 
    // Hist(lista);

    // delete [] lista;

    return 0;
}