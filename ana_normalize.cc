#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

void ana_normalize()
{
    //const char* input_flux = "flux/flux_mark.txt";
    const char* input_flux = "flux/flux_mark_for_integral.txt";
    std::ifstream infile;
    infile.open(input_flux);
    std::string line;
    Double_t integral;
    Double_t ene0=0.;
    Double_t prob0=0.;
    Double_t ene1,prob1;
    while (std::getline(infile, line))
    {
        std::stringstream ss(line);
        ss >> ene1 >> prob1;
	integral = (ene1-ene0) * (prob1+prob0)/2.;
	ene0=ene1;
	prob0=prob1;
    }
    std::cout<<"Integral = "<<integral<<std::endl;
    Double_t POT_per_second = 2.7 * TMath::Power(10,15);
    Double_t Scale = integral * POT_per_second; 
    std::cout<<"Integral = "<<Scale<<std::endl;
}
