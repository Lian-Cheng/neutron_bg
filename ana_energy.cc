#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#define ConvertToCSV 0
#define DRAW_ENE_BINS 0
#define XE 1
#define Cluster 1

void ana_energy()
{
    const int NRGBs = 5;
    const int NCont = 255;
    
    // Define the color stops (position on the color scale from 0.0 to 1.0).
    Double_t stops[NRGBs] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
    
    // Define the colors for each stop in RGB (Red, Green, Blue from 0.0 to 1.0).
    // The first color (at stop 0.00) is set to white (1.0, 1.0, 1.0).
    Double_t red[NRGBs]   = { 1.00, 0.00, 0.00, 0.50, 1.00 };
    Double_t green[NRGBs] = { 1.00, 0.00, 0.50, 1.00, 1.00 };
    Double_t blue[NRGBs]  = { 1.00, 0.50, 1.00, 1.00, 0.00 };

    // 4. Create the new palette and set it for the current ROOT style.
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    gStyle->SetOptStat(0);

    // Areas for the neutron sources
    Double_t dist1 = 82. + 1. * 2. + 15. * 2.;
    Double_t dist2 = 82. + 1. * 2. + 15. * 2.;
    Double_t dist3 = 111. + 1. * 2. + 15. * 2.;
    Double_t thick_5ch2 = 5.;
    Double_t thick_10ch2 = 10.;
    Double_t thick_15ch2 = 15.;
    Double_t thick_20ch2 = 20.;
    Double_t thick_25ch2 = 25.;
    Double_t thick_30ch2 = 30.;
    Double_t area_noShielding = 2. * TMath::Pi() * 82./4. + TMath::Pi() * 82. * (111. + 10. * 2.);
    Double_t area_5ch2 = 2 * ((dist1 + thick_5ch2) * (dist2 + thick_5ch2) + (dist2 + thick_5ch2) * (dist3 + thick_5ch2) + (dist3 + thick_5ch2) * (dist1 + thick_5ch2));
    Double_t area_10ch2 = 2 * ((dist1 + thick_10ch2) * (dist2 + thick_10ch2) + (dist2 + thick_10ch2) * (dist3 + thick_10ch2) + (dist3 + thick_10ch2) * (dist1 + thick_10ch2));
    Double_t area_15ch2 = 2 * ((dist1 + thick_15ch2) * (dist2 + thick_15ch2) + (dist2 + thick_15ch2) * (dist3 + thick_15ch2) + (dist3 + thick_15ch2) * (dist1 + thick_15ch2));
    Double_t area_20ch2 = 2 * ((dist1 + thick_20ch2) * (dist2 + thick_20ch2) + (dist2 + thick_20ch2) * (dist3 + thick_20ch2) + (dist3 + thick_20ch2) * (dist1 + thick_20ch2));
    Double_t area_25ch2 = 2 * ((dist1 + thick_25ch2) * (dist2 + thick_25ch2) + (dist2 + thick_25ch2) * (dist3 + thick_25ch2) + (dist3 + thick_25ch2) * (dist1 + thick_25ch2));
    Double_t area_30ch2 = 2 * ((dist1 + thick_30ch2) * (dist2 + thick_30ch2) + (dist2 + thick_30ch2) * (dist3 + thick_30ch2) + (dist3 + thick_30ch2) * (dist1 + thick_30ch2));
    Double_t norm_neutron = 1.09 * TMath::Power(10,-7) * 1.3 / 0.8 * 3. * 3600. * 3600.; //normalized neutrons/cm2/year

    //Reading result files
    std::vector<std::string> fileNames_noShielding;
    std::vector<std::string> fileNames_5ch2;
    std::vector<std::string> fileNames_10ch2;
    std::vector<std::string> fileNames_15ch2;
    std::vector<std::string> fileNames_20ch2;
    std::vector<std::string> fileNames_25ch2;
    std::vector<std::string> fileNames_30ch2;
    #if XE
    Double_t ene_bins[14] = {0.9,2.94451766,4.3987765,6.14103535,8.17129421,10.48955308,13.09581196,15.99007085,19.17232974,22.64258864,26.40084755,30.44710647,34.7813654,60.77266017};
        #if Cluster
    	    const std::filesystem::path dir_noShielding="results/results_xe_noShielding/";
    	    const std::filesystem::path dir_10ch2="results/results_xe_10ch2/";
    	    const std::filesystem::path dir_15ch2="results/results_xe_15ch2/";
    	    const std::filesystem::path dir_20ch2="results/results_xe_20ch2/";
    	    const std::filesystem::path dir_25ch2="results/results_xe_25ch2/";
    	    const std::filesystem::path dir_30ch2="results/results_xe_30ch2/";
        #else
    	    fileNames_noShielding = {"results/output_xe_noshielding1.txt","results/output_xe_noshielding2.txt"};
    	    fileNames_10ch2 = {"results/output_xe_10ch2.txt"};
    	    fileNames_15ch2 = {"results/output_xe_15ch2.txt"};
    	    fileNames_20ch2 = {"results/output_xe_20ch2.txt"};
    	    fileNames_25ch2 = {"results/output_xe_25ch2.txt"};
    	    fileNames_30ch2 = {"results/output_xe_30ch2.txt"};
        #endif
    #else
    Double_t ene_bins[18] = {0,8.17129421,10.48955308,13.09581196,15.99007085,19.17232974,22.64258864,26.40084755,30.44710647,34.7813654,39.40362434,44.31388329,49.51214224,54.9984012,60.77266017,73.18517814,86.74969615,143.29750943};
        #if Cluster
    	    const std::filesystem::path dir_noShielding="results/results_ar_noShielding/";
    	    const std::filesystem::path dir_10ch2="results/results_ar_10ch2/";
    	    const std::filesystem::path dir_15ch2="results/results_ar_15ch2/";
    	    const std::filesystem::path dir_20ch2="results/results_ar_20ch2/";
    	    const std::filesystem::path dir_25ch2="results/results_ar_25ch2/";
    	    const std::filesystem::path dir_30ch2="results/results_ar_30ch2/";
        #else
    	    fileNames_noShielding = {"results/output_ar_noshielding.txt"};
    	    fileNames_10ch2 = {"results/output_ar_10ch2.txt"};
    	    fileNames_15ch2 = {"results/output_ar_15ch2.txt"};
    	    fileNames_20ch2 = {"results/output_ar_20ch2.txt"};
    	    fileNames_25ch2 = {"results/output_ar_25ch2.txt"};
    	    fileNames_30ch2 = {"results/output_ar_30ch2_1.txt","results/output_ar_30ch2_2.txt"};
        #endif
    #endif
    #if Cluster
        for (auto const& entry : std::filesystem::directory_iterator{dir_noShielding}){
            if(entry.is_regular_file()){fileNames_noShielding.push_back(entry.path());}
        }
        for (auto const& entry : std::filesystem::directory_iterator{dir_10ch2}){
            if(entry.is_regular_file()){fileNames_10ch2.push_back(entry.path());}
        }
        for (auto const& entry : std::filesystem::directory_iterator{dir_15ch2}){
            if(entry.is_regular_file()){fileNames_15ch2.push_back(entry.path());}
        }
        for (auto const& entry : std::filesystem::directory_iterator{dir_20ch2}){
            if(entry.is_regular_file()){fileNames_20ch2.push_back(entry.path());}
        }
        for (auto const& entry : std::filesystem::directory_iterator{dir_25ch2}){
            if(entry.is_regular_file()){fileNames_25ch2.push_back(entry.path());}
        }
        for (auto const& entry : std::filesystem::directory_iterator{dir_30ch2}){
            if(entry.is_regular_file()){fileNames_30ch2.push_back(entry.path());}
        }
    #endif
    const char* output_filename = "energy_result.root";
    Int_t event_index;
    Float_t ene1,ene2,ene_event;
    std::string line;
    #if ConvertToCSV
    	#if XE
    	Double_t counter_noShielding[13]={0};
    	//Double_t counter_5ch2[13]={0};
    	Double_t counter_10ch2[13]={0};
    	Double_t counter_15ch2[13]={0};
    	Double_t counter_20ch2[13]={0};
    	Double_t counter_25ch2[13]={0};
    	Double_t counter_30ch2[13]={0};
    	#else
    	Double_t counter_noShielding[17]={0};
    	//Double_t counter_5ch2[17]={0};
    	Double_t counter_10ch2[17]={0};
    	Double_t counter_15ch2[17]={0};
    	Double_t counter_20ch2[17]={0};
    	Double_t counter_25ch2[17]={0};
    	Double_t counter_30ch2[17]={0};
	#endif
    #elif DRAW_ENE_BINS
	#if XE
    	TH1F* th1 = new TH1F("th1","Counts vs deposited energy;Energy [keV];Counts",13,ene_bins);
    	TH1F* th2 = new TH1F("th2","Counts vs deposited energy;Energy [keV];Counts",13,ene_bins);
    	TH1F* th3 = new TH1F("th3","Counts vs deposited energy;Energy [keV];Counts",13,ene_bins);
    	TH1F* th4 = new TH1F("th4","Counts vs deposited energy;Energy [keV];Counts",13,ene_bins);
    	TH1F* th5 = new TH1F("th5","Counts vs deposited energy;Energy [keV];Counts",13,ene_bins);
    	TH1F* th6 = new TH1F("th6","Counts vs deposited energy;Energy [keV];Counts",13,ene_bins);
	#else
    	TH1F* th1 = new TH1F("th1","Counts vs deposited energy;Energy [keV];Counts",17,ene_bins);
    	TH1F* th2 = new TH1F("th2","Counts vs deposited energy;Energy [keV];Counts",17,ene_bins);
    	TH1F* th3 = new TH1F("th3","Counts vs deposited energy;Energy [keV];Counts",17,ene_bins);
    	TH1F* th4 = new TH1F("th4","Counts vs deposited energy;Energy [keV];Counts",17,ene_bins);
    	TH1F* th5 = new TH1F("th5","Counts vs deposited energy;Energy [keV];Counts",17,ene_bins);
    	TH1F* th6 = new TH1F("th6","Counts vs deposited energy;Energy [keV];Counts",17,ene_bins);
	#endif
    #else
    TH1F* th1=new TH1F("th1","Counts vs deposited energy;Energy [keV];Counts",2000,1,2001);
    TH1F* th2=new TH1F("th2","Counts vs deposited energy;Energy [keV];Counts",2000,1,2001);
    TH1F* th3=new TH1F("th3","Counts vs deposited energy;Energy [keV];Counts",2000,1,2001);
    TH1F* th4=new TH1F("th4","Counts vs deposited energy;Energy [keV];Counts",2000,1,2001);
    TH1F* th5=new TH1F("th5","Counts vs deposited energy;Energy [keV];Counts",2000,1,2001);
    TH1F* th6=new TH1F("th6","Counts vs deposited energy;Energy [keV];Counts",2000,1,2001);
    #endif
    std::ifstream infile;
    
    ////////////////
    //No shielding//
    ////////////////
    ene_event=0;
    event_index=-1;
    for (string fileName : fileNames_noShielding)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> ene1 >> ene2;
            if(event_current != event_index && event_index!=-1)
            {
		#if ConvertToCSV
	        ene_event=ene_event*1000.;
			#if XE
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_noShielding[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_noShielding[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_noShielding[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_noShielding[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_noShielding[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_noShielding[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_noShielding[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_noShielding[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_noShielding[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_noShielding[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_noShielding[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_noShielding[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_noShielding[12]++;}
			#else
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_noShielding[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_noShielding[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_noShielding[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_noShielding[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_noShielding[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_noShielding[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_noShielding[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_noShielding[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_noShielding[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_noShielding[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_noShielding[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_noShielding[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_noShielding[12]++;}
			else if(ene_event>ene_bins[13] && ene_event<ene_bins[14]) {counter_noShielding[13]++;}
			else if(ene_event>ene_bins[14] && ene_event<ene_bins[15]) {counter_noShielding[14]++;}
			else if(ene_event>ene_bins[15] && ene_event<ene_bins[16]) {counter_noShielding[15]++;}
			else if(ene_event>ene_bins[16] && ene_event<ene_bins[17]) {counter_noShielding[16]++;}
			#endif
		#else
                th1->Fill(ene_event*1000.);
		#endif
                ene_event=0;
            }
            event_index=event_current;
            ene_event=ene_event+ene1-ene2;
        }
        infile.close();
    }
    
    /////////////
    //10 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    for (string fileName : fileNames_10ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> ene1 >> ene2;
            if(event_current != event_index && event_index!=-1)
            {
                #if ConvertToCSV
                ene_event=ene_event*1000.;
            		#if XE
            		if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_10ch2[0]++;}
            		else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_10ch2[1]++;}
            		else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_10ch2[2]++;}
            		else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_10ch2[3]++;}
            		else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_10ch2[4]++;}
            		else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_10ch2[5]++;}
            		else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_10ch2[6]++;}
            		else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_10ch2[7]++;}
            		else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_10ch2[8]++;}
            		else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_10ch2[9]++;}
            		else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_10ch2[10]++;}
            		else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_10ch2[11]++;}
            		else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_10ch2[12]++;}
            		#else
            		if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_10ch2[0]++;}
            		else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_10ch2[1]++;}
            		else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_10ch2[2]++;}
            		else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_10ch2[3]++;}
            		else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_10ch2[4]++;}
            		else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_10ch2[5]++;}
            		else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_10ch2[6]++;}
            		else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_10ch2[7]++;}
            		else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_10ch2[8]++;}
            		else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_10ch2[9]++;}
            		else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_10ch2[10]++;}
            		else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_10ch2[11]++;}
            		else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_10ch2[12]++;}
            		else if(ene_event>ene_bins[13] && ene_event<ene_bins[14]) {counter_10ch2[13]++;}
            		else if(ene_event>ene_bins[14] && ene_event<ene_bins[15]) {counter_10ch2[14]++;}
            		else if(ene_event>ene_bins[15] && ene_event<ene_bins[16]) {counter_10ch2[15]++;}
            		else if(ene_event>ene_bins[16] && ene_event<ene_bins[17]) {counter_10ch2[16]++;}
            		#endif
                #else
                th2->Fill(ene_event*1000.);
                #endif
                ene_event=0;
            }
            event_index=event_current;
            ene_event=ene_event+ene1-ene2;
        }
        infile.close();
    }
    
    /////////////
    //15 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    for (string fileName : fileNames_15ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> ene1 >> ene2;
            if(event_current == event_index){}
            if(event_current != event_index && event_index!=-1)
            {
		#if ConvertToCSV
	        ene_event=ene_event*1000.;
			#if XE
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_15ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_15ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_15ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_15ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_15ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_15ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_15ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_15ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_15ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_15ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_15ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_15ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_15ch2[12]++;}
			#else
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_15ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_15ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_15ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_15ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_15ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_15ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_15ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_15ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_15ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_15ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_15ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_15ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_15ch2[12]++;}
			else if(ene_event>ene_bins[13] && ene_event<ene_bins[14]) {counter_15ch2[13]++;}
			else if(ene_event>ene_bins[14] && ene_event<ene_bins[15]) {counter_15ch2[14]++;}
			else if(ene_event>ene_bins[15] && ene_event<ene_bins[16]) {counter_15ch2[15]++;}
			else if(ene_event>ene_bins[16] && ene_event<ene_bins[17]) {counter_15ch2[16]++;}
			#endif
		#else
                th3->Fill(ene_event*1000.);
                #endif
                ene_event=0;
            }
            event_index=event_current;
            ene_event=ene_event+ene1-ene2;
        }
        infile.close();
    }

    /////////////
    //20 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    for (string fileName : fileNames_20ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> ene1 >> ene2;
            if(event_current != event_index && event_index!=-1)
            {
		#if ConvertToCSV
	        ene_event=ene_event*1000.;
			#if XE
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_20ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_20ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_20ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_20ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_20ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_20ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_20ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_20ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_20ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_20ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_20ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_20ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_20ch2[12]++;}
			#else
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_20ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_20ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_20ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_20ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_20ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_20ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_20ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_20ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_20ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_20ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_20ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_20ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_20ch2[12]++;}
			else if(ene_event>ene_bins[13] && ene_event<ene_bins[14]) {counter_20ch2[13]++;}
			else if(ene_event>ene_bins[14] && ene_event<ene_bins[15]) {counter_20ch2[14]++;}
			else if(ene_event>ene_bins[15] && ene_event<ene_bins[16]) {counter_20ch2[15]++;}
			else if(ene_event>ene_bins[16] && ene_event<ene_bins[17]) {counter_20ch2[16]++;}
			#endif
		#else
                th4->Fill(ene_event*1000.);
                #endif
                ene_event=0;
            }
            event_index=event_current;
            ene_event=ene_event+ene1-ene2;
        }
        infile.close();
    }
    
    /////////////
    //25 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    for (string fileName : fileNames_25ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> ene1 >> ene2;
            if(event_current != event_index && event_index!=-1)
            {
		#if ConvertToCSV
	        ene_event=ene_event*1000.;
			#if XE
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_25ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_25ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_25ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_25ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_25ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_25ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_25ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_25ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_25ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_25ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_25ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_25ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_25ch2[12]++;}
			#else
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_25ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_25ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_25ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_25ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_25ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_25ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_25ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_25ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_25ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_25ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_25ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_25ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_25ch2[12]++;}
			else if(ene_event>ene_bins[13] && ene_event<ene_bins[14]) {counter_25ch2[13]++;}
			else if(ene_event>ene_bins[14] && ene_event<ene_bins[15]) {counter_25ch2[14]++;}
			else if(ene_event>ene_bins[15] && ene_event<ene_bins[16]) {counter_25ch2[15]++;}
			else if(ene_event>ene_bins[16] && ene_event<ene_bins[17]) {counter_25ch2[16]++;}
			#endif
		#else
                th5->Fill(ene_event*1000.);
                #endif
                ene_event=0;
            }
            event_index=event_current;
            ene_event=ene_event+ene1-ene2;
        }
        infile.close();
    }

    /////////////
    //30 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    for (string fileName : fileNames_30ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> ene1 >> ene2;
            if(event_current != event_index && event_index!=-1)
            {
		#if ConvertToCSV
	        ene_event=ene_event*1000.;
			#if XE
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_30ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_30ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_30ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_30ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_30ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_30ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_30ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_30ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_30ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_30ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_30ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_30ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_30ch2[12]++;}
			#else
			if(ene_event>ene_bins[0] && ene_event<ene_bins[1]) {counter_30ch2[0]++;}
			else if(ene_event>ene_bins[1] && ene_event<ene_bins[2]) {counter_30ch2[1]++;}
			else if(ene_event>ene_bins[2] && ene_event<ene_bins[3]) {counter_30ch2[2]++;}
			else if(ene_event>ene_bins[3] && ene_event<ene_bins[4]) {counter_30ch2[3]++;}
			else if(ene_event>ene_bins[4] && ene_event<ene_bins[5]) {counter_30ch2[4]++;}
			else if(ene_event>ene_bins[5] && ene_event<ene_bins[6]) {counter_30ch2[5]++;}
			else if(ene_event>ene_bins[6] && ene_event<ene_bins[7]) {counter_30ch2[6]++;}
			else if(ene_event>ene_bins[7] && ene_event<ene_bins[8]) {counter_30ch2[7]++;}
			else if(ene_event>ene_bins[8] && ene_event<ene_bins[9]) {counter_30ch2[8]++;}
			else if(ene_event>ene_bins[9] && ene_event<ene_bins[10]) {counter_30ch2[9]++;}
			else if(ene_event>ene_bins[10] && ene_event<ene_bins[11]) {counter_30ch2[10]++;}
			else if(ene_event>ene_bins[11] && ene_event<ene_bins[12]) {counter_30ch2[11]++;}
			else if(ene_event>ene_bins[12] && ene_event<ene_bins[13]) {counter_30ch2[12]++;}
			else if(ene_event>ene_bins[13] && ene_event<ene_bins[14]) {counter_30ch2[13]++;}
			else if(ene_event>ene_bins[14] && ene_event<ene_bins[15]) {counter_30ch2[14]++;}
			else if(ene_event>ene_bins[15] && ene_event<ene_bins[16]) {counter_30ch2[15]++;}
			else if(ene_event>ene_bins[16] && ene_event<ene_bins[17]) {counter_30ch2[16]++;}
			#endif
		#else
                th6->Fill(ene_event*1000.);
                #endif
                ene_event=0;
            }
            event_index=event_current;
            ene_event=ene_event+ene1-ene2;
        }
        infile.close();
    }
    
    ////////////////////
    //Drawing & output//
    ////////////////////
    #if ConvertToCSV
    	#if XE
    	    ofstream outFile("ene_xe_result.csv");
    	#else
    	    ofstream outFile("ene_ar_result.csv");
    	#endif
        outFile<<"E_min E_max Counts_NoShielding Counts_10_cm_CH2 Counts_15_cm_CH2 Counts_20_cm_CH2 Counts_25_cm_CH2 Counts_30_cm_CH2"<<endl;
    	#if XE
    	for(int it=0;it<13;it++)
    	#else
    	for(int it=0;it<17;it++)
        #endif
	{
        #if Cluster
	    counter_noShielding[it]=counter_noShielding[it]*norm_neutron*area_noShielding/TMath::Power(10,10);
	    counter_10ch2[it]=counter_10ch2[it]*norm_neutron*area_10ch2/TMath::Power(10,10);
	    counter_15ch2[it]=counter_15ch2[it]*norm_neutron*area_15ch2/TMath::Power(10,10);
	    counter_20ch2[it]=counter_20ch2[it]*norm_neutron*area_20ch2/TMath::Power(10,10);
	    counter_25ch2[it]=counter_25ch2[it]*norm_neutron*area_25ch2/TMath::Power(10,10);
	    counter_30ch2[it]=counter_30ch2[it]*norm_neutron*area_30ch2/TMath::Power(10,10);
    	    outFile<<ene_bins[it]<<" "<<ene_bins[it+1]<<" "<<counter_noShielding[it]<<" "<<counter_10ch2[it]<<" "<<counter_15ch2[it]<<" "<<counter_20ch2[it]<<" "<<counter_25ch2[it]<<" "<<counter_30ch2[it]<<endl;
        #else
	    counter_noShielding[it]=counter_noShielding[it]*norm_neutron*area_30ch2/TMath::Power(10,9);
	    counter_10ch2[it]=counter_10ch2[it]*norm_neutron*area_30ch2/TMath::Power(10,9);
	    counter_15ch2[it]=counter_15ch2[it]*norm_neutron*area_30ch2/TMath::Power(10,9);
	    counter_20ch2[it]=counter_20ch2[it]*norm_neutron*area_20ch2/TMath::Power(10,9);
	    counter_25ch2[it]=counter_25ch2[it]*norm_neutron*area_25ch2/TMath::Power(10,9);
	    counter_30ch2[it]=counter_30ch2[it]*norm_neutron*area_30ch2/TMath::Power(10,9);
    	    outFile<<ene_bins[it]<<" "<<ene_bins[it+1]<<" "<<counter_noShielding[it]<<" "<<counter_10ch2[it]<<" "<<counter_15ch2[it]<<" "<<counter_20ch2[it]<<" "<<counter_25ch2[it]<<" "<<counter_30ch2[it]<<endl;
        #endif
	}
    #else
    TCanvas* c1=new TCanvas();
        #if Cluster
        th1->Scale((norm_neutron*area_noShielding/TMath::Power(10,10)));
        th2->Scale((norm_neutron*area_10ch2/TMath::Power(10,10)));
        th3->Scale((norm_neutron*area_15ch2/TMath::Power(10,10)));
        th4->Scale((norm_neutron*area_20ch2/TMath::Power(10,10)));
        th5->Scale((norm_neutron*area_25ch2/TMath::Power(10,10)));
        th6->Scale((norm_neutron*area_30ch2/TMath::Power(10,10)));
        #else
        th1->Scale((norm_neutron*area_30ch2/TMath::Power(10,9)));
        th2->Scale((norm_neutron*area_30ch2/TMath::Power(10,9)));
        th3->Scale((norm_neutron*area_30ch2/TMath::Power(10,9)));
        th4->Scale((norm_neutron*area_20ch2/TMath::Power(10,9)));
        th5->Scale((norm_neutron*area_25ch2/TMath::Power(10,9)));
        th6->Scale((norm_neutron*area_30ch2/TMath::Power(10,9)));
        #endif
    th1->SetLineColor(1);
    th1->SetLineWidth(2);
    th2->SetLineColor(2);
    th2->SetLineWidth(2);
    //th3->SetLineColor(2);
    //th3->SetLineWidth(2);
    th4->SetLineColor(3);
    th4->SetLineWidth(2);
    //th5->SetLineColor(5);
    //th5->SetLineWidth(2);
    th6->SetLineColor(4);
    th6->SetLineWidth(2);
    th1->Draw();
    th1->GetYaxis()->SetRangeUser(1.,10.e8);
    th2->Draw("same");
    //th3->Draw("same");
    th4->Draw("same");
    //th5->Draw("same");
    th6->Draw("same");
    #endif
}
