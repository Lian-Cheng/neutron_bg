#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#define XE 1

using std::string;
struct event_struct
{
        Float_t ene;
        string name;
        string process;
        Float_t pos_x;
        Float_t pos_y;
        Float_t pos_z;
};

Int_t positions(std::vector<event_struct>& vec)
{
    float threshold = 10.;   // 1 cm
    float threshold2 = threshold * threshold;
    std::vector<bool> used(vec.size(), false);
    Int_t distinct = 0;
    for (size_t i = 0; i < vec.size(); i++)
    {
        if (used[i]) continue;
        distinct++;  // new cluster found
        used[i] = true;
        for (size_t j = i + 1; j < vec.size(); j++)
	{
            if (used[j]) continue;
            float dx = vec[i].pos_x - vec[j].pos_x;
            float dy = vec[i].pos_y - vec[j].pos_y;
            float dz = vec[i].pos_z - vec[j].pos_z;
            // Compare squared distance (avoid sqrt)
            if (dx*dx + dy*dy + dz*dz < threshold2){used[j] = true;}
        }
    }
    return distinct;
}

void ana_positions()
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

    ////////////////////////
    //Reading result files//
    ////////////////////////
    std::vector<std::string> fileNames_noShielding;
    std::vector<std::string> fileNames_10ch2;
    std::vector<std::string> fileNames_20ch2;
    std::vector<std::string> fileNames_30ch2;
    #if XE
    Double_t ene_bins[14] = {0.9,2.94451766,4.3987765,6.14103535,8.17129421,10.48955308,13.09581196,15.99007085,19.17232974,22.64258864,26.40084755,30.44710647,34.7813654,60.77266017};
    const std::filesystem::path dir_noShielding="results_positions/results_xe_noShielding/";
    const std::filesystem::path dir_10ch2="results_positions/results_xe_10ch2/";
    const std::filesystem::path dir_20ch2="results_positions/results_xe_20ch2/";
    const std::filesystem::path dir_30ch2="results_positions/results_xe_30ch2/";
    #else
    Double_t ene_bins[18] = {0,8.17129421,10.48955308,13.09581196,15.99007085,19.17232974,22.64258864,26.40084755,30.44710647,34.7813654,39.40362434,44.31388329,49.51214224,54.9984012,60.77266017,73.18517814,86.74969615,143.29750943};
    const std::filesystem::path dir_noShielding="results_positions/results_ar_noShielding/";
    const std::filesystem::path dir_10ch2="results_positions/results_ar_10ch2/";
    const std::filesystem::path dir_20ch2="results_positions/results_ar_20ch2/";
    const std::filesystem::path dir_30ch2="results_positions/results_ar_30ch2/";
    #endif
    for (auto const& entry : std::filesystem::directory_iterator{dir_noShielding}){
        if(entry.is_regular_file()){fileNames_noShielding.push_back(entry.path());}
    }
    //fileNames_noShielding.push_back("output.txt");
    for (auto const& entry : std::filesystem::directory_iterator{dir_10ch2}){
        if(entry.is_regular_file()){fileNames_10ch2.push_back(entry.path());}
    }
    for (auto const& entry : std::filesystem::directory_iterator{dir_20ch2}){
        if(entry.is_regular_file()){fileNames_20ch2.push_back(entry.path());}
    }
    for (auto const& entry : std::filesystem::directory_iterator{dir_30ch2}){
        if(entry.is_regular_file()){fileNames_30ch2.push_back(entry.path());}
    }

    ////////////////////////
    //Variables definition//
    ////////////////////////
    event_struct this_struct;
    Int_t event_index,track_id,parentID,mass_target;
    std::string name,process,cell;
    Float_t pre_Pos_x,pre_Pos_y,pre_Pos_z,post_Pos_x,post_Pos_y,post_Pos_z;
    Float_t ene,ene_event,globalTime;
    std::ifstream infile;
    std::string line;
    std::vector<event_struct> interactions;
    TH1F* th1_all=new TH1F("th1_all","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th1_1cm=new TH1F("th1_1cm","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th2_all=new TH1F("th2_all","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th2_1cm=new TH1F("th2_1cm","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th3_all=new TH1F("th3_all","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th3_1cm=new TH1F("th3_1cm","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th4_all=new TH1F("th4_all","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);
    TH1F* th4_1cm=new TH1F("th4_1cm","Counts vs Number of interactions per neutron; N ;Counts",10,0,10);

    ////////////////
    //No shielding//
    ////////////////
    ene_event=0;
    event_index=-1;
    interactions.clear();
    for (string fileName : fileNames_noShielding)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> track_id >> parentID >> this_struct.name >> this_struct.process >> cell >> this_struct.ene >> mass_target >> globalTime >> this_struct.pos_x >> this_struct.pos_y >> this_struct.pos_z;
            if(event_current != event_index && event_index!=-1)
            {
		    Int_t count = interactions.size();
	            //if (count>9) {std::cout<<fileName<<" "<<event_index<<" "<<interactions[10].name<<","<<interactions[10].process<<",pos_z= "<<interactions[10].pos_z<<std::endl;}
		    th1_all->Fill(count);
		    th1_1cm->Fill(positions(interactions));
		    interactions.clear();
		    std::vector<event_struct>().swap(interactions);
            }
	    if (this_struct.name == "neutron"){interactions.push_back(this_struct);}
            event_index=event_current;
            ene_event=ene_event+this_struct.ene;
	    this_struct = event_struct();
        }
        infile.close();
    }
    
    /////////////
    //10 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    interactions.clear();
    for (string fileName : fileNames_10ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> track_id >> parentID >> this_struct.name >> this_struct.process >> cell >> this_struct.ene >> mass_target >> globalTime >> this_struct.pos_x >> this_struct.pos_y >> this_struct.pos_z;
            if(event_current != event_index && event_index!=-1)
            {
		    Int_t count = interactions.size();
	            //if (count>9) {std::cout<<fileName<<" "<<event_index<<" "<<interactions[10].name<<","<<interactions[10].process<<",pos_z= "<<interactions[10].pos_z<<std::endl;}
		    th2_all->Fill(count);
		    th2_1cm->Fill(positions(interactions));
		    interactions.clear();
		    std::vector<event_struct>().swap(interactions);
            }
	    if (this_struct.name == "neutron"){interactions.push_back(this_struct);}
            event_index=event_current;
            ene_event=ene_event+this_struct.ene;
	    this_struct = event_struct();
        }
        infile.close();
    }
    
    /////////////
    //20 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    interactions.clear();
    for (string fileName : fileNames_20ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> track_id >> parentID >> this_struct.name >> this_struct.process >> cell >> this_struct.ene >> mass_target >> globalTime >> this_struct.pos_x >> this_struct.pos_y >> this_struct.pos_z;
            if(event_current != event_index && event_index!=-1)
            {
		    Int_t count = interactions.size();
	            //if (count>9) {std::cout<<fileName<<" "<<event_index<<" "<<interactions[10].name<<","<<interactions[10].process<<",pos_z= "<<interactions[10].pos_z<<std::endl;}
		    th3_all->Fill(count);
		    th3_1cm->Fill(positions(interactions));
		    interactions.clear();
		    std::vector<event_struct>().swap(interactions);
            }
	    if (this_struct.name == "neutron"){interactions.push_back(this_struct);}
            event_index=event_current;
            ene_event=ene_event+this_struct.ene;
	    this_struct = event_struct();
        }
        infile.close();
    }

    /////////////
    //30 cm CH2//
    /////////////
    ene_event=0;
    event_index=-1;
    interactions.clear();
    for (string fileName : fileNames_30ch2)
    {
        infile.open(fileName.c_str());
        while (std::getline(infile, line))
        {
            Int_t event_current;
            std::stringstream ss(line);
            ss >> event_current >> track_id >> parentID >> this_struct.name >> this_struct.process >> cell >> this_struct.ene >> mass_target >> globalTime >> this_struct.pos_x >> this_struct.pos_y >> this_struct.pos_z;
            if(event_current != event_index && event_index!=-1)
            {
		    Int_t count = interactions.size();
	            //if (count>9) {std::cout<<fileName<<" "<<event_index<<" "<<interactions[10].name<<","<<interactions[10].process<<",pos_z= "<<interactions[10].pos_z<<std::endl;}
		    th4_all->Fill(count);
		    th4_1cm->Fill(positions(interactions));
		    interactions.clear();
		    std::vector<event_struct>().swap(interactions);
            }
	    if (this_struct.name == "neutron"){interactions.push_back(this_struct);}
            event_index=event_current;
            ene_event=ene_event+this_struct.ene;
	    this_struct = event_struct();
        }
        infile.close();
    }
    
    ////////////////////
    //Drawing & output//
    ////////////////////
    TCanvas* c1=new TCanvas("c1","c1",1200,1200);
    c1->Divide(2,2);
    c1->cd(1);
    th1_all->Scale((norm_neutron*area_noShielding/TMath::Power(10,10)));
    th1_1cm->Scale((norm_neutron*area_noShielding/TMath::Power(10,10)));
    th1_all->GetYaxis()->SetRangeUser(0.,500.);
    th1_all->SetLineColor(1);
    th1_all->SetLineWidth(2);
    th1_1cm->SetLineColor(2);
    th1_1cm->SetLineWidth(2);
    th1_all->Draw();
    th1_1cm->Draw("same");
    c1->cd(2);
    th2_all->Scale((norm_neutron*area_10ch2/TMath::Power(10,10)));
    th2_1cm->Scale((norm_neutron*area_10ch2/TMath::Power(10,10)));
    th2_all->SetLineColor(1);
    th2_all->SetLineWidth(2);
    th2_1cm->SetLineColor(2);
    th2_1cm->SetLineWidth(2);
    th2_all->GetYaxis()->SetRangeUser(0.,500.);
    th2_all->Draw();
    th2_1cm->Draw("same");
    c1->cd(3);
    th3_all->Scale((norm_neutron*area_15ch2/TMath::Power(10,10)));
    th3_1cm->Scale((norm_neutron*area_10ch2/TMath::Power(10,10)));
    th3_all->SetLineColor(1);
    th3_all->SetLineWidth(2);
    th3_1cm->SetLineColor(2);
    th3_1cm->SetLineWidth(2);
    th3_all->GetYaxis()->SetRangeUser(0.,500.);
    th3_all->Draw();
    th3_1cm->Draw("same");
    c1->cd(4);
    th4_all->Scale((norm_neutron*area_20ch2/TMath::Power(10,10)));
    th4_1cm->Scale((norm_neutron*area_10ch2/TMath::Power(10,10)));
    th4_all->SetLineColor(1);
    th4_all->SetLineWidth(2);
    th4_1cm->SetLineColor(2);
    th4_1cm->SetLineWidth(2);
    th4_all->GetYaxis()->SetRangeUser(0.,500.);
    th4_all->Draw();
    th4_1cm->Draw("same");
}
