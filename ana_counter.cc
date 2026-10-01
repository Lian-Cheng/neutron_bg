#include <iostream>
#include <fstream>
#include <string>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

void ana_counter()
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

    const char* input_filename = "output.txt";
    const char* output_filename = "results/output_counter.root";
    Int_t event_index;
    Int_t surface;
    Float_t x;
    Float_t y;
    Float_t z;
    Float_t ene;
    Float_t dir_x;
    Float_t dir_y;
    Float_t dir_z;
    Float_t costheta;

    std::ifstream infile(input_filename);
    TFile* outfile = new TFile(output_filename, "RECREATE");

    TTree* tree = new TTree("data_tree", "Data from Text File");
    tree->Branch("event_index", &event_index, "event_index/I");
    tree->Branch("surface", &x, "surface/I");
    tree->Branch("x", &x, "x/F");
    tree->Branch("y", &y, "y/F");
    tree->Branch("z", &z, "z/F");
    tree->Branch("ene", &ene, "ene/F");
    tree->Branch("dir_x", &dir_x, "dir_x/F");
    tree->Branch("dir_y", &dir_y, "dir_y/F");
    tree->Branch("dir_z", &dir_z, "dir_z/F");

    //Energies
    TH1F* th1=new TH1F("th1","Counts vs Initial energy;Energy [MeV];Counts",20000,0,200);
    //Positions
    TH2F* th2_xy=new TH2F("th2_xy","Position distribution;X [mm];Y [mm]",3000,-1500,1500,3000,-1500,1500);
    TH2F* th2_xz=new TH2F("th2_xz","Position distribution;X [mm];Z [mm]",3000,-1500,1500,3600,-1800,1800);
    TH2F* th2_yz=new TH2F("th2_yz","Position distribution;Y [mm];Z [mm]",3000,-1500,1500,3600,-1800,1800);
    TH3F* th3_z=new TH3F("th3_z","Position distribution;X [mm];Y [mm];Z [mm]",300,-1500,1500,300,-1500,1500,360,-1800,1800);
    TH3F* th3_y=new TH3F("th3_y","Position distribution;X [mm];Y [mm];Z [mm]",300,-1500,1500,300,-1500,1500,360,-1800,1800);
    TH3F* th3_x=new TH3F("th3_x","Position distribution;X [mm];Y [mm];Z [mm]",300,-1500,1500,300,-1500,1500,360,-1800,1800);
    //Angles
    TH1F* th1_costheta=new TH1F("th1_costheta","Costheta distribution;Phi;Counts",100,-1,1);
    TH1F* th1_phi_xy=new TH1F("th1_phi_xy","Angular phi distribution;Phi;Counts",400,-2,2);
    TH1F* th1_phi_xz=new TH1F("th1_phi_xz","Angular phi distribution;Phi;Counts",400,-2,2);
    TH1F* th1_phi_yz=new TH1F("th1_phi_yz","Angular phi distribution;Phi;Counts",400,-2,2);
    TH1F* th1_theta_x=new TH1F("th1_theta_x","Angular theta distribution;Theta;Counts",400,-2,2);
    TH1F* th1_theta_y=new TH1F("th1_theta_y","Angular theta distribution;Theta;Counts",400,-2,2);
    TH1F* th1_theta_z=new TH1F("th1_theta_z","Angular theta distribution;Theta;Counts",400,-2,2);

    std::string line;
    Float_t src_x = 820. + 10. * 2. + 150. * 2. + 300. * 2.; //width direction
    Float_t src_y = 820. + 10. * 2. + 150. * 2. + 300. * 2.;
    Float_t src_z = 1110. + 100. * 2. + 10. * 2. + 150. * 2. + 300. * 2.;//length direction
    while (std::getline(infile, line))
    {
        std::stringstream ss(line);
        ss >> event_index >> ene >> surface >> x >> y >> z >> dir_x >> dir_y >> dir_z >> costheta;
	//if(event_index == 1000){break;}
	//cout<< event_index << " " << surface << " " << ene << " " << x << " " << y << " " << z << " " << dir_x << " " << dir_y << " " << dir_z << endl;
	th1->Fill(ene);
	Float_t theta;
	Float_t phi;
	if(surface == 0 || surface == 1)
	{
	    //th2_xy->Fill(x,y);
	    th3_z->Fill(x,y,z);
	    phi = std::atan(dir_y/dir_x);
	    th1_phi_xy->Fill(phi);
	    //theta = std::sqrt(1- TMath::Power(std::acos(-z/TMath::Abs(z)*dir_z),2));
	    theta = -z/TMath::Abs(z)*dir_z;
	    th1_theta_z->Fill(theta);
	    th1_costheta->Fill(costheta);
	}
	else if(surface == 2 || surface == 3)
	{
	    //th2_xz->Fill(x,z);
	    th3_x->Fill(x,y,z);
	    phi = std::atan(dir_z/dir_y);
	    th1_phi_yz->Fill(phi);
	    theta = -x/TMath::Abs(x)*dir_x;
	    th1_theta_x->Fill(theta);
	}
	else if(surface == 4 || surface == 5)
	{
	    //th2_yz->Fill(y,z);
	    th3_y->Fill(x,y,z);
	    phi = std::atan(dir_z/dir_x);
	    th1_phi_xz->Fill(phi);
	    theta = -y/TMath::Abs(y)*dir_y;
	    th1_theta_y->Fill(theta);
	}
	else
	{
	    cout<< event_index << " " << ene << " " << surface << " " << x << " " << y << " " << z << " " << dir_x << " " << dir_y << " " << dir_z << endl;
	}
        //tree->Fill();
    }
    th1->Write();
    //th3_x->Write();
    //th3_y->Write();
    //th3_z->Write();
    //th2_xy->Write();
    //th2_xz->Write();
    //th2_yz->Write();
    //th1_phi_xy->Write();
    //th1_phi_xz->Write();
    //th1_phi_yz->Write();
    //th1_theta_x->Write();
    //th1_theta_y->Write();
    //th1_theta_z->Write();
    th1_costheta->Write();
    infile.close();
    outfile->Close();
}
