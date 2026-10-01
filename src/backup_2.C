#include <n4-all.hh>
#include <G4ParticleGun.hh>
#include "G4EmStandardPhysics.hh"
#include <G4GenericMessenger.hh>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <G4PrimaryParticle.hh>
#include <G4Material.hh>
#include "G4NistManager.hh"
#include <G4String.hh>
#include <G4SystemOfUnits.hh>   // physical units such as `m` for metre
#include <G4Event.hh>           // needed to inject primary particles into an event
#include <G4Box.hh>             // for creating shapes in the geometry
#include <G4SubtractionSolid.hh>// for boolean solids
#include <G4UnionSolid.hh>      //
#include <G4MultiUnion.hh>      //
#include <G4Sphere.hh>          // for creating shapes in the geometry
#include <FTFP_BERT.hh>         // our choice of physics list
#include <G4RandomDirection.hh> // for launching particles in random directions
#include <Randomize.hh>
#include <G4ThreeVector.hh>
#include <cstdlib>
#define VISUAL 1
#define SHIELDING 1
#define XE 0
#define DIR_CENTERED 0

std::ofstream output;
std::ifstream inFile;
std::vector<double> flux;
std::vector<double> energy;

std::map<std::string, G4double> tempInitialEnergies;
std::map<std::string, G4double> tempInitialTimes;
std::map<std::string, G4String> tempParticleNames;
std::map<std::string, G4String> tempInitialVolumes;
std::map<G4int, G4double> secondaryInitialEnergies;
std::map<G4int, G4String> secondaryParticleNames;
std::map<std::string, G4double> reducingEnergies;
std::map<std::string, G4String> reducingMothers;

const G4double thick_ch2 = 50. * mm;
const G4double thick_pb = 150 * mm;
const G4double thick_gap = 10 * mm;
const G4double MeV_to_amu = 1/931.5;
const G4double pi = M_PI;

std::vector<std::string> pendingIdentifiers;  // identifiers created in this step
std::map<std::string, int> identifierToManualTrackID;  // temp mapping

G4int manualTrackID;

struct my
{
	G4double msg_pressure{20 * bar};
	unsigned int ch2_counter = 0;
	unsigned int lead_counter = 0;
	unsigned int vessel_counter = 0;
	unsigned int drift_counter = 0;
	bool flag_ch2=false;
	bool flag_lead=false;
	bool flag_vessel=false;
	bool flag_drift=false;
};

auto my_generator(const my& my)
{
	return [&](G4Event* event)
	{
		std::random_device rd;
		std::mt19937 gen(rd());
		////////////////////////////////////////////////
		//Random distribution for positions and angles//
		////////////////////////////////////////////////
		//Surface of the Ch2 shielding
		std::uniform_real_distribution<G4double> dist0(0.0, 1.0);
		#if SHIELDING
			std::uniform_real_distribution<G4double> dist1(0.0, 1.0);
			std::uniform_real_distribution<G4double> dist2(0.0, 1.0);
			G4double src_x = 820. + 10. * 2. + 150. * 2. + thick_ch2 * 2.; //width direction
			G4double src_y = 820. + 10. * 2. + 150. * 2. + thick_ch2 * 2.;
			G4double src_z = 1110. + 100. * 2. + 10. * 2. + 150. * 2. + thick_ch2 * 2.;//length direction
			G4double area_1 = src_x * src_y;
			G4double area_2 = src_y * src_z;
			G4double area_3 = src_z * src_x;
			G4double area_tot = 2. * (area_1 + area_2 + area_3);//total surface
		#else
			std::uniform_real_distribution<G4double> dist1(0.0, 1.0);
			std::uniform_real_distribution<G4double> alpha(0.0, 2. * M_PI);
			G4double diameter_xy = 820.;
			G4double src_z = 1110. + 100. * 2.;//length direction
			G4double area_1 = M_PI * diameter_xy * diameter_xy / 4.;
			G4double area_2 = M_PI * diameter_xy * src_z;
			G4double area_tot = 2. * area_1 + area_2;//total surface
		#endif
		G4double src_a = dist0(gen) * area_tot;
		std::uniform_real_distribution<G4double> dist_cosTheta(0.0, 1.);
		std::uniform_real_distribution<G4double> dist_phi(0.0, 2 * M_PI);
		G4double phi = dist_phi(gen);
		G4double cosPhi = std::cos(phi);
		G4double sinPhi = std::sin(phi);
		G4double cosTheta = dist_cosTheta(gen);
		G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
		#if DIR_CENTERED
			G4int dir_centered = 1;
		#else
			G4int dir_centered = 0;
		#endif

		////////////////////////////////////////////////////////////////
		//Position and momentum direction distribution of neutron flux//
		////////////////////////////////////////////////////////////////
		G4double x,y,z;
		G4double x_dir,y_dir,z_dir;
		short surface=0;
		#if SHIELDING
		if(src_a <= area_1)
		{
			surface=0;
			x = (dist1(gen)-0.5)*src_x;
			y = (dist2(gen)-0.5)*src_y;
			z = - src_z/2.;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = sinTheta * cosPhi;
				y_dir = sinTheta * sinPhi;
				z_dir = cosTheta;
			}
		}
		else if(src_a <= area_1*2.)
		{
			surface=1;
			x = (dist1(gen)-0.5)*src_x;
			y = (dist2(gen)-0.5)*src_y;
			z = src_z/2.;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = sinTheta * cosPhi;
				y_dir = sinTheta * sinPhi;
				z_dir = - cosTheta;
			}
		}
		else if(src_a <= area_1*2. + area_2)
		{
			surface=2;
			x = src_x/2.;
			y = (dist1(gen)-0.5)*src_y;
			z = (dist2(gen)-0.5)*src_z;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = - cosTheta;
				y_dir = sinTheta * cosPhi;
				z_dir = sinTheta * sinPhi;
			}
		}
		else if(src_a <= area_tot - area_3*2.)
		{
			surface=3;
			x = - src_x/2.;
			y = (dist1(gen)-0.5)*src_y;
			z = (dist2(gen)-0.5)*src_z;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = cosTheta;
				y_dir = sinTheta * cosPhi;
				z_dir = sinTheta * sinPhi;
			}
		}
		else if(src_a <= area_tot - area_3)
		{
			surface=4;
			x = (dist1(gen)-0.5)*src_x;
			y = src_y/2.;
			z = (dist2(gen)-0.5)*src_z;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = sinTheta * sinPhi;
				y_dir = - cosTheta;
				z_dir = sinTheta * cosPhi;
			}
		}
		else
		{
			surface=5;
			x = (dist1(gen)-0.5)*src_x;
			y = - src_y/2.;
			z = (dist2(gen)-0.5)*src_z;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = sinTheta * sinPhi;
				y_dir = cosTheta;
				z_dir = sinTheta * cosPhi;
			}
		}
		#else
		G4double fvalue = alpha(gen);
		if(src_a <= area_1)
		{
			surface=0;
			x = dist1(gen)*diameter_xy/2. * std::cos(fvalue);
			y = dist1(gen)*diameter_xy/2. * std::sin(fvalue);
			z = - src_z/2.;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = cosTheta * cosPhi;
				y_dir = cosTheta * sinPhi;
				z_dir = sinTheta;
			}
		}
		else if(src_a <= area_1*2.)
		{
			surface=1;
			x = dist1(gen)*diameter_xy/2. * std::cos(fvalue);
			y = dist1(gen)*diameter_xy/2. * std::sin(fvalue);
			z = src_z/2.;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				x_dir = cosTheta * cosPhi;
				y_dir = cosTheta * sinPhi;
				z_dir = - sinTheta;
			}
		}
		else
		{
			surface=2;
			x = diameter_xy/2. * std::cos(fvalue);
			y = diameter_xy/2. * std::sin(fvalue);
			z = (dist1(gen)-0.5)*src_z;
			if(dir_centered)
			{
				G4double m_dir = std::sqrt(x*x + y*y + z*z);
				x_dir = - x/m_dir;
				y_dir = - y/m_dir;
				z_dir = - z/m_dir;
			}
			else
			{
				///////////////////////////////////////////
				//the normal inward angle for x,y,z:     //
				//u=(-cos(alpha),-sin(alpha),0)          //
				//new coordinates, (p,q,u):              //
				//p=(sin(alpha),-cos(alpha))             //
				//q=(0,0,1)                              //
				//the half sphere centered on u:         //
				//m= cos(phi)*u                          //
				//n= sin(phi)*cos(theta)*p               //
				//t= sin(phi)*sin(theta)*q               //
				//then the random momentum direction is: //
				//x= sin(theta)*cos(phi)*sin(alpha)      //
				//   - cos(theta)*cos(alpha)             //
				//y= -sin(theta)*cos(phi)*cos(alpha)     //
				//   - cos(theta)*cos(alpha)             //
				//z= sin(theta)*sin(phi)                 //
				///////////////////////////////////////////
				x_dir = sinTheta * cosPhi * std::sin(fvalue) - cosTheta * std::cos(fvalue);
				y_dir = - sinTheta * cosPhi * std::cos(fvalue) - cosTheta * std::sin(fvalue);
				z_dir = sinTheta * sinPhi;
			}
		}
		#endif

		///////////////////////////////////////
		//Energy distribution of neutron flux//
		///////////////////////////////////////
		//To generate primary neutrons with energies according to the flux
		std::uniform_real_distribution<> dist_flux(0.0, 1.0);
		if (!inFile.is_open())
		{
			G4cerr << "Error : impossible to open the flux file" << G4endl;
			return;
		}
		std::string line;
		while (std::getline(inFile, line))
		{
			std::istringstream iss(line);
			double energy_value;
			double flux_value;
			
			if (!(iss >> energy_value >> flux_value))
			{
				G4cerr << "Error at the line : " << line << G4endl;
				continue;
			}
			//G4cout <<"flux value is :" << flux_value <<  "energy value is :" << energy_value << G4endl;
			energy.emplace_back(energy_value);
			flux.emplace_back(flux_value);
		}

		std::vector<G4double> cdf;
		G4double somme = 0.0;
		cdf.push_back(0.0);
		for (G4int it = 0; it < flux.size() - 1; it++)
		{
			G4double trapezoidal = 0.5 * (flux[it] + flux[it+1]) * (energy[it+1] - energy[it]);
			somme += trapezoidal;
			cdf.push_back(somme);
		}

		G4double p_flux = dist_flux(gen) * somme;
		auto it = std::lower_bound(cdf.begin(), cdf.end(), p_flux);
		size_t index = std::distance(cdf.begin(), it);
		if(index == 0){index=1;}
		G4double f1 = flux[index-1];
		G4double f2 = flux[index];
		G4double c1 = cdf[index-1];
		G4double c2 = cdf[index];
		G4double e1 = energy[index-1];
		G4double e2 = energy[index];
		G4double diff_f = f2 - f1;
		G4double diff_c = c2 - c1;
		G4double diff_e = e2 - e1;
		G4double sampled_energy;
		//sampled_energy =1;
		//sampled_energy = e1 + (p_flux - c1) * diff_e / diff_c; //linear extrapolation
		// --- Quadratic Interpolation ---
		if (std::abs(diff_f) > 1.0e-12)
		{
			G4double factor_a = 0.5 * diff_f * diff_e;
			G4double factor_b = diff_e * f1;
			G4double factor_c = c1 - p_flux;
			//discriminant of the quadratic equation
			G4double discriminant = factor_b * factor_b - 4 * factor_a * factor_c;
			if (discriminant < 0)
			{
				G4cerr << "Warning: Negative discriminant in quadratic sampling. Falling back to linear interpolation." << G4endl;
				sampled_energy = e1 + diff_e * (p_flux - c1) / diff_c; 
			}
			else
			{
			    G4double sqrt_discriminant = std::sqrt(discriminant);
			    // The fractional distance of the energy in the bin (0 to 1)
			    G4double frac_e = (- factor_b + sqrt_discriminant) / 2.0 / factor_a;
			    // Final sampled energy
			    sampled_energy = e1 + diff_e * frac_e;
			}
		} 
		else
		{
			// If flux is constant within the bin e1-e2, the CDF is linear in the bin, 
			// and the quadratic solution simplifies to linear interpolation.
			sampled_energy = e1 + diff_e * (p_flux - c1) / diff_c;
		}

		//////////////////////////////
		//Generation of neutron flux//
		//////////////////////////////
		auto neutron = n4::find_particle("neutron");
		auto particle_gun = std::make_shared<G4ParticleGun>(neutron, 1); // 1 particle
		G4double mass = neutron -> GetPDGMass();
		//std::cout << "e="<<sampled_energy*MeV<<" p="<<p_flux <<" c1="<<c1<<" c2="<<c2<<" e1="<<e1 <<" e2="<<e2 << std::endl;
		particle_gun -> SetParticleEnergy(sampled_energy*MeV);
		particle_gun -> SetParticlePosition({x, y, z});
		particle_gun -> SetParticleMomentumDirection({x_dir, y_dir ,z_dir});
		particle_gun -> GeneratePrimaryVertex(event);
		//output << n_event << "   " << sampled_energy << " " << surface << " " << x << "  " << y << "  " << z << "  " << x_dir << "  " << y_dir << "  " << z_dir << "  " << cosTheta <<std::endl;
	};
}

n4::sensitive_detector* sensitive_detector(unsigned& n_event)
{
	auto process_hits = [&] (G4Step* step){return true;};// See https://jacg.github.io/nain4/explanation/process-hits-return-value.html
	return new n4::sensitive_detector{"Detector", process_hits};
}

n4::actions* create_actions(my& my, unsigned& n_event, int* ene_count)
{
	auto my_run_action_end = [&] (const G4Run*){std::cout<<"Neutrons in ch2: "<<my.ch2_counter<<"; in pb: "<<my.lead_counter<<"; in the vessel: "<<my.vessel_counter<<"; in the drift vessel: "<<my.drift_counter<<std::endl;};
	auto my_event_action_begin = [&] (const G4Event* event)
	{
		// Reset variables to zero for the new event.
		my.flag_ch2 = false;
		my.flag_lead = false;
		my.flag_vessel = false;
		my.flag_drift = false;
		//G4ThreeVector event_initial = event -> GetPrimaryVertex() -> GetPosition();
		manualTrackID = 1;
		tempInitialEnergies.clear();
		tempParticleNames.clear();
		tempInitialTimes.clear();
		tempInitialVolumes.clear();
		reducingEnergies.clear();
		reducingMothers.clear();
	};
	auto my_event_action_end = [&] (const G4Event*){n_event++;};
	auto per_step = [&] (const G4Step* step)
	{
		auto current_vol = step -> GetPreStepPoint() -> GetTouchableHandle() -> GetVolume() -> GetLogicalVolume();
		auto sd = current_vol -> GetSensitiveDetector();
		auto current_vol_name = current_vol -> GetName();
		auto en_temp = step -> GetTotalEnergyDeposit();
		auto pre           = step -> GetPreStepPoint();
		auto post           = step -> GetPostStepPoint();
		auto track = step -> GetTrack();
		auto name = track -> GetDefinition() -> GetParticleName();
		auto trk_id = track -> GetTrackID();
		if(name == "neutron")
		{
			if(0){}
			#if SHIELDING
			else if(my.flag_ch2 == false && current_vol_name == "ch2_hollow")
			{
				my.flag_ch2 = true;
				my.ch2_counter++;
			}
			else if(my.flag_lead == false && current_vol_name == "pb_hollow")
			{
				my.flag_lead = true;
				my.lead_counter++;
			}
			#endif
			else if(my.flag_vessel == false && (current_vol_name == "vessel" || current_vol_name == "closure_vessel1" || current_vol_name == "closure_vessel2"))
			{
				my.flag_vessel = true;
				my.vessel_counter++;
			}
			else if(my.flag_drift == false && current_vol_name == "gasdrift")
			{
				my.flag_drift = true;
				my.drift_counter++;
			}
		}
		if(trk_id > manualTrackID)
		{
			G4int skipped = trk_id - manualTrackID;
			std::sort(pendingIdentifiers.begin(), pendingIdentifiers.end(),
				[](const std::string& a, const std::string& b) {
					int a_num = std::stoi(a.substr(a.find('_') + 1));
					int b_num = std::stoi(b.substr(b.find('_') + 1));
					return a_num > b_num;  // descending order
				});
			for (auto& id : pendingIdentifiers)
			{
				identifierToManualTrackID[id] += skipped; // increment all pending identifiers by the number of skipped track IDs
				std::string newID; // = std::to_string(track -> GetParentID()) + "_" + std::to_string(identifierToManualTrackID[id]);
				std::string first = id.substr(0, id.find('_')); // get the parent ID part of the identifier
				int second = std::stoi(id.substr(id.find('_') + 1)); // get the manual track ID part of the identifier
				int newSecond = second + skipped; // increment the manual track ID part by the number of skipped track IDs
				newID = first + "_" + std::to_string(newSecond); // create the new
				std::string oldID = id;
				// Update the maps with the new identifier
				if(tempInitialEnergies.count(oldID))
				{
					tempInitialEnergies[newID] = tempInitialEnergies[oldID];
					tempInitialEnergies.erase(oldID);
				}
				if(tempParticleNames.count(oldID))
				{
					tempParticleNames[newID] = tempParticleNames[oldID];
					tempParticleNames.erase(oldID);
				}
				if(tempInitialTimes.count(oldID))
				{
					tempInitialTimes[newID] = tempInitialTimes[oldID];
					tempInitialTimes.erase(oldID);
				}
				if(tempInitialVolumes.count(oldID))
				{
					tempInitialVolumes[newID] = tempInitialVolumes[oldID];
					tempInitialVolumes.erase(oldID);
				}
			}
			manualTrackID += skipped; // increment the manual track ID by the number of skipped track IDs
		}
		if(trk_id == manualTrackID)
		{
			pendingIdentifiers.clear(); // clear pending identifiers when everything is aligned
			identifierToManualTrackID.clear(); // clear the mapping of identifiers to manual track IDs
		}

		auto parentID = track -> GetParentID();
		auto globalTime = track -> GetGlobalTime();
		auto initial_k_energy = pre  -> GetKineticEnergy();
		auto pos = pre -> GetPosition();
		auto pos_initial = track -> GetVertexPosition(); // gets where it was converted, if done so
		auto pos_interaction = post -> GetPosition();
		auto cell = pre -> GetTouchable() -> GetVolume() -> GetName();
		auto process = post -> GetProcessDefinedStep() -> GetProcessName();
		const G4VProcess* preProcess = pre->GetProcessDefinedStep();
		std::string preProcessName = (preProcess) ? preProcess->GetProcessName() : "LackOfProcess";
		auto k_energy_left = post -> GetKineticEnergy();
		G4ThreeVector p_i = pre -> GetMomentum();
		G4ThreeVector p_f = post -> GetMomentum();
		G4double E_i = pre->GetTotalEnergy(); // gets kinetic and rest mass energy
		G4double E_f = post->GetTotalEnergy(); // E^2 = p^2 + m^2
		G4ThreeVector p_T = p_i - p_f; // Recoil nucleus momentum , because p_i = p_f + p_T
		G4double pT2 = p_T.mag2(); // Squared magnitude of recoil momentum
		// Compute the estimated target mass using relativistic energy conservation
		// E_i + M = E_f + sqrt(pT^2 + M^2)
		G4double numerator = std::pow((E_i - E_f), 2) - pT2;
		G4double denominator = 2 * (E_i - E_f);
		G4double M = numerator / denominator * MeV_to_amu * -1; // relativistic approach, works for all mass target particles
		int target = std::round(M);
		if(sd != nullptr)
		{
			if(process == "hadElastic" && name == "neutron"){output << n_event << " " << trk_id << " " << parentID << " " << name << " " << process << " " << cell << " " << initial_k_energy-k_energy_left << " " << target << " " << globalTime << " " << pos_interaction.x() << " " << pos_interaction.y() << " " << pos_interaction.z() << "\n";}
			if(process != "Transportation" && name == "gamma"){output << n_event << " " << trk_id << " " << parentID << " " << name << " " << process << " " << cell << " " << en_temp << " " << target << " " << globalTime << " " << pos_interaction.x() << " " << pos_interaction.y() << " " << pos_interaction.z() << "\n";}// this accounts for the small amount of energy ( O(1-100s eV) ) tossed into a local deposition by Geant4 energy/distance cut
		}
		//  need to track all secondaries in order to keep manualTrackID in line with true trackID
		auto secondaries = step -> GetSecondaryInCurrentStep();
		if(secondaries -> size() > 0)
		{
			// G4cout << "Secondaries: " << secondaries -> size() << G4endl;
			for(int i=0; i<secondaries -> size(); i++)
			{
				auto secondary = secondaries -> at(i);
				auto particle = secondary -> GetDefinition() -> GetParticleName();
				auto energy_sec = secondary -> GetKineticEnergy();
				auto sec_parentID = secondary -> GetParentID();
				manualTrackID++;
				std::string identifier = std::to_string(sec_parentID) + "_" + std::to_string(manualTrackID);
				pendingIdentifiers.push_back(identifier);
				identifierToManualTrackID[identifier] = manualTrackID;
				if(process == "eBrem" && particle == "gamma")
				{ // i think this only accounst for 1 gamma produced by eBrem, not all of them
					reducingEnergies[identifier] = energy_sec; // store the energy of the gamma produced by eBrem
					std::string mother = std::to_string(parentID) + "_" + std::to_string(trk_id); // since "mother" is not unique for each gamma, would need to store identifiers in an array
					// and then loop through them in the track end-logic to subtract more than one gamma energy
					reducingMothers[mother] = identifier; // store the parent ID of the particle that produced the gamma
				}
				// "Scintillation" is a process that produces the decay products when using the IonTable
				if(process == "neutronInelastic" || process == "neutronCapture" || name == "gamma" || process == "Scintillation")
				{
					tempInitialEnergies[identifier] = energy_sec;
					tempParticleNames[identifier] = particle;
					tempInitialTimes[identifier] = globalTime;
					tempInitialVolumes[identifier] = cell;
				}
			}
		}

		std::string identifier = std::to_string(parentID) + "_" + std::to_string(trk_id);
		// std::cout << "Identifier: " << identifier << ", manualTrackID: " << manualTrackID << ", pendingIdentifiers size: " << pendingIdentifiers.size() << ", is ID saved? " << std::to_string(tempInitialEnergies.find(identifier) != tempInitialEnergies.end()) << std::endl;
		if ((tempInitialEnergies.find(identifier) != tempInitialEnergies.end() && track -> GetTrackStatus() == fStopAndKill && name != "neutron" && name != "gamma") || (tempInitialEnergies.find(identifier) != tempInitialEnergies.end() && process == "Transportation" && name != "neutron" && name != "gamma"))
		{
			auto sec_energy = tempInitialEnergies[identifier];
			auto sec_particle = tempParticleNames[identifier];
			auto sec_time = tempInitialTimes[identifier];
			auto sec_cell = tempInitialVolumes[identifier];

			if(reducingMothers.find(identifier) != reducingMothers.end())
			{
			    // if (n_event == 7416) {
			    //     G4cout << "Reducing energy for particle " << name << " with track ID " << trk_id << " and parent ID " << parentID << " and identifier " << identifier << " by gamma energy produced by eBrem: " << reducingEnergies[reducingMothers[identifier]] << G4endl;
			    // }
			    sec_energy -= reducingEnergies[reducingMothers[identifier]]; // reduce the energy by the gamma produced by eBrem
			    reducingEnergies.erase(reducingMothers[identifier]); // remove the gamma energy so it doesn't affect other particles
			    reducingMothers.erase(identifier);
			}
			auto sec_final_energy = track -> GetKineticEnergy();
			G4String sec_process; // = pre -> GetProcessDefinedStep() -> GetProcessName();
			// if you see "Transportation" as the process for something like an e-, then it left the volume before the end of its track, not that it mistakenly
			// left energy via that process
			if (process == "Transportation" || preProcessName == "LackOfProcess" ){sec_process = process;}
			else{sec_process = pre -> GetProcessDefinedStep() -> GetProcessName();}
			if (sec_energy != sec_final_energy && sd != nullptr){output << n_event << " " << trk_id << " " << parentID << " " << name << " " << sec_process << " " << sec_cell << " " << sec_energy-sec_final_energy << " " << "0" << " " << sec_time << " " << pos_interaction.x() << " " << pos_interaction.y() << " " << pos_interaction.z() << "\n";}
			if (sec_energy-sec_final_energy < 0){std::cout << "Warning: negative energy deposition detected! event " << n_event << ", energy " << sec_energy-sec_final_energy << " MeV in " << sec_cell << " for particle " << name << " with track ID " << trk_id << " and parent ID " << parentID << " and identifier " << identifier << std::endl;}
			if (process == "Transportation")
			{
				tempInitialVolumes[identifier] = cell;
				tempInitialEnergies[identifier] = sec_final_energy;
				tempInitialTimes[identifier] = globalTime;
			}
			else
			{
				tempInitialEnergies.erase(identifier);
				tempParticleNames.erase(identifier);
				tempInitialTimes.erase(identifier);
				tempInitialVolumes.erase(identifier);
			}
		}
	};

	return (new n4::        actions{my_generator(my)})
		-> set((new n4::   run_action{}) -> end(my_run_action_end))
		-> set((new n4::   event_action{}) -> begin(my_event_action_begin)
														-> end(my_event_action_end))
		-> set( new n4::stepping_action{per_step});
}

auto my_geometry(const my& my, unsigned& n_event)
{
	G4double vessel_cylinder_length = 1110 * mm;
	G4double vessel_diam_ext = 665.2 * mm;
	G4double vessel_diam_int = 640 * mm;
	G4double closure_vessel_diam_ext = 820 * mm;
	G4double closure_vessel_length = 100 * mm;
	G4double closure_vessel_diam_int = 560 * mm;
	//auto closure_vessel_z = vessel_cylinder_length/2 + closure_vessel_length/2;
	G4double closure_vessel_z = 475 * mm;
	G4double disc_cath_diam = 600.4 * mm;
	G4double disc_cath_length = 4 * mm;
	G4double disc_cath_z = 0 * mm;
	G4double cont_EL_diam = 601.4 * mm;
	G4double cont_EL_length = 10.5 * mm;
	G4double cont_EL_z = disc_cath_z + 12.5 * mm + cont_EL_length/2;
	G4double EL_diam = 525.4 * mm;
	G4double EL_length = 0.5 * mm;
	G4double EL_z = cont_EL_z - 4.75 * mm;
	//
	G4double tube_teflon_diam_int = 570.618 * mm;
	G4double tube_teflon_diam_ext = tube_teflon_diam_int + 2 * 3.453 * mm;
	G4double tube_teflon_length = 400 * mm;
	G4double tube_teflon_z = cont_EL_z + tube_teflon_length/2 + 5 * mm;
	//
	G4double gasdrift_diam =  tube_teflon_diam_int;
	G4double gasdrift_length = 2 * tube_teflon_length + 4 * 5 * mm;
	G4double gasdrift_z = disc_cath_z;
	//
	G4double ring_diam_int = 582 * mm;
	G4double ring_diam_ext = 588 * mm;
	G4double ring_length = 7.01 * mm;
	G4double ring_z = disc_cath_z + 418.595 * mm + ring_length/2;

	G4double world_rad = gasdrift_diam/2 + 180 * cm;

	auto water  = n4::material("G4_WATER");
	auto air    = n4::material("G4_AIR");
  
	auto vaccum = n4::material("G4_Galactic");

	auto steel  = n4::material("G4_STAINLESS-STEEL");
	auto copper = n4::material("G4_Cu");
	auto teflon = n4::material("G4_TEFLON");

	auto lead = n4::material("G4_Pb");
	auto ch2 = n4::material("G4_POLYETHYLENE");
 
	G4double pressure = 20*bar; 
	G4double temperature = 293.*kelvin; 
	G4double R =  8.3145*joule/(mole*kelvin);
  
	G4Material* gas = nullptr;
	G4double M,density;
	G4Element *Ar,*Xe;
	#if XE
		G4String gas_sel = "Xe";
		M = 131.29*g/mole;
		density = (pressure * M) / (R * temperature);
		Xe = new G4Element("Xenon","Xe", 54., M);
		gas = new G4Material("Xenon", density, 1, kStateGas, temperature, pressure); 
		gas->AddElement(Xe, 1);
		//std::cout << "M is : " << M << " Xe is : " << Xe << " density is : " << density << " gas is : " << gas << std::endl;
	#else
		G4String gas_sel = "Ar";
		M =  39.9480*g/mole;
		density = (pressure * M) / (R * temperature);
		Ar = new G4Element("Argon","Ar", 18., M);
		gas = new G4Material("Argon", density, 1, kStateGas, temperature, pressure); 
		gas->AddElement(Ar, 1);
		//std::cout << "M is : " << M << " Ar is : " << Ar << " density is : " << density << " gas is : " << gas << std::endl;
	#endif
	auto world  = n4::sphere("world").r(world_rad).volume(air);

	////////////////////
	//Shielding layers//
	////////////////////
	G4double leng_pb = vessel_cylinder_length + closure_vessel_length * 2 + thick_gap * 2;
	G4double wid_pb = closure_vessel_diam_ext + thick_gap * 2;
	G4double leng_ch2 = leng_pb + thick_pb * 2;
	G4double wid_ch2 = wid_pb + thick_pb * 2;
	auto pb_inner_solid = new G4Box("pb_inner", wid_pb/2, wid_pb/2, leng_pb/2);
	auto pb_outer_solid = new G4Box("pb_outer", wid_pb/2 + thick_pb, wid_pb/2 + thick_pb, leng_pb/2 + thick_pb);
	G4SubtractionSolid* pb_hollow = new G4SubtractionSolid("pb_hollow", pb_outer_solid, pb_inner_solid);
	auto pb_logical = new G4LogicalVolume(pb_hollow, lead, "pb_hollow");
	#if SHIELDING
	new G4PVPlacement(nullptr, {}, pb_logical, "pb_hollow", world, false, 0,true);
	#endif
	auto ch2_inner_solid = new G4Box("ch2_inner", wid_ch2/2, wid_ch2/2, leng_ch2/2);
	auto ch2_outer_solid = new G4Box("ch2_outer", wid_ch2/2 + thick_ch2, wid_ch2/2 + thick_ch2, leng_ch2/2 + thick_ch2);
	G4SubtractionSolid* ch2_hollow = new G4SubtractionSolid("ch2_hollow", ch2_outer_solid, ch2_inner_solid);
	auto ch2_logical = new G4LogicalVolume(ch2_hollow, ch2, "ch2_hollow");
	#if SHIELDING
	new G4PVPlacement(nullptr, {}, ch2_logical, "ch2_hollow", world, false, 0,true);
	#endif

	////////////
	//Detector//
	////////////
	auto vessel = n4::tubs  ("vessel" ).r_inner(vessel_diam_int/2).r_delta(vessel_diam_ext/2 - vessel_diam_int/2).z(vessel_cylinder_length).place(steel).at(0, 0, disc_cath_z).in(world).now();
	auto closure_vessel1 = n4::tubs  ("closure_vessel1" ).r(closure_vessel_diam_ext/2).z(closure_vessel_length).place(steel).at(0, 0, closure_vessel_z).in(world).now();
	auto closure_vessel2 = n4::tubs  ("closure_vessel2" ).r(closure_vessel_diam_ext/2).z(closure_vessel_length).place(steel).at(0, 0, - closure_vessel_z).in(world).now();

	auto disc_cathode = n4::tubs  ("disc_cathode" ).r(disc_cath_diam/2).z(disc_cath_length).place(steel).at(0, 0, disc_cath_z).in(world).now();
  
	auto EL1 = n4::tubs  ("EL1" ).r(EL_diam/2).z(EL_length).place(steel).at(0, 0, EL_z).in(world).now();
	auto EL2 = n4::tubs  ("EL2" ).r(EL_diam/2).z(EL_length).place(steel).at(0, 0, - EL_z).in(world).now();

	auto cont_EL1 = n4::tubs  ("cont_EL1" ).r_inner(EL_diam/2).r_delta(cont_EL_diam/2 - EL_diam/2).z(cont_EL_length).place(steel).at(0, 0, cont_EL_z).in(world).now();
	auto cont_EL2 = n4::tubs  ("cont_EL2" ).r_inner(EL_diam/2).r_delta(cont_EL_diam/2 - EL_diam/2).z(cont_EL_length).place(steel).at(0, 0, - cont_EL_z).in(world).now();

	auto tube_teflon1 = n4::tubs  ("tube_teflon1" ).r_inner(tube_teflon_diam_int/2).r_delta(tube_teflon_diam_ext/2 - tube_teflon_diam_int/2).z(tube_teflon_length).place(teflon).at(0, 0, tube_teflon_z).in(world).now();
	auto tube_teflon2 = n4::tubs  ("tube_teflon2" ).r_inner(tube_teflon_diam_int/2).r_delta(tube_teflon_diam_ext/2 - tube_teflon_diam_int/2).z(tube_teflon_length).place(teflon).at(0, 0,  - tube_teflon_z).in(world).now();

	for (G4int i = 0; i < 27; i++)
	{
		G4double ring_delta = 7.99 * mm + ring_length;
		G4double ring_space = ring_delta * i;
		auto ring1 = n4::tubs  ("ring1" ).r_inner(ring_diam_int/2).r_delta(ring_diam_ext/2 - ring_diam_int/2).z(ring_length).place(copper).at(0, 0, ring_z - ring_space).in(world).now();
		auto ring2 = n4::tubs  ("ring2" ).r_inner(ring_diam_int/2).r_delta(ring_diam_ext/2 - ring_diam_int/2).z(ring_length).place(copper).at(0, 0, -(ring_z - ring_space)).in(world).now();
	}
	auto gasdrift = n4::tubs  ("gasdrift" ).r(gasdrift_diam/2).z(gasdrift_length).sensitive(sensitive_detector(n_event)).place(gas).in(world).now();

	return n4::place(world).now();
}

int main(int argc, char* argv[])
{
	unsigned n_event = 0;
	#if VISUAL
	int argc_first_three=argc;
	#else
	int argc_first_three=3;
	if(argc>3)
	{
   	std::ostringstream filename_stream;
   	filename_stream<<"results/output_"<<argv[3]<<".txt";
   	std::string filename = filename_stream.str();
   	output.open (filename);
	}
	else
	{
		output.open ("output.txt");
	}
	#endif
	//inFile.open("flux/flux_mark.txt");
	inFile.open("flux/flux_mark_uniform.txt");
	my my;
	G4int physics_verbosity = 0;
	int ene_count[2000] = {0};
	// The trailing slash after '/my_geometry' is CRUCIAL: without it, the
	// messenger violates the principle of least surprise.
	auto messenger = new G4GenericMessenger{nullptr, "/my/", "docs: bla bla bla"};
	messenger -> DeclarePropertyWithUnit("msg_pressure"      , "bar"  , my.msg_pressure  );
	/*messenger -> DeclarePropertyWithUnit("bubble_radius"     , "m"  , my.bubble_radius  );
	messenger -> DeclarePropertyWithUnit("socket_rot"        , "deg", my.socket_rot     );
	messenger -> DeclarePropertyWithUnit("particle_energy"   , "keV", my.particle_energy);
	messenger -> DeclareProperty        ("particle"          ,        my.particle_name  );
	messenger -> DeclareProperty        ("particle_direction",        my.particle_dir   );
	messenger -> DeclareProperty        ("physics_verbosity" ,        physics_verbosity );*/
	n4::run_manager::create()
		.ui("neutron_bg", argc_first_three, argv)
		.macro_path("macs")
		.apply_early_macro("early-hard-wired.mac")
		.physics<FTFP_BERT>(physics_verbosity)
		.geometry([&] { return my_geometry(my,n_event); })
		.actions(create_actions(my, n_event,ene_count))
		.run();
	std::cout << "Number of events is : " << n_event << std::endl;
	// Important! physics list has to be set before the generator!
}
