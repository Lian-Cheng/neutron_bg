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
#include <G4Tubs.hh>             // for creating shapes in the geometry
#include <G4SubtractionSolid.hh>// for boolean solids
#include <G4UnionSolid.hh>      //
#include <G4MultiUnion.hh>      //
#include <G4Sphere.hh>          // for creating shapes in the geometry
#include <FTFP_BERT.hh>         // our choice of physics list
#include <G4RandomDirection.hh> // for launching particles in random directions
#include <Randomize.hh>
#include <G4ThreeVector.hh>
#include <G4OpticalSurface.hh>
#include <G4LogicalSkinSurface.hh>
#include <cstdlib>
#define VISUAL 1
#define SHIELDING 0
#define XE 0
#define DIR_CENTERED 0

using namespace CLHEP;

thread_local std::ofstream output;
thread_local std::ifstream inFile;
thread_local std::vector<double> flux;
thread_local std::vector<double> energy;

thread_local std::map<std::string, G4double> tempInitialEnergies;
thread_local std::map<std::string, G4double> tempInitialTimes;
thread_local std::map<std::string, G4String> tempParticleNames;
thread_local std::map<std::string, G4String> tempInitialVolumes;
thread_local std::map<G4int, G4double> secondaryInitialEnergies;
thread_local std::map<G4int, G4String> secondaryParticleNames;
thread_local std::map<std::string, G4double> reducingEnergies;
thread_local std::map<std::string, G4String> reducingMothers;

const G4double MeV_to_amu = 1/931.5;
const G4double thick_ch2 = 100. * mm;
const G4double thick_pb = 150. * mm;
const G4double thick_gap = 10. * mm;
const G4double length_detector = 1436. * mm;
const G4double y_top_detector = 459. * mm;
const G4double y_bot_detector = -1498. * mm;
const G4double x_detector = 900. * mm;

thread_local std::vector<std::string> pendingIdentifiers;  // identifiers created in this step
thread_local std::map<std::string, int> identifierToManualTrackID;  // temp mapping

thread_local G4int manualTrackID;

thread_local struct my
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
		//////////////////////////////////////////////////
		// Random distribution for positions and angles //
		//////////////////////////////////////////////////
		//Surface of the Ch2 shielding
		std::uniform_real_distribution<G4double> dist0(0.0, 1.0);
		#if SHIELDING
			std::uniform_real_distribution<G4double> dist1(0.0, 1.0);
			std::uniform_real_distribution<G4double> dist2(0.0, 1.0);
			G4double src_x = x_detector + thick_gap * 2. + thick_pb * 2. + thick_ch2 * 2.;// width direction
			G4double src_y = (y_top_detector - y_bot_detector) + thick_gap * 2. + thick_pb * 2. + thick_ch2 * 2.;
			G4double src_z = length_detector + thick_gap * 2. + thick_pb * 2. + thick_ch2 * 2.;// length direction
			G4double area_1 = src_x * src_y;
			G4double area_2 = src_y * src_z;
			G4double area_3 = src_z * src_x;
			G4double area_tot = 2. * (area_1 + area_2 + area_3);//total surface
		#else
			std::uniform_real_distribution<G4double> dist1(0.0, 1.0);
			std::uniform_real_distribution<G4double> alpha(0.0, 2. * M_PI);
			G4double diameter_xy = x_detector;
			G4double src_z = length_detector;//length direction
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

		//////////////////////////////////////////////////////////////////
		// Position and momentum direction distribution of neutron flux //
		//////////////////////////////////////////////////////////////////
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

		/////////////////////////////////////////
		// Energy distribution of neutron flux //
		/////////////////////////////////////////
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

		////////////////////////////////
		// Generation of neutron flux //
		////////////////////////////////
		auto neutron = n4::find_particle("neutron");
		auto particle_gun = std::make_shared<G4ParticleGun>(neutron, 1); // 1 particle
		//G4double mass = neutron -> GetPDGMass();
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
	/////////////////////////
	// Material definition //
	/////////////////////////
	auto air = n4::material("G4_AIR");
	//auto water = n4::material("G4_WATER");
	auto vacuum = n4::material("G4_Galactic");
	auto steel = n4::material("G4_STAINLESS-STEEL");
	auto copper = n4::material("G4_Cu");
	auto teflon = n4::material("G4_TEFLON");
	auto lead = n4::material("G4_Pb");
	auto ch2 = n4::material("G4_POLYETHYLENE");
	G4NistManager* nist = G4NistManager::Instance();
	G4Material* silica = G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
	G4Material* aluminum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
	G4Material* sapphire = G4NistManager::Instance()->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
	// Kovar
	G4Material* Fe = nist->FindOrBuildMaterial("G4_Fe");
	G4Material* Ni = nist->FindOrBuildMaterial("G4_Ni");
	G4Material* Co = nist->FindOrBuildMaterial("G4_Co");
	G4Material* Kovar = new G4Material("Kovar", 8.35 * g/cm3, 3);
	Kovar->AddMaterial(Fe, 54.0 * perCent);
	Kovar->AddMaterial(Ni, 29.0 * perCent);
	Kovar->AddMaterial(Co, 17.0 * perCent);
	// TPB
	G4Element* H = nist->FindOrBuildElement("H");
	G4Element* C = nist->FindOrBuildElement("C");
	G4Material* tpb = new G4Material("TPB", 1*g/cm3, 2, kStateSolid);
	tpb->AddElement(H, 22);
	tpb->AddElement(C, 28);

	////////////////////////
	// TPC gas definition //
	////////////////////////
	G4Material* gas = nullptr;
	G4double pressure = 20*bar; 
	G4double temperature = 293.*kelvin; 
	G4double R =  8.3145*joule/(mole*kelvin);
	G4double M_gas,density;
	G4Element *Ar,*Xe;
	#if XE
		G4String gas_sel = "Xe";
		M_gas = 131.29*g/mole;
		density = (pressure * M_gas) / (R * temperature);
		Xe = new G4Element("Xenon","Xe", 54., M_gas);
		gas = new G4Material("Xenon", density, 1, kStateGas, temperature, pressure); 
		gas->AddElement(Xe, 1);
		//std::cout << "M_gas is : " << M_gas << " Xe is : " << Xe << " density is : " << density << " gas is : " << gas << std::endl;
	#else
		G4String gas_sel = "Ar";
		M_gas =  39.9480*g/mole;
		density = (pressure * M_gas) / (R * temperature);
		Ar = new G4Element("Argon","Ar", 18., M_gas);
		gas = new G4Material("Argon", density, 1, kStateGas, temperature, pressure); 
		gas->AddElement(Ar, 1);
		//std::cout << "M_gas is : " << M_gas << " Ar is : " << Ar << " density is : " << density << " gas is : " << gas << std::endl;
	#endif

	////////////////
	// Geometries //
	////////////////
	#if SHIELDING
	G4double world_rad = length_detector/2. + thick_gap/2. + thick_pb/2. + thick_ch2/2. + 130. * cm;
	#else
	G4double world_rad = length_detector/2. + 10. * cm;
	#endif
	auto world = n4::sphere("world").r(world_rad).volume(air);

	//////////////////////
	// Shielding layers //
	//////////////////////
	#if SHIELDING
	// Dimensions
	G4double x_air = x_detector + thick_gap * 2;
	G4double y_air = (y_top_detector - y_bot_detector) + thick_gap * 2;
	G4double z_air = length_detector + thick_gap * 2;
	G4double x_pb = x_air + thick_pb * 2;
	G4double y_pb = y_air + thick_pb * 2;
	G4double z_pb = z_air + thick_pb * 2;
	G4double x_ch2 = x_pb + thick_ch2 * 2;
	G4double y_ch2 = y_pb + thick_ch2 * 2;
	G4double z_ch2 = z_pb + thick_ch2 * 2;
	// Polyethylene
	auto ch2_solid = new G4Box("ch2_layer_solid", x_ch2/2., y_ch2/2., z_ch2/2.);
	auto ch2_logical = new G4LogicalVolume(ch2_solid, ch2, "ch2_layer_logic");
	new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, 0.0), ch2_logical, "ch2_layer_logic", world, false, 0,true);
	// Lead
	auto pb_solid = new G4Box("pb_layer_solid", x_pb/2., y_pb/2., z_pb/2.);
	auto pb_logical = new G4LogicalVolume(pb_solid, lead, "pb_layer_logic");
	new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, 0.0), pb_logical, "pb_layer_logic", ch2_logical, false, 0,true);
	// Air
	auto air_solid = new G4Box("air_layer_solid", x_air/2., y_air/2., z_air/2.);
	auto air_logical = new G4LogicalVolume(air_solid, air, "air_layer_logic");
	new G4PVPlacement(nullptr, G4ThreeVector(0.0, 0.0, 0.0), air_logical, "air_layer_logic", pb_logical, false, 0,true);
	#endif

	////////////
	// Vessel //
	////////////
	// Dimensions
	G4double cylinder_length = 850. * mm;
	G4double cylinder_outer_diameter = 664. * mm;
	G4double cylinder_inner_diameter = 640. * mm;
	G4double flange_thick = 100. * mm;
	G4double flange_inner_thick = 45. * mm;
	G4double flange_outer_diameter = 820. * mm;
	G4double flange_inner_diameter = 605. * mm;
	G4double endcap_base = 85.94 * mm; // 35 mm for the cylinder and 50.94 mm for the approximation of the 55.73 deg round part
	G4double endcap_inner_rad = 532. * mm; //=(665.2/2 - 12.6)
	G4double endcap_theta = 37.2 * deg; //(38.7 * deg)- 1.5º visfit
	G4double endcap_thick = 12.6 * mm;
	G4double endcap_outer_rad = endcap_inner_rad + endcap_thick;
	// Cylinder solid
	auto cylinder_solid = new G4Tubs("vessel_cylinder_solid", 0., cylinder_outer_diameter/2., cylinder_length/2. + flange_thick + endcap_base, 0.,twopi);
	// Flange solid
	G4Tubs* flange_solid = new G4Tubs("vessel_flange", cylinder_outer_diameter/2., flange_outer_diameter/2., flange_thick/2., 0., twopi);
	// Sphere solid
	G4Sphere* endcap1_solid = new G4Sphere("vessel_endcap1",
						0. * cm,  endcap_outer_rad,   //radius
						0. * deg, 360. * deg,       // phi
						0. * deg, endcap_theta);   // theta
	G4Sphere* endcap2_solid = new G4Sphere("vessel_endcap2",
						0. * cm,  endcap_outer_rad,     //radius
						0. * deg, 360. * deg,        // phi
						180. * deg - endcap_theta, endcap_theta); // theta
	// Union
	G4double endcap_in_z_width = 96.48 * mm;// inner sphere to flange (96,28)-visfit (1.8mm)
	G4double endcap_z_pos = (cylinder_length/2. + flange_thick + endcap_base) - (endcap_inner_rad - endcap_in_z_width);
	G4ThreeVector endcap1_pos(0, 0, endcap_z_pos);
	G4ThreeVector endcap2_pos(0, 0, -1. * endcap_z_pos);
	G4double flange_z_pos = cylinder_length/2. + flange_thick/2.;
	G4ThreeVector flange1_pos(0, 0, flange_z_pos);
	G4ThreeVector flange2_pos(0, 0, -1. * flange_z_pos);
	// Body + 2 flanges
	G4UnionSolid* vessel_solid = new G4UnionSolid("vessel_assembly_solid", cylinder_solid, flange_solid, 0,flange1_pos);
	vessel_solid = new G4UnionSolid("vessel_assembly_solid", vessel_solid, flange_solid, 0, flange2_pos);
	// Body + 2 flanges + 2 endcaps
	vessel_solid = new G4UnionSolid("vessel_assembly_solid", vessel_solid, endcap1_solid, 0, endcap1_pos);
	vessel_solid = new G4UnionSolid("vessel_assembly_solid", vessel_solid, endcap2_solid, 0, endcap2_pos);
	G4LogicalVolume* vessel_logic = new G4LogicalVolume(vessel_solid,	steel, "vessel_logic");
	#if SHIELDING
	G4double y_shield_transit = (y_top_detector + y_bot_detector)/2;
	G4ThreeVector trans_shield(0., -y_shield_transit, 0.);
	new G4PVPlacement(nullptr, trans_shield, vessel_logic, "vessel_logic", air_logical, false, 0, true);
	#else
	new G4PVPlacement(nullptr, G4ThreeVector(0., 0., 0.), vessel_logic, "vessel_logic", world, false, 0, true);
	#endif

	///////////////
	// Gas solid //
	///////////////
	// Gas tube
	G4double gas_tube_length = 960. * mm;
	G4Tubs* gas_tube_solid = new G4Tubs("gas_tube_solid", 0., cylinder_inner_diameter/2., gas_tube_length/2., 0., twopi);
	G4Tubs* gas_flange_solid = new G4Tubs("gas_flange_solid", 0., flange_inner_diameter/2., flange_inner_thick/2., 0., twopi);
	G4Tubs* gas_base_solid = new G4Tubs("gas_base_solid", 0., cylinder_inner_diameter/2., endcap_base/2., 0., twopi);
	// Gas endcaps
	G4double endcap_theta2 = 35.2 * deg; //(38.7 * deg)- 1.5º visfit
	G4Sphere* gas_endcap1_solid = new G4Sphere("vessel_gas_endcap1",
						0. * cm,  endcap_inner_rad,   //radius
						0. * deg, 360. * deg,       // phi
						0. * deg, endcap_theta2);   // theta
	G4Sphere* gas_endcap2_solid = new G4Sphere("vessel_gas_endcap2",
						0. * cm,  endcap_inner_rad,    //radius
						0. * deg, 360. * deg,        // phi
						180. * deg - endcap_theta2, endcap_theta2); // theta
	// Gas union
	G4double endcap_base_z_pos = cylinder_length/2. + flange_thick + endcap_base/2.;
	G4ThreeVector endcap_base1_pos(0, 0, endcap_base_z_pos);
	G4ThreeVector endcap_base2_pos(0, 0, -1. * endcap_base_z_pos);
	G4double flange_inner_z_pos = cylinder_length/2. + flange_thick - flange_inner_thick/2.;
	G4ThreeVector flange_inner1_pos(0, 0, flange_inner_z_pos);
	G4ThreeVector flange_inner2_pos(0, 0, -1. * flange_inner_z_pos);
	// Body
	G4UnionSolid* vessel_gas_solid = new G4UnionSolid("vessel_assembly_gas", gas_tube_solid, gas_flange_solid, 0, flange_inner1_pos);
	vessel_gas_solid = new G4UnionSolid("vessel_assembly_gas", vessel_gas_solid, gas_flange_solid, 0, flange_inner2_pos);
	vessel_gas_solid = new G4UnionSolid("vessel_assembly_gas", vessel_gas_solid, gas_base_solid, 0, endcap_base1_pos);
	vessel_gas_solid = new G4UnionSolid("vessel_assembly_gas", vessel_gas_solid, gas_base_solid, 0, endcap_base2_pos);
	// Body gas + 2 endcaps gas
	vessel_gas_solid = new G4UnionSolid("vessel_assembly_gas", vessel_gas_solid, gas_endcap1_solid, 0, endcap1_pos);
	vessel_gas_solid = new G4UnionSolid("vessel_assembly_gas", vessel_gas_solid, gas_endcap2_solid, 0, endcap2_pos);
	G4LogicalVolume* gas_logic = new G4LogicalVolume(vessel_gas_solid, gas, "gas_logic");
	new G4PVPlacement(nullptr, G4ThreeVector(0.,0.,0.), gas_logic, "gas_logic", vessel_logic, false, 0, true);

	/////////////
	// Cathode //
	/////////////
	G4double cathode_diam = 600.4 * mm;
	G4double cathode_thick = 4 * mm;
	auto cathode = n4::tubs("cathode").r(cathode_diam/2).z(cathode_thick).place(steel).at(0., 0., 0.).in(gas_logic).now();

	////////////
	// Teflon //
	////////////
	G4double teflon_inner_diam = 567.8 * mm; // equivalent area to the surface;
	G4double teflon_outer_diam = teflon_inner_diam + 2 * 5. * mm;
	G4double teflon_length = 370. * mm;
	G4double teflon_z_pos = cathode_thick/2 + teflon_length/2 + 5.2 * mm;
	auto plate_teflon = n4::tubs("teflon").r_inner(teflon_inner_diam/2).r_delta(teflon_outer_diam/2 - teflon_inner_diam/2).z(teflon_length).place(teflon).in(gas_logic);
	plate_teflon.clone().at(0, 0, teflon_z_pos).now();
	plate_teflon.clone().at(0, 0, -teflon_z_pos).now();

	////////////////
	// Field cage //
	////////////////
	G4double ring_inner_diam = 582 * mm;
	G4double ring_outer_diam = 588 * mm;
	G4double ring_length = 8. * mm;
	G4double ring_gap = 7 * mm;
	G4double ring_offset_z = cathode_thick/2 + (5.3 + 1) * mm + ring_length/2;
	auto ring_FC = n4::tubs("ring").r_inner(ring_inner_diam/2).r_delta(ring_outer_diam/2 - ring_inner_diam/2).z(ring_length).place(copper).in(gas_logic);
	for (G4int i = 0; i < 25; i++)
	{
		G4double ring_z_pos = ring_offset_z + (ring_gap + ring_length)* i;
		ring_FC.clone().at(0, 0, ring_z_pos).now();
		ring_FC.clone().at(0, 0, -ring_z_pos).now();
	}

	////////
	// EL //
	////////
	G4double support_EL_diam = 601.4 * mm;
	G4double EL_diam = 541.4 * mm;
	G4double support_anode_diam = 632.4 * mm;
	G4double anode_diam = 573.4 * mm;
	G4double support_thick = 10.0 * mm;
	G4double electrode_thick = 0.5 * mm;
	G4double gasdrift_length1 = 380.3 * mm;
	G4double support_EL_z_pos = gasdrift_length1 + cathode_thick/2 + support_thick/2;
	G4double EL_z_pos = gasdrift_length1 + cathode_thick/2 + support_thick - electrode_thick/2;
	G4double anode_z_pos = support_EL_z_pos + support_thick/2 + 10. * mm + electrode_thick/2;
	G4double support_anode_z_pos = support_EL_z_pos + support_thick/2 + 10. * mm + support_thick/2;
	auto ring_support_EL = n4::tubs("support_EL").r_inner(EL_diam/2).r_delta(support_EL_diam/2 - EL_diam/2).z(support_thick).place(steel).in(gas_logic);
	ring_support_EL.clone().at(0, 0, support_EL_z_pos).now();
	ring_support_EL.clone().at(0, 0, -support_EL_z_pos).now();
	auto plate_EL = n4::tubs("EL").r(EL_diam/2).z(electrode_thick).place(steel).in(gas_logic);
	plate_EL.clone().at(0, 0, EL_z_pos).now();
	plate_EL.clone().at(0, 0, -EL_z_pos).now();
	auto ring_support_anode = n4::tubs("support_anode").r_inner(anode_diam/2).r_delta(support_anode_diam/2 - anode_diam/2).z(support_thick).place(steel).in(gas_logic);
	ring_support_anode.clone().at(0, 0, support_anode_z_pos).now();
	ring_support_anode.clone().at(0, 0, -support_anode_z_pos).now();
	auto plate_anode = n4::tubs("anode").r(anode_diam/2).z(electrode_thick).place(steel).in(gas_logic);
	plate_anode.clone().at(0, 0, anode_z_pos).now();
	plate_anode.clone().at(0, 0, -anode_z_pos).now();

	//////////////////
	// Energy plane //
	//////////////////
	G4int num_PMTs_each = 19;
	G4double pos_x_PMTs[19] = {-118.5, 0, 118.5, -177.75, -59.25, 59.25, 177.75, -237, -118.5, 0, 118.5, 237, -177.75, -59.25, 59.25, 177.75, -118.5, 0, 118.5};
	G4double pos_y_PMTs[19] = {-205.248, -205.248, -205.248, -102.624, -102.624, -102.624, -102.624, 0, 0, 0, 0, 0, 102.624, 102.624, 102.624, 102.624, 205.248, 205.248, 205.248};
	// Plate
	G4double copper_plate_thickn_ = 38 * mm;
	G4double copper_plate_diam_ = 630. * mm;
	G4Tubs* copper_plate1_origin_solid = new G4Tubs("COPPER_PLATE1_ORIGIN", 0., copper_plate_diam_/2., copper_plate_thickn_/2., 0., twopi);
	G4Tubs* copper_plate2_origin_solid = new G4Tubs("COPPER_PLATE2_ORIGIN", 0., copper_plate_diam_/2., copper_plate_thickn_/2., 0., twopi);
	// Huts
	G4double hut_base_diam = 116. * mm;
	G4double hut_base_length = 8. * mm;
	G4double hut_cylinder_diam = 78. * mm; // hut external diam
	G4double hut_cylinder_length = 130. * mm; // hut external diam
	G4Tubs* hut_base_solid = new G4Tubs("hut_base", 0., hut_base_diam/2., hut_base_length/2., 0., twopi);
	G4Tubs* hut_cylinder_solid = new G4Tubs("hut_cylinder", 0., hut_cylinder_diam/2., hut_cylinder_length/2., 0., twopi);
	G4double hut_cylinder_z_pos = hut_cylinder_length/2. + hut_base_length/2.;
	G4ThreeVector hut_cylinder_pos(0., 0., hut_cylinder_z_pos);
	G4UnionSolid* hut1_solid = new G4UnionSolid("hut1_assembly", hut_base_solid, hut_cylinder_solid, 0, hut_cylinder_pos);
	G4UnionSolid* hut2_solid = new G4UnionSolid("hut2_assembly", hut_base_solid, hut_cylinder_solid, 0, -hut_cylinder_pos);
	// Huts on the plane
	G4double plate_hut_z_pos = hut_base_length/2. + copper_plate_thickn_/2.;
	G4UnionSolid* copper_plate1_hut_solid = new G4UnionSolid("EP_COPPER_PLATE1_hut0", copper_plate1_origin_solid, hut1_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], plate_hut_z_pos));
	G4UnionSolid* copper_plate2_hut_solid = new G4UnionSolid("EP_COPPER_PLATE2_hut0", copper_plate2_origin_solid, hut2_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], -plate_hut_z_pos));
	for(G4int pos_i=1;pos_i<num_PMTs_each;pos_i++)
	{
		copper_plate1_hut_solid = new G4UnionSolid("EP_COPPER_PLATE1", copper_plate1_hut_solid, hut1_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], plate_hut_z_pos));
		copper_plate2_hut_solid = new G4UnionSolid("EP_COPPER_PLATE2", copper_plate2_hut_solid, hut2_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], -plate_hut_z_pos));
	}
	// PMT cover
	G4double gap_win_pmt = 0.1 * mm;
	G4double cover_base_length = gap_win_pmt + 22.8 * mm; // 19 + 3.8
	G4double cover_base_diam = 118. * mm;
	G4double cover_win_length = gap_win_pmt + 3.8 * mm; // 8 - (19 - 14.8)
	G4double cover_win_diam = 96.3 * mm;
	G4double cover_rubber_length = 3.6 * mm;
	G4double cover_rubber_inner_diam = 80. * mm;
	G4Tubs* cover_base_solid = new G4Tubs("cover_base", 0., cover_base_diam/2., cover_base_length/2., 0., twopi);
	G4Tubs* cover_win_solid = new G4Tubs("cover_win", 0., cover_win_diam/2., cover_win_length/2., 0., twopi);
	G4Tubs* cover_rubber_solid = new G4Tubs("cover_rubber", cover_rubber_inner_diam/2., cover_win_diam/2., cover_rubber_length/2., 0., twopi);
	// Covers on the plane
	G4double cover_base_z_pos = cover_base_length/2. + copper_plate_thickn_/2.;
	G4double cover_win_z_pos = cover_win_length/2. + cover_base_length + copper_plate_thickn_/2.;
	G4double cover_rubber_z_pos = cover_rubber_length/2. + cover_win_length + cover_base_length + copper_plate_thickn_/2.;
	G4UnionSolid* copper_plate1_cover_solid = new G4UnionSolid("EP_COPPER_PLATE1_cover0", copper_plate1_hut_solid, cover_base_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], -cover_base_z_pos));
	copper_plate1_cover_solid = new G4UnionSolid("EP_COPPER_PLATE1_cover0", copper_plate1_cover_solid, cover_win_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], -cover_win_z_pos));
	copper_plate1_cover_solid = new G4UnionSolid("EP_COPPER_PLATE1_cover0", copper_plate1_cover_solid, cover_rubber_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], -cover_rubber_z_pos));
	G4UnionSolid* copper_plate2_cover_solid = new G4UnionSolid("EP_COPPER_PLATE2_cover0", copper_plate2_hut_solid, cover_base_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], cover_base_z_pos));
	copper_plate2_cover_solid = new G4UnionSolid("EP_COPPER_PLATE2_cover0", copper_plate2_cover_solid, cover_win_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], cover_win_z_pos));
	copper_plate2_cover_solid = new G4UnionSolid("EP_COPPER_PLATE2_cover0", copper_plate2_cover_solid, cover_rubber_solid, 0, G4ThreeVector(pos_x_PMTs[0], pos_y_PMTs[0], cover_rubber_z_pos));
	for(G4int pos_i=1;pos_i<num_PMTs_each;pos_i++)
	{
		copper_plate1_cover_solid = new G4UnionSolid("EP_COPPER_PLATE1_cover0", copper_plate1_cover_solid, cover_base_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], -cover_base_z_pos));
		copper_plate1_cover_solid = new G4UnionSolid("EP_COPPER_PLATE1_cover0", copper_plate1_cover_solid, cover_win_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], -cover_win_z_pos));
		copper_plate1_cover_solid = new G4UnionSolid("EP_COPPER_PLATE1_cover0", copper_plate1_cover_solid, cover_rubber_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], -cover_rubber_z_pos));
		copper_plate2_cover_solid = new G4UnionSolid("EP_COPPER_PLATE2_cover0", copper_plate2_cover_solid, cover_base_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], cover_base_z_pos));
		copper_plate2_cover_solid = new G4UnionSolid("EP_COPPER_PLATE2_cover0", copper_plate2_cover_solid, cover_win_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], cover_win_z_pos));
		copper_plate2_cover_solid = new G4UnionSolid("EP_COPPER_PLATE2_cover0", copper_plate2_cover_solid, cover_rubber_solid, 0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], cover_rubber_z_pos));
	}
	// Placement
	G4double plate_z_pos = 441.5 * mm + copper_plate_thickn_/2.;
	G4LogicalVolume* energy_plane1_logic = new G4LogicalVolume(copper_plate1_cover_solid, steel, "energy plane");
	G4LogicalVolume* energy_plane2_logic = new G4LogicalVolume(copper_plate2_cover_solid, steel, "energy plane");
	new G4PVPlacement(0, G4ThreeVector(0., 0., plate_z_pos), energy_plane1_logic, "energy plane", gas_logic, false, 0, true);
	new G4PVPlacement(0, G4ThreeVector(0., 0., -plate_z_pos), energy_plane2_logic, "energy plane", gas_logic, false, 0, true);

	//////////////////
	// Vacuum holes //
	//////////////////
	G4double hole_front_diam = 80. * mm;
	G4double hole_offset = 14.8 * mm;
	G4double hole_front_depth = 39. * mm;
	G4double hole_rear_diam = 62. * mm;
	G4double hole_rear_depth = 122. * mm + copper_plate_thickn_ - (hole_front_depth - hole_offset);
	G4double hole_win_thick = 8.1 * mm;
	G4double hole_win_diam = 89.8 * mm;
	G4Tubs* hole_front_solid = new G4Tubs("hole_front", 0., hole_front_diam/2., hole_front_depth/2., 0., twopi);
	G4Tubs* hole_rear_solid = new G4Tubs("hole_rear", 0., hole_rear_diam/2., hole_rear_depth/2., 0., twopi);
	G4Tubs* hole_win_solid = new G4Tubs("hole_win", 0., hole_win_diam/2., hole_win_thick/2., 0., twopi);
	G4double hole_rear_z_pos = hole_front_depth/2. + hole_rear_depth/2.;
	G4double hole_win_z_pos = hole_front_depth/2. + hole_win_thick/2.;
	G4UnionSolid* hole1_solid = new G4UnionSolid("HOLE1", hole_front_solid, hole_rear_solid, 0, G4ThreeVector(0., 0., hole_rear_z_pos));
	hole1_solid = new G4UnionSolid("HOLE1", hole1_solid, hole_win_solid, 0, G4ThreeVector(0., 0., -hole_win_z_pos));
	G4UnionSolid* hole2_solid = new G4UnionSolid("HOLE2", hole_front_solid, hole_rear_solid, 0, G4ThreeVector(0., 0., -hole_rear_z_pos));
	hole2_solid = new G4UnionSolid("HOLE2", hole2_solid, hole_win_solid, 0, G4ThreeVector(0., 0., hole_win_z_pos));
	G4LogicalVolume* hole1_logic = new G4LogicalVolume(hole1_solid, vacuum, "EP_HOLE1");
	G4LogicalVolume* hole2_logic = new G4LogicalVolume(hole2_solid, vacuum, "EP_HOLE2");
	// Placement
	G4double plate_hole_z_pos = (hole_front_depth - copper_plate_thickn_)/2. -  hole_offset;
	for(G4int pos_i=0;pos_i<num_PMTs_each;pos_i++)
	{
		new G4PVPlacement(0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], plate_hole_z_pos), hole1_logic, "hole_vacuum", energy_plane1_logic, false, 0, true);
		new G4PVPlacement(0, G4ThreeVector(pos_x_PMTs[pos_i], pos_y_PMTs[pos_i], -plate_hole_z_pos), hole2_logic, "hole_vacuum", energy_plane2_logic, false, 0, true);
	}

	/////////////////////
	// Sapphire window //
	/////////////////////
	G4double sapphire_window_thick = 8. * mm;
	G4double optical_pad_thick = 0.0 * mm;
	G4double tpb_thick = 5.4 * micrometer;
	G4Tubs* sapphire_window_solid = new G4Tubs("SAPPHIRE_WINDOW", 0., hole_win_diam/2., (sapphire_window_thick + tpb_thick)/2., 0., twopi);
	G4LogicalVolume* sapphire_window_logic = new G4LogicalVolume(sapphire_window_solid, sapphire, "SAPPHIRE_WINDOW");
	G4double sapphire_win_posz = -hole_front_depth/2. - hole_win_thick + (sapphire_window_thick + tpb_thick)/2.;
	G4VPhysicalVolume* sapphire_window1_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., sapphire_win_posz), sapphire_window_logic, "SAPPHIRE_WINDOW1", hole1_logic, false, 0, true);
	G4VPhysicalVolume* sapphire_window2_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., -sapphire_win_posz), sapphire_window_logic, "SAPPHIRE_WINDOW2", hole2_logic, false, 0, true);
	// TPB coating on sapphire window
	G4Tubs* tpb_solid = new G4Tubs("SAPPHIRE_WDW_TPB", 0., hole_win_diam/2, tpb_thick/2., 0., twopi);
	G4LogicalVolume* tpb_logic = new G4LogicalVolume(tpb_solid, tpb, "SAPPHIRE_WDW_TPB");
	G4double tpb_posz = - sapphire_window_thick/2. + tpb_thick/2.;
	new G4PVPlacement(0, G4ThreeVector(0., 0., tpb_posz), tpb_logic, "SAPPHIRE_WDW_TPB", sapphire_window_logic, false, 0, false);
	// Optical surface on TPB to model roughness
	G4OpticalSurface* tpb_surf = new G4OpticalSurface("tpb_sapphire_surf", glisur, ground, dielectric_dielectric, .01);
	new G4LogicalSkinSurface("tpb_sapphire_surf", tpb_logic, tpb_surf);
	// /// Optical pad ///
	// G4Tubs* optical_pad_solid = new G4Tubs("OPTICAL_PAD", 0., hole_diam_front_/2., optical_pad_thick/2., 0., twopi);
	// G4LogicalVolume* optical_pad_logic = new G4LogicalVolume(optical_pad_solid, optical_coupler, "OPTICAL_PAD");
	// G4double pad_posz = -vacuum_front_length/2. + sapphire_window_thickn_ + tpb_thick + optical_pad_thick/2.;
	// G4VPhysicalVolume* optical_pad_phys = new G4PVPlacement(0, G4ThreeVector(0., 0., pad_posz), optical_pad_logic, "OPTICAL_PAD", hole_logic, false, 0, false);

	//////////
	// PMTs //
	//////////
	// PMT body
	G4double front_body_diam = 76. * mm;
	G4double front_body_length = 38. * mm;
	G4double rear_body_diam = 53. * mm;
	G4double rear_body_length = 76. * mm;
	G4Tubs* front_body_solid = new G4Tubs("FRONT_BODY", 0., front_body_diam/2., front_body_length/2., 0., twopi);
	G4Tubs* rear_body_solid = new G4Tubs("REAR_BODY", 0., rear_body_diam/2., rear_body_length/2., 0., twopi);
	// Union of the two volumes of the phototube body
	G4double z_transl = -front_body_length/2. - rear_body_length/2.;
	G4ThreeVector transl(0., 0., z_transl);
	G4UnionSolid* pmt_solid = new G4UnionSolid("PMT_R11410",front_body_solid,rear_body_solid,0,transl);
	G4LogicalVolume* pmt_logic = new G4LogicalVolume(pmt_solid, Kovar, "PMT_R11410");
	// PMT gas
	G4double body_thickness = .5 * mm;
	G4double front_body_gas_diam = front_body_diam - 2. * body_thickness;
	G4double front_body_gas_length = front_body_length - body_thickness;
	G4Tubs* front_body_gas_solid = new G4Tubs("FRONT_BODY_GAS", 0., front_body_gas_diam/2., front_body_gas_length/2., 0., twopi);
	G4double rear_body_gas_diam = rear_body_diam - 2. * body_thickness;
	G4double rear_body_gas_length = rear_body_length;
	G4Tubs* rear_body_gas_solid = new G4Tubs("REAR_BODY_GAS", 0., rear_body_gas_diam/2., rear_body_gas_length/2., 0., twopi);
	// Union of the two volumes of the phototube body gas
	G4double z_gas_transl = -front_body_gas_length/2. - rear_body_gas_length/2.;
	G4ThreeVector gas_transl(0., 0., z_gas_transl);
	G4UnionSolid* pmt_gas_solid = new G4UnionSolid("PMT_GAS", front_body_gas_solid, rear_body_gas_solid, 0, gas_transl);
	G4LogicalVolume* pmt_gas_logic = new G4LogicalVolume(pmt_gas_solid, vacuum, "PMT_GAS");
	// PMT window
	G4double window_diam_ = front_body_gas_diam;
	G4double window_thickness_ = 2. * mm;
	G4Tubs* window_solid = new G4Tubs("PMT_WINDOW", 0, window_diam_/2., window_thickness_/2., 0., twopi);
	G4LogicalVolume* window_logic = new G4LogicalVolume(window_solid, silica, "PMT_WINDOW");
	// PMT photocathode
	G4double photocathode_diam_ = 64. * mm;
	G4double photocathode_thickness_ = .1 * mm;
	G4Tubs* photocathode_solid = new G4Tubs("PMT_PHOTOCATHODE", 0, photocathode_diam_/2., photocathode_thickness_/2., 0., twopi);
	G4LogicalVolume* photocathode_logic = new G4LogicalVolume(photocathode_solid, aluminum, "PMT_PHOTOCATHODE");
	// Placement
	G4RotationMatrix* rot = new G4RotationMatrix();
	rot->rotateX(pi);
	G4ThreeVector pmt_posz = G4ThreeVector(0., 0., 0);
	G4double pmt_gas_posz = body_thickness/2.;
	G4double window_posz = front_body_gas_length/2. - window_thickness_/2.;
	G4double photocathode_posz = window_posz - window_thickness_/2. - photocathode_thickness_/2.;
	new G4PVPlacement(0, G4ThreeVector(0., 0., pmt_gas_posz), pmt_gas_logic, "PMT_GAS", pmt_logic, false, 0, true);
	new G4PVPlacement(0, G4ThreeVector(0., 0., window_posz), window_logic, "PMT_WINDOW", pmt_gas_logic, false, 0, true);
	new G4PVPlacement(0, G4ThreeVector(0., 0., photocathode_posz), photocathode_logic, "PMT_PHOTOCATHODE", pmt_gas_logic, false, 0, true);
	new G4PVPlacement(rot, pmt_posz, pmt_logic, "PMT", hole1_logic, false, 0, true);
	new G4PVPlacement(0, pmt_posz, pmt_logic, "PMT", hole2_logic, false, 0, true);
	// Optical properties
	const G4int entries = 57;
	G4double ENERGIES[entries] =
	{2.06640321682*eV, 2.10142700016*eV, 2.13765850016*eV, 2.17516128086*eV,
	2.21400344659*eV, 2.25425805471*eV, 2.29600357425*eV, 2.3393243964*eV,
	2.38431140402*eV, 2.43106260802*eV, 2.47968386018*eV, 2.53028965325*eV,
	2.58300402103*eV, 2.63796155339*eV, 2.69530854368*eV, 2.75520428909*eV,
	2.81782256839*eV, 2.8833533258*eV,  2.95200459546*eV, 3.02400470754*eV,
	3.09960482523*eV, 3.17908187203*eV, 3.2627419213*eV,  3.35092413538*eV,
	3.44400536137*eV, 3.54240551455*eV, 3.64659391204*eV, 3.75709675786*eV,
	3.87450603154*eV, 3.99949009707*eV, 4.13280643364*eV, 4.27531700032*eV,
	4.42800689319*eV, 4.59200714849*eV, 4.76862280805*eV, 4.95936772037*eV,
	5.16600804205*eV, 5.39061708736*eV, 5.63564513678*eV, 5.90400919092*eV,
	6.19920965046*eV, 6.35816374406*eV, 6.52548384259*eV, 6.70184827077*eV,
	6.88801072274*eV, 7.08481102910*eV, 7.29318782407*eV, 7.31470165246*eV,
	7.33634278161*eV, 7.35811234476*eV, 7.38001148864*eV, 7.40204137369*eV,
	7.42420317421*eV, 7.44649807863*eV, 7.46892728971*eV, 7.49149202473*eV,
	7.51419351571*eV };
	G4double EFFICIENCY[entries] =
	{ 0.0530, 0.0625, 0.0720, 0.0850,
	0.1050, 0.1190, 0.1335, 0.1550,
	0.1770, 0.1970, 0.2100, 0.2200,
	0.2300,	0.2430, 0.2580, 0.2770,
	0.2920,	0.3050, 0.3150, 0.3270,
	0.3320, 0.3400, 0.3480, 0.3500,
	0.3530,	0.3600, 0.3680, 0.3650,
	0.3640, 0.3640, 0.3560, 0.3420,
	0.3280, 0.3180, 0.3050, 0.2980,
	0.2920,	0.2900, 0.2920, 0.2945,
	0.3100,	0.3280, 0.3560, 0.3880,
	0.3920,	0.3900, 0.4040, 0.3930,
	0.3700,	0.3500, 0.3300, 0.3150,
	0.2950,	0.2750, 0.2550, 0.2450,
	0.2400 };
	G4double REFLECTIVITY[entries] =
	{ 0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0.,
	0., 0., 0., 0., 0. };
	G4MaterialPropertiesTable* phcath_mpt = new G4MaterialPropertiesTable();
	phcath_mpt->AddProperty("EFFICIENCY", ENERGIES, EFFICIENCY, entries);
	phcath_mpt->AddProperty("REFLECTIVITY", ENERGIES, REFLECTIVITY, entries);
	G4OpticalSurface* pmt_opt_surf = new G4OpticalSurface("PHOTOCATHODE", unified, polished, dielectric_metal);
	pmt_opt_surf->SetMaterialPropertiesTable(phcath_mpt);
	new G4LogicalSkinSurface("PMT_PHOTOCATHODE", photocathode_logic, pmt_opt_surf);

	///////////////////
	// Active volume //
	///////////////////
	G4double gasdrift_length2 = 9.5 * mm;
	G4Tubs* drift1_solid = new G4Tubs("drift1", 0., teflon_inner_diam/2., gasdrift_length1/2., 0., twopi);
	G4Tubs* drift2_solid = new G4Tubs("drift2", 0., EL_diam/2., gasdrift_length2/2., 0., twopi);
	G4double drift_z_pos = gasdrift_length1/2. + gasdrift_length2/2.;
	G4UnionSolid* drift_solid = new G4UnionSolid("drift", drift1_solid, drift2_solid, 0, G4ThreeVector(0., 0., drift_z_pos));
	G4double gasdrift_z_pos = cathode_thick/2. + gasdrift_length1/2.;
	G4LogicalVolume* drift_logic = new G4LogicalVolume(drift_solid, gas, "gasdrift");
	auto gasSD = sensitive_detector(n_event);
	G4SDManager::GetSDMpointer()->AddNewDetector(gasSD);
	drift_logic->SetSensitiveDetector(gasSD);

	new G4PVPlacement(0, G4ThreeVector(0., 0., gasdrift_z_pos), drift_logic, "gasdrift", gas_logic, false, 0, true);
	new G4PVPlacement(rot, G4ThreeVector(0., 0., -gasdrift_z_pos), drift_logic, "gasdrift", gas_logic, false, 1, true);

	//////////////////
	// VISIBILITIES //
	//////////////////
	#if VISUAL
	G4Colour transGrey(0.5, 0.5, 0.5, 0.9);
	G4Colour transGrey2(0.5, 0.5, 0.5, 0.4);
	G4Colour transYellow(1.0, 1.0, 0.0, 0.4);
	G4Colour transRed(1.0, 0.0, 0.0, 0.4);
	G4Colour transCyan(0.0, 1.0, 1.0, 0.6);
	G4Colour transGreen(0.0, 1.0, 0.0, 0.7);
	G4Colour transBlue(0.0, 0.0, 1.0, 0.9);
	G4Colour transPink(1.0, 0.75, 0.8, 0.7);
	// Vessel color
	G4VisAttributes* vessel_col = new G4VisAttributes(transGrey);
	vessel_col->SetForceSolid(true);
	vessel_logic->SetVisAttributes(vessel_col);
	// Gas color
	G4VisAttributes* gas_col = new G4VisAttributes(transGrey);
	gas_col->SetForceSolid(true);
	gas_logic->SetVisAttributes(gas_col);
	// Energy plane color
	G4VisAttributes* energy_plane_col = new G4VisAttributes(transRed);
	energy_plane_col->SetForceSolid(true);
	energy_plane1_logic->SetVisAttributes(energy_plane_col);
	energy_plane2_logic->SetVisAttributes(energy_plane_col);
	// Hole color
	hole1_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
	hole2_logic->SetVisAttributes(G4VisAttributes::GetInvisible());
	// Hole-win color
	G4VisAttributes* sapphire_col = new G4VisAttributes(transPink);
	sapphire_col->SetForceSolid(true);
	sapphire_window_logic->SetVisAttributes(sapphire_col);
	// PMT color
	G4VisAttributes* pmt_col = new G4VisAttributes(transYellow);
	pmt_col->SetForceSolid(true);
	pmt_logic->SetVisAttributes(pmt_col);
	// PMT-gas color
	G4VisAttributes* pmt_gas_col = new G4VisAttributes(transCyan);
	pmt_gas_col->SetForceSolid(true);
	pmt_gas_logic->SetVisAttributes(pmt_gas_col);
	// PMT-window color
	G4VisAttributes* pmt_win_col = new G4VisAttributes(transGreen);
	pmt_win_col->SetForceSolid(true);
	window_logic->SetVisAttributes(pmt_win_col);
	// PMT-photocathode color
	G4VisAttributes* pmt_photo_col = new G4VisAttributes(transBlue);
	pmt_photo_col->SetForceSolid(true);
	photocathode_logic->SetVisAttributes(pmt_photo_col);
	// drift color
	G4VisAttributes* drift_col = new G4VisAttributes(transGreen);
	drift_col->SetForceSolid(true);
	drift_logic->SetVisAttributes(drift_col);
		#if SHIELDING
		G4VisAttributes* ch2_col = new G4VisAttributes(transGreen);
		ch2_col->SetForceSolid(true);
		ch2_logical->SetVisAttributes(ch2_col);
		G4VisAttributes* pb_col = new G4VisAttributes(transRed);
		pb_col->SetForceSolid(true);
		pb_logical->SetVisAttributes(pb_col);
		G4VisAttributes* air_col = new G4VisAttributes(transYellow);
		air_col->SetForceSolid(true);
		air_logical->SetVisAttributes(air_col);
		#endif
	#endif
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
