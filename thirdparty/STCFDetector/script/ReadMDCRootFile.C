#include "DDG4/Geant4Particle.h"
//#include "DDG4/Geant4Data.h"
#include "/Users/danny/software/stcf/software/DD4hep/source/DD4hep/examples/DDG4_MySensDet/src/MyTrackerHit.h"

void ReadMDCRootFile()
{
  gSystem->Load("libDDG4Plugins.dylib");
  gSystem->Load("libDDG4.dylib");
  gSystem->Load("/Users/danny/software/stcf/software/DD4hep-work/lib/libClientTests.dylib");
  //gSystem->Load("/Users/danny/software/stcf/software/DD4hep-work/lib/libDDG4_MySensDet.dylib");
  //gSystem->Load("libDDG4_MySensDet.dylib");
  typedef ROOT::Math::XYZVector Position;
  //typedef SomeExperiment::MyTrackerHit Hit; // class in initial example
  typedef dd4hep::sim::Geant4Tracker::Hit Hit;

  TFile *file = new TFile("MDC_2018-12-27_00-36.root");
  if (file==0) return;
  TTree* tree = (TTree*)file->Get("EVENT");
  

  // access to data
  vector<dd4hep::sim::Geant4Particle*> *mc_v = new vector<dd4hep::sim::Geant4Particle*>;
  tree->SetBranchAddress("MCParticles", &mc_v);
  //EcalBarrelHits : vector<dd4hep::sim::Geant4Calorimeter::Hit*> 
  //vector<dd4hep::sim::Geant4Calorimeter::Hit*> *ecalBarrel_v = new vector<dd4hep::sim::Geant4Calorimeter::Hit*>;
  //tree->SetBranchAddress("EcalBarrelHits", &ecalBarrel_v);
  // SiliconUpperHits : vector<SomeExperiment::MyTrackerHit*> 
  //vector<Hit*> *hit_v = new vector<Hit*>;
  //tree->SetBranchAddress("SiliconUpperHits", &hit_v);
  //dd4hep::sim::Geant4Tracker::Hit, Branch: SiliconDownHits : vector<dd4hep::sim::Geant4Tracker::Hit*>
  //ECALPMTBARREL : vector<dd4hep::sim::Geant4Tracker::Hit*> 
  vector<Hit*> *hit_v = new vector<Hit*>;
  //tree->SetBranchAddress("SiliconUpperHits", &hit_v);
  //MDCCollection : vector<dd4hep::sim::Geant4Tracker::Hit*> 
  tree->SetBranchAddress("MDCCollection", &hit_v);
  
  int nentries = tree->GetEntries();
  cout<<"The number of events is "<<nentries<<endl;

  // test info of each event
  for (int ievt=0; ievt<nentries; ievt++ ) {
    tree->GetEntry(ievt);

	// read truth informaiton
	int np = mc_v->size();
	cout<<"eventID:"<<ievt<<", the number of mc particles is "<< mc_v->size()<<endl;
    // list particle in each event
////cout<<"\t";
////for (int ip=0; ip<np; ip++) {
////  int pdgcode = mc_v->at(ip)->pdgID;
////  cout<<" "<<pdgcode;
////}
////cout<<endl;

	// read echo information
	int nhits = hit_v->size();
	cout<<"n hits "<<nhits<<endl;
	//continue;
	int nshow = nhits>10?10:nhits;
	for (int ihit=0; ihit<nshow; ihit++) {
	  Hit* hit = hit_v->at(ihit);
	  long long int cellID = hit->cellID;
	  int g4ID = hit->g4ID;
	  // typedef ROOT::Math::XYZVector Position
	  Position pos = hit->position;
      double ed = hit->energyDeposit;
	  cout<<dec<<"cellID="<<cellID<<", g4ID="<<g4ID<<"; pos=("<<pos.X()<<", "<<pos.Y()<<", "<<pos.Z()<<"); ed="<<ed<<endl;
	  int sysid = cellID%(1<<10);
	  int sidid = cellID/(1<<10)%(1<<2);
	  int modid = cellID/(1<<12)%(1<<6);
	  int ladid = cellID/(1<<18);
	  cout<<"cellID system:side:module:ladder="<<(bitset<sizeof(long long int)*8>)cellID<<endl;//<<sysid<<":"<<sidid<<":"<<modid<<":"<<ladid<<endl;
	  cout<<"cellID system:side:module:ladder="<<sysid<<":"<<sidid<<":"<<modid<<":"<<ladid<<endl;
	}
  }
}
