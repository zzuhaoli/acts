#include "DDG4/Geant4Config.h"
#include "DDG4/Geant4TestActions.h"
#include "CLHEP/Units/SystemOfUnits.h"
#include "TSystem.h"
//#include "DDG4/Geant4TrackHandler.h"
#include <iostream>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::sim;
using namespace dd4hep::sim::Test;
using namespace dd4hep::sim::Setup;

#if defined(__MAKECINT__)
#pragma link C++ class Geant4RunActionSequence;
#pragma link C++ class Geant4EventActionSequence;
#pragma link C++ class Geant4SteppingActionSequence;
#pragma link C++ class Geant4StackingActionSequence;
#pragma link C++ class Geant4GeneratorActionSequence;
#pragma link C++ class Geant4Action;
#pragma link C++ class Geant4Kernel;
#endif

//SensitiveSeq::handled_type* setupDetector(Kernel& kernel, const std::string& name) {
Geant4SensDetActionSequence* setupDetector(Kernel& kernel, const std::string& name) {
	SensitiveSeq sd = SensitiveSeq(kernel,name);
	Sensitive  sens = Sensitive(kernel,"Geant4TestSensitive/"+name+"Handler",name);
	sd->adopt(sens);
	sens = Sensitive(kernel,"Geant4TestSensitive/"+name+"Monitor",name);
	sd->adopt(sens);
	return sd;
}

void exampleAClick() {
	//Geant4Kernel& kernel = Geant4Kernel::instance(LCDD::getInstance());
	Geant4Kernel& kernel = Geant4Kernel::instance(Detector::getInstance());
	//kernel.loadGeometry("file:../DD4hep.trunk/DDExamples/CLICSiD/compact/compact.xml");
	string installPath = gSystem->Getenv("DD4hepExamplesINSTALL");
	kernel.loadGeometry("file:"+installPath+"/examples/CLICSiD/compact/compact.xml");
	//kernel.loadXML("DDG4_field.xml");
	kernel.loadXML((installPath+"/examples/LHeD/compact/DDG4_field.xml").c_str());

	GenAction gun(kernel,"Geant4ParticleGun/Gun");
	gun["energy"] = 0.5*CLHEP::GeV;
	gun["particle"] = "e-";
	gun["multiplicity"] = 1;
	kernel.generatorAction().adopt(gun);

	//Action run_init(kernel,"Geant4TestRunAction/RunInit");
	RunAction run_init(kernel,"Geant4TestRunAction/RunInit");
	run_init["Property_int"] = 12345;
	//kernel.runAction().callAtBegin (run_init.get(),&Geant4TestRunAction::begin);
	//kernel.runAction().callAtBegin (run_init.get(),&Geant4TestRunAction::begin);
    kernel.runAction().adopt(run_init);
	kernel.eventAction().callAtBegin(run_init.get(),&Geant4TestRunAction::beginEvent);
	kernel.eventAction().callAtEnd (run_init.get(),&Geant4TestRunAction::endEvent);

	//Action evt_1(kernel,"Geant4TestEventAction/UserEvent_1");
	EventAction evt_1(kernel,"Geant4TestEventAction/UserEvent_1");
	evt_1["Property_int"] = 12345; // Set properties
	evt_1["Property_string"] = "Events";
	kernel.eventAction().adopt(evt_1);
	EventAction evt_2(kernel,"Geant4TestEventAction/UserEvent_2");
	kernel.eventAction().adopt(evt_2);
  
    Phase p;
	p = kernel.addPhase<const G4Run*>("BeginRun");
	p->add(evt_2.get(),&Geant4TestEventAction::beginRun);
	kernel.runAction().callAtBegin(p.get(),&Geant4ActionPhase::call<const G4Run*>);

	p = kernel.addPhase<const G4Run*>("EndRun");
	p->add(evt_2.get(),&Geant4TestEventAction::endRun);
	kernel.runAction().callAtEnd(p.get(),&Geant4ActionPhase::call<const G4Run*>);


////kernel.runAction().callAtBegin(evt_2.get(),&Geant4TestEventAction::begin);
////kernel.runAction().callAtEnd (evt_2.get(),&Geant4TestEventAction::end);

	setupDetector(kernel,"SiVertexBarrel");
	setupDetector(kernel,"SiVertexEndcap");
	// .... more subdetectors here .....
	setupDetector(kernel,"LumiCal");
	setupDetector(kernel,"BeamCal");

	kernel.configure();
	kernel.initialize();
	kernel.run();
	std::cout << "Successfully executed application .... " << std::endl;
	kernel.terminate();
}
