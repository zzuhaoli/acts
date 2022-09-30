#
#
import os, sys, time, DDG4, dd4hep
from DDG4 import OutputLevel as Output
from SystemOfUnits import *
#
#
"""

   dd4hep simulation example setup using the python configuration

   @author  M.Frank
   @version 1.0

"""
def run():
  kernel = DDG4.Kernel()
  #install_dir = os.environ['DD4hepExamplesINSTALL']
  #install_dir = "/Users/danny/software/stcf/software/DD4hep/source/DD4hep"
  #kernel.loadGeometry("file:"+install_dir+"/examples/ClientTests/compact/SiliconBlock.xml")
  #kernel.loadGeometry("file:./SiliconBlock.xml")
  kernel.loadGeometry("file:../compact/detectorECal.xml")
  DDG4.importConstants(kernel.detectorDescription(),debug=False)
  # =======================================================================================
  # ===> This is actually the ONLY difference to ClientTests/scripts/SiliconBlock.py
  # =======================================================================================
  #geant4 = DDG4.Geant4(kernel,tracker='MyTrackerSDAction')
  # add by Dong Liu, save mc information?
  # use default tracker action
  geant4 = DDG4.Geant4(kernel, calo='Geant4CalorimeterAction')

  geant4.printDetectors()
  kernel.NumEvents = 50
  kernel.UI = ''

  # Configure field
  field = geant4.setupTrackingField(prt=True)
  # Configure Event actions
  prt = DDG4.EventAction(kernel,'Geant4ParticlePrint/ParticlePrint')
  prt.OutputLevel = Output.WARNING
  prt.OutputType  = 3 # Print both: table and tree
  kernel.eventAction().adopt(prt)

  # Configure I/O
  evt_root = geant4.setupROOTOutput('RootOutput','ECal_'+time.strftime('%Y-%m-%d_%H-%M'),mc_truth=False)
  # Setup particle gun
  #gun = geant4.setupGun("Gun",particle='mu-',energy=5*GeV,multiplicity=1,Standalone=True,position=(0,0,0))
  gun = geant4.setupGun("Gun",particle='e-',energy=5*GeV,multiplicity=1,Standalone=True,position=(0,0,0))
  
  # add by Dong Liu, save mc information?
  #logging.info("#  ....and handle the simulation particles.")
  part = DDG4.GeneratorAction(kernel,"Geant4ParticleHandler/ParticleHandler")
  kernel.generatorAction().adopt(part)
  
  geant4.setupTracker('STCFBEMC')
  #geant4.setupCalorimeter('STCFBEMC')
  #geant4.setupCalorimeter('STCFEEMC')
  #geant4.setupTracker('SiliconBlockUpper')
  #geant4.setupTracker('SiliconBlockDown')
  
  ## Now build the physics list:
  #phys = kernel.physicsList()
  #phys.extends = 'QGSP_BERT'
  #phys.enableUI()
  #phys.dump()

  #phys = geant4.setupPhysics('QGSP_BERT')
  #ph = DDG4.PhysicsList(kernel,'Geant4PhysicsList/Myphysics')
  #ph.addParticleConstructor('G4Geantino')
  #ph.addParticleConstructor('G4BosonConstructor')
  #ph.enableUI()
  #phys.adopt(ph)
  #phys.dump()

  phys = geant4.setupPhysics('QGSP_BERT') 
  #ph = geant4.addPhysics('Geant4PhysicsList/Myphysics')  
  part = geant4.addPhysics('Geant4ExtraParticles/ExtraParticles') 
  part.pdgfile = '/Users/danny/software/stcf/software/DD4hep/source/DD4hep/DDG4/examples/particle.tbl'
  rg = geant4.addPhysics('Geant4DefaultRangeCut/GlobalRangeCut')
  rg.RangeCut = 0.7*mm
  phys.dump()

  # run
  kernel.configure()
  kernel.initialize()
  kernel.run()
  kernel.terminate()

if __name__ == "__main__":
  run()
