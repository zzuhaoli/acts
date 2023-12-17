// This file is part of the Acts project.
//
// Copyright (C) 2017-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/IReader.hpp"
#include "ActsExamples/Framework/IService.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Propagator/MaterialInteractor.hpp>
#include <Acts/Utilities/Logger.hpp>
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <mutex>
#include <vector>
#include <string>

class TTreeReader;

namespace ActsExamples
{
class RootHoughCanTracksReader : public IReader { 
  public: 
    struct Config
    {
      std::string outputcantrackparameters;
      std::string treeName = "events";
      std::string filePath;

      /// Constant term of the Loc0 resolution.
      double sigmaLoc0 = 100 * Acts::UnitConstants::um;
      /// Constant term of the Loc1 resolution.
      double sigmaLoc1 = 400 * Acts::UnitConstants::um;
      /// Time resolution.
      double sigmaT0 = 1400 * Acts::UnitConstants::s;
      /// Phi angular resolution.
      double sigmaPhi = 0.2 * Acts::UnitConstants::degree;
      /// Theta angular resolution.
      double sigmaTheta = 0.2 * Acts::UnitConstants::degree;
      /// Relative momentum resolution.
      double sigmaQOverP = 0.1/1 * Acts::UnitConstants::GeV;
      /// Inflate the initial covariance matrix
      std::array<double, 6> initialVarInflation = {1., 1., 1., 1., 1., 1.};
      std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
    };
     RootHoughCanTracksReader(const Config& config,Acts::Logging::Level level);
    
     ~RootHoughCanTracksReader();
     
     std::string name() const final override{return "RootHoughCanTracksReader";}
    
     ProcessCode read(const ActsExamples::AlgorithmContext& context) final override;
     
     /// Return the available events range.
     std::pair<size_t, size_t> availableEvents() const final override;     

     const Config& config() const{return m_cfg;}
     
   private:
     std::unique_ptr<const Acts::Logger> m_logger;
     
     const Acts::Logger& logger() const {return *m_logger;} 

     Config m_cfg;
     
     std::mutex m_read_mutex;

     size_t m_evtCounter =0;
     size_t m_events; 
     TTreeReader* m_treeReader = nullptr;
     TTreeReaderArray<int>* PDG = nullptr;
     TTreeReaderArray<int>* trackID = nullptr;
     TTreeReaderArray<int>* Charge = nullptr;
     TTreeReaderArray<Acts::ActsScalar>* initialPositionX = nullptr;
     TTreeReaderArray<Acts::ActsScalar>* initialPositionY = nullptr;
     TTreeReaderArray<Acts::ActsScalar>* initialPositionZ = nullptr;
     TTreeReaderArray<Acts::ActsScalar>* initialMomentumX = nullptr;
     TTreeReaderArray<Acts::ActsScalar>* initialMomentumY = nullptr;
     TTreeReaderArray<Acts::ActsScalar>* initialMomentumZ = nullptr;
     TTreeReaderArray<double>* chiSquare = nullptr;

};//class RootHoughCanTracksReader

}//namespace ActsExamples



