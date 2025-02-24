// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Root/RootHoughCanTracksReader.hpp"

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
//#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/WriterT.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

ActsExamples::RootHoughCanTracksReader::RootHoughCanTracksReader(
    const ActsExamples::RootHoughCanTracksReader::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_logger{Acts::getDefaultLogger(name(), level)},
      m_cfg(config),
      m_events(0),
      m_treeReader(nullptr) {
  if (m_cfg.outputcantrackparameters.empty()) {
    throw std::invalid_argument("Missing candidate tracks parameters");
  }
  if (m_cfg.filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }

  TFile* file = TFile::Open(m_cfg.filePath.c_str(), "READ");

  m_treeReader = new TTreeReader(m_cfg.treeName.c_str(), file);

  PDG = new TTreeReaderArray<int>(*m_treeReader, "CanTrackCol.pdgID");
  trackID = new TTreeReaderArray<int>(*m_treeReader, "CanTrackCol.trackID");
  Charge = new TTreeReaderArray<int>(*m_treeReader, "CanTrackCol.charge");
  initialPositionX = new TTreeReaderArray<Acts::ActsScalar>(
      *m_treeReader, "CanTrackCol.initialPos.x");
  initialPositionY = new TTreeReaderArray<Acts::ActsScalar>(
      *m_treeReader, "CanTrackCol.initialPos.y");
  initialPositionZ = new TTreeReaderArray<Acts::ActsScalar>(
      *m_treeReader, "CanTrackCol.initialPos.z");
  initialMomentumX = new TTreeReaderArray<Acts::ActsScalar>(
      *m_treeReader, "CanTrackCol.initialMom.x");
  initialMomentumY = new TTreeReaderArray<Acts::ActsScalar>(
      *m_treeReader, "CanTrackCol.initialMom.y");
  initialMomentumZ = new TTreeReaderArray<Acts::ActsScalar>(
      *m_treeReader, "CanTrackCol.initialMom.z");
  chiSquare =
      new TTreeReaderArray<double>(*m_treeReader, "CanTrackCol.chiSquare");

  ACTS_DEBUG("Adding File " << m_cfg.filePath << " to tree '" << m_cfg.treeName
                            << "'.");

  m_events = m_treeReader->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

std::pair<size_t, size_t>
ActsExamples::RootHoughCanTracksReader::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::RootHoughCanTracksReader::~RootHoughCanTracksReader() {
  delete PDG;
  delete trackID;
  delete Charge;
  delete initialPositionX;
  delete initialPositionY;
  delete initialPositionZ;
  delete initialMomentumX;
  delete initialMomentumY;
  delete initialMomentumZ;
  delete chiSquare;
}

ActsExamples::ProcessCode ActsExamples::RootHoughCanTracksReader::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read Hough candidate tracks parameters.");
  auto rng = m_cfg.randomNumbers->spawnGenerator(context);
  std::normal_distribution<double> stdNormal(0.0, 1.0);
  TrackParametersContainer trackparameters;
  if (m_treeReader != nullptr && context.eventNumber < m_events) {
    std::lock_guard<std::mutex> lock(m_read_mutex);
    if (m_treeReader->Next()) {
      std::cout << "Reading event" << m_evtCounter++ << std::endl;
      if (m_evtCounter == 94386) {
        context.eventStore.add(m_cfg.outputcantrackparameters,
                               std::move(trackparameters));
        return ActsExamples::ProcessCode::SUCCESS;
      }
      for (size_t i = 0; i < trackID->GetSize(); i++) {
        Acts::Vector3 inipos((*initialPositionX)[i] * Acts::UnitConstants::mm,
                             (*initialPositionY)[i] * Acts::UnitConstants::mm,
                             (*initialPositionZ)[i] * Acts::UnitConstants::mm);
        auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(inipos);

        Acts::Vector3 inimom(
            (*initialMomentumX)[i] / 1000 * Acts::UnitConstants::GeV,
            (*initialMomentumY)[i] / 1000 * Acts::UnitConstants::GeV,
            (*initialMomentumZ)[i] / 1000 * Acts::UnitConstants::GeV);

        // caulate phi and theta
        auto p = std::sqrt(inimom.x() * inimom.x() + inimom.y() * inimom.y() +
                           inimom.z() * inimom.z());
        auto theta = std::acos(inimom.z() / p) * 180 / M_PI *
                     Acts::UnitConstants::degree;
        auto phi = std::atan2(inimom.y(), inimom.x()) * 180 / M_PI *
                   Acts::UnitConstants::degree;

        auto q = (*Charge)[i];
        std::cout << "P is:" << p << std::endl;
        std::cout << "Q is:" << q << std::endl;
        //////////////////////////////////
        const double sigmaLoc0 = m_cfg.sigmaLoc0;
        const double sigmaLoc1 = m_cfg.sigmaLoc1;
        const double sigmaP = m_cfg.sigmaQOverP * (p * p);
        // var(q/p) = (d(1/p)/dp)² * var(p) = (-1/p²)² * var(p)
        const double sigmaQOverP = m_cfg.sigmaQOverP;
        // shortcuts for other resolutions
        const double sigmaT0 = m_cfg.sigmaT0;
        const double sigmaPhi = m_cfg.sigmaPhi;
        const double sigmaTheta = m_cfg.sigmaTheta;
        /////////////////////////////////////
        Acts::BoundVector params = Acts::BoundVector::Zero();
        // no needing smear Loc0 and Loc1
        params[Acts::eBoundLoc0] = 0;
        params[Acts::eBoundLoc1] = 0;
        // set time=0
        params[Acts::eBoundTime] = 0;
        // no needing smear Phi,Theta,P
        const auto [newPhi, newTheta] = Acts::detail::normalizePhiTheta(
            phi,theta);
        params[Acts::eBoundPhi] = newPhi;
        params[Acts::eBoundTheta] = newTheta;
        // compute smeared absolute momentum vector
        const double newP = std::max(0.0, p);
        params[Acts::eBoundQOverP] = (q != 0) ? (q / newP) : (1 / newP);
        std::cout << "q_p:" << params[Acts::eBoundQOverP] << std::endl;
        // build the track covariance matrix using the smearing sigmas//
        // Acts::BoundSquareMatrix cov = Acts::BoundSquareMatrix::Zero();
        Acts::BoundSymMatrix cov = Acts::BoundSymMatrix::Zero();
        cov(Acts::eBoundLoc0, Acts::eBoundLoc0) =
            m_cfg.initialVarInflation[Acts::eBoundLoc0] * sigmaLoc0 * sigmaLoc0;
        cov(Acts::eBoundLoc1, Acts::eBoundLoc1) =
            m_cfg.initialVarInflation[Acts::eBoundLoc1] * sigmaLoc1 * sigmaLoc1;
        cov(Acts::eBoundTime, Acts::eBoundTime) =
            m_cfg.initialVarInflation[Acts::eBoundTime] * sigmaT0 * sigmaT0;
        cov(Acts::eBoundPhi, Acts::eBoundPhi) =
            m_cfg.initialVarInflation[Acts::eBoundPhi] * sigmaPhi * sigmaPhi;
        cov(Acts::eBoundTheta, Acts::eBoundTheta) =
            m_cfg.initialVarInflation[Acts::eBoundTheta] * sigmaTheta *
            sigmaTheta;
        cov(Acts::eBoundQOverP, Acts::eBoundQOverP) =
            m_cfg.initialVarInflation[Acts::eBoundQOverP] * sigmaQOverP *
            sigmaQOverP;

        trackparameters.emplace_back(perigee, params, q, cov);

      }  //读取每个事例的track
    }    //循环if(m_treeReader->Next())
    // Write the tracks parameters collections to the EventStore
    context.eventStore.add(m_cfg.outputcantrackparameters,
                           std::move(trackparameters));
  } else {
    ACTS_WARNING("Could not read in event.");
  }  // if判断if(m_treeReader !=nullptr && context.eventNumber < m_events)

  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
