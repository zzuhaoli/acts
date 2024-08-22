// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#include "ActsExamples/Io/Root/RootSTCFMeasurementReader_BKG.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <iostream>

#include <TChain.h>
#include <TFile.h>
#include <TMath.h>

ActsExamples::RootSTCFMeasurementReader_BKG::RootSTCFMeasurementReader_BKG(
    const ActsExamples::RootSTCFMeasurementReader_BKG::Config& config,
    Acts::Logging::Level level)
    : ActsExamples::IReader(),
      m_logger{Acts::getDefaultLogger(name(), level)},
      m_cfg(config),
      m_events(0),
      m_particle_treeReader(nullptr),
      m_simhit_treeReader(nullptr) {
  if (m_cfg.v_filePath.empty()) {
    throw std::invalid_argument("Missing input filename");
  }
  if (m_cfg.outputParticles.empty()) {
    throw std::invalid_argument(
        "Missing simulated particles output collection");
  }
  if (m_cfg.outputSimHits.empty()) {
    throw std::invalid_argument("Missing simulated hits output collection");
  }
  if (m_cfg.outputMeasurements.empty()) {
    throw std::invalid_argument("Missing measurements output collection");
  }
  if (m_cfg.outputSourceLinks.empty()) {
    throw std::invalid_argument("Missing source links output collection");
  }
  if (m_cfg.outputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-particles map output collection");
  }
  if (m_cfg.outputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument(
        "Missing hit-to-simulated-hits map output collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }

  TFile* MCParticle_file = TFile::Open(m_cfg.v_filePath[0].c_str(), "READ");
  TFile* SimHitBKG_file = TFile::Open(m_cfg.v_filePath[1].c_str(), "READ");

  m_particle_treeReader =
      new TTreeReader(m_cfg.v_treeName[0].c_str(), MCParticle_file);
  m_simhit_treeReader =
      new TTreeReader(m_cfg.v_treeName[1].c_str(), SimHitBKG_file);

  // Set the branches

  particlePDG =
      new TTreeReaderArray<int>(*m_particle_treeReader, "MCParticleCol.PDG");
  particleCharge = new TTreeReaderArray<float>(*m_particle_treeReader,
                                               "MCParticleCol.charge");
  particleMass = new TTreeReaderArray<Acts::ActsScalar>(*m_particle_treeReader,
                                                        "MCParticleCol.mass");
  particleVertexX = new TTreeReaderArray<Acts::ActsScalar>(
      *m_particle_treeReader, "MCParticleCol.vertex.x");
  particleVertexY = new TTreeReaderArray<Acts::ActsScalar>(
      *m_particle_treeReader, "MCParticleCol.vertex.y");
  particleVertexZ = new TTreeReaderArray<Acts::ActsScalar>(
      *m_particle_treeReader, "MCParticleCol.vertex.z");
  particleMomentumX = new TTreeReaderArray<float>(*m_particle_treeReader,
                                                  "MCParticleCol.momentum.x");
  particleMomentumY = new TTreeReaderArray<float>(*m_particle_treeReader,
                                                  "MCParticleCol.momentum.y");
  particleMomentumZ = new TTreeReaderArray<float>(*m_particle_treeReader,
                                                  "MCParticleCol.momentum.z");
  particleTrackID = new TTreeReaderArray<int>(*m_particle_treeReader,
                                              "MCParticleCol.trackID");
  // particleTime  = new
  // TTreeReaderArray<int>(*m_treeReader,"MCParticleCol.time");

  // The ITKhitcol_mcParticleIndex means the "i"th particle in the particle
  // collection above
  ITKtype = new TTreeReaderArray<int>(*m_simhit_treeReader, "ITKhitcol_type");
  ITKlayerID =
      new TTreeReaderArray<int>(*m_simhit_treeReader, "ITKhitcol_layerId");
  ITKparentID =
      new TTreeReaderArray<int>(*m_simhit_treeReader, "ITKhitcol_parentId");
  ITKparticleId = new TTreeReaderArray<int>(*m_simhit_treeReader,
                                            "ITKhitcol_mcParticleIndex");
  ITKpositionX =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "ITKhitcol_positionX");
  ITKpositionY =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "ITKhitcol_positionY");
  ITKpositionZ =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "ITKhitcol_positionZ");
  ITKmomentumX =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "ITKhitcol_momentumX");
  ITKmomentumY =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "ITKhitcol_momentumY");
  ITKmomentumZ =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "ITKhitcol_momentumZ");
  ITKtime = new TTreeReaderArray<Acts::ActsScalar>(*m_simhit_treeReader,
                                                   "ITKhitcol_time");
  ITKmass =
      new TTreeReaderArray<double>(*m_simhit_treeReader, "ITKhitcol_mass");

  // The MDChitcol_mcParticleIndex means the "i"th particle in the particle
  // collection above
  MDCtype = new TTreeReaderArray<int>(*m_simhit_treeReader, "MDChitcol_type");
  MDCcellID =
      new TTreeReaderArray<int>(*m_simhit_treeReader, "MDChitcol_cellId");
  MDClayerID =
      new TTreeReaderArray<int>(*m_simhit_treeReader, "MDChitcol_layerId");
  MDCparentID =
      new TTreeReaderArray<int>(*m_simhit_treeReader, "MDChitcol_parentId");
  MDCparticleId = new TTreeReaderArray<int>(*m_simhit_treeReader,
                                            "MDChitcol_mcParticleIndex");
  MDCpositionX =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_positionX");
  MDCpositionY =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_positionY");
  MDCpositionZ =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_positionZ");
  MDCmomentumX =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_momentumX");
  MDCmomentumY =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_momentumY");
  MDCmomentumZ =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_momentumZ");
  MDCdriftDistance = new TTreeReaderArray<Acts::ActsScalar>(
      *m_simhit_treeReader, "MDChitcol_driftdistance");
  MDCwirePoint1X =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_wire1X");
  MDCwirePoint1Y =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_wire1Y");
  MDCwirePoint1Z =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_wire1Z");
  MDCwirePoint2X =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_wire2X");
  MDCwirePoint2Y =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_wire2Y");
  MDCwirePoint2Z =
      new TTreeReaderArray<float>(*m_simhit_treeReader, "MDChitcol_wire2Z");
  MDCtime = new TTreeReaderArray<Acts::ActsScalar>(*m_simhit_treeReader,
                                                   "MDChitcol_time");
  MDCmass =
      new TTreeReaderArray<double>(*m_simhit_treeReader, "MDChitcol_mass");

  ACTS_DEBUG("Adding File " << m_cfg.v_filePath[0] << " to tree '"
                            << m_cfg.v_treeName[0] << "'.");
  ACTS_DEBUG("Adding File " << m_cfg.v_filePath[1] << " to tree '"
                            << m_cfg.v_treeName[1] << "'.");
  m_events = m_simhit_treeReader->GetEntries();
  ACTS_DEBUG("The full chain has " << m_events << " entries.");
}

std::pair<size_t, size_t>
ActsExamples::RootSTCFMeasurementReader_BKG::availableEvents() const {
  return {0u, m_events};
}

ActsExamples::RootSTCFMeasurementReader_BKG::~RootSTCFMeasurementReader_BKG() {
  delete particlePDG;
  delete particleCharge;
  delete particleMass;
  delete particleTrackID;
  delete particleVertexX;
  delete particleVertexY;
  delete particleVertexZ;
  delete particleMomentumX;
  delete particleMomentumY;
  delete particleMomentumZ;

  delete ITKtype;
  delete ITKlayerID;
  delete ITKparentID;
  delete ITKparticleId;
  delete ITKpositionX;
  delete ITKpositionY;
  delete ITKpositionZ;
  delete ITKmomentumX;
  delete ITKmomentumY;
  delete ITKmomentumZ;
  delete ITKtime;
  delete ITKmass;

  delete MDCtype;
  delete MDCcellID;
  delete MDClayerID;
  delete MDCparentID;
  delete MDCparticleId;
  delete MDCpositionX;
  delete MDCpositionY;
  delete MDCpositionZ;
  delete MDCmomentumX;
  delete MDCmomentumY;
  delete MDCmomentumZ;
  delete MDCdriftDistance;
  delete MDCwirePoint1X;
  delete MDCwirePoint1Y;
  delete MDCwirePoint1Z;
  delete MDCwirePoint2X;
  delete MDCwirePoint2Y;
  delete MDCwirePoint2Z;
  delete MDCtime;
  delete MDCmass;
}

ActsExamples::ProcessCode ActsExamples::RootSTCFMeasurementReader_BKG::read(
    const ActsExamples::AlgorithmContext& context) {
  ACTS_DEBUG("Trying to read STCF tracker measurements.");

  // setup random number generator and standard gaussian
  auto rng = m_cfg.randomNumbers->spawnGenerator(context);
  std::normal_distribution<double> stdNormal(0.0, 1.0);

  // Prepare output containers
  // need list here for stable addresses
  SimParticleContainer particles;
  SimParticleContainer::sequence_type unordered_particles;

  SimHitContainer simHits;
  SimHitContainer::sequence_type unordered_hits;

  std::list<IndexSourceLink> sourceLinkStorage;
  IndexSourceLinkContainer sourceLinks;
  MeasurementContainer measurements;
  IndexMultimap<ActsFatras::Barcode> measurementParticlesMap;
  IndexMultimap<Index> measurementSimHitsMap;

  // read in the fitted track parameters and particles
  if (m_particle_treeReader != nullptr && m_simhit_treeReader != nullptr &&
      context.eventNumber < m_events) {
    // lock the mutex
    std::lock_guard<std::mutex> lock(m_read_mutex);
    std::map<int, std::vector<int>> particleITKLayer;
    std::map<int, int> particleITKAllHitIdx;
    std::map<int, int> particleITKSigHitIdx;
    std::map<int, int> particleHitIdx;
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // I. read sim hits
    /////////////////////////////////////////////////////////////////////////////////////////////////
    if (m_simhit_treeReader->Next()) {
      std::cout << "Reading event " << m_evtCounter++ << std::endl;
      // int nParticles = 0;
      int nITKHits = 0;
      int nMDCHits = 0;
      // The index of the sim hit among all the sim hits from this particle
      particleITKAllHitIdx.clear();
      particleITKSigHitIdx.clear();
      particleITKLayer.clear();
      particleHitIdx.clear();

      /////////////////////////////////////////////////////////////////////////////////////////////////
      // 1) Reading ITK hits
      /////////////////////////////////////////////////////////////////////////////////////////////////
      //for (size_t i = 0; i < ITKpositionX->GetSize(); ++i) {
      for (size_t i = 0; i < ITKtype->GetSize(); ++i) {
        int parentID = (*ITKparentID)[i];
        if (parentID != 0){
          continue;
	}
        int type = (*ITKtype)[i];
	// std::cout<<"ITK type is :"<<type<<std::endl;
        int particleId = (*ITKparticleId)[i];
        //bool isNoiseHit = ((*ITKmomentumX)[i] ==0 and (*ITKmomentumY)[i]==0 and (*ITKmomentumZ)[i] == 0); 
        bool isNoiseHit = (type == -1); 
        if(m_cfg.ignoreNoiseHits and isNoiseHit) {
	  continue;
	}
        nITKHits++;
        Acts::Vector3 pos((*ITKpositionX)[i], (*ITKpositionY)[i],
                          (*ITKpositionZ)[i]);
        Acts::Vector3 mom((*ITKmomentumX)[i], (*ITKmomentumY)[i],
                          (*ITKmomentumZ)[i]);

        // std::cout<<"r = " << std::hypot((*ITKpositionX)[i],
        // (*ITKpositionY)[i]) <<", phi = " << std::atan2((*ITKpositionY)[i],
        // (*ITKpositionX)[i]) << std::endl;
        int layerID = (*ITKlayerID)[i];
        Acts::GeometryIdentifier geoId = Acts::GeometryIdentifier()
                                             .setVolume(m_volumeIDs[0])
                                             .setLayer(2 * (layerID + 1))
                                             .setSensitive(1);
        const Acts::Surface* surfacePtr =
            m_cfg.trackingGeometry->findSurface(geoId);

        Acts::Vector3 posUpdated = pos;
        if (not isNoiseHit) {
          auto intersection = surfacePtr->intersect(context.geoContext, pos,
                                                    mom.normalized(), true);
          posUpdated = intersection.intersection.position;
        }

        auto cylinderSurface =
            dynamic_cast<const Acts::CylinderSurface*>(surfacePtr);
        if (cylinderSurface == nullptr) {
          std::cout << "Cast ITK surface to cylinder surface failed "
                    << std::endl;
        } else {
          auto bounds = cylinderSurface->bounds();
          auto values = bounds.values();
          // std::cout<<"ITK surface r = "<< values[0] << std::endl;
        }

        ActsFatras::Hit::Vector4 pos4{
            posUpdated.x() * Acts::UnitConstants::mm,
            posUpdated.y() * Acts::UnitConstants::mm,
            posUpdated.z() * Acts::UnitConstants::mm,
            (*ITKtime)[i] * Acts::UnitConstants::ns,
        };
        auto energy =
            std::sqrt(mom.x() * mom.x() + mom.y() * mom.y() +
                      mom.z() * mom.z() + (*ITKmass)[i] * (*ITKmass)[i]);
        ActsFatras::Hit::Vector4 mom4{
            mom.x() / 1000 * Acts::UnitConstants::GeV,
            mom.y() / 1000 * Acts::UnitConstants::GeV,
            mom.z() / 1000 * Acts::UnitConstants::GeV,
            energy / 1000. * Acts::UnitConstants::GeV,
        };
        ActsFatras::Hit::Vector4 delta4{
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
        };

        // The last parameter could be incorrect!!!
        ActsFatras::Hit hit(geoId, particleId, pos4, mom4, mom4 + delta4,
                            particleHitIdx[particleId]);
        unordered_hits.push_back(std::move(hit));

        particleHitIdx[particleId]++;
        particleITKAllHitIdx[particleId]++;
        if (not isNoiseHit) {
          particleITKSigHitIdx[particleId]++;
        }
        particleITKLayer[particleId].push_back(layerID);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////
      // 2) Reading MDC hits
      /////////////////////////////////////////////////////////////////////////////////////////////////
      for (size_t i = 0; i < MDCcellID->GetSize(); ++i) {
        int parentID = (*MDCparentID)[i];
        if (parentID != 0){
          continue;
	}
        int type = (*MDCtype)[i];
        // std::cout<<"MDC type is :"<<type<<std::endl;
        int particleId = (*MDCparticleId)[i];
	//bool isNoiseHit = ((*MDCmomentumX)[i]==0 and (*MDCmomentumY)[i]==0 and (*MDCmomentumZ)[i]==0);
	bool isNoiseHit = (type == -1);
        if(m_cfg.ignoreNoiseHits and isNoiseHit) {
	  continue;
	}
        nMDCHits++;

        Acts::Vector3 pos((*MDCpositionX)[i], (*MDCpositionY)[i],
                          (*MDCpositionZ)[i]);
        Acts::Vector3 mom((*MDCmomentumX)[i], (*MDCmomentumY)[i],
                          (*MDCmomentumZ)[i]);

        // double dang = 2*M_PI/m_MDCnCells[(*MDClayerID)[i]];
        int layerID = (*MDClayerID)[i];
        int cellID = (*MDCcellID)[i] - m_MDCsCells[(*MDClayerID)[i]];

        Acts::GeometryIdentifier geoId = Acts::GeometryIdentifier()
                                             .setVolume(m_volumeIDs[1])
                                             .setLayer((layerID + 1) * 2)
                                             .setSensitive(cellID + 1);
        const Acts::Surface* surfacePtr =
            m_cfg.trackingGeometry->findSurface(geoId);

        Acts::Vector3 posUpdated = pos;
        if (not isNoiseHit) {
          // Only update the position if the MDC hit is not a noise
          auto intersection = surfacePtr->intersect(context.geoContext, pos,
                                                    mom.normalized(), true);
          posUpdated = intersection.intersection.position;
          auto lpResult = surfacePtr->globalToLocal(
              context.geoContext, posUpdated, mom.normalized());
          if (not lpResult.ok()) {
            ACTS_FATAL("Global to local transformation did not succeed.");
            return ProcessCode::ABORT;
          }
          // auto lPosition = lpResult.value();
          // std::cout<<"MDC " << geoId <<" has drift distance = " <<
          // (*MDCdriftDistance)[i] << ", local Pos x = " << lPosition[0] <<",
          // localPos.y() " << lPosition[1] << std::endl;
        }


        ActsFatras::Hit::Vector4 pos4{
            posUpdated.x() * Acts::UnitConstants::mm,
            posUpdated.y() * Acts::UnitConstants::mm,
            posUpdated.z() * Acts::UnitConstants::mm,
            (*MDCtime)[i] * Acts::UnitConstants::ns,
        };
        auto energy =
            std::sqrt(mom.x() * mom.x() + mom.y() * mom.y() +
                      mom.z() * mom.z() + (*MDCmass)[i] * (*MDCmass)[i]);
        ActsFatras::Hit::Vector4 mom4{
            mom.x() / 1000 * Acts::UnitConstants::GeV,
            mom.y() / 1000 * Acts::UnitConstants::GeV,
            mom.z() / 1000 * Acts::UnitConstants::GeV,
            energy / 1000. * Acts::UnitConstants::GeV,
        };
        ActsFatras::Hit::Vector4 delta4{
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
            0 * Acts::UnitConstants::GeV,
            // 0 * Acts::UnitConstants::GeV,
            (*MDCdriftDistance)[i] * Acts::UnitConstants::mm,
        };

        // The last parameter is not simHitIdx?
        ActsFatras::Hit hit(geoId, particleId, pos4, mom4, mom4 + delta4,
                            particleHitIdx[particleId]);
        unordered_hits.push_back(std::move(hit));
        particleHitIdx[particleId]++;
      }

      simHits.insert(unordered_hits.begin(), unordered_hits.end());

      sourceLinks.reserve(simHits.size());
      measurements.reserve(simHits.size());
      measurementParticlesMap.reserve(simHits.size());
      measurementSimHitsMap.reserve(simHits.size());

      /////////////////////////////////////////////////////////////////////////////////////////////////
      // 3) Smear truth hits to emulate measurements
      /////////////////////////////////////////////////////////////////////////////////////////////////
      ACTS_DEBUG("Starting loop over modules ...");
      for (auto simHitsGroup : groupByModule(simHits)) {
        // Manual pair unpacking instead of using
        //   auto [moduleGeoId, moduleSimHits] : ...
        // otherwise clang on macos complains that it is unable to capture the
        // local binding in the lambda used for visiting the smearer below.
        Acts::GeometryIdentifier moduleGeoId = simHitsGroup.first;
        const auto& moduleSimHits = simHitsGroup.second;

        const Acts::Surface* surfacePtr =
            m_cfg.trackingGeometry->findSurface(moduleGeoId);

        if (surfacePtr == nullptr) {
          // this is either an invalid geometry id or a misconfigured smearer
          // setup; both cases can not be handled and should be fatal.
          ACTS_ERROR("Could not find surface " << moduleGeoId
                                               << " for configured smearer");
          return ProcessCode::ABORT;
        }

        for (auto h = moduleSimHits.begin(); h != moduleSimHits.end(); ++h) {
          const auto& simHit = *h;
          const auto simHitIdx = simHits.index_of(h);
          auto pos = simHit.position();
          auto dir = simHit.unitDirection();
          auto particleId = simHit.particleId();

          // This is used to decide whether this hit is a noise (noise hit have
          // zero momentum)
          auto momBefore = simHit.momentum4Before();
          bool isNoise = false;
          if (momBefore[0] == 0 and momBefore[1] == 0 and momBefore[2] == 0) {
            isNoise = true;
          }

          Index measurementIdx = measurements.size();
          sourceLinkStorage.emplace_back(moduleGeoId, measurementIdx);
          IndexSourceLink& sourceLink = sourceLinkStorage.back();
          sourceLinks.insert(sourceLinks.end(), sourceLink);

          if (surfacePtr->type() == Acts::Surface::SurfaceType::Cylinder) {
            std::array<Acts::BoundIndices, 2> indices = {Acts::eBoundLoc0,
                                                         Acts::eBoundLoc1};

            Acts::ActsSymMatrix<2> cov = Acts::ActsSymMatrix<2>::Identity();
            cov(0, 0) = 0.1 * 0.1;
            cov(1, 1) = 0.4 * 0.4;

            // Don't smear noise
            if (not isNoise) {
              Acts::ActsVector<2> par{m_ITKRadius[moduleGeoId.layer() / 2 - 1] *
                                              Acts::VectorHelpers::phi(pos) +
                                          0.1 * stdNormal(rng),
                                      pos.z() + 0.4 * stdNormal(rng)};
              measurements.emplace_back(
                  Acts::Measurement<Acts::BoundIndices, 2>(sourceLink, indices,
                                                           par, cov));
            } else {
              Acts::ActsVector<2> par{m_ITKRadius[moduleGeoId.layer() / 2 - 1] *
                                          Acts::VectorHelpers::phi(pos),
                                      pos.z()};
              measurements.emplace_back(
                  Acts::Measurement<Acts::BoundIndices, 2>(sourceLink, indices,
                                                           par, cov));
            }
          } else if (surfacePtr->type() == Acts::Surface::SurfaceType::Straw) {
            // The truth deposit energy
            double driftDistance = simHit.depositedEnergy();
            std::array<Acts::BoundIndices, 1> indices = {Acts::eBoundLoc0};

            // Gives the drift distance a sign if it's not noise
            if (not isNoise) {
              auto lpResult =
                  surfacePtr->globalToLocal(context.geoContext, pos, dir);
              if (not lpResult.ok()) {
                ACTS_FATAL("Global to local transformation did not succeed.");
                return ProcessCode::ABORT;
              }
              auto lPosition = lpResult.value();
              // auto driftDistance = std::copysign(simHit.fourPosition()[3],
              // lPosition[0]);
              driftDistance =
                  std::copysign(simHit.depositedEnergy(), lPosition[0]);
            }

            double sigma = 0.125;
            // std::uniform_real_distribution<double> uniform(0,1);
            // if(uniform(rng)<0.782){
            //   sigma= 0.09;
            // }  else {
            //   sigma = 0.24;
            // }

            Acts::ActsSymMatrix<1> cov = Acts::ActsSymMatrix<1>::Identity();
            cov(0, 0) = sigma * sigma;
            if (not isNoise) {
              Acts::ActsVector<1> par{driftDistance + sigma * stdNormal(rng)};
              measurements.emplace_back(
                  Acts::Measurement<Acts::BoundIndices, 1>(sourceLink, indices,
                                                           par, cov));
            } else {
              Acts::ActsVector<1> par{driftDistance};
              measurements.emplace_back(
                  Acts::Measurement<Acts::BoundIndices, 1>(sourceLink, indices,
                                                           par, cov));
            }
          } else {
            ACTS_ERROR("The surface type must be Cylinder or Straw");
            return ProcessCode::ABORT;
          }

          measurementParticlesMap.emplace_hint(measurementParticlesMap.end(),
                                               measurementIdx, particleId);
          measurementSimHitsMap.emplace_hint(measurementSimHitsMap.end(),
                                             measurementIdx, simHitIdx);
        }
      }
      for (size_t i = 0; i < measurements.size(); i++) {
        const auto sourceLink_ = sourceLinks.nth(i);
        auto geoId_ = sourceLink_->get().geometryId();
      }
      std::cout << "nITKHits = " << nITKHits << ", nMDCHits = " << nMDCHits
                << std::endl;
    } // Finish reading sim hits and transform them into measurements
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////
    // II. now read particles
    /////////////////////////////////////////////////////////////////////////////////////////////////
    int nParticles = 0;
    if (m_particle_treeReader->Next()) {
      // std::cout << "Reading event " << m_evtCounter++ << std::endl;
      // Reading truth particles
      for (size_t i = 0; i < particlePDG->GetSize(); ++i) {
        if (particleITKSigHitIdx[i] < 3) {
          continue;
        };
        nParticles++;
        int oncelayer = 0;
        int secondlayer = 0;
        int thirdlayer = 0;
        int particleITKnum = 0;
        auto particleITKlayer = particleITKLayer[i];
        for (size_t j = 0; j < particleITKlayer.size(); j++) {
          if (particleITKlayer[j] == 0) {
            oncelayer++;
          } else if (particleITKlayer[j] == 1) {
            secondlayer++;
          } else {
            thirdlayer++;
          }
        }
        if (oncelayer != 0 && secondlayer != 0 && thirdlayer != 0) {
          particleITKnum = 3;
        }
        if (oncelayer == 0 && secondlayer != 0 && thirdlayer != 0) {
          particleITKnum = 2;
        }
        if (oncelayer == 0 && secondlayer == 0 && thirdlayer != 0) {
          particleITKnum = 1;
        }

        Acts::Vector3 pos((*particleVertexX)[i], (*particleVertexY)[i],
                          (*particleVertexZ)[i]);
        Acts::Vector3 mom((*particleMomentumX)[i], (*particleMomentumY)[i],
                          (*particleMomentumZ)[i]);
        double charge = 0;
        if ((*particlePDG)[i] == 13 or (*particlePDG)[i] == 11 or
            (*particlePDG)[i] == -211 or (*particlePDG)[i] == -2212) {
          charge = -1;
        } else if ((*particlePDG)[i] == -13 or (*particlePDG)[i] == -11 or
                   (*particlePDG)[i] == 211 or (*particlePDG)[i] == 2212) {
          charge = 1;
        }
        ActsFatras::Particle particle(
            ActsFatras::Barcode(i), Acts::PdgParticle((*particlePDG)[i]),
            charge * Acts::UnitConstants::e,
            (*particleMass)[i] / 1000 * Acts::UnitConstants::GeV);
        // 1 means "Undefined"
        particle.setITKHits(particleITKnum);  //////******set  particle ITKHits
        particle.setProcess(static_cast<ActsFatras::ProcessType>(1));
        particle.setPosition4(pos.x() * Acts::UnitConstants::mm,
                              pos.y() * Acts::UnitConstants::mm,
                              pos.z() * Acts::UnitConstants::mm,
                              0 * Acts::UnitConstants::ns);
        //// Only used for direction; normalization/units do not matter
        particle.setDirection(mom.x(), mom.y(), mom.z());  // in GeV
        particle.setAbsoluteMomentum(std::sqrt(mom.x() * mom.x() +
                                               mom.y() * mom.y() +
                                               mom.z() * mom.z()) *
                                     Acts::UnitConstants::GeV);
        unordered_particles.push_back(std::move(particle));
      }

      for (const auto& [particleIndex, nHits] : particleHitIdx) {
        std::cout << "particle " << particleIndex << " has " << nHits
                  << " nHits"
                  << ", nITKSigHits = " << particleITKSigHitIdx[particleIndex]
                  << std::endl;
      }
      std::cout << "nParticles = " << nParticles << std::endl;
    }
    
    // Write the collections to the EventStore
    context.eventStore.add(m_cfg.outputSourceLinks, std::move(sourceLinks));
    context.eventStore.add(m_cfg.outputSourceLinks + "__storage",
                           std::move(sourceLinkStorage));
    context.eventStore.add(m_cfg.outputMeasurements, std::move(measurements));
    context.eventStore.add(m_cfg.outputMeasurementParticlesMap,
                           std::move(measurementParticlesMap));
    context.eventStore.add(m_cfg.outputMeasurementSimHitsMap,
                           std::move(measurementSimHitsMap));

    particles.insert(unordered_particles.begin(), unordered_particles.end());
    context.eventStore.add(m_cfg.outputParticles, std::move(particles));

    context.eventStore.add(m_cfg.outputSimHits, std::move(simHits));
  } else {
    ACTS_WARNING("Could not read in event.");
  }
  // Return success flag
  return ActsExamples::ProcessCode::SUCCESS;
}
