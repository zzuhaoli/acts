// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"

#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/Utilities/Range.hpp"

#include <algorithm>
#include <stdexcept>
#include <vector>
#include <map>

using namespace ActsExamples;

TruthTrackFinder::TruthTrackFinder(const Config& config,
                                   Acts::Logging::Level level)
    : BareAlgorithm("TruthTrackFinder", level), m_cfg(config) {
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input truth particles collection");
  }
  if (m_cfg.inputMeasurementParticlesMap.empty()) {
    throw std::invalid_argument("Missing input hit-particles map collection");
  }
  ///////////////////////////////// 
  if (m_cfg.inputSimulatedHits.empty()) {
    throw std::invalid_argument("Missing input simulated hits collection");
  }
  if (m_cfg.inputMeasurementSimHitsMap.empty()) {
    throw std::invalid_argument("Missing input measurement sim hits map");
  }
  ///////////////////////////////// 
  if (m_cfg.outputProtoTracks.empty()) {
    throw std::invalid_argument("Missing output proto tracks collection");
  }
}

ProcessCode TruthTrackFinder::execute(const AlgorithmContext& ctx) const {
  using HitParticlesMap = IndexMultimap<ActsFatras::Barcode>;
  using HitSimHitsMap = IndexMultimap<Index>;

  // prepare input collections
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);
  const auto& hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputMeasurementParticlesMap);
  
  const auto& simHits =
      ctx.eventStore.get<SimHitContainer>(m_cfg.inputSimulatedHits);
  const auto& simHitsMap =
      ctx.eventStore.get<HitSimHitsMap>(m_cfg.inputMeasurementSimHitsMap);
  
  // compute particle_id -> {hit_id...} map from the
  // hit_id -> {particle_id...} map on the fly.
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);

  // prepare output collection
  ProtoTrackContainer tracks;
  tracks.reserve(particles.size());

  ACTS_VERBOSE("Create prototracks for " << particles.size() << " particles");
  for (const auto& particle : particles) {
    // find the corresponding hits for this particle
    const auto& hits =
        makeRange(particleHitsMap.equal_range(particle.particleId()));
    ACTS_VERBOSE(" - Prototrack from " << hits.size() << " hits");
 

    // fill hit indices to create the proto track
    ProtoTrack track;
    if(m_cfg.removeHitsFromLoops){
    
      std::map<const double, std::pair<const Index, Acts::GeometryIdentifier>> trackHitList;
      for (const auto& hit : hits) {
        const auto measurementIndex = hit.second;
        const auto simHitIndex = simHitsMap.find(measurementIndex)->second;
        const auto simHit = simHits.nth(simHitIndex);
        const auto geoId = simHit->geometryId();
        const auto simHitTime = simHit->time();
        trackHitList.insert(std::make_pair(simHitTime, std::make_pair(measurementIndex, geoId)));
      }
      std::cout<<"There are " << trackHitList.size() << " sim hits for the particle " << std::endl;

      // remove the hits from looping tracks from the list
      std::map<const Acts::GeometryIdentifier, const Index> trackHitFilteredList;
      // The hit in trackHitList are sorted in time 
      for(const auto& [time, measurement] : trackHitList) {
         const auto measurementIndex = measurement.first;
         const auto geoId = measurement.second;
           // only store hit with surface that is not intersected yet 
	   if(trackHitFilteredList.find(geoId)!=trackHitFilteredList.end()){
              break; 
           } else {
	      trackHitFilteredList.emplace(geoId, measurementIndex);
	   }
      }
      std::cout<<"After removing hits from loops, there are " << trackHitFilteredList.size() << " sim hits for the particle " << std::endl;

      track.reserve(trackHitFilteredList.size());
      for (const auto& [geoId, measurementIndex] : trackHitFilteredList) {
        track.emplace_back(measurementIndex);
      }
    
    } else{  
      track.reserve(hits.size());
      for (const auto& hit : hits) {
        const auto simHitIndex = simHitsMap.find(hit.second)->second;
        const auto simHit = simHits.nth(simHitIndex);
        const auto geoId = simHit->geometryId();
	track.emplace_back(hit.second);
      }
    }
    
    // add proto track to the output collection
    tracks.emplace_back(std::move(track));
  }

  ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(tracks));
  return ProcessCode::SUCCESS;
}
