// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>

ActsExamples::TrackFittingAlgorithm::TrackFittingAlgorithm(
    Config config, Acts::Logging::Level level)
    : ActsExamples::BareAlgorithm("TrackFittingAlgorithm", level),
      m_cfg(std::move(config)) {
  if (m_cfg.inputMeasurements.empty()) {
    throw std::invalid_argument("Missing input measurement collection");
  }
  if (m_cfg.inputSourceLinks.empty()) {
    throw std::invalid_argument("Missing input source links collection");
  }
  if (m_cfg.inputProtoTracks.empty()) {
    throw std::invalid_argument("Missing input proto tracks collection");
  }
  if (m_cfg.inputInitialTrackParameters.empty()) {
    throw std::invalid_argument(
        "Missing input initial track parameters collection");
  }
  if (not m_cfg.trackingGeometry) {
    throw std::invalid_argument("Missing tracking geometry");
  }
  if (m_cfg.outputTrajectories.empty()) {
    throw std::invalid_argument("Missing output trajectories collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TrackFittingAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Read input data
  const auto& measurements =
      ctx.eventStore.get<MeasurementContainer>(m_cfg.inputMeasurements);
  const auto& sourceLinks =
      ctx.eventStore.get<IndexSourceLinkContainer>(m_cfg.inputSourceLinks);
  const auto& protoTracks =
      ctx.eventStore.get<ProtoTrackContainer>(m_cfg.inputProtoTracks);
  const auto& initialParameters = ctx.eventStore.get<TrackParametersContainer>(
      m_cfg.inputInitialTrackParameters);

  // Consistency cross checks
  if (protoTracks.size() != initialParameters.size()) {
    ACTS_FATAL("Inconsistent number of proto tracks and parameters");
    return ProcessCode::ABORT;
  }

  // Prepare the output data with MultiTrajectory
  TrajectoriesContainer trajectories;
  trajectories.reserve(protoTracks.size());

  // Construct a perigee surface as the target surface
  auto pSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(
      Acts::Vector3{0., 0., 0.});

  // Measurement calibrator must be instantiated here, because we need the
  // measurements to construct it. The other extensions are hold by the
  // fit-function-object
  ActsExamples::MeasurementCalibrator calibrator(measurements);
  auto pOptions = Acts::PropagatorPlainOptions();
  pOptions.maxSteps = 10000;

  GeneralFitterOptions options{ctx.geoContext,   ctx.magFieldContext,
                               ctx.calibContext, calibrator,
                               &(*pSurface),     Acts::LoggerWrapper{logger()},
                               pOptions};

  // Perform the fit for each input track
  std::vector<std::reference_wrapper<const IndexSourceLink>> trackSourceLinks;
  std::vector<const Acts::Surface*> surfSequence;
  for (std::size_t itrack = 0; itrack < protoTracks.size(); ++itrack) {
    // Check if you are not in picking mode
    if (m_cfg.pickTrack > -1 and m_cfg.pickTrack != static_cast<int>(itrack)) {
      continue;
    }

    // The list of hits and the initial start parameters
    const auto& protoTrack = protoTracks[itrack];
    const auto& initialParams = initialParameters[itrack];

    // We can have empty tracks which must give empty fit results so the number
    // of entries in input and output containers matches.
    if (protoTrack.empty()) {
      trajectories.push_back(Trajectories());
      ACTS_WARNING("Empty track " << itrack << " found.");
      continue;
    }

    ACTS_VERBOSE("Initial parameters: "
                 << initialParams.fourPosition(ctx.geoContext).transpose()
                 << " -> " << initialParams.unitDirection().transpose());

    // Clear & reserve the right size
    trackSourceLinks.clear();
    trackSourceLinks.reserve(protoTrack.size());

    // Fill the source links via their indices from the container
    for (auto hitIndex : protoTrack) {
      if (auto it = sourceLinks.nth(hitIndex); it != sourceLinks.end()) {
        const IndexSourceLink& sourceLink = *it;
        auto geoId = sourceLink.geometryId();
        trackSourceLinks.push_back(std::cref(sourceLink));
        surfSequence.push_back(m_cfg.trackingGeometry->findSurface(geoId));
      } else {
        ACTS_FATAL("Proto track " << itrack << " contains invalid hit index"
                                  << hitIndex);
        return ProcessCode::ABORT;
      }
    }

    ACTS_DEBUG("Invoke fitter");
    auto result =
        fitTrack(trackSourceLinks, initialParams, options, surfSequence);

    if (result.ok()) {
      // Get the fit output object
      auto& fitOutput = result.value();
      // The track entry indices container. One element here.
      std::vector<Acts::MultiTrajectoryTraits::IndexType> trackTips;
      trackTips.reserve(1);
      trackTips.emplace_back(fitOutput.lastMeasurementIndex);
      // The fitted parameters container. One element (at most) here.
      Trajectories::IndexedParameters indexedParams;
      if (fitOutput.fittedParameters) {
        const auto& params = fitOutput.fittedParameters.value();
        ACTS_VERBOSE("Fitted paramemeters for track " << itrack);
        ACTS_VERBOSE("  " << params.parameters().transpose());
        // Push the fitted parameters to the container
        indexedParams.emplace(fitOutput.lastMeasurementIndex, params);
      } else {
        ACTS_DEBUG("No fitted paramemeters for track " << itrack);
      }
      // store the result
      trajectories.emplace_back(fitOutput.fittedStates, std::move(trackTips),
                                std::move(indexedParams));
    } else {
      ACTS_WARNING("Fit failed for track "
                   << itrack << " with error: " << result.error() << ", "
                   << result.error().message());
      // Fit failed. Add an empty result so the output container has
      // the same number of entries as the input.
      trajectories.push_back(Trajectories());
    }
  }

  ctx.eventStore.add(m_cfg.outputTrajectories, std::move(trajectories));
  return ActsExamples::ProcessCode::SUCCESS;
}
