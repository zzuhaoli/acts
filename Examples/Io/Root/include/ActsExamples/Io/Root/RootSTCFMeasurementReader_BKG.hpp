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
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"

#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TTreeReaderValue.h>

#include <mutex>
#include <vector>

class TTreeReader;
//class TTreeReaderArray;

namespace ActsExamples {

/// @class RootSTCFMeasurementReader
///
/// @brief Reads in TrackParameter information from a root file
/// and fills it into a Acts::BoundTrackParameter format
class RootSTCFMeasurementReader_BKG : public IReader {
 public:
  /// @brief The nested configuration struct
  struct Config {
    /// Output simulated (truth) hits collection.
    std::string outputSimHits;
    /// Which particle collection to read into.
    std::string outputParticles;
    /// Output source links collection.
    std::string outputSourceLinks = "sourcelinks_";
    /// Output measurements collection.
    std::string outputMeasurements = "measurements_";
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap = "measurement_particles_map";
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap = "measurement_simhits_map";
    /// Tracking geometry required to get approximated sim hits on surface.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// Random numbers service.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
   
    std::vector<std::string> v_treeName = {"events","info_withBKG"}; ///< name of the input tree 
    //std::string treeName = "events";  ///< name of the input tree
    //std::string filePath;                   ///< The name of the input file
    std::vector<std::string> v_filePath;   ///< The name of the input file
    /// Whether the events are ordered or not
    bool orderedEvents = true;
  };

  /// Constructor
  /// @param config The Configuration struct
  /// @param level The log level
  RootSTCFMeasurementReader_BKG(const Config& config, Acts::Logging::Level level);

  /// Destructor
  ~RootSTCFMeasurementReader_BKG();

  /// Framework name() method
  std::string name() const final override { return "RootSTCFMeasurementReader_BKG"; }

  /// Return the available events range.
  std::pair<size_t, size_t> availableEvents() const final override;

  /// Read out data from the input stream
  ///
  /// @param context The algorithm context
  ProcessCode read(
      const ActsExamples::AlgorithmContext& context) final override;

  /// Readonly access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// The logger
  std::unique_ptr<const Acts::Logger> m_logger;

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// The config class
  Config m_cfg;

  /// mutex used to protect multi-threaded reads
  std::mutex m_read_mutex;

  /// The number of events
  size_t m_events = 0;

  std::array<int, 48> m_MDCnCells = {128, 128, 128, 128, 128, 128, 160, 160, 160, 160, 160, 160, 192, 192, 192, 192, 192, 192, 224, 224, 224, 224, 224, 224, 256, 256, 256, 256, 256, 256, 288, 288, 288, 288, 288, 288, 320, 320, 320, 320, 320, 320, 352, 352, 352, 352, 352, 352};

  std::array<int, 48> m_MDCsCells = { 0,128,256,384,512,640,
  768,928,1088,1248,1408,1568,
  1728,1920,2112,2304,2496,2688,
  2880,3104,3328,3552,3776,4000,
  4224,4480,4736,4992,5248,5504,
  5760,6048,6336,6624,6912,7200,
  7488,7808,8128,8448,8768,9088,
  9408,9760,10112,10464,10816,11168};

  std::array<int, 2> m_volumeIDs = {9,11};
  //std::array<Acts::ActsScalar, 3> m_ITKRadius = {65.2276,115.228,165.228};
  std::array<Acts::ActsScalar, 3> m_ITKRadius = {65.1125,115.113,165.112};
//  std::map<int, std::pair<Acts::ActsScalar, Acts::ActsScalar>> m_ITKVariances = {
//	  {0, {}}, 
//	  {1, {}},
//	  {2, {}},
//  };

  size_t m_evtCounter =0;
  //std::map<int, std::vector<int>> m_particleITKLayer;
  /// The input tree name
  TTreeReader* m_particle_treeReader = nullptr;
  TTreeReader* m_simhit_treeReader = nullptr;

  TTreeReaderArray<int>* particlePDG = nullptr;
  TTreeReaderArray<float>* particleCharge = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleMass = nullptr;
  TTreeReaderArray<int>* particleTrackID = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleVertexX = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleVertexY = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* particleVertexZ = nullptr;
  TTreeReaderArray<float>* particleMomentumX = nullptr;
  TTreeReaderArray<float>* particleMomentumY = nullptr;
  TTreeReaderArray<float>* particleMomentumZ = nullptr;

  
  TTreeReaderArray<int>* ITKtype = nullptr;
  TTreeReaderArray<int>* ITKlayerID = nullptr;
  TTreeReaderArray<int>* ITKparentID = nullptr;
  TTreeReaderArray<int>* ITKparticleId = nullptr;
  TTreeReaderArray<float>* ITKpositionX = nullptr;
  TTreeReaderArray<float>* ITKpositionY = nullptr;
  TTreeReaderArray<float>* ITKpositionZ = nullptr;
  TTreeReaderArray<float>* ITKmomentumX = nullptr;
  TTreeReaderArray<float>* ITKmomentumY = nullptr;
  TTreeReaderArray<float>* ITKmomentumZ = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* ITKtime = nullptr;
  TTreeReaderArray<double>* ITKmass = nullptr;
  
  TTreeReaderArray<int>* MDCtype = nullptr;
  TTreeReaderArray<int>* MDCcellID = nullptr;
  TTreeReaderArray<int>* MDClayerID = nullptr;
  TTreeReaderArray<int>* MDCparentID = nullptr;
  TTreeReaderArray<int>* MDCparticleId = nullptr;
  TTreeReaderArray<float>* MDCpositionX = nullptr;
  TTreeReaderArray<float>* MDCpositionY = nullptr;
  TTreeReaderArray<float>* MDCpositionZ = nullptr;
  TTreeReaderArray<float>* MDCmomentumX = nullptr;
  TTreeReaderArray<float>* MDCmomentumY = nullptr;
  TTreeReaderArray<float>* MDCmomentumZ = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* MDCdriftDistance = nullptr;
  TTreeReaderArray<float>* MDCwirePoint1X = nullptr;
  TTreeReaderArray<float>* MDCwirePoint1Y = nullptr;
  TTreeReaderArray<float>* MDCwirePoint1Z = nullptr;
  TTreeReaderArray<float>* MDCwirePoint2X = nullptr;
  TTreeReaderArray<float>* MDCwirePoint2Y = nullptr;
  TTreeReaderArray<float>* MDCwirePoint2Z = nullptr;
  TTreeReaderArray<Acts::ActsScalar>* MDCtime = nullptr;
  TTreeReaderArray<double>* MDCmass = nullptr;
};

}  // namespace ActsExamples