// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Csv/CsvMultiTrajectoryWriter.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsReader.hpp"
#include "ActsExamples/Io/Csv/CsvOptionsWriter.hpp"
#include "ActsExamples/Io/Csv/CsvParticleReader.hpp"
#include "ActsExamples/Io/Csv/CsvSimHitReader.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootSTCFMeasurementReader.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectorySummaryWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Reconstruction/ReconstructionBase.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/SpacePointMakerOptions.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"
#include "ActsExamples/TrackFinding/TrackParamsEstimationAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingOptions.hpp"
#include "ActsExamples/TruthTracking/TruthSeedSelector.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Helpers.hpp"
#include "ActsExamples/Utilities/Options.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Definitions/Units.hpp>

#include <memory>
#include <boost/filesystem.hpp>

#include "ActsExamples/Io/Root/RootHoughCanTracksReader.hpp"

using namespace Acts::UnitLiterals;
using namespace ActsExamples;
using namespace boost::filesystem;
using namespace std::placeholders;

void addRecCKFOptions(ActsExamples::Options::Description& desc) {
  using namespace ActsExamples;
  using boost::program_options::bool_switch;
  using boost::program_options::value;

  auto opt = desc.add_options();
  opt("ckf-truth-smeared-seeds", bool_switch(),
      "Use track parameters smeared from truth particles for steering CKF");
  opt("ckf-truth-estimated-seeds", bool_switch(),
      "Use track parameters estimated from truth tracks for steering CKF");
  opt("perf-eta-range",
      value<ActsExamples::Options::Reals<3>>()->default_value(
          {{35, -1.75, 1.75}}),
      "Eta bins, min and max for plotting the performance, must be "
      "of form i:j:k.");
  opt("perf-pt-range",
      value<ActsExamples::Options::Reals<3>>()->default_value({{20, 0.0, 2.0}}),
      "pT bins, min and max for plotting the performance, must be "
      "of form i:j.");
 opt("ckf-prop-steps", value<int>()->default_value(10000),
                    "The max propagation steps during ckf.");
 opt("seed-sigma-scattering", value<double>()->default_value(200),
                    "The seeding sigma scattering.");
 opt("seed-rad-length-per-seed", value<double>()->default_value(0.1),
                    "The seeding radLengthPerSeed.");
 opt("seed-max-ptscattering", value<double>()->default_value(10),
                    "The seeding maxPtScattering.");
 opt("seed-impact-max", value<double>()->default_value(10),
                    "The seeding2impactMax.");
 opt("seed-deltar-max", value<double>()->default_value(80),
                    "The seeding deltaRMax.");
 opt("seed-deltar-min", value<double>()->default_value(20),
                    "The seeding deltaRMin.");
 opt("seed-cottheta-max", value<double>()->default_value(2.74),
                    "The seeding cotThetaMax.");
 opt("seed-max-seeds", value<int>()->default_value(2),
                    "The maximum number of seeds per event");
 opt("ckf-prop-tolerance", value<double>()->default_value(0.0001),
                    "The stepper tolerance during ckf.");
 opt("ckf-prop-mass", value<double>()->default_value(139.57018),
                    "The particle mass during ckf.");
}

int main(int argc, char* argv[]) {
  auto detector = std::make_shared<ActsExamples::TGeoDetector>();

  // setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  Options::addSequencerOptions(desc);
  Options::addRandomNumbersOptions(desc);
  Options::addGeometryOptions(desc);
  Options::addMaterialOptions(desc);
  Options::addInputOptions(desc);
  Options::addOutputOptions(desc,
                            OutputFormat::Csv | OutputFormat::DirectoryOnly);
  detector->addOptions(desc);
  Options::addMagneticFieldOptions(desc);
  // Options::addFittingOptions(desc);
  Options::addTrackFindingOptions(desc);
  addRecCKFOptions(desc);
  Options::addParticleSmearingOptions(desc);
  Options::addSpacePointMakerOptions(desc);
  Options::addCsvWriterOptions(desc);

  auto vm = Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

  Sequencer sequencer(Options::readSequencerConfig(vm));

  // Read some standard options
  auto logLevel = Options::readLogLevel(vm);
  auto inputDir = vm["input-dir"].as<std::string>();
  //*****************************************************************
  auto inputfiles = vm["input-files"].template as<std::vector<std::string>>();
  if (inputfiles.size() > 2) {
    throw std::invalid_argument("Sorry, only two input files are possible here");
  } else if (inputfiles.size() < 2) {
    throw std::invalid_argument("Sorry, please provide two input files");
  }
  //*******************************************************************
  auto outputDir = ensureWritableDirectory(vm["output-dir"].as<std::string>());
  auto rnd = std::make_shared<ActsExamples::RandomNumbers>(
      Options::readRandomNumbersConfig(vm));
  bool truthSmearedSeeded = vm["ckf-truth-smeared-seeds"].template as<bool>();
  bool truthEstimatedSeeded =
      vm["ckf-truth-estimated-seeds"].template as<bool>();
  auto etaRange = vm["perf-eta-range"].template as<Options::Reals<3>>();
  auto ptRange = vm["perf-pt-range"].template as<Options::Reals<3>>();

  // Setup detector geometry
  auto geometry = Geometry::build(vm, *detector);
  auto trackingGeometry = geometry.first;
  // Add context decorators
  for (auto cdr : geometry.second) {
    sequencer.addContextDecorator(cdr);
  }
  // Setup the magnetic field
  auto magneticField = Options::readMagneticField(vm);
  //****************************************************
  // Read seed after houghtransform from CSV files
  RootHoughCanTracksReader::Config RootHoughCanTracksReaderCfg;
  RootHoughCanTracksReaderCfg.outputcantrackparameters="initialtracks_parameters";
  RootHoughCanTracksReaderCfg.filePath=inputDir + "/" + inputfiles.at(1);
  RootHoughCanTracksReaderCfg.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<RootHoughCanTracksReader>(
      RootHoughCanTracksReaderCfg, logLevel));
  //*****************************************************
  // Read truth hits from CSV files
  RootSTCFMeasurementReader::Config STCFMeasurementReaderCfg;
  STCFMeasurementReaderCfg.filePath = inputDir + "/" + inputfiles.front();
  STCFMeasurementReaderCfg.outputSimHits = "hits";
  STCFMeasurementReaderCfg.outputParticles = "particles";
  STCFMeasurementReaderCfg.trackingGeometry = trackingGeometry;
  STCFMeasurementReaderCfg.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<RootSTCFMeasurementReader>(
      STCFMeasurementReaderCfg, logLevel));
  
  //****************************************************
  
  // Run the particle selection
  // The pre-selection will select truth particles satisfying provided criteria
  // from all particles read in by particle reader for further processing. It
  // has no impact on the truth hits read-in by the cluster reader.
  TruthSeedSelector::Config particleSelectorCfg;
  particleSelectorCfg.inputParticles = STCFMeasurementReaderCfg.outputParticles;
  particleSelectorCfg.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  particleSelectorCfg.outputParticles = "particles_selected";
  particleSelectorCfg.ptMin = 40_MeV;
  particleSelectorCfg.nHitsMin = 5;
  sequencer.addAlgorithm(
      std::make_shared<TruthSeedSelector>(particleSelectorCfg, logLevel));

  // The selected particles
  const auto& inputParticles = particleSelectorCfg.outputParticles;

  // Setup the track finding algorithm with CKF
  // It takes all the source links created from truth hit smearing, seeds from
  // truth particle smearing and source link selection config
  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
  trackFindingCfg.inputMeasurements =
      STCFMeasurementReaderCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = STCFMeasurementReaderCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = RootHoughCanTracksReaderCfg.outputcantrackparameters;//initialtracks_parameters from houghtransform
  trackFindingCfg.outputTrajectories = "trajectories";
  trackFindingCfg.outputTrackParameters = "parameters";
  trackFindingCfg.computeSharedHits = true;
  trackFindingCfg.maxPropSteps = vm["ckf-prop-steps"].template as<int>();
  trackFindingCfg.tolerance = vm["ckf-prop-tolerance"].template as<double>();
  trackFindingCfg.mass = vm["ckf-prop-mass"].template as<double>();
  trackFindingCfg.findTracks = TrackFindingAlgorithm::makeTrackFinderFunction(
      trackingGeometry, magneticField);
  sequencer.addAlgorithm(
      std::make_shared<TrackFindingAlgorithm>(trackFindingCfg, logLevel));

  std::cout << "Added CKF " << std::endl;

/*
  // write track states from CKF
  RootTrajectoryStatesWriter::Config trackStatesWriter;
  trackStatesWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. This could be avoided when a seperate track
  // selection algorithm is used.
  trackStatesWriter.inputParticles = STCFMeasurementReaderCfg.outputParticles;
  trackStatesWriter.inputSimHits = STCFMeasurementReaderCfg.outputSimHits;
  trackStatesWriter.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  trackStatesWriter.inputMeasurementSimHitsMap =
      STCFMeasurementReaderCfg.outputMeasurementSimHitsMap;
  trackStatesWriter.filePath = outputDir + "/trackstates_ckf.root";
  trackStatesWriter.treeName = "trackstates";
  sequencer.addWriter(std::make_shared<RootTrajectoryStatesWriter>(
      trackStatesWriter, logLevel));
*/

  // write track summary from CKF
  RootTrajectorySummaryWriter::Config trackSummaryWriter;
  trackSummaryWriter.inputTrajectories = trackFindingCfg.outputTrajectories;
  // @note The full particles collection is used here to avoid lots of warnings
  // since the unselected CKF track might have a majority particle not in the
  // filtered particle collection. This could be avoided when a seperate track
  // selection algorithm is used.
  trackSummaryWriter.inputParticles = inputParticles;//after particleSelect
  trackSummaryWriter.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  trackSummaryWriter.filePath = outputDir + "/tracksummary_ckf_houghseed.root";
  trackSummaryWriter.treeName = "tracksummary_houghseed";
  sequencer.addWriter(std::make_shared<RootTrajectorySummaryWriter>(
      trackSummaryWriter, logLevel));

  // Write CKF performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = STCFMeasurementReaderCfg.outputParticles;
  perfWriterCfg.inputTrajectories = trackFindingCfg.outputTrajectories;
  perfWriterCfg.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  // The bottom seed could be the first, second or third hits on the truth
  // track?
  perfWriterCfg.nMeasurementsMin = particleSelectorCfg.nHitsMin;
  perfWriterCfg.ptMin = 40_MeV;
  perfWriterCfg.effPlotToolConfig.varBinning["Eta"] =
      PlotHelpers::Binning("#eta", etaRange[0], etaRange[1], etaRange[2]);
  perfWriterCfg.effPlotToolConfig.varBinning["Pt"] =
      PlotHelpers::Binning("pT [GeV/c]", ptRange[0], ptRange[1], ptRange[2]);
  perfWriterCfg.duplicationPlotToolConfig.varBinning["Eta"] =
      PlotHelpers::Binning("#eta", etaRange[0], etaRange[1], etaRange[2]);
  perfWriterCfg.duplicationPlotToolConfig.varBinning["Pt"] =
      PlotHelpers::Binning("pT [GeV/c]", ptRange[0], ptRange[1], ptRange[2]);
  perfWriterCfg.fakeRatePlotToolConfig.varBinning["Eta"] =
      PlotHelpers::Binning("#eta", etaRange[0], etaRange[1], etaRange[2]);
  perfWriterCfg.fakeRatePlotToolConfig.varBinning["Pt"] =
      PlotHelpers::Binning("pT [GeV/c]", ptRange[0], ptRange[1], ptRange[2]);
  perfWriterCfg.trackSummaryPlotToolConfig.varBinning["Eta"] =
      PlotHelpers::Binning("#eta", etaRange[0], etaRange[1], etaRange[2]);
  perfWriterCfg.trackSummaryPlotToolConfig.varBinning["Pt"] =
      PlotHelpers::Binning("pT [GeV/c]", ptRange[0], ptRange[1], ptRange[2]);
  perfWriterCfg.trackSummaryPlotToolConfig.varBinning["Num"] =
      PlotHelpers::Binning("N", 60, -0.5, 59.5);
  perfWriterCfg.filePath = outputDir + "/performance_ckf_houghseed.root";
  sequencer.addWriter(
      std::make_shared<CKFPerformanceWriter>(perfWriterCfg, logLevel));

  if (vm["output-csv"].template as<bool>()) {
    // Write the CKF track as Csv
    CsvMultiTrajectoryWriter::Config trackWriterCsvConfig;
    trackWriterCsvConfig.inputTrajectories = trackFindingCfg.outputTrajectories;
    trackWriterCsvConfig.outputDir = outputDir;
    trackWriterCsvConfig.inputMeasurementParticlesMap =
        STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
    sequencer.addWriter(std::make_shared<CsvMultiTrajectoryWriter>(
        trackWriterCsvConfig, logLevel));
  }

  return sequencer.run();
}





