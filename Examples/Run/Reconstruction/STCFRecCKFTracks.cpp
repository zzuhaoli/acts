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
  auto inputfiles = vm["input-files"].template as<std::vector<std::string>>();
  if (inputfiles.size() > 1) {
    throw std::invalid_argument("Sorry, only one input file is possible here");
  } else if (inputfiles.size() < 1) {
    throw std::invalid_argument("Sorry, please provide one input file");
  }
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

  // Read truth hits from CSV files
  RootSTCFMeasurementReader::Config STCFMeasurementReaderCfg;
  STCFMeasurementReaderCfg.filePath = inputDir + "/" + inputfiles.front();
  STCFMeasurementReaderCfg.outputSimHits = "hits";
  STCFMeasurementReaderCfg.outputParticles = "particles";
  STCFMeasurementReaderCfg.trackingGeometry = trackingGeometry;
  STCFMeasurementReaderCfg.randomNumbers = rnd;
  sequencer.addReader(std::make_shared<RootSTCFMeasurementReader>(
      STCFMeasurementReaderCfg, logLevel));

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

  // Create starting parameters from either particle smearing or combined seed
  // finding and track parameters estimation
  std::string outputTrackParameters;
  if (truthSmearedSeeded) {
    // Run the particle smearing
    auto particleSmearingCfg =
        setupParticleSmearing(vm, sequencer, rnd, inputParticles);
    outputTrackParameters = particleSmearingCfg.outputTrackParameters;
  } else {
    // Create space points
    SpacePointMaker::Config spCfg = Options::readSpacePointMakerConfig(vm);
    spCfg.inputSourceLinks = STCFMeasurementReaderCfg.outputSourceLinks;
    spCfg.inputMeasurements = STCFMeasurementReaderCfg.outputMeasurements;
    spCfg.outputSpacePoints = "spacepoints";
    spCfg.trackingGeometry = trackingGeometry;
    sequencer.addAlgorithm(std::make_shared<SpacePointMaker>(spCfg, logLevel));

    // Create seeds (i.e. proto tracks) using either truth track finding or seed
    // finding algorithm
    std::string inputProtoTracks = "";
    std::string inputSeeds = "";
    if (truthEstimatedSeeded) {
      // Truth track finding algorithm
      TruthTrackFinder::Config trackFinderCfg;
      trackFinderCfg.inputParticles = inputParticles;
      trackFinderCfg.inputMeasurementParticlesMap =
          STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
      trackFinderCfg.inputMeasurementSimHitsMap =
          STCFMeasurementReaderCfg.outputMeasurementSimHitsMap;
      trackFinderCfg.inputSimulatedHits =
          STCFMeasurementReaderCfg.outputSimHits;

      trackFinderCfg.outputProtoTracks = "prototracks";
      sequencer.addAlgorithm(
          std::make_shared<TruthTrackFinder>(trackFinderCfg, logLevel));
      inputProtoTracks = trackFinderCfg.outputProtoTracks;
    } else {
      // Seeding algorithm
      SeedingAlgorithm::Config seedingCfg;
      seedingCfg.inputSpacePoints = {
          spCfg.outputSpacePoints,
      };
      seedingCfg.outputSeeds = "seeds";
      seedingCfg.outputProtoTracks = "prototracks";

      seedingCfg.gridConfig.rMax = 200._mm;
      seedingCfg.seedFinderConfig.rMax = seedingCfg.gridConfig.rMax;

      seedingCfg.seedFilterConfig.deltaRMin =
          1._mm * vm["seed-deltar-min"].template as<double>();
      seedingCfg.seedFinderConfig.deltaRMin =
          seedingCfg.seedFilterConfig.deltaRMin;
      // seedingCfg.seedFinderConfig.deltaZMax = 20_mm;

      seedingCfg.gridConfig.deltaRMax =
          1._mm * vm["seed-deltar-max"].template as<double>();
      seedingCfg.seedFinderConfig.deltaRMax = seedingCfg.gridConfig.deltaRMax;

      seedingCfg.seedFinderConfig.maxPtScattering =
          vm["seed-max-ptscattering"].template as<double>() *
          Acts::UnitConstants::GeV;

      seedingCfg.seedFinderConfig.collisionRegionMin = -250._mm;
      seedingCfg.seedFinderConfig.collisionRegionMax = 250._mm;

      seedingCfg.gridConfig.zMin = -500._mm;
      seedingCfg.gridConfig.zMax = 500._mm;
      seedingCfg.seedFinderConfig.zMin = seedingCfg.gridConfig.zMin;
      seedingCfg.seedFinderConfig.zMax = seedingCfg.gridConfig.zMax;

      seedingCfg.seedFilterConfig.maxSeedsPerSpM = 1;
      seedingCfg.seedFinderConfig.maxSeedsPerSpM =
          seedingCfg.seedFilterConfig.maxSeedsPerSpM;

      seedingCfg.gridConfig.cotThetaMax =
          vm["seed-cottheta-max"].template as<double>();  // 1.75 eta
      seedingCfg.seedFinderConfig.cotThetaMax =
          seedingCfg.gridConfig.cotThetaMax;

      seedingCfg.seedFinderConfig.sigmaScattering =
          vm["seed-sigma-scattering"].template as<double>();
      // seedingCfg.seedFinderConfig.sigmaScattering = 10;
      // seedingCfg.seedFinderConfig.sigmaScattering = 2; // 2 for 75 and 100
      // mev
      seedingCfg.seedFinderConfig.radLengthPerSeed =
          vm["seed-rad-length-per-seed"].template as<double>();
      seedingCfg.maxSeeds = vm["seed-max-seeds"].template as<int>();

      seedingCfg.gridConfig.minPt = 40._MeV;
      seedingCfg.seedFinderConfig.minPt = seedingCfg.gridConfig.minPt;

      seedingCfg.gridConfig.bFieldInZ = 1.0_T;
      seedingCfg.seedFinderConfig.bFieldInZ = seedingCfg.gridConfig.bFieldInZ;

      seedingCfg.seedFinderConfig.beamPos = {0_mm, 0_mm};

      // seedingCfg.seedFinderConfig.impactMax = 3._mm;
      seedingCfg.seedFinderConfig.impactMax =
          vm["seed-impact-max"].template as<double>() * 1._mm;

      sequencer.addAlgorithm(
          std::make_shared<SeedingAlgorithm>(seedingCfg, logLevel));
      inputProtoTracks = seedingCfg.outputProtoTracks;
      inputSeeds = seedingCfg.outputSeeds;
    }

    // write track finding/seeding performance
    TrackFinderPerformanceWriter::Config tfPerfCfg;
    tfPerfCfg.inputProtoTracks = inputProtoTracks;
    // using selected particles
    tfPerfCfg.inputParticles = inputParticles;
    tfPerfCfg.inputMeasurementParticlesMap =
        STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
    // tfPerfCfg.filePath = outputDir + "/muon_performance_seeding_trees.root";
    tfPerfCfg.filePath = outputDir + "/performance_seeding_trees20w.root";
    sequencer.addWriter(
        std::make_shared<TrackFinderPerformanceWriter>(tfPerfCfg, logLevel));

    SeedingPerformanceWriter::Config seedPerfCfg;
    seedPerfCfg.inputProtoTracks = inputProtoTracks;
    seedPerfCfg.inputParticles = inputParticles;
    seedPerfCfg.inputMeasurementParticlesMap =
        STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
    // seedPerfCfg.filePath = outputDir +
    // "/muon_performance_seeding_hists.root";
    seedPerfCfg.filePath = outputDir + "/performance_seeding_hists20w.root";
    seedPerfCfg.effPlotToolConfig.varBinning["Eta"] =
        PlotHelpers::Binning("#eta", etaRange[0], etaRange[1], etaRange[2]);
    seedPerfCfg.effPlotToolConfig.varBinning["Pt"] =
        PlotHelpers::Binning("pT [GeV/c]", ptRange[0], ptRange[1], ptRange[2]);
    seedPerfCfg.duplicationPlotToolConfig.varBinning["Eta"] =
        PlotHelpers::Binning("#eta", etaRange[0], etaRange[1], etaRange[2]);
    seedPerfCfg.duplicationPlotToolConfig.varBinning["Pt"] =
        PlotHelpers::Binning("pT [GeV/c]", ptRange[0], ptRange[1], ptRange[2]);
    sequencer.addWriter(
        std::make_shared<SeedingPerformanceWriter>(seedPerfCfg, logLevel));

    // Algorithm estimating track parameter from seed
    TrackParamsEstimationAlgorithm::Config paramsEstimationCfg;
    paramsEstimationCfg.inputSeeds = inputSeeds;
    paramsEstimationCfg.inputProtoTracks = inputProtoTracks;
    paramsEstimationCfg.inputSpacePoints = {
        spCfg.outputSpacePoints,
    };
    paramsEstimationCfg.inputSourceLinks =
        STCFMeasurementReaderCfg.outputSourceLinks;
    paramsEstimationCfg.outputTrackParameters = "estimatedparameters";
    paramsEstimationCfg.outputProtoTracks = "prototracks_estimated";
    paramsEstimationCfg.trackingGeometry = trackingGeometry;
    paramsEstimationCfg.magneticField = magneticField;
    paramsEstimationCfg.bFieldMin = 0.1_T;
    paramsEstimationCfg.deltaRMax = 100._mm;
    paramsEstimationCfg.deltaRMin = 30._mm;
    paramsEstimationCfg.sigmaLoc0 = 100._um;
    paramsEstimationCfg.sigmaLoc1 = 400._um;
    paramsEstimationCfg.sigmaPhi = 0.2_degree;
    paramsEstimationCfg.sigmaTheta = 0.2_degree;
    paramsEstimationCfg.sigmaQOverP = 0.1 / 1._GeV;
    paramsEstimationCfg.sigmaT0 = 1400._s;
    paramsEstimationCfg.initialVarInflation =
        vm["ckf-initial-variance-inflation"].template as<Options::Reals<6>>();

    sequencer.addAlgorithm(std::make_shared<TrackParamsEstimationAlgorithm>(
        paramsEstimationCfg, logLevel));

    outputTrackParameters = paramsEstimationCfg.outputTrackParameters;
  }

  // Setup the track finding algorithm with CKF
  // It takes all the source links created from truth hit smearing, seeds from
  // truth particle smearing and source link selection config
  auto trackFindingCfg = Options::readTrackFindingConfig(vm);
  trackFindingCfg.inputMeasurements =
      STCFMeasurementReaderCfg.outputMeasurements;
  trackFindingCfg.inputSourceLinks = STCFMeasurementReaderCfg.outputSourceLinks;
  trackFindingCfg.inputInitialTrackParameters = outputTrackParameters;
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

  /*生成的trackststes_ckf文件
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
  trackSummaryWriter.inputParticles = STCFMeasurementReaderCfg.outputParticles;
  trackSummaryWriter.inputMeasurementParticlesMap =
      STCFMeasurementReaderCfg.outputMeasurementParticlesMap;
  // trackSummaryWriter.filePath = outputDir + "/muon_tracksummary_ckf.root";
  trackSummaryWriter.filePath = outputDir + "/tracksummary_ckf20w.root";
  trackSummaryWriter.treeName = "tracksummary";
  sequencer.addWriter(std::make_shared<RootTrajectorySummaryWriter>(
      trackSummaryWriter, logLevel));

  // Write CKF performance data
  CKFPerformanceWriter::Config perfWriterCfg;
  perfWriterCfg.inputParticles = inputParticles;  // after select
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
  // perfWriterCfg.filePath = outputDir + "/muon_performance_ckf.root";
  perfWriterCfg.filePath = outputDir + "/performance_ckf20w.root";
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
