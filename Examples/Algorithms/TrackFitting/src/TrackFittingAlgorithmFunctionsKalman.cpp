// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/MagneticField/MagneticField.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

namespace {

using Updater = Acts::GainMatrixUpdater;
using Smoother = Acts::GainMatrixSmoother;
using Stepper = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
using Fitter = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;
using DirectPropagator = Acts::Propagator<Stepper, Acts::DirectNavigator>;
using DirectFitter =
    Acts::KalmanFitter<DirectPropagator, Acts::VectorMultiTrajectory>;

struct SimpleReverseFilteringLogic {
  double momentumThreshold;

  bool doBackwardFiltering(
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy
          trackState) const {
    auto momentum = fabs(1 / trackState.filtered()[Acts::eBoundQOverP]);
    return (momentum <= momentumThreshold);
  }
};

struct SimpleOutlierFinder {
  double chi2Cut;
  bool isOutlier(
      Acts::MultiTrajectory<Acts::VectorMultiTrajectory>::ConstTrackStateProxy
          state) const {
    if (chi2Cut <= 0) {
      return false;
    }

    // can't determine an outlier w/o a measurement or predicted parameters
    if (not state.hasCalibrated() or not state.hasPredicted()) {
      return false;
    }
    // const size_t measdim = state.calibratedSize();

    // auto geoId = state.referenceSurface().geometryId();

    auto H = state.effectiveProjector();
    auto res = state.effectiveCalibrated() - H * state.predicted();
    auto predictedCovariance = state.predictedCovariance();
    auto calibratedCovariance = state.effectiveCalibratedCovariance();

    auto cov = calibratedCovariance + H * predictedCovariance * H.transpose();

    auto chi2 = (res.transpose() * cov.inverse() * res).eval()(0, 0);
    if (chi2 > chi2Cut) {
      return true;
    }

    return false;
  }
};

template <typename TrackFitterFunktion>
auto makeKfOptions(
    const TrackFitterFunktion& f,
    ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions options) {
  Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> extensions;
  extensions.updater.connect<
      &Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
      &f.kfUpdater);
  extensions.smoother.connect<
      &Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
      &f.kfSmoother);
  extensions.reverseFilteringLogic
      .connect<&SimpleReverseFilteringLogic::doBackwardFiltering>(
          &f.reverseFilteringLogic);
  extensions.outlierFinder.connect<&SimpleOutlierFinder::isOutlier>(
      &f.outlierFinder);

  Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
      options.geoContext, options.magFieldContext, options.calibrationContext,
      extensions, options.logger, options.propOptions,
      &(*options.referenceSurface));

  kfOptions.multipleScattering = f.multipleScattering;
  kfOptions.energyLoss = f.energyLoss;
  kfOptions.freeToBoundCorrection = f.freeToBoundCorrection;

  return kfOptions;
}

struct TrackFitterFunctionImpl
    : public ActsExamples::TrackFittingAlgorithm::TrackFitterFunction {
  Fitter trackFitter;

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  SimpleReverseFilteringLogic reverseFilteringLogic;
  SimpleOutlierFinder outlierFinder;

  bool multipleScattering;
  bool energyLoss;
  Acts::FreeToBoundCorrection freeToBoundCorrection;

  TrackFitterFunctionImpl(Fitter&& f) : trackFitter(std::move(f)) {}

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsExamples::IndexSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options)
      const override {
    auto kfOptions = makeKfOptions(*this, options);
    kfOptions.extensions.calibrator
        .connect<&ActsExamples::MeasurementCalibrator::calibrate>(
            &options.calibrator.get());
    return trackFitter.fit(sourceLinks.begin(), sourceLinks.end(),
                           initialParameters, kfOptions);
  };
};

struct DirectedFitterFunctionImpl
    : public ActsExamples::TrackFittingAlgorithm::DirectedTrackFitterFunction {
  DirectFitter fitter;

  Acts::GainMatrixUpdater kfUpdater;
  Acts::GainMatrixSmoother kfSmoother;
  SimpleReverseFilteringLogic reverseFilteringLogic;
  SimpleOutlierFinder outlierFinder;

  bool multipleScattering;
  bool energyLoss;
  Acts::FreeToBoundCorrection freeToBoundCorrection;

  DirectedFitterFunctionImpl(DirectFitter&& f) : fitter(std::move(f)) {}

  ActsExamples::TrackFittingAlgorithm::TrackFitterResult operator()(
      const std::vector<std::reference_wrapper<
          const ActsExamples::IndexSourceLink>>& sourceLinks,
      const ActsExamples::TrackParameters& initialParameters,
      const ActsExamples::TrackFittingAlgorithm::GeneralFitterOptions& options,
      const std::vector<const Acts::Surface*>& sSequence) const override {
    auto kfOptions = makeKfOptions(*this, options);
    kfOptions.extensions.calibrator
        .connect<&ActsExamples::MeasurementCalibrator::calibrate>(
            &options.calibrator.get());
    return fitter.fit(sourceLinks.begin(), sourceLinks.end(), initialParameters,
                      kfOptions, sSequence);
  };
};

}  // namespace

std::shared_ptr<ActsExamples::TrackFittingAlgorithm::TrackFitterFunction>
ActsExamples::TrackFittingAlgorithm::makeKalmanFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold, double chi2Cut,
    Acts::FreeToBoundCorrection freeToBoundCorrection) {
  Stepper stepper(std::move(magneticField));
  Acts::Navigator::Config cfg{trackingGeometry};
  cfg.resolvePassive = false;
  cfg.resolveMaterial = true;
  cfg.resolveSensitive = true;
  Acts::Navigator navigator(cfg);
  Propagator propagator(std::move(stepper), std::move(navigator));
  Fitter trackFitter(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<TrackFitterFunctionImpl>(std::move(trackFitter));
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->outlierFinder.chi2Cut = chi2Cut;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}

std::shared_ptr<
    ActsExamples::TrackFittingAlgorithm::DirectedTrackFitterFunction>
ActsExamples::TrackFittingAlgorithm::makeKalmanFitterFunction(
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    bool multipleScattering, bool energyLoss,
    double reverseFilteringMomThreshold, double chi2Cut,
    Acts::FreeToBoundCorrection freeToBoundCorrection) {
  // construct all components for the fitter
  Stepper stepper(std::move(magneticField));
  Acts::DirectNavigator navigator;
  DirectPropagator propagator(std::move(stepper), navigator);
  DirectFitter fitter(std::move(propagator));

  // build the fitter functions. owns the fitter object.
  auto fitterFunction =
      std::make_shared<DirectedFitterFunctionImpl>(std::move(fitter));
  fitterFunction->multipleScattering = multipleScattering;
  fitterFunction->energyLoss = energyLoss;
  fitterFunction->reverseFilteringLogic.momentumThreshold =
      reverseFilteringMomThreshold;
  fitterFunction->outlierFinder.chi2Cut = chi2Cut;
  fitterFunction->freeToBoundCorrection = freeToBoundCorrection;

  return fitterFunction;
}
