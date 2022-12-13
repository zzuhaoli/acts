// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

double MeasurementSelector::calculateChi2(
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Measurement fullCalibrated,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::MeasurementCovariance fullCalibratedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Parameters predicted,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Covariance predictedCovariance,
    TrackStateTraits<MultiTrajectoryTraits::MeasurementSizeMax,
                     false>::Projector projector,
    unsigned int calibratedSize, const Surface& surface) const {
  return visit_measurement(
      fullCalibrated, fullCalibratedCovariance, calibratedSize,
      [&](const auto calibrated, const auto calibratedCovariance) -> double {
        constexpr size_t kMeasurementSize =
            decltype(calibrated)::RowsAtCompileTime;

        using ParametersVector = ActsVector<kMeasurementSize>;

        // Take the projector (measurement mapping function)
        const auto H =
            projector.template topLeftCorner<kMeasurementSize, eBoundSize>()
                .eval();

        // Get the residuals
        ParametersVector res;
        //ParametersVector res_;
	auto updatedCalibrated = calibrated;
        if(surface.type() == Surface::SurfaceType::Straw){
	   auto localParameters = H * predicted;
          updatedCalibrated[eBoundLoc0] = std::copysign(calibrated[eBoundLoc0], localParameters[eBoundLoc0]); 
	}	
        res = updatedCalibrated - H * predicted;
       
       	//res = calibrated - H * predicted;

        //auto chi2 = (res.transpose() *
        //        ((calibratedCovariance +
        //          H * predictedCovariance * H.transpose()))
        //            .inverse() *
        //        res)
        //    .eval()(0, 0);

        //auto chi2_ = (res_.transpose() *
        //        ((calibratedCovariance +
        //          H * predictedCovariance * H.transpose()))
        //            .inverse() *
        //        res_)
        //    .eval()(0, 0);
     
        //if(chi2_ > chi2)	
        //std::cout<<"chi2 =  "<<  chi2 <<", chi2_ " << chi2_ << std::endl; 	

	//std::cout<<"Surface " << surface.geometryId() <<", calibrated = " << calibrated.transpose();
        //std::cout<<"calibratedCovariance = " << calibratedCovariance << std::endl;
        //std::cout<<", H*predicted = " << (H*predicted).transpose();
        //std::cout<<", res = " << res.transpose() << std::endl;

        // Get the chi2
        return (res.transpose() *
                ((calibratedCovariance +
                  H * predictedCovariance * H.transpose()))
                    .inverse() *
                res)
            .eval()(0, 0);
      });
}

}  // namespace Acts
