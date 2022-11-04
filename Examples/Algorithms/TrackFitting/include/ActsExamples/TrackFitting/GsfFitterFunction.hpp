// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"

namespace ActsExamples {

std::shared_ptr<TrackFittingAlgorithm::TrackFitterFunction>
makeGsfFitterFunction(
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
    std::size_t maxComponents, bool abortOnError,
    bool disableAllMaterialHandling);

}  // namespace ActsExamples