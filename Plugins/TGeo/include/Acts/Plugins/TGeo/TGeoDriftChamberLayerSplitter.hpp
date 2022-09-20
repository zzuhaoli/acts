// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/TGeo/ITGeoDetectorElementSplitter.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <vector>

class TGeoNode;

namespace Acts {

class TGeoDetectorElement;

/// @brief TGeoDriftChamberLayerSplitter
///
/// Split Cylinder and disks into submodules
class TGeoDriftChamberLayerSplitter : public ITGeoDetectorElementSplitter {
 public:
  /// Nested configuration struct
  struct Config {
    /// Number of segments in phi
    int phiSegments = -1;

    float unitScalar = 10;
  };

  /// Constructor
  ///
  /// @param cfg the configuration struct
  /// @param logger the logging object
  TGeoDriftChamberLayerSplitter(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "TGeoDriftChamberLayerSplitter", Acts::Logging::INFO));

  virtual ~TGeoDriftChamberLayerSplitter() = default;

  /// Take a geometry context and TGeoElement and split it into sub elements
  ///
  /// @param gctx is a geometry context object
  /// @param tgde is a TGeoDetectorElement that is eventually split
  ///
  /// @note If no split is performed the unsplit detector element is returned
  ///
  /// @return a vector of TGeoDetectorElement objects
  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>> split(
      const GeometryContext& gctx,
      std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const;

 private:
  Config m_cfg;

  /// Private access to the logger
  const Acts::Logger& logger() const { return *m_logger; }

  /// Logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace Acts
