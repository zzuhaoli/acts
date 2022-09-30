// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoURwellLayerSplitter.hpp"

#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RadialBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Surfaces/TrapezoidBounds.hpp"

#include "TGeoTube.h"

Acts::TGeoURwellLayerSplitter::TGeoURwellLayerSplitter(
    const TGeoURwellLayerSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoURwellLayerSplitter::split(
    const GeometryContext& /*gctx*/,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  // Thickness
  auto tgIdentifier = tgde->identifier();
  const TGeoNode& tgNode = tgde->tgeoNode();
  std::cout << "TGeoURwellLayerSplitter splitting detElement "
            << tgNode.GetName() << std::endl;

  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
      tgDetectorElements = {};

  // The daughters of the 'base' on this layer
  auto daugthers = tgNode.GetVolume()->GetNodes();
  TIter iObj(daugthers);
  while (TObject* obj = iObj()) {
    TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
    if (node != nullptr) {
      TGeoTube* daughter =
          dynamic_cast<TGeoTube*>(node->GetVolume()->GetShape());
      if (daughter == nullptr) {
        std::cout << "cast bad (base daughters is always a tube. This should "
                     "not happen)"
                  << std::endl;
      } else {
        std::string nodeName = node->GetName();
        // make a new detElement using a tube having the same minR, maxR and
        // halfZ as the uRwell top layer
        if (nodeName.find("topStrip_") != std::string::npos) {
          std::cout << tgNode.GetName() << " has daughter with  " << nodeName
                    << std::endl;
          ActsScalar parameters[5];
          daughter->GetBoundingCylinder(parameters);
          ActsScalar minR = daughter->GetRmin() * m_cfg.unitScalar;
          ActsScalar maxR = daughter->GetRmax() * m_cfg.unitScalar;
          ActsScalar deltaR = maxR - minR;
          // Maybe we should use the base thickness?
          ActsScalar thickness = deltaR;
          ActsScalar medR = 0.5 * (minR + maxR);
          ActsScalar halfZ = daughter->GetDz() * m_cfg.unitScalar;
          double halfPhi = M_PI;
          double avgPhi = 0.;
          std::cout << "cast good: base top strip has rmin = " << minR
                    << ", rmax = " << maxR << " and dphi "
                    << std::abs(parameters[2] - parameters[3])
                    << " and thickness " << thickness << " and halfZ " << halfZ
                    << std::endl;

          // make a cylinderBounds
          if (halfZ > deltaR) {
            auto bounds = std::make_shared<Acts::CylinderBounds>(
                medR, halfZ, halfPhi, avgPhi);

            // Create a new detector element with shifted medR
            // This assumes the base tube has Identity transform
            auto tgDetectorElement =
                std::make_shared<Acts::TGeoDetectorElement>(
                    tgIdentifier, *node, Transform3::Identity(), bounds,
                    thickness);

            tgDetectorElements.push_back(tgDetectorElement);
            break;
          }
        }
      }
    }
  }

  std::cout << "Found " << tgDetectorElements.size()
            << " topStrip detElement on this layer" << std::endl;

  if (not tgDetectorElements.empty()) {
    return tgDetectorElements;
  }

  return {tgde};
}
