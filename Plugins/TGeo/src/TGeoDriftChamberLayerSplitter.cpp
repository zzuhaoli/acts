// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoDriftChamberLayerSplitter.hpp"

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

Acts::TGeoDriftChamberLayerSplitter::TGeoDriftChamberLayerSplitter(
    const TGeoDriftChamberLayerSplitter::Config& cfg,
    std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {}

std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
Acts::TGeoDriftChamberLayerSplitter::split(
    const GeometryContext& /*gctx*/,
    std::shared_ptr<const Acts::TGeoDetectorElement> tgde) const {
  // Thickness
  auto tgIdentifier = tgde->identifier();
  const TGeoNode& tgNode = tgde->tgeoNode();
  ACTS_DEBUG("TGeoDriftChamberLayerSplitter splitting detElement "
             << tgNode.GetName());

  std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
      tgDetectorElements = {};

  // get the number of sense wires in this layer
  int nSenseWires = 0;
  // The cells on this layer
  auto daugthers = tgNode.GetVolume()->GetNodes();
  TIter iObj(daugthers);
  while (TObject* obj = iObj()) {
    TGeoNode* node = dynamic_cast<TGeoNode*>(obj);
    if (node != nullptr) {
      // This is a cell. Check its size
      TGeoTube* cell = dynamic_cast<TGeoTube*>(node->GetVolume()->GetShape());
      if (cell == nullptr) {
        ACTS_WARNING(
            "cast bad (drift cell is always a tube. This is not supposed to "
            "happen)");
      } else {
        ActsScalar parameters[5];
        cell->GetBoundingCylinder(parameters);
        ActsScalar minR = cell->GetRmin() * m_cfg.unitScalar;  // in cm
        ActsScalar maxR = cell->GetRmax() * m_cfg.unitScalar;  // in cm
        ActsScalar deltaPhi =
            std::abs(parameters[2] - parameters[3]) * M_PI / 180;
        // calcuate the approximate lineBounds r
        ActsScalar r = std::hypot(maxR * std::cos(deltaPhi) - (minR + maxR) / 2,
                                  maxR * std::sin(deltaPhi));
        ActsScalar thickness = maxR - minR;
        ACTS_DEBUG("cast good: drift cell has minR = "
                   << minR << ", maxR = " << maxR << " and deltaPhi "
                   << deltaPhi << " (deg) and thickness " << thickness);

        // The relative position of cell w.r.t. layer
        const TGeoMatrix* cmatrix = node->GetMatrix();

        auto grandDaugthers = node->GetVolume()->GetNodes();
        TIter jObj(grandDaugthers);
        while (TObject* dobj = jObj()) {
          TGeoNode* dNode = dynamic_cast<TGeoNode*>(dobj);
          if (dNode != nullptr) {
            std::string dNodeName = dNode->GetName();
            if (dNodeName.find("sense_wire") != std::string::npos) {
              ACTS_DEBUG("Found sense wire" << dNodeName
                                            << " in the drift cell");
              // create a Line surface
              nSenseWires++;
              TGeoTube* wire =
                  dynamic_cast<TGeoTube*>(dNode->GetVolume()->GetShape());
              if (wire == nullptr) {
                ACTS_WARNING("Failed to cast to TGeoTube");
              } else {
                ActsScalar halfZ = wire->GetDz() * m_cfg.unitScalar;
                ACTS_DEBUG("half length of the cell: " << halfZ);

                const TGeoMatrix* wmatrix = dNode->GetMatrix();
                // Is this correct?
                TGeoHMatrix transform =
                    TGeoCombiTrans(*cmatrix) * TGeoCombiTrans(*wmatrix);

                // Get the placement and orientation in respect to its mother
                const Double_t* rotation = transform.GetRotationMatrix();
                const Double_t* translation = transform.GetTranslation();
                // Create a eigen transform
                Vector3 t(translation[0] * m_cfg.unitScalar,
                          translation[1] * m_cfg.unitScalar,
                          translation[2] * m_cfg.unitScalar);
                Vector3 cx(rotation[0], rotation[3], rotation[6]);
                Vector3 cy(rotation[1], rotation[4], rotation[7]);
                Vector3 cz(rotation[2], rotation[5], rotation[8]);
                auto etrf = TGeoPrimitivesHelper::makeTransform(cx, cy, cz, t);

                // make a lineBounds
                auto tgWire = std::make_shared<Acts::LineBounds>(r*1.5, halfZ);

                // Create a new detector element per split
                auto tgDetectorElement =
                    std::make_shared<Acts::TGeoDetectorElement>(
                        tgIdentifier, *node, etrf, tgWire, thickness);

                tgDetectorElements.push_back(tgDetectorElement);
              }
            }
          }
        }
      }
    }
  }

  ACTS_DEBUG("Found " << nSenseWires << " sense wires on this layer");

  if (not tgDetectorElements.empty()) {
    return tgDetectorElements;
  }

  return {tgde};
}
