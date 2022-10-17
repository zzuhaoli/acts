// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/TGeo/TGeoLayerBuilder.hpp"

#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/ProtoLayer.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/TGeo/ITGeoDetectorElementSplitter.hpp"
#include "Acts/Plugins/TGeo/TGeoDetectorElement.hpp"
#include "Acts/Plugins/TGeo/TGeoParser.hpp"
#include "Acts/Plugins/TGeo/TGeoPrimitivesHelper.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"

#include <stdio.h>

#include "TGeoManager.h"
#include "TGeoMatrix.h"

Acts::TGeoLayerBuilder::TGeoLayerBuilder(
    const Acts::TGeoLayerBuilder::Config& config,
    std::unique_ptr<const Logger> logger)
    : m_cfg(), m_logger(std::move(logger)) {
  setConfiguration(config);
}

Acts::TGeoLayerBuilder::~TGeoLayerBuilder() = default;

void Acts::TGeoLayerBuilder::setConfiguration(
    const Acts::TGeoLayerBuilder::Config& config) {
  m_cfg = config;
}

void Acts::TGeoLayerBuilder::setLogger(
    std::unique_ptr<const Logger> newLogger) {
  m_logger = std::move(newLogger);
}

const Acts::LayerVector Acts::TGeoLayerBuilder::negativeLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector nVector;
  mutableThis->buildLayers(gctx, nVector, -1);
  return nVector;
}

const Acts::LayerVector Acts::TGeoLayerBuilder::centralLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector cVector;
  mutableThis->buildLayers(gctx, cVector, 0);
  return cVector;
}

const Acts::LayerVector Acts::TGeoLayerBuilder::positiveLayers(
    const GeometryContext& gctx) const {
  // @todo Remove this hack once the m_elementStore mess is sorted out
  auto mutableThis = const_cast<TGeoLayerBuilder*>(this);
  LayerVector pVector;
  mutableThis->buildLayers(gctx, pVector, 1);
  return pVector;
}

// build layers for a volume?
void Acts::TGeoLayerBuilder::buildLayers(const GeometryContext& gctx,
                                         LayerVector& layers, int type) {
  // Bail out if you have no gGeoManager
  if (gGeoManager == nullptr) {
    ACTS_WARNING("No gGeoManager found - bailing out.");
    return;
  }

  std::cout << "TGeoLayerBuilder::buildLayers called " << std::endl;

  using LayerSurfaceVector = std::vector<std::shared_ptr<const Surface>>;
  LayerSurfaceVector layerSurfaces;

  std::vector<LayerConfig> layerConfigs = m_cfg.layerConfigurations[type + 1];
  std::string layerType = m_layerTypes[type + 1];
  std::cout << "layerType " << layerType << std::endl;

  // Appropriate screen output
  std::string addonOutput = m_cfg.layerSplitToleranceR[type + 1] > 0.
                                ? std::string(", splitting in r")
                                : std::string("");
  addonOutput += m_cfg.layerSplitToleranceZ[type + 1] > 0.
                     ? std::string(", splitting in z")
                     : std::string("");
  addonOutput += std::string(".");

  // Screen output of the configuration
  ACTS_DEBUG(layerType << " layers : found " << layerConfigs.size()
                       << " configuration(s)" + addonOutput);

  std::cout << layerType << " layers : found " << layerConfigs.size()
            << " configuration(s)" + addonOutput << std::endl;

  // Helper function to fill the layer
  auto fillLayer = [&](const LayerSurfaceVector lSurfaces,
                       const LayerConfig& lCfg,
                       unsigned int pl_id = 0) -> void {
    int nb0 = 0, nt0 = 0;
    bool is_autobinning = ((lCfg.binning0.size() == 1) and
                           (std::get<int>(lCfg.binning0.at(0)) <= 0));
    if (is_autobinning) {
      std::cout << "is_autobinning" << std::endl;
    }

    if (!is_autobinning and std::get<int>(lCfg.binning0.at(pl_id)) <= 0) {
      throw std::invalid_argument(
          "Incorrect binning configuration found for loc0 protolayer #" +
          std::to_string(pl_id) +
          ". Layer is autobinned: No mixed binning (manual and auto) for loc0 "
          "possible between layers in a single subvolume. Quitting");
    }
    if (is_autobinning) {
      // Set binning by hand if nb0 > 0 and nb1 > 0
      nb0 = std::get<int>(lCfg.binning0.at(0));
      // Read the binning type
      nt0 = std::get<BinningType>(lCfg.binning0.at(0));
    } else if (pl_id < lCfg.binning0.size()) {
      // Set binning by hand if nb0 > 0 and nb1 > 0
      nb0 = std::get<int>(lCfg.binning0.at(pl_id));
    }

    int nb1 = 0, nt1 = 0;
    is_autobinning = (lCfg.binning1.size() == 1) and
                     (std::get<int>(lCfg.binning1.at(0)) <= 0);
    if (!is_autobinning and std::get<int>(lCfg.binning1.at(pl_id)) <= 0) {
      throw std::invalid_argument(
          "Incorrect binning configuration found for loc1 protolayer #" +
          std::to_string(pl_id) +
          ". Layer is autobinned: No mixed binning (manual and auto) for loc1 "
          "possible between layers in a single subvolume. Quitting");
    }
    if (is_autobinning) {
      // Set binning by hand if nb0 > 0 and nb1 > 0
      nb1 = std::get<int>(lCfg.binning1.at(0));
      // For a binning type
      nt1 = std::get<BinningType>(lCfg.binning1.at(0));
    } else if (pl_id < lCfg.binning1.size()) {
      // Set binning by hand if nb0 > 0 and nb1 > 0
      nb1 = std::get<int>(lCfg.binning1.at(pl_id));
    }

    std::cout << "type = " << type << std::endl;

    if (type == 0) {
      ProtoLayer pl(gctx, lSurfaces);
      ACTS_DEBUG("- creating CylinderLayer with "
                 << lSurfaces.size() << " surfaces at r = " << pl.medium(binR));

      std::cout << " creating CylinderLayer with " << lSurfaces.size()
                << " surfaces at r = " << pl.medium(binR) << std::endl;
      std::cout << "lCfg.envelope.first = " << lCfg.envelope.first
                << " and lCfg.envelope.second " << lCfg.envelope.second
                << std::endl;

      pl.envelope[Acts::binR] = {lCfg.envelope.first, lCfg.envelope.second};
      pl.envelope[Acts::binZ] = {lCfg.envelope.second, lCfg.envelope.second};
      std::cout << "nb0 = " << nb0 << ", nb1 = " << nb1 << std::endl;
      std::cout << "nt0 = " << nb0 << ", nt1 = " << nb1 << std::endl;
      if (nb0 >= 0 and nb1 >= 0) {
        layers.push_back(
            m_cfg.layerCreator->cylinderLayer(gctx, lSurfaces, nb0, nb1, pl));
      } else {
        layers.push_back(
            m_cfg.layerCreator->cylinderLayer(gctx, lSurfaces, nt0, nt1, pl));
      }
    } else {
      ProtoLayer pl(gctx, lSurfaces);
      ACTS_DEBUG("- creating DiscLayer with "
                 << lSurfaces.size() << " surfaces at z = " << pl.medium(binZ));

      pl.envelope[Acts::binR] = {lCfg.envelope.first, lCfg.envelope.second};
      pl.envelope[Acts::binZ] = {lCfg.envelope.second, lCfg.envelope.second};
      if (nb0 >= 0 and nb1 >= 0) {
        layers.push_back(
            m_cfg.layerCreator->discLayer(gctx, lSurfaces, nb0, nb1, pl));
      } else {
        layers.push_back(
            m_cfg.layerCreator->discLayer(gctx, lSurfaces, nt0, nt1, pl));
      }
    }
  };

  for (auto layerCfg : layerConfigs) {
    ACTS_DEBUG("- layer configuration found for layer " << layerCfg.volumeName
                                                        << " with sensors ");
    for (auto& sensor : layerCfg.sensorNames) {
      ACTS_DEBUG("  - sensor: " << sensor);
    }
    if (not layerCfg.parseRanges.empty()) {
      for (const auto& pRange : layerCfg.parseRanges) {
        ACTS_DEBUG("- layer parsing restricted in "
                   << binningValueNames()[pRange.first] << " to ["
                   << pRange.second.first << "/" << pRange.second.second
                   << "].");
      }
    }
    if (not layerCfg.splitConfigs.empty()) {
      for (const auto& sConfig : layerCfg.splitConfigs) {
        ACTS_DEBUG("- layer splitting attempt in "
                   << binningValueNames()[sConfig.first] << " with tolerance "
                   << sConfig.second << ".");
      }
    }

    // Either pick the configured volume or take the top level volume
    TGeoVolume* tVolume =
        gGeoManager->FindVolumeFast(layerCfg.volumeName.c_str());
    if (tVolume == nullptr) {
      tVolume = gGeoManager->GetTopVolume();
      ACTS_DEBUG("- search volume is TGeo top volume");
    } else {
      ACTS_DEBUG("- setting search volume to " << tVolume->GetName());
    }

    if (tVolume != nullptr) {
      TGeoParser::Options tgpOptions;
      tgpOptions.volumeNames = {layerCfg.volumeName};
      tgpOptions.targetNames = layerCfg.sensorNames;
      tgpOptions.parseRanges = layerCfg.parseRanges;
      tgpOptions.unit = m_cfg.unit;
      TGeoParser::State tgpState;
      tgpState.volume = tVolume;

      ACTS_DEBUG("- applying  " << layerCfg.parseRanges.size()
                                << " search restrictions.");
      for (const auto& prange : layerCfg.parseRanges) {
        ACTS_VERBOSE(" - range " << binningValueNames()[prange.first]
                                 << " within [ " << prange.second.first << ", "
                                 << prange.second.second << "]");
      }

      std::cout << "Start select nodes for volume " << tVolume->GetName()
                << std::endl;
      TGeoParser::select(tgpState, tgpOptions);

      ACTS_DEBUG("- number of selected nodes found : "
                 << tgpState.selectedNodes.size());
      std::cout << "- number of selected nodes found : "
                << tgpState.selectedNodes.size() << std::endl;

      // This assume each selectedNode is one potential layer
      std::vector<ProtoLayer> protoLayers;

      // Loop over the nodes to split for more nodes using splitter
      for (auto& snode : tgpState.selectedNodes) {
        auto identifier =
            m_cfg.identifierProvider != nullptr
                ? m_cfg.identifierProvider->identify(gctx, *snode.node)
                : Identifier();

        // create the elementary detElement
        auto tgElement =
            m_cfg.elementFactory(identifier, *snode.node, *snode.transform,
                                 layerCfg.localAxes, m_cfg.unit, nullptr);

        const TGeoNode& tgNode = tgElement->tgeoNode();
        ActsScalar tgThickness = tgElement->thickness();  // turn cm into mm
        std::cout << "original detElement based on node " << tgNode.GetName()
                  << " has thickness " << tgThickness << std::endl;

        if (m_cfg.detectorElementSplitter == nullptr) {
          std::cout << " Empty splitter" << std::endl;
        }
        std::vector<std::shared_ptr<const Acts::TGeoDetectorElement>>
            tgElements =
                (m_cfg.detectorElementSplitter == nullptr)
                    ? std::vector<std::shared_ptr<
                          const Acts::TGeoDetectorElement>>{tgElement}
                    : m_cfg.detectorElementSplitter->split(gctx, tgElement);

        std::cout << tgElements.size() << " detElements are created for node "
                  << snode.node->GetName() << std::endl;

        LayerSurfaceVector thisLayerSurfaces;
        for (auto tge : tgElements) {
          m_elementStore.push_back(tge);
          layerSurfaces.push_back(tge->surface().getSharedPtr());
          thisLayerSurfaces.push_back(tge->surface().getSharedPtr());
        }
        protoLayers.push_back(ProtoLayer(gctx, thisLayerSurfaces));
      }

      ACTS_DEBUG("- created TGeoDetectorElements : " << layerSurfaces.size());
      std::cout << "- created a total of TGeoDetectorElements : "
                << layerSurfaces.size() << " for volume " << layerCfg.volumeName
                << std::endl;

      if (m_cfg.protoLayerHelper != nullptr and
          not layerCfg.splitConfigs.empty()) {
        //@TODO
        // auto protoLayers = m_cfg.protoLayerHelper->protoLayers(
        //     gctx, unpack_shared_vector(layerSurfaces),
        //     layerCfg.splitConfigs);
        ACTS_DEBUG("- splitting into " << protoLayers.size() << " layers.");
        std::cout << "The TGeoDetectorElements are then splitted into "
                  << protoLayers.size() << " layers." << std::endl;

        // Number of options mismatch and has not been configured for
        // auto-binning
        const bool is_loc0_n_config =
            layerCfg.binning0.size() == protoLayers.size();
        const bool is_loc0_autobinning =
            (layerCfg.binning0.size() == 1) and
            (std::get<int>(layerCfg.binning0.at(0)) <= 0);
        const bool is_loc1_n_config =
            layerCfg.binning1.size() == protoLayers.size();
        const bool is_loc1_autobinning =
            (layerCfg.binning1.size() == 1) and
            (std::get<int>(layerCfg.binning1.at(0)) <= 0);
        if ((!is_loc0_n_config and !is_loc0_autobinning) or
            (!is_loc1_n_config and !is_loc1_autobinning)) {
          throw std::invalid_argument(
              "Incorrect binning configuration found: Number of configurations "
              "does not match number of protolayers in subvolume " +
              layerCfg.volumeName + ". Quitting.");
        }
        unsigned int layer_id = 0;
        for (auto& pLayer : protoLayers) {
          layerSurfaces.clear();

          for (const auto& lsurface : pLayer.surfaces()) {
            layerSurfaces.push_back(lsurface->getSharedPtr());
          }
          fillLayer(layerSurfaces, layerCfg, layer_id);
          layer_id++;
        }
      } else {
        fillLayer(layerSurfaces, layerCfg);
        std::cout << "The TGeoDetectorElements are put into one layer."
                  << std::endl;
      }
    }
  }  // end loop of negative subvol, middle subvol, positive subvol

  std::vector<std::pair<const std::string, Acts::BinningOption>> binning = {
      {"binPhi", Acts::closed}, {"binZ", Acts::open}};
  for (auto& layer : layers) {
    // if(layer->trackingVolume()!=nullptr){
    //   auto volName = layer->trackingVolume()->volumeName();
    //   std::cout<<"volName = " << volName << std::endl;
    // }
    // auto geoID = layer->geometryId();
    // std::cout<<"layer has geoID " << geoID << std::endl;
    size_t nSensitiveSurfaces = layer->surfaceArray()->surfaces().size();
    std::cout << "layer has " << nSensitiveSurfaces
              << " sensitive surfaces with layer thickness "
              << layer->thickness() << std::endl;
    for(const auto& ss : layer->surfaceArray()->surfaces()){
      if(ss->surfaceMaterial()){
        std::cout<<"sensitive surface: "<< ss->geometryId() <<" has material " << std::endl;
      } 
    } 
    std::vector<std::pair<std::string, int>> surfaceBins = {{"binPhi", 1},
                                                            {"binZ", 20}};
    // use the same binning for representing surface and approach surfaces
    Layer* castLayer = const_cast<Layer*>(&(*layer));
    addCylinderLayerProtoMaterial(*castLayer, binning, surfaceBins, surfaceBins,
                                  surfaceBins);

    if (layer->surfaceRepresentation().surfaceMaterial() != nullptr) {
      std::cout << "layer representative surface has material " << std::endl;
    }
    auto aDescriptor = layer->approachDescriptor();
    if (aDescriptor != nullptr) {
      for (const auto& surface : aDescriptor->containedSurfaces()) {
        const CylinderBounds* cbounds =
            dynamic_cast<const CylinderBounds*>(&surface->bounds());
        std::cout << "approach cylinder surface at r =  "
                  << cbounds->values()[0] << std::endl;
        if (surface->surfaceMaterial() != nullptr) {
          std::cout << "approach surface has material " << std::endl;
        }
      }
    }
  }

  return;
}

std::shared_ptr<Acts::TGeoDetectorElement>
Acts::TGeoLayerBuilder::defaultElementFactory(
    const Identifier& identifier, const TGeoNode& tGeoNode,
    const TGeoMatrix& tGeoMatrix, const std::string& axes, double scalor,
    std::shared_ptr<const Acts::ISurfaceMaterial> material) {
  return std::make_shared<TGeoDetectorElement>(identifier, tGeoNode, tGeoMatrix,
                                               axes, scalor, material);
}

std::shared_ptr<Acts::ProtoSurfaceMaterial>
Acts::TGeoLayerBuilder::createProtoMaterial(
    const std::vector<std::pair<const std::string, Acts::BinningOption>>&
        binning,
    const std::vector<std::pair<std::string, int>>& surfaceBins) const {
  // Create the bin utility
  Acts::BinUtility bu;
  // Loop over the bins
  for (auto& bin : binning) {
    // finding the iterator position to determine the binning value
    auto bit = std::find(Acts::binningValueNames().begin(),
                         Acts::binningValueNames().end(), bin.first);
    size_t indx = std::distance(Acts::binningValueNames().begin(), bit);
    Acts::BinningValue bval = Acts::BinningValue(indx);
    Acts::BinningOption bopt = bin.second;
    double min = 0.;
    double max = 0.;
    if (bopt == Acts::closed) {
      min = -M_PI;
      max = M_PI;
    }
    auto binIt = std::find_if(surfaceBins.begin(), surfaceBins.end(),
                              [=](const std::pair<std::string, int>& value) {
                                return (value.first == bin.first);
                              });
    if (binIt != surfaceBins.end()) {
      int bins = (*binIt).second;
      std::cout << "bins =  " << bins << std::endl;
      if (bins >= 1) {
        bu += Acts::BinUtility(bins, min, max, bopt, bval);
      }
    }
  }

  return std::make_shared<Acts::ProtoSurfaceMaterial>(bu);
}

void Acts::TGeoLayerBuilder::addCylinderLayerProtoMaterial(
    Layer& layer,
    const std::vector<std::pair<const std::string, Acts::BinningOption>>&
        binning,
    const std::vector<std::pair<std::string, int>>& representSurfaceBins,
    const std::vector<std::pair<std::string, int>>& innerApproachSurfaceBins,
    const std::vector<std::pair<std::string, int>>& outerApproachSurfaceBins)
    const {
  //Surface* representationSurface =
  //    const_cast<Surface*>(&(layer.surfaceRepresentation()));
  //representationSurface->assignSurfaceMaterial(
  //    createProtoMaterial(binning, representSurfaceBins));

  auto aDescriptor = layer.approachDescriptor();
  if (aDescriptor != nullptr and aDescriptor->containedSurfaces().size() >= 2) {
    // Add the inner and outer approach surface
    const std::vector<const Surface*>& aSurfaces =
        aDescriptor->containedSurfaces();
    
    Surface* innerApproachSurface = const_cast<Surface*>(aSurfaces[0]);
    innerApproachSurface->assignSurfaceMaterial(
        createProtoMaterial(binning, innerApproachSurfaceBins));
    
    //Surface* outerApproachSurface = const_cast<Surface*>(aSurfaces[1]);
    //outerApproachSurface->assignSurfaceMaterial(
    //    createProtoMaterial(binning, outerApproachSurfaceBins));
  }
}
