// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <array>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TFile.h>

#include "CommonUtils.h"
#include "TreeReader.h"
#include "SetStyleATLAS.hpp"

/// This script/function reads all the reconstructed tracks from e.g.
/// 'tracksummary_ckf.root' and the (possibly selected) truth particles from
/// e.g. 'track_finder_particles.root' (which contains the info of 'nHits'), and
/// defines the efficiency, fake rate and duplicaiton rate. It aims to make
/// custom definition and tuning of the reconstruction performance easier.
/// Multiple files for the reconstructed tracks are allowed.
/// 
/// NB: It's very likely that fiducal cuts are already imposed on the truth
/// particles. Please check the selection criteria in the truth fitting example
/// which writes out the 'track_finder_particles.root'. For instance, if the
/// truth particles are already required to have pT > 1 GeV, it does not make
/// sense to have ptMin = 0.5 GeV here.
///
void defineReconstructionPerformance(
    const std::vector<std::string>& inputSimParticleFileNames =
        {
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_90deg/performance_seeding_trees.root",
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_30deg/performance_seeding_trees.root",
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_90deg/performance_seeding_trees.root",
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_30deg/performance_seeding_trees.root",
	"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/2mu2pi_test/performance_seeding_trees.root",
	},
    const std::vector<std::string>& inputTrackSummaryFileNames =
        {
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_90deg/tracksummary_ckf.root",
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_30deg/tracksummary_ckf.root",
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_90deg/tracksummary_ckf.root",
	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_30deg/tracksummary_ckf.root",
	"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/2mu2pi_test/tracksummary_ckf.root",
	},
    const std::vector<std::string>& trackSummaryFileLegends =
        {
	//"muon (#theta=90 deg)",
	//"muon (#theta=30 deg)",
	//"pion (#theta=90 deg)",
	//"pion (#theta=30 deg)",
	"pion in #psi(2S)->#pi^{+}#pi^{-}J/#psi(#rightarrow#mu^{+}#mu^{-})",
	},
	std::vector<int> colors={
          872,
          866,
          854,
          896,
	},
	std::vector<int> markers ={24, 26, 20, 22},
    const std::string& simParticleTreeName = "track_finder_particles",
    const std::string& trackSummaryTreeName = "tracksummary",
    unsigned int nHitsMin = 5, unsigned int nMeasurementsMin = 5,
    double ptMin = 0.05, double absEtaMax=2, double truthMatchProbMin = 0.5) {
  gStyle->SetOptFit(0011);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.05, "xy");
  gStyle->SetLabelSize(0.05, "xy");
  gStyle->SetTitleOffset(1., "x");
  gStyle->SetTitleOffset(1.5, "y");
  gStyle->SetNdivisions(505, "y");

  //std::vector<double> ptRanges = {40, 0.05, 0.45};
  std::vector<double> ptRanges = {20, 0.05, 0.45};
  //std::vector<double> ptRanges = {19, 0.05, 1.95};
  //std::vector<double> pts={0.05, 0.15, 0.3, 0.45};

  // Check the inputs are valid
  if (inputSimParticleFileNames.size()!=inputTrackSummaryFileNames.size() or inputTrackSummaryFileNames.size() != trackSummaryFileLegends.size()) {
    throw std::invalid_argument(
        "Please specify the legends you want to show for all the track files");
  }

  // Open the files for reading
  //TFile* particleFile = TFile::Open(inputSimParticleFileName.c_str(), "read");
  
  // The number of track files to read
  unsigned int nTrackFiles = inputTrackSummaryFileNames.size();
  std::vector<TFile*> particleFiles;
  std::vector<TFile*> trackFiles;
  particleFiles.reserve(nTrackFiles);
  trackFiles.reserve(nTrackFiles);
  for (const auto& fileName : inputSimParticleFileNames) {
    particleFiles.push_back(TFile::Open(fileName.c_str(), "read"));
  }
  for (const auto& fileName : inputTrackSummaryFileNames) {
    trackFiles.push_back(TFile::Open(fileName.c_str(), "read"));
  }

  // Define variables for tree reading (turn on the events sorting since we have more than one root files to read)
  //ParticleReader pReader(
  //    (TTree*)particleFile->Get(simParticleTreeName.c_str()), true);
  std::vector<ParticleReader> pReaders;
  std::vector<TrackSummaryReader> tReaders;
  pReaders.reserve(nTrackFiles);
  tReaders.reserve(nTrackFiles);
  for (const auto& particleFile : particleFiles) {
    pReaders.emplace_back((TTree*)particleFile->Get(simParticleTreeName.c_str()), true);
  }
  for (const auto& trackFile : trackFiles) {
    tReaders.emplace_back((TTree*)trackFile->Get(trackSummaryTreeName.c_str()), true);
  }

  std::vector<size_t> nEvents;
  nEvents.reserve(nTrackFiles);
  for (const auto& tReader : tReaders) {
    size_t entries = tReader.tree->GetEntries();
    nEvents.push_back(entries);
  }

  // Define the efficiency plots
  std::vector<TEfficiency*> trackEff_vs_eta;
  std::vector<TEfficiency*> fakeRate_vs_eta;
  std::vector<TEfficiency*> duplicateRate_vs_eta;
  std::vector<TEfficiency*> trackEff_vs_pt;
  std::vector<TEfficiency*> fakeRate_vs_pt;
  std::vector<TEfficiency*> duplicateRate_vs_pt;

  for (int i = 0; i < nTrackFiles; ++i) {
    trackEff_vs_eta.push_back(new TEfficiency(
        Form("trackeff_vs_eta_%i", i), ";Truth #theta;Efficiency", 40, -4, 4));
    fakeRate_vs_eta.push_back(new TEfficiency(
        Form("fakerate_vs_eta_%i", i), ";#theta;fake rate", 40, -4, 4));
    duplicateRate_vs_eta.push_back(new TEfficiency(
        Form("duplicaterate_vs_eta_%i", i), ";#theta;Duplicate rate", 40, -4, 4));
    trackEff_vs_pt.push_back(new TEfficiency(
        Form("trackeff_vs_pt_%i", i), ";Truth pt [GeV];Efficiency", ptRanges[0], ptRanges[1], ptRanges[2]));
    fakeRate_vs_pt.push_back(new TEfficiency(
        Form("fakerate_vs_pt_%i", i), ";pt [GeV];fake rate", ptRanges[0], ptRanges[1], ptRanges[2]));
    duplicateRate_vs_pt.push_back(new TEfficiency(
        Form("duplicaterate_vs_pt_%i", i), ";pt [GeV];Duplicate rate", ptRanges[0], ptRanges[1], ptRanges[2]));
  }

  // Set styles
  for (int i = 0; i < nTrackFiles; ++i) {
    auto color = colors[i];
    auto marker  = markers[i];
    setEffStyle(trackEff_vs_eta[i], color, marker);
    setEffStyle(fakeRate_vs_eta[i], color, marker);
    setEffStyle(duplicateRate_vs_eta[i], color, marker);
    setEffStyle(trackEff_vs_pt[i], color, marker);
    setEffStyle(fakeRate_vs_pt[i], color, marker);
    setEffStyle(duplicateRate_vs_pt[i], color, marker);
  }

  // Loop over the track files
  for (unsigned int ifile = 0; ifile < nTrackFiles; ++ifile) {
    std::cout << "Processing track file: " << inputTrackSummaryFileNames[ifile]
              << std::endl;
  
    // The particles in each event
    std::map<unsigned int, std::vector<ParticleInfo>> particles;


    // The container from track-particle matching info (Flushed per event)
    std::map<uint64_t, std::vector<RecoTrackInfo>> matchedParticles;

    // Loop over the events to fill plots
    for (size_t i = 0; i < nEvents[ifile]; ++i) {
      if (i % 10 == 0) {
        std::cout << "Processed events: " << i << std::endl;
      }

      // Get the tracks
      tReaders[ifile].getEntry(i);

      // Get the particles (do nothing if the particles for this event already
      // read)
      auto it = particles.find(i);
      if (it == particles.end()) {
        particles.emplace(i, pReaders[ifile].getParticles(i));
      }

      // Loop over the tracks
      // The fake rate is defined as the ratio of selected truth-matched tracks
      // over all selected tracks
      for (size_t j = 0; j < tReaders[ifile].nStates->size(); ++j) {
        bool hasFittedParameters = tReaders[ifile].hasFittedParams->at(j);
        auto nMeasurements = tReaders[ifile].nMeasurements->at(j);
        auto nOutliers = tReaders[ifile].nOutliers->at(j);
        auto nHoles = tReaders[ifile].nHoles->at(j);
        auto theta = tReaders[ifile].eTHETA_fit->at(j);
        auto qop = tReaders[ifile].eQOP_fit->at(j);
        auto pt = std::abs(1 / qop) * std::sin(theta);
        auto eta = std::atanh(std::cos(theta));
        auto nMajorityHits = tReaders[ifile].nMajorityHits->at(j);
        auto majorityParticleId = tReaders[ifile].majorityParticleId->at(j);

        // Select the track, e.g. you might also want to add cuts on the
        // nOutliers, nHoles
        if ((!hasFittedParameters) or nMeasurements < nMeasurementsMin or
            pt < ptMin) {
          continue;
        }

        // Fill the fake rate plots
        if (nMajorityHits * 1. / nMeasurements >= truthMatchProbMin) {
          matchedParticles[majorityParticleId].push_back(
              {eta, pt, nMajorityHits, nMeasurements});
          fakeRate_vs_eta[ifile]->Fill(false, theta);
          fakeRate_vs_pt[ifile]->Fill(false, pt);
        } else {
          fakeRate_vs_eta[ifile]->Fill(true, theta);
          fakeRate_vs_pt[ifile]->Fill(true, pt);
        }
      }  // end of all tracks

      // Loop over all selected and truth-matched tracks
      // The duplicate rate is defined as the ratio of duplicate tracks among
      // all the selected truth-matched tracks (only one track is 'real'; others
      // are 'duplicated')
      for (auto& [id, matchedTracks] : matchedParticles) {
        // Sort all tracks matched to this particle according to majority prob
        // and track quality
        std::sort(matchedTracks.begin(), matchedTracks.end(),
                  [](const RecoTrackInfo& lhs, const RecoTrackInfo& rhs) {
                    if (lhs.nMajorityHits > rhs.nMajorityHits) {
                      return true;
                    }
                    if (lhs.nMajorityHits < rhs.nMajorityHits) {
                      return false;
                    }
                    if (lhs.nMeasurements > rhs.nMeasurements) {
                      return true;
                    }
                    return false;
                  });
        // Fill the duplication rate plots
        for (size_t k = 0; k < matchedTracks.size(); ++k) {
          auto eta = matchedTracks[k].eta;
          auto pt = matchedTracks[k].pt;
          double theta = std::atan(std::exp(-eta))*2;
          if (k == 0) {
            duplicateRate_vs_eta[ifile]->Fill(false, theta);
            duplicateRate_vs_pt[ifile]->Fill(false, pt);
          } else {
            duplicateRate_vs_eta[ifile]->Fill(true, theta);
            duplicateRate_vs_pt[ifile]->Fill(true, pt);
          }
        }
      }  // end of all selected truth-matched tracks

      // Loop over all truth particles in this event
      // The effiency is define as the ratio of selected particles that have
      // been matched with reco
      if(particles[i].size()==0){
         std::cout<<"Skip event " << i << " as there is no truth particle in this event" << std::endl; 
      } 
      for (const auto& particle : particles[i]) {
        auto nHits = particle.nHits;
        auto eta = particle.eta;
        auto pt = particle.pt;
        if (abs(particle.particlePdg) != 211){
          continue;	
	}	
        if (abs(particle.eta) > absEtaMax){
          continue;	
	}	
	if (nHits < nHitsMin or pt < ptMin ) {
          continue;
        }
        uint64_t id = particle.particleId;

        double theta = std::atan(std::exp(-eta))*2;
        // Fill the efficiency plots
        auto ip = matchedParticles.find(id);
        if (ip != matchedParticles.end()) {
          trackEff_vs_eta[ifile]->Fill(true, theta);
          trackEff_vs_pt[ifile]->Fill(true, pt);
        } else {
          trackEff_vs_eta[ifile]->Fill(false, theta);
          trackEff_vs_pt[ifile]->Fill(false, pt);
        }
      }  // end of all particles

      matchedParticles.clear();
    }  // end of all events
  }    // end of all track files

  std::cout << "All good. Now plotting..." << std::endl;

  // The legends
  std::vector<TLegend*> legs;
  for (int i = 0; i < 6; ++i) {
    TLegend* legend = new TLegend(0.4, 0.8, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legs.push_back(legend);
  }
  // Add entry for the legends
  for (int i = 0; i < nTrackFiles; ++i) {
    legs[0]->AddEntry(trackEff_vs_eta[i],
                      Form("%s", trackSummaryFileLegends[i].c_str()), "lp");
    legs[1]->AddEntry(fakeRate_vs_eta[i],
                      Form("%s", trackSummaryFileLegends[i].c_str()), "lp");
    legs[2]->AddEntry(duplicateRate_vs_eta[i],
                      Form("%s", trackSummaryFileLegends[i].c_str()), "lp");
    legs[3]->AddEntry(trackEff_vs_pt[i],
                      Form("%s", trackSummaryFileLegends[i].c_str()), "lp");
    legs[4]->AddEntry(fakeRate_vs_pt[i],
                      Form("%s", trackSummaryFileLegends[i].c_str()), "lp");
    legs[5]->AddEntry(duplicateRate_vs_pt[i],
                      Form("%s", trackSummaryFileLegends[i].c_str()), "lp");
  }

  // Now draw the plots
  TCanvas* c1 = new TCanvas("recoPerf", " ", 1500, 800);
  c1->Divide(3, 2);

  float scaleRangeMax = 1.1;
  for (int i = 0; i < nTrackFiles; ++i) {
    std::string mode = (i == 0) ? "" : "same";
    c1->cd(1);
    trackEff_vs_eta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[0]->Draw();
    }
    adaptEffRange(trackEff_vs_eta[i], 1, scaleRangeMax);

    c1->cd(2);
    fakeRate_vs_eta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[1]->Draw();
    }
    adaptEffRange(fakeRate_vs_eta[i], 1, scaleRangeMax);

    c1->cd(3);
    duplicateRate_vs_eta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[2]->Draw();
    }
    adaptEffRange(duplicateRate_vs_eta[i], 1, scaleRangeMax);

    c1->cd(4);
    trackEff_vs_pt[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[3]->Draw();
    }
    adaptEffRange(trackEff_vs_pt[i], 1, scaleRangeMax);

    c1->cd(5);
    fakeRate_vs_pt[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[4]->Draw();
    }
    adaptEffRange(fakeRate_vs_pt[i], 1, scaleRangeMax);

    c1->cd(6);
    duplicateRate_vs_pt[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[5]->Draw();
    }
    adaptEffRange(duplicateRate_vs_pt[i], 1, scaleRangeMax);
  }

  c1->Update();
}
