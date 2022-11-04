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







void defineReconstructionPerformance_v1(
//    const std::vector<std::string>& inputSimParticleFileNames =
//        {
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_90deg/performance_seeding_trees.root",
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_30deg/performance_seeding_trees.root",
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_90deg/performance_seeding_trees.root",
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_30deg/performance_seeding_trees.root",
//	"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/2mu2pi_test/performance_seeding_trees.root",
//	},
//    const std::vector<std::string>& inputTrackSummaryFileNames =
//        {
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_90deg/tracksummary_ckf.root",
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/muon_30deg/tracksummary_ckf.root",
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_90deg/tracksummary_ckf.root",
//	//"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/pion_30deg/tracksummary_ckf.root",
//	"/home/xiaocong/Software/Oscar/acts/build/bin/data/reco_STCF_OscarSim_reco_seeds/landauLoss/2mu2pi_test/tracksummary_ckf.root",
//	},
//    const std::vector<std::string>& trackSummaryFileLegends =
//        {
//	//"muon (#theta=90 deg)",
//	//"muon (#theta=30 deg)",
//	//"pion (#theta=90 deg)",
//	//"pion (#theta=30 deg)",
//	"muon in #psi(2S)->#pi^{+}#pi^{-}J/#psi(#rightarrow#mu^{+}#mu^{-})",
//	},
//

       std::string inputPath = "/home/xiaocong/Software/Oscar/acts/RunSpace/scan/v1.1.testLandauLowPt",
       std::string particle = "mu-",
       //std::vector<std::string> degs = {"90","60","30"},
       std::vector<std::string> degs = {"90"},
       //in GeV
       bool usePt = true,
       //std::vector<double> ps = {0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8},
       std::vector<double> ps = {0.05, 0.075, 0.1, 0.125, 0.15},

       std::map<std::string, std::string> tags = {
        {"90","muon (#theta=90 deg)"},
        {"60","muon (#theta=60 deg)"},
        {"30","muon (#theta=30 deg)"},
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
    unsigned int nHitsMin = 5, unsigned int nMeasurementsMin = 5, unsigned int nOutliersMax = 100,
    double ptMin = 0.03, double absEtaMax=2, double truthMatchProbMin = 0.5, int absPdgId=13) {
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
  //std::vector<double> ptRanges = {20, 0.05, 0.45};
  std::vector<double> ptRanges = {19, 0.05, 1.95};
  //std::vector<double> pts={0.05, 0.15, 0.3, 0.45};


  std::vector<double> ptbins;
  double plow;
  double pup;
  for(int j=0; j < ps.size(); ++j){
    if(j==0)  {
            plow = ps[0] - (ps[1]-ps[0])/2;
            pup = ps[0] + (ps[1]-ps[0])/2;
    } else{
            plow = pup;
            pup = ps[j]*2 - plow;
    }

    ptbins.push_back(plow);
  }
  ptbins.push_back(pup);


  //std::map<std::string, std::vector<TFile*>> inputSimParticleFileNames;
  //std::map<std::string, std::vector<TFile*>> inputTrackSummaryFileNames;
  //std::map<std::string, std::vector<ParticleReader>> pReaders;
  //std::map<std::string, std::vector<TrackSummaryReader>> tReaders;
  std::vector<TFile*> inputSimParticleFileNames;
  std::vector<TFile*> inputTrackSummaryFileNames;

  for(int i=0; i<degs.size(); ++i){
    for(int j=0; j<ps.size(); ++j){

      std::string psstring = std::to_string(static_cast<int>(ps[j]*1000));
      std::string ptUnit = (usePt)? "Pt":"";
      std::string inTrackSummaryFile = inputPath + "/" + particle + "/absThetaDeg_" + degs[i] + "_" + degs[i] + "_momentum" + ptUnit + "Mev_" + psstring + "_" + psstring + "/tracksummary_ckf.root";
      std::cout<<"Reading track summary file " << inTrackSummaryFile << std::endl;
     
      inputTrackSummaryFileNames.push_back(TFile::Open(inTrackSummaryFile.c_str(), "read")); 
    }
  }

   for(int i=0; i<degs.size(); ++i){
    for(int j=0; j<ps.size(); ++j){

      std::string psstring = std::to_string(static_cast<int>(ps[j]*1000));
      std::string ptUnit = (usePt)? "Pt":"";
      std::string inSimParticleFile = inputPath + "/" + particle + "/absThetaDeg_" + degs[i] + "_" + degs[i] + "_momentum" + ptUnit + "Mev_" + psstring + "_" + psstring + "/performance_seeding_trees.root";

      std::cout<<"Reading sim particle file " << inSimParticleFile << std::endl;
      inputSimParticleFileNames.push_back(TFile::Open(inSimParticleFile.c_str(), "read")); 
    }
  }

  std::vector<ParticleReader> pReaders;
  pReaders.reserve(inputSimParticleFileNames.size()); 
  std::vector<TrackSummaryReader> tReaders;
  tReaders.reserve(inputSimParticleFileNames.size()); 

  for (int i=0; i<inputSimParticleFileNames.size(); ++i){
    auto file = inputSimParticleFileNames[i]; 
    pReaders.emplace_back((TTree*)file->Get(simParticleTreeName.c_str()), true);
  }
  for (int i=0; i<inputTrackSummaryFileNames.size(); ++i){
    auto file = inputTrackSummaryFileNames[i]; 
    tReaders.emplace_back((TTree*)file->Get(trackSummaryTreeName.c_str()), true);
  }



/*
  // Check the inputs are valid
  if (inputSimParticleFileNames.size()!=inputTrackSummaryFileNames.size() or inputTrackSummaryFileNames.size() != trackSummaryFileLegends.size()) {
    throw std::invalid_argument(
        "Please specify the legends you want to show for all the track files");
  }

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
*/

  // Define the efficiency plots
  std::vector<TEfficiency*> trackEff_vs_theta;
  std::vector<TEfficiency*> fakeRate_vs_theta;
  std::vector<TEfficiency*> duplicateRate_vs_theta;
  std::vector<TEfficiency*> trackEff_vs_pt;
  std::vector<TEfficiency*> fakeRate_vs_pt;
  std::vector<TEfficiency*> duplicateRate_vs_pt;

  for (int i = 0; i < degs.size(); ++i) {
    trackEff_vs_theta.push_back(new TEfficiency(
        Form("trackeff_vs_theta_%i", i), ";Truth #theta;Efficiency", 40, -4, 4));
    fakeRate_vs_theta.push_back(new TEfficiency(
        Form("fakerate_vs_theta_%i", i), ";#theta;fake rate", 40, -4, 4));
    duplicateRate_vs_theta.push_back(new TEfficiency(
        Form("duplicaterate_vs_theta_%i", i), ";#theta;Duplicate rate", 40, -4, 4));
    trackEff_vs_pt.push_back(new TEfficiency(
        Form("trackeff_vs_pt_%i", i), ";Truth pt [GeV];Efficiency", ptbins.size()-1, ptbins.data()));
    fakeRate_vs_pt.push_back(new TEfficiency(
        Form("fakerate_vs_pt_%i", i), ";pt [GeV];fake rate", ptbins.size()-1, ptbins.data()));
    duplicateRate_vs_pt.push_back(new TEfficiency(
        Form("duplicaterate_vs_pt_%i", i), ";pt [GeV];Duplicate rate", ptbins.size()-1, ptbins.data()));
  }

  // Set styles
  for (int i = 0; i < degs.size(); ++i) {
    auto color = colors[i];
    auto marker  = markers[i];
    setEffStyle(trackEff_vs_theta[i], color, marker);
    setEffStyle(fakeRate_vs_theta[i], color, marker);
    setEffStyle(duplicateRate_vs_theta[i], color, marker);
    setEffStyle(trackEff_vs_pt[i], color, marker);
    setEffStyle(fakeRate_vs_pt[i], color, marker);
    setEffStyle(duplicateRate_vs_pt[i], color, marker);
  }



  // Loop over the track files
  for (unsigned int iFile = 0; iFile < tReaders.size();++iFile) {
     int i = iFile/ps.size(); 
   //for(int i=0; i<degs.size(); ++i){
  //  for(int j=0; j<ps.size(); ++j){
  //    int iFile = i*ps.size() + j;
      //std::string deg = degs[i];
      //std::cout<<"deg =  " << deg << std::endl; 
      std::cout << "Processing track file: " << inputTrackSummaryFileNames[iFile]->GetName()
                << std::endl;
  
      // The particles in each event
      std::map<unsigned int, std::vector<ParticleInfo>> particles;


      // The container from track-particle matching info (Flushed per event)
      std::map<uint64_t, std::vector<RecoTrackInfo>> matchedParticles;

      //inputTrackSummaryFileNames[iFile]->cd();

      // Loop over the events to fill plots
      for (size_t k = 0; k < tReaders[iFile].tree->GetEntries(); ++k) {
        //if (k % 10 == 0) {
        //  std::cout << "Processed events: " << k << std::endl;
        //}

        // Get the tracks
        tReaders[iFile].tree->GetEntry(k);

        // Get the particles (do nothing if the particles for this event already
        // read)
        auto it = particles.find(k);
        if (it == particles.end()) {
          particles.emplace(k, pReaders[iFile].getParticles(k));
        }

        // Loop over the tracks
        // The fake rate is defined as the ratio of selected truth-matched tracks
        // over all selected tracks
        for (size_t l = 0; l < tReaders[iFile].nStates->size(); ++l) {
          bool hasFittedParameters = tReaders[iFile].hasFittedParams->at(l);
          auto nMeasurements = tReaders[iFile].nMeasurements->at(l);
          auto nOutliers = tReaders[iFile].nOutliers->at(l);
          auto nHoles = tReaders[iFile].nHoles->at(l);
          auto theta = tReaders[iFile].eTHETA_fit->at(l);
          auto qop = tReaders[iFile].eQOP_fit->at(l);
          auto pt = std::abs(1 / qop) * std::sin(theta);
          auto eta = std::atanh(std::cos(theta));
          auto nMajorityHits = tReaders[iFile].nMajorityHits->at(l);
          auto majorityParticleId = tReaders[iFile].majorityParticleId->at(l);

          // Select the track, e.g. you might also want to add cuts on the
          // nOutliers, nHoles
          if ((!hasFittedParameters) or nMeasurements < nMeasurementsMin or
               nOutliers > nOutliersMax or pt < ptMin) {
            continue;
          }

    //      // Fill the fake rate plots
          if (nMajorityHits * 1. / nMeasurements >= truthMatchProbMin) {
            matchedParticles[majorityParticleId].push_back(
                {eta, pt, nMajorityHits, nMeasurements});
            fakeRate_vs_theta[i]->Fill(false, theta);
            fakeRate_vs_pt[i]->Fill(false, pt);
          } else {
            fakeRate_vs_theta[i]->Fill(true, theta);
            fakeRate_vs_pt[i]->Fill(true, pt);
          }
        }  // end of all tracks in this event


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
          for (size_t m = 0; m < matchedTracks.size(); ++m) {
            auto eta = matchedTracks[m].eta;
            auto pt = matchedTracks[m].pt;
            double theta = std::atan(std::exp(-eta))*2;
            if (m == 0) {
              duplicateRate_vs_theta[i]->Fill(false, theta);
              duplicateRate_vs_pt[i]->Fill(false, pt);
            } else {
              duplicateRate_vs_theta[i]->Fill(true, theta);
              duplicateRate_vs_pt[i]->Fill(true, pt);
            }
          }
        }  // end of all selected truth-matched tracks in this event

        // Loop over all truth particles in this event
        // The effiency is define as the ratio of selected particles that have
        // been matched with reco
        if(particles[k].size()==0){
           std::cout<<"Skip event " << i << " as there is no truth particle in this event" << std::endl; 
        } 
        for (const auto& particle : particles[k]) {
          auto nHits = particle.nHits;
          auto eta = particle.eta;
          auto pt = particle.pt;
          if (absPdgId!=999 and abs(particle.particlePdg) != absPdgId){
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
            trackEff_vs_theta[i]->Fill(true, theta);
            trackEff_vs_pt[i]->Fill(true, pt);
          } else {
            trackEff_vs_theta[i]->Fill(false, theta);
            trackEff_vs_pt[i]->Fill(false, pt);
          }
        }  // end of all particles

        matchedParticles.clear();
      }  // end of all events
    }    // end of all files 

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
  for (int i = 0; i < degs.size(); ++i) {
    legs[0]->AddEntry(trackEff_vs_theta[i],
                      Form("%s", tags[degs[i]].c_str()), "lp");
    legs[1]->AddEntry(fakeRate_vs_theta[i],
                      Form("%s", tags[degs[i]].c_str()), "lp");
    legs[2]->AddEntry(duplicateRate_vs_theta[i],
                      Form("%s", tags[degs[i]].c_str()), "lp");
    legs[3]->AddEntry(trackEff_vs_pt[i],
                      Form("%s", tags[degs[i]].c_str()), "lp");
    legs[4]->AddEntry(fakeRate_vs_pt[i],
                      Form("%s", tags[degs[i]].c_str()), "lp");
    legs[5]->AddEntry(duplicateRate_vs_pt[i],
                      Form("%s", tags[degs[i]].c_str()), "lp");
  }

  // Now draw the plots
  TCanvas* c1 = new TCanvas("recoPerf", " ", 1500, 800);
  c1->Divide(3, 2);

  auto nTrackFiles = degs.size();
  float scaleRangeMax = 1.1;
  for (int i = 0; i < nTrackFiles; ++i) {
    std::string mode = (i == 0) ? "" : "same";
    c1->cd(1);
    trackEff_vs_theta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[0]->Draw();
    }
    adaptEffRange(trackEff_vs_theta[i], 1, scaleRangeMax);

    c1->cd(2);
    fakeRate_vs_theta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[1]->Draw();
    }
    adaptEffRange(fakeRate_vs_theta[i], 1, scaleRangeMax);

    c1->cd(3);
    duplicateRate_vs_theta[i]->Draw(mode.c_str());
    if (i == nTrackFiles - 1) {
      legs[2]->Draw();
    }
    adaptEffRange(duplicateRate_vs_theta[i], 1, scaleRangeMax);

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
