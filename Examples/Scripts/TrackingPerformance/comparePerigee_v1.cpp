// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <bitset>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>
#include <TVectorF.h>

//#include "CommonUtils.h"
#include "TreeReader.h"

using namespace ROOT;

template <typename T>
void setThisHistStyle(T* hist, short color = 1, short marker = 20,
                      float xTitleSize = 0.05, float yTitleSize = 0.05,
                      float xLabelSize = 0.05, float yLabelSize = 0.05,
                      float xTitleOffset = 1.2, float yTitleOffset = 1.2,
                      int nDiv = 510) {
  hist->GetXaxis()->SetTitleSize(xTitleSize);
  hist->GetYaxis()->SetTitleSize(yTitleSize);
  hist->GetXaxis()->SetLabelSize(xLabelSize);
  hist->GetYaxis()->SetLabelSize(yLabelSize);
  hist->GetXaxis()->SetTitleOffset(xTitleOffset);
  hist->GetYaxis()->SetTitleOffset(yTitleOffset);
  hist->GetXaxis()->SetNdivisions(510);
  hist->GetYaxis()->SetNdivisions(nDiv);
  hist->SetMarkerStyle(marker);
  hist->SetMarkerSize(1.0);
  hist->SetLineWidth(2);
  hist->SetTitle("");
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
}

void myText(Double_t x, Double_t y, Color_t color, float font,
            const char* text) {
  Double_t tsize = 0.05;
  TLatex l;  //
  l.SetTextAlign(12);
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextSize(font);
  l.SetTextColor(color);
  l.DrawLatex(x, y, text);
}

void setYRange(TH1F* hist, double yMinScale = 0.5, double yMaxScale = 1.5) {
  int binmin = hist->GetMinimumBin();
  int binmax = hist->GetMaximumBin();
  double binMinContent = hist->GetBinContent(binmin);
  double binMaxContent = hist->GetBinContent(binmax);
  hist->GetYaxis()->SetRangeUser(binMinContent * yMinScale,
                                 binMaxContent * yMaxScale);
}

void anaHisto(TH1D* inputHist, int j, TH1F* meanHist, TH1F* widthHist,
              bool fit = false, double scale = 1) {
  // evaluate mean and width via the Gauss fit
  assert(inputHist != nullptr);
  std::cout << "anaHisto for " << inputHist->GetName() << std::endl;
  if (inputHist->GetEntries() > 0) {
    if (fit) {
      TFitResultPtr r = inputHist->Fit("gaus", "QS0");
      if (r.Get() and ((r->Status() % 1000) == 0)) {
        // fill the mean and width into 'j'th bin of the meanHist and widthHist,
        // respectively
        meanHist->SetBinContent(j, r->Parameter(1));
        meanHist->SetBinError(j, r->ParError(1));
        widthHist->SetBinContent(j, r->Parameter(2) * scale);
        widthHist->SetBinError(j, r->ParError(2));
      }
    } else {
      meanHist->SetBinContent(j, inputHist->GetMean());
      meanHist->SetBinError(j, inputHist->GetMeanError());
      // meanHist->SetBinError(j, inputHist->GetRMS());

      widthHist->SetBinContent(j, inputHist->GetRMS() * scale);
      widthHist->SetBinError(j, inputHist->GetRMSError());
    }
  }
}

void comparePerigee_v1(std::string inFile1 =
                           "/home/xiaocong/Software/Oscar/acts/build/bin/data/"
                           "reco_STCF_OscarSim_reco_seeds/landauLoss/"
                           "muon_90deg/tracksummary_ckf.root",
                       std::string inFile2 =
                           "/home/xiaocong/Software/Oscar/acts/build/bin/data/"
                           "reco_STCF_OscarSim_reco_seeds/landauLoss/"
                           "muon_30deg/tracksummary_ckf.root",
                       std::string inFile3 =
                           "/home/xiaocong/Software/Oscar/acts/build/bin/data/"
                           "reco_STCF_OscarSim_reco_seeds/landauLoss/"
                           "pion_90deg/tracksummary_ckf.root",
                       std::string inFile4 =
                           "/home/xiaocong/Software/Oscar/acts/build/bin/data/"
                           "reco_STCF_OscarSim_reco_seeds/landauLoss/"
                           "pion_30deg/tracksummary_ckf.root",
                       double absEtaMin = 0, double absEtaMax = 1.75,
                       double ptMin = 0.05, double ptMax = 1.8,
                       bool saveAs = false, bool showEta = false,
                       bool showPt = true, bool fit = true,
                       bool plotResidual = true,
                       // plotType 0: mean,   1:width,   2: mean and width
                       int plotType = 1, bool plotResidualRatio = false,
                       bool absEta = true, bool variablePtBin = true,
                       std::vector<std::string> legs =
                           {
                               "muon (#theta=90 deg)",
                               "muon (#theta=30 deg)",
                               "pion (#theta=90 deg)",
                               "pion (#theta=30 deg)",
                           },
                       std::vector<int> colors = {872, 866, 854, 896},
                       std::vector<int> markers = {24, 26, 20, 22}) {
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetNdivisions(510, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  if (plotType == 2 or not plotResidual)
    plotResidualRatio = false;

  if (plotType != 0 and plotType != 1 and plotType != 2) {
    throw std::invalid_argument("The plotType must be 0, 1 or 2!");
  }

  int etaMinT = static_cast<int>(absEtaMin * 10);
  int etaMaxT = static_cast<int>(absEtaMax * 10);
  int ptMinT = static_cast<int>(ptMin);
  int ptMaxT = static_cast<int>(ptMax);
  std::string etaTag = "_etaMin_" + std::to_string(etaMinT) + "_etaMax_" +
                       std::to_string(etaMaxT);
  std::string ptTag =
      "_ptMin_" + std::to_string(ptMinT) + "_ptMax_" + std::to_string(ptMaxT);

  std::string path = etaTag + ptTag;

  gSystem->Exec(Form("mkdir %s", path.c_str()));

  std::string ratioTag = "";
  if (plotResidualRatio and plotType == 1 and plotResidual) {
    ratioTag = "_ratio";
  }

  std::string detector = "ODD_ATLASField_Fatras";
  // std::string detector = "Generic_ATLASField";
  std::string variable = plotResidual ? "res" : "pull";
  detector += "_" + variable;

  double xTitleSize = 0.05;
  double yTitleSize = 0.05;
  double xLabelSize = 0.05;
  double yLabelSize = 0.05;
  double xTitleOffset = 1.2;
  double yTitleOffset = plotResidual ? 1.9 : 1.7;

  double topXTitleSize = 0.075;
  double topYTitleSize = 0.075;
  double topXLabelSize = 0.075;
  double topYLabelSize = 0.075;
  double topXTitleOffset = 1.1;
  double topYTitleOffset = 0.95;

  double botXTitleSize = 0.1;
  double botYTitleSize = 0.1;
  double botXLabelSize = 0.1;
  double botYLabelSize = 0.1;
  double botXTitleOffset = 1.3;
  double botYTitleOffset = 0.7;

  // Load the tree chain
  TChain* treeChain1 = new TChain("tracksummary");
  TChain* treeChain2 = new TChain("tracksummary");
  TChain* treeChain3 = new TChain("tracksummary");
  TChain* treeChain4 = new TChain("tracksummary");

  std::string inpath = "";
  inFile1 = inpath + inFile1;
  inFile2 = inpath + inFile2;
  inFile3 = inpath + inFile3;
  inFile4 = inpath + inFile4;

  treeChain1->Add(inFile1.c_str());
  treeChain2->Add(inFile2.c_str());
  treeChain3->Add(inFile3.c_str());
  treeChain4->Add(inFile4.c_str());

  if (treeChain1->GetEntries() == 0) {
    std::cout << "[x] No entries found ... " << std::endl;
    return -1;
  }

  //=======================================================================
  std::map<std::string, TH2F*> ref1_vs_eta;
  std::map<std::string, TH2F*> ref1_vs_pt;
  std::map<std::string, TH2F*> ref1_vs_d0;
  std::map<std::string, TH2F*> ref1_vs_z0;

  std::map<std::string, TH1F*> ref1mean_vs_eta;
  std::map<std::string, TH1F*> ref1mean_vs_pt;
  std::map<std::string, TH1F*> ref1mean_vs_d0;
  std::map<std::string, TH1F*> ref1mean_vs_z0;
  std::map<std::string, TH1F*> ref1width_vs_eta;
  std::map<std::string, TH1F*> ref1width_vs_pt;
  std::map<std::string, TH1F*> ref1width_vs_d0;
  std::map<std::string, TH1F*> ref1width_vs_z0;

  std::map<std::string, TH1F*> ref1meanwidth_vs_eta;
  std::map<std::string, TH1F*> ref1meanwidth_vs_pt;

  std::map<std::string, TH2F*> ref2_vs_eta;
  std::map<std::string, TH2F*> ref2_vs_pt;
  std::map<std::string, TH2F*> ref2_vs_d0;
  std::map<std::string, TH2F*> ref2_vs_z0;

  std::map<std::string, TH1F*> ref2mean_vs_eta;
  std::map<std::string, TH1F*> ref2mean_vs_pt;
  std::map<std::string, TH1F*> ref2mean_vs_d0;
  std::map<std::string, TH1F*> ref2mean_vs_z0;
  std::map<std::string, TH1F*> ref2width_vs_eta;
  std::map<std::string, TH1F*> ref2width_vs_pt;
  std::map<std::string, TH1F*> ref2width_vs_d0;
  std::map<std::string, TH1F*> ref2width_vs_z0;

  std::map<std::string, TH1F*> ref2meanwidth_vs_eta;
  std::map<std::string, TH1F*> ref2meanwidth_vs_pt;

  //=======================================================================
  std::map<std::string, TH2F*> val1_vs_eta;
  std::map<std::string, TH2F*> val1_vs_pt;
  std::map<std::string, TH2F*> val1_vs_d0;
  std::map<std::string, TH2F*> val1_vs_z0;

  std::map<std::string, TH1F*> val1mean_vs_eta;
  std::map<std::string, TH1F*> val1mean_vs_pt;
  std::map<std::string, TH1F*> val1mean_vs_d0;
  std::map<std::string, TH1F*> val1mean_vs_z0;
  std::map<std::string, TH1F*> val1width_vs_eta;
  std::map<std::string, TH1F*> val1width_vs_pt;
  std::map<std::string, TH1F*> val1width_vs_d0;
  std::map<std::string, TH1F*> val1width_vs_z0;

  std::map<std::string, TH1F*> val1meanwidth_vs_eta;
  std::map<std::string, TH1F*> val1meanwidth_vs_pt;

  std::map<std::string, TH2F*> val2_vs_eta;
  std::map<std::string, TH2F*> val2_vs_pt;
  std::map<std::string, TH2F*> val2_vs_d0;
  std::map<std::string, TH2F*> val2_vs_z0;

  std::map<std::string, TH1F*> val2mean_vs_eta;
  std::map<std::string, TH1F*> val2mean_vs_pt;
  std::map<std::string, TH1F*> val2mean_vs_d0;
  std::map<std::string, TH1F*> val2mean_vs_z0;
  std::map<std::string, TH1F*> val2width_vs_eta;
  std::map<std::string, TH1F*> val2width_vs_pt;
  std::map<std::string, TH1F*> val2width_vs_d0;
  std::map<std::string, TH1F*> val2width_vs_z0;

  std::map<std::string, TH1F*> val2meanwidth_vs_eta;
  std::map<std::string, TH1F*> val2meanwidth_vs_pt;

  //=======================================================================

  double etaRange = 3.;
  int etaBins = 12;
  // std::vector<double> pts = {0.05, 0.1, 0.15, 0.2, 0.4, 0.6,
  // 0.8, 1.0, 1.5, 2.0};
  std::vector<double> pts = {0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0};
  std::vector<double> ps;
  std::vector<double> ps_;
  for (int i = 0; i < pts.size() - 1; ++i) {
    ps.push_back((pts[i] + pts[i + 1]) / 2);
    ps_.push_back((pts[i] + pts[i + 1]) / 2 / sin(30.0 / 180 * M_PI));
  }

  std::pair<double, double> pullRange = {-5, 5};
  std::vector<std::pair<double, double>> resRanges = {
      {-1, 1}, {-1, 1}, {-0.02, 0.02}, {-0.02, 0.02}, {-0.1, 0.1}, {-3.5, 3.5},
  };

  std::vector<std::string> names = {"l0", "l1", "phi", "theta", "qop", "t"};
  std::vector<std::string> pullTitles = {
      "(#frac{d_{0}^{fit} - d_{0}^{truth}}{#sigma(d_{0})})",
      "(#frac{z_{0}^{fit} - z_{0}^{truth}}{#sigma(z_{0})})",
      "(#frac{#phi^{fit} - #phi^{truth}}{#sigma(#phi)})",
      "(#frac{#theta^{fit} - #theta^{truth}}{#sigma(#theta)})",
      "(#frac{(q/p)^{fit} - (q/p)^{truth}}{#sigma(q/p)})",
      "(#frac{t^{fit} - t^{truth}}{#sigma(t)})",
  };

  std::vector<std::string> resTitles = {
      "(d_{0}^{fit} - d_{0}^{truth}) [mm]",
      "(z_{0}^{fit} - z_{0}^{truth}) [mm]",
      "(#phi^{fit} - #phi^{truth}) [rad]",
      "(#theta^{fit} - #theta^{truth}) [rad]",
      //"((q/p)^{fit} - (q/p)^{truth}) [GeV^{-1}]",
      "((q/p)^{fit} - (q/p)^{truth})xp^{truth}",
      "(t^{fit} - t^{truth}) [ns]",
  };

  std::vector<std::string> resRatioTitles = {
      "(d_{0}^{fit} - d_{0}^{truth})", "(z_{0}^{fit} - z_{0}^{truth})",
      "(#phi^{fit} - #phi^{truth})",   "(#theta^{fit} - #theta^{truth})",
      "((q/p)^{fit} - (q/p)^{truth})", "(t^{fit} - t^{truth})",
  };

  std::vector<std::pair<double, double>> yRange_resmean_vs_eta;
  std::vector<std::pair<double, double>> yRange_resmean_vs_pt;
  std::vector<std::pair<double, double>> yRange_reswidth_vs_eta;
  std::vector<std::pair<double, double>> yRange_reswidth_vs_pt;
  std::vector<std::pair<double, double>> yRange_pullmean_vs_eta;
  std::vector<std::pair<double, double>> yRange_pullmean_vs_pt;
  std::vector<std::pair<double, double>> yRange_pullwidth_vs_eta;
  std::vector<std::pair<double, double>> yRange_pullwidth_vs_pt;

  if (1) {
    // y axis range of resmean_vs_eta plots
    yRange_resmean_vs_eta = {
        //{-0.1, 0.15},     {-0.2, 0.45},  {-0.0003, 0.0006},
        {-0.1, 0.1},      {-0.2, 0.3},   {-0.0003, 0.0006},
        {-0.0006, 0.001}, {-2e-4, 3e-4}, {-0.10, 0.2},
    };
    // y axis range of resmean_vs_pt plots
    yRange_resmean_vs_pt = {
        {-0.03, 0.03},
        {-0.02, 0.03},
        {-0.15e-3, 0.3e-3},
        {-2e-5, 0.5e-4},
        //{-1.5e-4, 1.5e-4},
        {-0.005, 0.005},
        {-0.15, 0.2},
    };

    //================================================================================
    // y axis range of reswidth_vs_eta plots
    yRange_reswidth_vs_eta = {
        //{-0.04, 0.2},     {-0.1, 0.6},  {-0.0002, 0.001},
        {0.0, 0.2},       {0.0, 0.4},    {-0.0002, 0.001},
        {-0.0002, 0.001}, {0.0, 0.0015}, {0.95, 1.1},
    };
    // y axis range of reswidth_vs_pt plots
    yRange_reswidth_vs_pt = {
        {0.0, 0.5}, {0.25, 0.6}, {0.0, 0.015},
        {0, 0.015}, {0, 0.02},   {0.95, 1.1},
    };

    //================================================================================
    // y axis range of pullmean_vs_eta plots
    yRange_pullmean_vs_eta = {
        {-2, 2}, {-3, 5}, {-1, 3}, {-3, 5}, {-0.2, 0.6}, {-0.1, 0.15},
    };
    // y axis range of pullmean_vs_pt plots
    yRange_pullmean_vs_pt = {
        {-0.8, 0.5}, {-0.2, 0.4}, {-0.3, 1.0},
        {-0.3, 0.4}, {-0.3, 0.5}, {-0.15, 0.2},
    };

    //================================================================================
    // y axis range of pullwidth_vs_eta plots
    yRange_pullwidth_vs_eta = {
        //{0.5, 3.}, {0.5, 4}, {0.5, 3}, {0.5, 4}, {0.9, 1.4}, {0.9, 1.2},
        {0.5, 2.5}, {0.5, 3}, {0.5, 2.5}, {0.5, 3.5}, {0.9, 1.3}, {0.9, 1.2},
    };
    // y axis range of pullwidth_vs_pt plots
    yRange_pullwidth_vs_pt = {
        //{0.9, 1.5},
        {0.8, 2.0}, {0.5, 3.5},  // 0.8, 2.5
        {0.8, 2.1}, {0.8, 3.0}, {0.85, 1.3}, {0.9, 1.2},
    };
    //================================================================================
  } else {
    // y axis range of resmean_vs_eta plots
    yRange_resmean_vs_eta = {
        {-0.05, 0.05},     {-0.1, 0.2},   {-0.00015, 0.0003},
        {-0.0002, 0.0003}, {-1e-4, 2e-4}, {-0.10, 0.2},
    };
    // y axis range of resmean_vs_pt plots
    yRange_resmean_vs_pt = {
        {-0.01, 0.01},
        {-0.015, 0.02},
        {-0.07e-3, 0.15e-3},
        {-2e-5, 0.5e-4},
        //{-1.5e-4, 1.5e-4},
        {-0.003, 0.004},
        {-0.15, 0.2},
    };

    //================================================================================
    // y axis range of reswidth_vs_eta plots
    yRange_reswidth_vs_eta = {
        {-0.02, 0.1},      {-0.1, 0.4},   {-0.0002, 0.001},
        {-0.0002, 0.0006}, {0.0, 0.0015}, {0.95, 1.1},
    };
    // y axis range of reswidth_vs_pt plots
    yRange_reswidth_vs_pt = {
        {0.0, 0.14},        {0.0, 0.25},    {-0.0002, 0.0015},
        {-0.00015, 0.0012}, {-0.003, 0.02}, {0.95, 1.1},
    };

    //================================================================================
    // y axis range of pullmean_vs_eta plots
    yRange_pullmean_vs_eta = {
        {-1.5, 1.5}, {-2, 3}, {-0.5, 1.5}, {-2, 3}, {-0.1, 0.3}, {-0.05, 0.1},
    };
    // y axis range of pullmean_vs_pt plots
    yRange_pullmean_vs_pt = {
        {-0.4, 0.4}, {-0.1, 0.15}, {-0.2, 0.5},
        {-0.1, 0.2}, {-0.3, 0.3},  {-0.1, 0.1},
    };

    //================================================================================
    // y axis range of pullwidth_vs_eta plots
    yRange_pullwidth_vs_eta = {
        {0.5, 2.5}, {0.6, 3.0}, {0.8, 2.2}, {0.5, 3.5}, {0.9, 1.2}, {0.95, 1.1},
    };
    // y axis range of pullwidth_vs_pt plots
    yRange_pullwidth_vs_pt = {
        //{0.9, 1.5},
        {0.9, 1.4}, {0.8, 2.5},  // 0.8, 2.5
        {0.8, 1.7}, {0.8, 2.5}, {0.9, 1.2}, {0.95, 1.1},
    };
  }

  auto timeInns = [](double timeInmm) -> double {
    return timeInmm / 299792458000.0 * 1e+9;
  };

  // book hists
  for (int i = 0; i < names.size(); ++i) {
    const auto& name = names[i];

    // The 2D hist y range
    std::pair<double, double> yRange = pullRange;
    if (plotResidual) {
      yRange = resRanges[i];
    }

    // ref1 plots
    ref1_vs_eta.emplace(name, new TH2F(Form("ref1_%s_vs_eta", name.c_str()), "",
                                       etaBins, -1.0 * etaRange, etaRange, 100,
                                       yRange.first, yRange.second));
    if (variablePtBin) {
      ref1_vs_pt.emplace(name, new TH2F(Form("ref1_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.data(), 100,
                                        yRange.first, yRange.second));
    } else {
      ref1_vs_pt.emplace(name, new TH2F(Form("ref1_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.front(), pts.back(),
                                        100, yRange.first, yRange.second));
    }
    ref1_vs_d0.emplace(
        name, new TH2F(Form("ref1_%s_vs_d0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    ref1_vs_z0.emplace(
        name, new TH2F(Form("ref1_%s_vs_z0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));

    ref1mean_vs_eta.emplace(
        name, new TH1F(Form("ref1mean_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      ref1mean_vs_pt.emplace(
          name, new TH1F(Form("ref1mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      ref1mean_vs_pt.emplace(
          name, new TH1F(Form("ref1mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    ref1mean_vs_d0.emplace(
        name,
        new TH1F(Form("ref1mean_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    ref1mean_vs_z0.emplace(
        name,
        new TH1F(Form("ref1mean_%s_vs_z0", name.c_str()), "", 100, -100, 100));
    ref1width_vs_eta.emplace(
        name, new TH1F(Form("ref1width_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      ref1width_vs_pt.emplace(
          name, new TH1F(Form("ref1width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      ref1width_vs_pt.emplace(
          name, new TH1F(Form("ref1width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    ref1width_vs_d0.emplace(
        name,
        new TH1F(Form("ref1width_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    ref1width_vs_z0.emplace(
        name,
        new TH1F(Form("ref1width_%s_vs_z0", name.c_str()), "", 100, -100, 100));

    ////////////////////////////
    // ref2 plots
    ////////////////////////////
    ref2_vs_eta.emplace(name, new TH2F(Form("ref2_%s_vs_eta", name.c_str()), "",
                                       etaBins, -1.0 * etaRange, etaRange, 100,
                                       yRange.first, yRange.second));
    if (variablePtBin) {
      ref2_vs_pt.emplace(name, new TH2F(Form("ref2_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.data(), 100,
                                        yRange.first, yRange.second));
    } else {
      ref2_vs_pt.emplace(name, new TH2F(Form("ref2_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.front(), pts.back(),
                                        100, yRange.first, yRange.second));
    }
    ref2_vs_d0.emplace(
        name, new TH2F(Form("ref2_%s_vs_d0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    ref2_vs_z0.emplace(
        name, new TH2F(Form("ref2_%s_vs_z0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    ref2mean_vs_eta.emplace(
        name, new TH1F(Form("ref2mean_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      ref2mean_vs_pt.emplace(
          name, new TH1F(Form("ref2mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      ref2mean_vs_pt.emplace(
          name, new TH1F(Form("ref2mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    ref2mean_vs_d0.emplace(
        name,
        new TH1F(Form("ref2mean_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    ref2mean_vs_z0.emplace(
        name,
        new TH1F(Form("ref2mean_%s_vs_z0", name.c_str()), "", 100, -100, 100));
    ref2width_vs_eta.emplace(
        name, new TH1F(Form("ref2width_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      ref2width_vs_pt.emplace(
          name, new TH1F(Form("ref2width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      ref2width_vs_pt.emplace(
          name, new TH1F(Form("ref2width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    ref2width_vs_d0.emplace(
        name,
        new TH1F(Form("ref2width_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    ref2width_vs_z0.emplace(
        name,
        new TH1F(Form("ref2width_%s_vs_z0", name.c_str()), "", 100, -100, 100));

    ////////////////////////////
    // val1 plots
    ////////////////////////////
    val1_vs_eta.emplace(name, new TH2F(Form("val1_%s_vs_eta", name.c_str()), "",
                                       etaBins, -1.0 * etaRange, etaRange, 100,
                                       yRange.first, yRange.second));
    if (variablePtBin) {
      val1_vs_pt.emplace(name, new TH2F(Form("val1_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.data(), 100,
                                        yRange.first, yRange.second));
    } else {
      val1_vs_pt.emplace(name, new TH2F(Form("val1_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.front(), pts.back(),
                                        100, yRange.first, yRange.second));
    }
    val1_vs_d0.emplace(
        name, new TH2F(Form("val1_%s_vs_d0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    val1_vs_z0.emplace(
        name, new TH2F(Form("val1_%s_vs_z0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    val1mean_vs_eta.emplace(
        name, new TH1F(Form("val1mean_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      val1mean_vs_pt.emplace(
          name, new TH1F(Form("val1mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      val1mean_vs_pt.emplace(
          name, new TH1F(Form("val1mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    val1mean_vs_d0.emplace(
        name,
        new TH1F(Form("val1mean_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    val1mean_vs_z0.emplace(
        name,
        new TH1F(Form("val1mean_%s_vs_z0", name.c_str()), "", 100, -100, 100));
    val1width_vs_eta.emplace(
        name, new TH1F(Form("val1width_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      val1width_vs_pt.emplace(
          name, new TH1F(Form("val1width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      val1width_vs_pt.emplace(
          name, new TH1F(Form("val1width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    val1width_vs_d0.emplace(
        name,
        new TH1F(Form("val1width_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    val1width_vs_z0.emplace(
        name,
        new TH1F(Form("val1width_%s_vs_z0", name.c_str()), "", 100, -100, 100));

    ////////////////////////////
    // val2 plots
    ////////////////////////////
    val2_vs_eta.emplace(name, new TH2F(Form("val2_%s_vs_eta", name.c_str()), "",
                                       etaBins, -1.0 * etaRange, etaRange, 100,
                                       yRange.first, yRange.second));
    if (variablePtBin) {
      val2_vs_pt.emplace(name, new TH2F(Form("val2_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.data(), 100,
                                        yRange.first, yRange.second));
    } else {
      val2_vs_pt.emplace(name, new TH2F(Form("val2_%s_vs_pt", name.c_str()), "",
                                        pts.size() - 1, pts.front(), pts.back(),
                                        100, yRange.first, yRange.second));
    }
    val2_vs_d0.emplace(
        name, new TH2F(Form("val2_%s_vs_d0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    val2_vs_z0.emplace(
        name, new TH2F(Form("val2_%s_vs_z0", name.c_str()), "", 100, -100, 100,
                       100, yRange.first, yRange.second));
    val2mean_vs_eta.emplace(
        name, new TH1F(Form("val2mean_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      val2mean_vs_pt.emplace(
          name, new TH1F(Form("val2mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      val2mean_vs_pt.emplace(
          name, new TH1F(Form("val2mean_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    val2mean_vs_d0.emplace(
        name,
        new TH1F(Form("val2mean_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    val2mean_vs_z0.emplace(
        name,
        new TH1F(Form("val2mean_%s_vs_z0", name.c_str()), "", 100, -100, 100));
    val2width_vs_eta.emplace(
        name, new TH1F(Form("val2width_%s_vs_eta", name.c_str()), "", etaBins,
                       -1.0 * etaRange, etaRange));
    if (variablePtBin) {
      val2width_vs_pt.emplace(
          name, new TH1F(Form("val2width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.data()));
    } else {
      val2width_vs_pt.emplace(
          name, new TH1F(Form("val2width_%s_vs_pt", name.c_str()), "",
                         pts.size() - 1, pts.front(), pts.back()));
    }
    val2width_vs_d0.emplace(
        name,
        new TH1F(Form("val2width_%s_vs_d0", name.c_str()), "", 100, -100, 100));
    val2width_vs_z0.emplace(
        name,
        new TH1F(Form("val2width_%s_vs_z0", name.c_str()), "", 100, -100, 100));
  }

  TrackSummaryReader tracks1(treeChain1, true);
  TrackSummaryReader tracks2(treeChain2, true);
  TrackSummaryReader tracks3(treeChain3, true);
  TrackSummaryReader tracks4(treeChain4, true);

  unsigned long entries1 = 10000;
  unsigned long entries2 = 10000;
  unsigned long entries3 = 10000;
  unsigned long entries4 = 10000;

  // Tree1
  std::cout << "entries1 = " << entries1 << std::endl;
  std::cout << "entries2 = " << entries2 << std::endl;
  std::cout << "entries3 = " << entries3 << std::endl;
  std::cout << "entries4 = " << entries4 << std::endl;

  //  for(const auto& [name, h2]: ref1_vs_eta){
  //     h2->StatOverflows(kTRUE);
  //  }
  //  for(const auto& [name, h2]: ref2_vs_eta){
  //     h2->StatOverflows(kTRUE);
  //  }
  //  for(const auto& [name, h2]: val1_vs_eta){
  //     h2->StatOverflows(kTRUE);
  //  }
  //  for(const auto& [name, h2]: val2_vs_eta){
  //     h2->StatOverflows(kTRUE);
  //  }
  //
  //  for(const auto& [name, h2]: ref1_vs_pt){
  //     h2->StatOverflows(kTRUE);
  //  }
  //  for(const auto& [name, h2]: ref2_vs_pt){
  //     h2->StatOverflows(kTRUE);
  //  }
  //  for(const auto& [name, h2]: val1_vs_pt){
  //     h2->StatOverflows(kTRUE);
  //  }
  //  for(const auto& [name, h2]: val2_vs_pt){
  //     h2->StatOverflows(kTRUE);
  //  }

  // Lambda to fill the plots
  auto fillHists = [&](unsigned long entries, const TrackSummaryReader& tracks,
                       std::map<std::string, TH2F*>& etaHists,
                       std::map<std::string, TH2F*>& ptHists) {
    for (unsigned long ie = 0; ie < entries; ++ie) {
      // Make sure you have the entry
      tracks.tree->GetEntry(ie);
      size_t nTracks = tracks.hasFittedParams->size();
      for (size_t it = 0; it < nTracks; ++it) {
        if (tracks.hasFittedParams->at(it)) {
          bool withinPt =
              tracks.t_pT->at(it) >= ptMin and tracks.t_pT->at(it) < ptMax;
          bool withinEta =
              absEta ? (std::abs(tracks.t_eta->at(it)) >= absEtaMin and
                        std::abs(tracks.t_eta->at(it)) < absEtaMax)
                     : (tracks.t_eta->at(it) >= absEtaMin and
                        tracks.t_eta->at(it) < absEtaMax);
          if (plotResidual) {
            if (withinPt) {
              etaHists["l0"]->Fill(tracks.t_eta->at(it),
                                   tracks.res_eLOC0_fit->at(it));
              etaHists["l1"]->Fill(tracks.t_eta->at(it),
                                   tracks.res_eLOC1_fit->at(it));
              etaHists["phi"]->Fill(tracks.t_eta->at(it),
                                    tracks.res_ePHI_fit->at(it));
              etaHists["theta"]->Fill(tracks.t_eta->at(it),
                                      tracks.res_eTHETA_fit->at(it));
              etaHists["qop"]->Fill(tracks.t_eta->at(it),
                                    tracks.res_eQOP_fit->at(it));
              // tracks.res_eQOP_fit->at(it)*tracks.t_p->at(it));
              etaHists["t"]->Fill(tracks.t_eta->at(it),
                                  timeInns(tracks.res_eT_fit->at(it)));
            }

            // if (withinEta) {
            if (1) {
              ptHists["l0"]->Fill(tracks.t_pT->at(it),
                                  tracks.res_eLOC0_fit->at(it));
              ptHists["l1"]->Fill(tracks.t_pT->at(it),
                                  tracks.res_eLOC1_fit->at(it));
              ptHists["phi"]->Fill(tracks.t_pT->at(it),
                                   tracks.res_ePHI_fit->at(it));
              ptHists["theta"]->Fill(tracks.t_pT->at(it),
                                     tracks.res_eTHETA_fit->at(it));
              ptHists["qop"]->Fill(tracks.t_pT->at(it),
                                   tracks.res_eQOP_fit->at(it));
              // tracks.res_eQOP_fit->at(it)*tracks.t_p->at(it));
              ptHists["t"]->Fill(tracks.t_pT->at(it),
                                 timeInns(tracks.res_eT_fit->at(it)));
            }
          } else {
            if (withinPt) {
              etaHists["l0"]->Fill(tracks.t_eta->at(it),
                                   tracks.pull_eLOC0_fit->at(it));
              etaHists["l1"]->Fill(tracks.t_eta->at(it),
                                   tracks.pull_eLOC1_fit->at(it));
              etaHists["phi"]->Fill(tracks.t_eta->at(it),
                                    tracks.pull_ePHI_fit->at(it));
              etaHists["theta"]->Fill(tracks.t_eta->at(it),
                                      tracks.pull_eTHETA_fit->at(it));
              etaHists["qop"]->Fill(tracks.t_eta->at(it),
                                    tracks.pull_eQOP_fit->at(it));
              etaHists["t"]->Fill(tracks.t_eta->at(it),
                                  tracks.pull_eT_fit->at(it));
            }

            // if (withinEta) {
            if (1) {
              ptHists["l0"]->Fill(tracks.t_pT->at(it),
                                  tracks.pull_eLOC0_fit->at(it));
              ptHists["l1"]->Fill(tracks.t_pT->at(it),
                                  tracks.pull_eLOC1_fit->at(it));
              ptHists["phi"]->Fill(tracks.t_pT->at(it),
                                   tracks.pull_ePHI_fit->at(it));
              ptHists["theta"]->Fill(tracks.t_pT->at(it),
                                     tracks.pull_eTHETA_fit->at(it));
              ptHists["qop"]->Fill(tracks.t_pT->at(it),
                                   tracks.pull_eQOP_fit->at(it));
              ptHists["t"]->Fill(tracks.t_pT->at(it),
                                 tracks.pull_eT_fit->at(it));
            }
          }
        }
      }
    }
  };

  fillHists(entries1, tracks1, ref1_vs_eta, ref1_vs_pt);
  fillHists(entries2, tracks2, val1_vs_eta, val1_vs_pt);
  fillHists(entries3, tracks3, ref2_vs_eta, ref2_vs_pt);
  fillHists(entries4, tracks4, val2_vs_eta, val2_vs_pt);

  std::cout << "Done " << std::endl;

  // ana the hists
  for (const auto& name : names) {
    // mean and width vs. eta
    for (int j = 1; j <= ref1mean_vs_eta["l0"]->GetNbinsX(); ++j) {
      TH1D* temp = ref1_vs_eta[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_eta_Histo", j), j, j);
      anaHisto(temp, j, ref1mean_vs_eta.at(name), ref1width_vs_eta.at(name),
               fit);
    }
    for (int j = 1; j <= val1mean_vs_eta["l0"]->GetNbinsX(); ++j) {
      TH1D* temp = val1_vs_eta[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_eta_Histo", j), j, j);
      anaHisto(temp, j, val1mean_vs_eta.at(name), val1width_vs_eta.at(name),
               fit);
    }
    for (int j = 1; j <= ref2mean_vs_eta["l0"]->GetNbinsX(); ++j) {
      TH1D* temp = ref2_vs_eta[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_eta_Histo", j), j, j);
      anaHisto(temp, j, ref2mean_vs_eta.at(name), ref2width_vs_eta.at(name),
               fit);
    }
    for (int j = 1; j <= val2mean_vs_eta["l0"]->GetNbinsX(); ++j) {
      TH1D* temp = val2_vs_eta[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_eta_Histo", j), j, j);
      anaHisto(temp, j, val2mean_vs_eta.at(name), val2width_vs_eta.at(name),
               fit);
    }

    // mean and width vs. pt
    std::cout << "There are " << ref1mean_vs_pt["l0"]->GetNbinsX() << " pt bins"
              << std::endl;
    for (int j = 1; j <= ref1mean_vs_pt["l0"]->GetNbinsX(); ++j) {
      double scale = 1;
      if (name == "qop")
        scale = ps[j - 1];
      TH1D* temp = ref1_vs_pt[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_pt_Histo", j), j, j);
      anaHisto(temp, j, ref1mean_vs_pt.at(name), ref1width_vs_pt.at(name), fit,
               scale);
    }
    for (int j = 1; j <= val1mean_vs_pt["l0"]->GetNbinsX(); ++j) {
      double scale = 1;
      if (name == "qop")
        scale = ps_[j - 1];
      TH1D* temp = val1_vs_pt[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_pt_Histo", j), j, j);
      anaHisto(temp, j, val1mean_vs_pt.at(name), val1width_vs_pt.at(name), fit,
               scale);
    }
    for (int j = 1; j <= ref2mean_vs_pt["l0"]->GetNbinsX(); ++j) {
      double scale = 1;
      if (name == "qop")
        scale = ps[j - 1];
      TH1D* temp = ref2_vs_pt[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_pt_Histo", j), j, j);
      anaHisto(temp, j, ref2mean_vs_pt.at(name), ref2width_vs_pt.at(name), fit,
               scale);
    }
    for (int j = 1; j <= val2mean_vs_pt["l0"]->GetNbinsX(); ++j) {
      double scale = 1;
      if (name == "qop")
        scale = ps_[j - 1];
      TH1D* temp = val2_vs_pt[name]->ProjectionY(
          Form("%s_projy_bin%d", "Pull_vs_pt_Histo", j), j, j);
      anaHisto(temp, j, val2mean_vs_pt.at(name), val2width_vs_pt.at(name), fit,
               scale);
    }
  }

  std::cout << "Start plotting" << std::endl;

  // plotting
  std::vector<TCanvas*> cs_eta;
  std::vector<TCanvas*> cs_pt;
  std::vector<TLegend*> legs_eta;
  std::vector<TLegend*> legs_pt;
  std::vector<TLine*> lines_0_eta;
  std::vector<TLine*> lines_0_pt;
  std::vector<TLine*> lines_1_eta;
  std::vector<TLine*> lines_1_pt;

  std::vector<TPad*> tpads_eta;
  std::vector<TPad*> tpads_pt;
  std::vector<TPad*> bpads_eta;
  std::vector<TPad*> bpads_pt;

  for (const auto& name : names) {
    if (showEta) {
      if (plotType == 2) {
        cs_eta.push_back(
            new TCanvas(Form("c_%s_vs_eta", name.c_str()), "", 500, 450));
        tpads_eta.push_back(new TPad(Form("tpad_%s_vs_eta", name.c_str()), "",
                                     0, 0.46, 1, 1.0));
        bpads_eta.push_back(new TPad(Form("bpad_%s_vs_eta", name.c_str()), "",
                                     0, 0.05, 1, 0.45));
      } else {
        cs_eta.push_back(
            new TCanvas(Form("c_%s_vs_eta", name.c_str()), "", 500, 450));
      }
    }
    if (showPt) {
      if (plotType == 2) {
        cs_pt.push_back(
            new TCanvas(Form("c_%s_vs_pt", name.c_str()), "", 500, 450));
        tpads_pt.push_back(
            new TPad(Form("tpad_%s_vs_pt", name.c_str()), "", 0, 0.46, 1, 1.0));
        bpads_pt.push_back(new TPad(Form("bpad_%s_vs_pt", name.c_str()), "", 0,
                                    0.05, 1, 0.45));
      } else {
        cs_pt.push_back(
            new TCanvas(Form("c_%s_vs_pt", name.c_str()), "", 500, 450));
      }
    }

    if (plotType == 2) {
      legs_eta.push_back(new TLegend(0.6, 0.6, 0.9, 0.85));
      legs_pt.push_back(new TLegend(0.6, 0.6, 0.9, 0.85));
    } else {
      legs_eta.push_back(new TLegend(0.6, 0.6, 0.9, 0.85));
      legs_pt.push_back(new TLegend(0.6, 0.6, 0.9, 0.85));
    }

    lines_1_eta.push_back(
        new TLine(ref1mean_vs_eta["l0"]->GetXaxis()->GetXmin(), 1.,
                  ref1mean_vs_eta["l0"]->GetXaxis()->GetXmax(), 1.));
    lines_1_pt.push_back(
        new TLine(ref1mean_vs_pt["l0"]->GetXaxis()->GetXmin(), 1.,
                  ref1mean_vs_pt["l0"]->GetXaxis()->GetXmax(), 1.));
    lines_0_eta.push_back(
        new TLine(ref1mean_vs_eta["l0"]->GetXaxis()->GetXmin(), 0.,
                  ref1mean_vs_eta["l0"]->GetXaxis()->GetXmax(), 0.));
    lines_0_pt.push_back(
        new TLine(ref1mean_vs_pt["l0"]->GetXaxis()->GetXmin(), 0.,
                  ref1mean_vs_pt["l0"]->GetXaxis()->GetXmax(), 0.));
  }

  for (int i = 0; i < names.size(); i++) {
    legs_eta[i]->SetLineStyle(0);
    legs_eta[i]->SetBorderSize(0);
    legs_eta[i]->SetFillStyle(0);
    lines_0_eta[i]->SetLineWidth(2);
    lines_0_eta[i]->SetLineStyle(kDashed);
    lines_1_eta[i]->SetLineWidth(2);
    lines_1_eta[i]->SetLineStyle(kDashed);
    legs_pt[i]->SetLineStyle(0);
    legs_pt[i]->SetBorderSize(0);
    legs_pt[i]->SetFillStyle(0);
    lines_0_pt[i]->SetLineWidth(2);
    lines_0_pt[i]->SetLineStyle(kDashed);
    lines_1_pt[i]->SetLineWidth(2);
    lines_1_pt[i]->SetLineStyle(kDashed);

    if (not tpads_eta.empty()) {
      tpads_eta[i]->SetRightMargin(0.05);
      tpads_eta[i]->SetBottomMargin(0.005);
      bpads_eta[i]->SetTopMargin(0.005);
      bpads_eta[i]->SetRightMargin(0.05);
      bpads_eta[i]->SetBottomMargin(0.3);
    }

    if (not tpads_pt.empty()) {
      tpads_pt[i]->SetRightMargin(0.05);
      tpads_pt[i]->SetBottomMargin(0.005);
      bpads_pt[i]->SetTopMargin(0.005);
      bpads_pt[i]->SetRightMargin(0.05);
      bpads_pt[i]->SetBottomMargin(0.3);
    }
  }

  std::cout << "will do plotHist " << std::endl;

  // Labmda for plotting
  auto plotHists = [&](std::map<std::string, TH1F*>& ref1Hists,
                       std::map<std::string, TH1F*>& val1Hists,
                       std::map<std::string, TH1F*>& ref2Hists,
                       std::map<std::string, TH1F*>& val2Hists, bool plotEta,
                       bool plotMean = 0, int labelSize = 0,
                       double scaleYMin = 1., double scaleYMax = 1,
                       bool addEntry = true) {
    for (int i = 0; i < names.size(); i++) {
      if (plotEta) {
        ref1Hists[names[i]]->GetXaxis()->SetTitle("Truth #eta");
        val1Hists[names[i]]->GetXaxis()->SetTitle("Truth #eta");
        ref2Hists[names[i]]->GetXaxis()->SetTitle("Truth #eta");
        val2Hists[names[i]]->GetXaxis()->SetTitle("Truth #eta");
      } else {
        ref1Hists[names[i]]->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        val1Hists[names[i]]->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        ref2Hists[names[i]]->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
        val2Hists[names[i]]->GetXaxis()->SetTitle("Truth p_{T} [GeV]");
      }

      std::string resTitle;
      if (plotResidualRatio) {
        resTitle = resRatioTitles[i];
      } else {
        resTitle = resTitles[i];
      }
      std::string leg = plotResidual ? resTitle : pullTitles[i];
      if (fit) {
        if (not plotMean) {
          ref1Hists[names[i]]->GetYaxis()->SetTitle(
              Form("#sigma%s", leg.c_str()));
          val1Hists[names[i]]->GetYaxis()->SetTitle(
              Form("#sigma%s", leg.c_str()));
          ref2Hists[names[i]]->GetYaxis()->SetTitle(
              Form("#sigma%s", leg.c_str()));
          val2Hists[names[i]]->GetYaxis()->SetTitle(
              Form("#sigma%s", leg.c_str()));
        } else if (plotMean) {
          ref1Hists[names[i]]->GetYaxis()->SetTitle(Form("#mu%s", leg.c_str()));
          val1Hists[names[i]]->GetYaxis()->SetTitle(Form("#mu%s", leg.c_str()));
          ref2Hists[names[i]]->GetYaxis()->SetTitle(Form("#mu%s", leg.c_str()));
          val2Hists[names[i]]->GetYaxis()->SetTitle(Form("#mu%s", leg.c_str()));
        }
      } else {
        if (not plotMean) {
          ref1Hists[names[i]]->GetYaxis()->SetTitle(
              Form("RMS of %s", leg.c_str()));
          val1Hists[names[i]]->GetYaxis()->SetTitle(
              Form("RMS of %s", leg.c_str()));
          ref2Hists[names[i]]->GetYaxis()->SetTitle(
              Form("RMS of %s", leg.c_str()));
          val2Hists[names[i]]->GetYaxis()->SetTitle(
              Form("RMS of %s", leg.c_str()));
        } else if (plotMean) {
          ref1Hists[names[i]]->GetYaxis()->SetTitle(
              Form("Mean of %s", leg.c_str()));
          val1Hists[names[i]]->GetYaxis()->SetTitle(
              Form("Mean of %s", leg.c_str()));
          ref2Hists[names[i]]->GetYaxis()->SetTitle(
              Form("Mean of %s", leg.c_str()));
          val2Hists[names[i]]->GetYaxis()->SetTitle(
              Form("Mean of %s", leg.c_str()));
        }
      }

      if (plotResidual) {
        // residual
        if (plotEta) {
          std::pair<double, double> range = (not plotMean)
                                                ? yRange_reswidth_vs_eta[i]
                                                : yRange_resmean_vs_eta[i];
          ref1Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
          ref2Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
        } else {
          std::pair<double, double> range = (not plotMean)
                                                ? yRange_reswidth_vs_pt[i]
                                                : yRange_resmean_vs_pt[i];
          ref1Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
          ref2Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
        }
      } else {
        // pull
        if (plotEta) {
          std::pair<double, double> range = (not plotMean)
                                                ? yRange_pullwidth_vs_eta[i]
                                                : yRange_pullmean_vs_eta[i];
          ref1Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
          ref2Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
        } else {
          std::pair<double, double> range = (not plotMean)
                                                ? yRange_pullwidth_vs_pt[i]
                                                : yRange_pullmean_vs_pt[i];
          ref1Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
          ref2Hists[names[i]]->GetYaxis()->SetRangeUser(
              range.first * scaleYMin, range.second * scaleYMax);
        }
      }

      if (labelSize == 0) {
        setThisHistStyle(ref1Hists[names[i]], colors[0], markers[0], xTitleSize,
                         yTitleSize, xLabelSize, yLabelSize, xTitleOffset,
                         yTitleOffset);
        setThisHistStyle(val1Hists[names[i]], colors[1], markers[1], xTitleSize,
                         yTitleSize, xLabelSize, yLabelSize, xTitleOffset,
                         yTitleOffset);
        setThisHistStyle(ref2Hists[names[i]], colors[2], markers[2], xTitleSize,
                         yTitleSize, xLabelSize, yLabelSize, xTitleOffset,
                         yTitleOffset);
        setThisHistStyle(val2Hists[names[i]], colors[3], markers[3], xTitleSize,
                         yTitleSize, xLabelSize, yLabelSize, xTitleOffset,
                         yTitleOffset);

      } else if (labelSize == 1) {
        setThisHistStyle(ref1Hists[names[i]], colors[0], markers[0],
                         topXTitleSize, topYTitleSize, topXLabelSize,
                         topYLabelSize, topXTitleOffset, topYTitleOffset, 505);
        setThisHistStyle(val1Hists[names[i]], colors[1], markers[1],
                         topXTitleSize, topYTitleSize, topXLabelSize,
                         topYLabelSize, topXTitleOffset, topYTitleOffset, 505);
        setThisHistStyle(ref2Hists[names[i]], colors[2], markers[2],
                         topXTitleSize, topYTitleSize, topXLabelSize,
                         topYLabelSize, topXTitleOffset, topYTitleOffset, 505);
        setThisHistStyle(val2Hists[names[i]], colors[3], markers[3],
                         topXTitleSize, topYTitleSize, topXLabelSize,
                         topYLabelSize, topXTitleOffset, topYTitleOffset, 505);

      } else if (labelSize == 2) {
        setThisHistStyle(ref1Hists[names[i]], colors[0], markers[0],
                         botXTitleSize, botYTitleSize, botXLabelSize,
                         botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
        setThisHistStyle(val1Hists[names[i]], colors[1], markers[1],
                         botXTitleSize, botYTitleSize, botXLabelSize,
                         botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
        setThisHistStyle(ref2Hists[names[i]], colors[2], markers[2],
                         botXTitleSize, botYTitleSize, botXLabelSize,
                         botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
        setThisHistStyle(val2Hists[names[i]], colors[3], markers[3],
                         botXTitleSize, botYTitleSize, botXLabelSize,
                         botYLabelSize, botXTitleOffset, botYTitleOffset, 505);
      }

      if (plotEta) {
        cs_eta[i]->cd();
        if (plotType == 2 and plotMean) {
          tpads_eta[i]->Draw();
          tpads_eta[i]->cd();
        }
        if (plotType == 2 and not plotMean) {
          bpads_eta[i]->Draw();
          bpads_eta[i]->cd();
        }
      } else {
        cs_pt[i]->cd();
        if (plotType == 2 and plotMean) {
          tpads_pt[i]->Draw();
          tpads_pt[i]->cd();
        }
        if (plotType == 2 and not plotMean) {
          bpads_pt[i]->Draw();
          bpads_pt[i]->cd();
        }
      }

      if (plotResidual and not plotMean and plotResidualRatio) {
        val1Hists[names[i]]->Divide(ref1Hists[names[i]]);
        val2Hists[names[i]]->Divide(ref2Hists[names[i]]);
        val1Hists[names[i]]->GetYaxis()->SetRangeUser(0.1, 1.6);
        val2Hists[names[i]]->GetYaxis()->SetRangeUser(0.1, 1.6);
        val1Hists[names[i]]->Draw("esame");
        val2Hists[names[i]]->Draw("esame");
        val1Hists[names[i]]->SetMarkerColor(1);
        val2Hists[names[i]]->SetMarkerColor(1);
        val1Hists[names[i]]->SetLineColor(1);
        val2Hists[names[i]]->SetLineColor(1);
        val1Hists[names[i]]->SetMarkerStyle(24);
        val2Hists[names[i]]->SetMarkerStyle(20);
        val1Hists[names[i]]->GetYaxis()->SetTitle(
            Form("Ratio of RMS of %s", leg.c_str()));
        val2Hists[names[i]]->GetYaxis()->SetTitle(
            Form("Ratio of RMS of %s", leg.c_str()));
      } else {
        ref1Hists[names[i]]->Draw("e");
        ref2Hists[names[i]]->Draw("esame");
        val1Hists[names[i]]->Draw("esame");
        val2Hists[names[i]]->Draw("esame");
      }

      if (addEntry) {
        if (plotEta) {
          if (plotResidualRatio) {
            legs_eta[i]->AddEntry(val1Hists[names[i]], "val1 tag", "lep");
            legs_eta[i]->AddEntry(val2Hists[names[i]], "val2 tag", "lep");
          } else {
            legs_eta[i]->AddEntry(ref2Hists[names[i]], legs[2].c_str(), "lep");
            legs_eta[i]->AddEntry(val2Hists[names[i]], legs[3].c_str(), "lep");
          }
          legs_eta[i]->Draw();

        } else {
          if (plotResidualRatio) {
            legs_pt[i]->AddEntry(val1Hists[names[i]], "val1 tag", "lep");
            legs_pt[i]->AddEntry(val2Hists[names[i]], "val2 tag", "lep");
          } else {
            legs_pt[i]->AddEntry(ref1Hists[names[i]], legs[0].c_str(), "lep");
            legs_pt[i]->AddEntry(val1Hists[names[i]], legs[1].c_str(), "lep");
            legs_pt[i]->AddEntry(ref2Hists[names[i]], legs[2].c_str(), "lep");
            legs_pt[i]->AddEntry(val2Hists[names[i]], legs[3].c_str(), "lep");
          }
          legs_pt[i]->Draw();
        }
      }

      std::cout << "draw lines" << std::endl;
      if (plotMean) {
        if (plotEta)
          lines_0_eta[i]->Draw();
        else
          lines_0_pt[i]->Draw();
      } else if (not plotResidual) {
        if (plotEta)
          lines_1_eta[i]->Draw();
        else
          lines_1_pt[i]->Draw();
      }
      std::cout << "draw lines" << std::endl;

      // labelSize==1 means plotting the top panel
      if (saveAs and not(plotType == 2 and labelSize == 1)) {
        std::string meritType;
        if (plotType == 0) {
          meritType = fit ? "mu" : "Mean";
        } else if (plotType == 1) {
          meritType = fit ? "width" : "RMS";
        } else if (plotType == 2) {
          meritType = fit ? "muwidth" : "MeanRMS";
        }

        if (plotEta) {
          cs_eta[i]->Update();
          cs_eta[i]->SaveAs(Form("%s/%s_%s%s_vs_eta_%s.pdf", path.c_str(),
                                 detector.c_str(), meritType.c_str(),
                                 ratioTag.c_str(), names[i].c_str()));
        } else {
          cs_pt[i]->Update();
          cs_pt[i]->SaveAs(Form("%s/%s_%s%s_vs_pt_%s.pdf", path.c_str(),
                                detector.c_str(), meritType.c_str(),
                                ratioTag.c_str(), names[i].c_str()));
        }
      }

    }  // end of all variables
  };

  double scaleYMin_res_eta = 1.9;
  double scaleYMax_res_eta = 0.5;
  double scaleYMin_res_pt = 1.5;
  double scaleYMax_res_pt = 0.6;

  double scaleYMin_pull_eta = 1.;
  double scaleYMax_pull_eta = 1.;
  double scaleYMin_pull_pt = 1.;
  double scaleYMax_pull_pt = 1.;

  double scaleYMin_eta = plotResidual ? scaleYMin_res_eta : scaleYMin_pull_eta;
  double scaleYMax_eta = plotResidual ? scaleYMax_res_eta : scaleYMax_pull_eta;
  double scaleYMin_pt = plotResidual ? scaleYMin_res_pt : scaleYMin_pull_pt;
  double scaleYMax_pt = plotResidual ? scaleYMax_res_pt : scaleYMax_pull_pt;

  // parameters of plotHists in addition to the hists
  // bool plotEta, bool plotMean = true/false, int labelSize = 0/1/2, double
  // scaleYMin = 1., double scaleYMax = 1, bool addEntry = true
  if (plotType == 0) {
    if (showEta)
      plotHists(ref1mean_vs_eta, val1mean_vs_eta, ref2mean_vs_eta,
                val2mean_vs_eta, true, true);
    if (showPt)
      plotHists(ref1mean_vs_pt, val1mean_vs_pt, ref2mean_vs_pt, val2mean_vs_pt,
                false, true);
  } else if (plotType == 1) {
    if (showEta)
      plotHists(ref1width_vs_eta, val1width_vs_eta, ref2width_vs_eta,
                val2width_vs_eta, true, false);
    if (showPt)
      plotHists(ref1width_vs_pt, val1width_vs_pt, ref2width_vs_pt,
                val2width_vs_pt, false, false);
  } else {
    if (showEta) {
      plotHists(ref1mean_vs_eta, val1mean_vs_eta, ref2mean_vs_eta,
                val2mean_vs_eta, true, true, 1, 0.9, 1.1);
      plotHists(ref1width_vs_eta, val1width_vs_eta, ref2width_vs_eta,
                val2width_vs_eta, true, false, 2, scaleYMin_eta, scaleYMax_eta,
                false);
    }
    if (showPt) {
      plotHists(ref1mean_vs_pt, val1mean_vs_pt, ref2mean_vs_pt, val2mean_vs_pt,
                false, true, 1, 0.9, 1.1);
      plotHists(ref1width_vs_pt, val1width_vs_pt, ref2width_vs_pt,
                val2width_vs_pt, false, false, 2, scaleYMin_pt, scaleYMax_pt,
                false);
    }
  }

  /*
  // pull width vs eta
  for(int i=0; i<names.size(); i++){
  ref1width_vs_eta[names[i]]->GetXaxis()->SetTitle("Truth #eta");

  //std::string leg = (names[i]=="phi" or names[i] == "theta")? "#"+names[i]:
  names[i]; std::string leg = plotResidual? resTitles[i]:pullTitles[i]; if(fit){
  ref1width_vs_eta[names[i]]->GetYaxis()->SetTitle(Form("Fitted #sigma of %s",
  leg.c_str())); }else{
  ref1width_vs_eta[names[i]]->GetYaxis()->SetTitle(Form("RMS of %s",
  leg.c_str()));
  }

  if(plotResidual){
  ref1width_vs_eta[names[i]]->GetYaxis()->SetRangeUser(yRange_reswidth_vs_eta[i].first,
  yRange_reswidth_vs_eta[i].second); } else {
  ref1width_vs_eta[names[i]]->GetYaxis()->SetRangeUser(yRange_pullwidth_vs_eta[i].first,
  yRange_pullwidth_vs_eta[i].second);
  }
  //setYRange(ref1width_vs_eta[names[i]], 0.2, 2.5);
  //setYRange(val1width_vs_eta[names[i]], 0.2, 2.5);


  setThisHistStyle(ref1width_vs_eta[names[i]], 64, 24, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
  setThisHistStyle(val1width_vs_eta[names[i]], 95, 26, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
  setThisHistStyle(ref2width_vs_eta[names[i]], 64, 20, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
  setThisHistStyle(val2width_vs_eta[names[i]], 95, 22, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);

  cs_eta[i]->cd();
  ref1width_vs_eta[names[i]]->Draw("e");
  val1width_vs_eta[names[i]]->Draw("esame");
  ref2width_vs_eta[names[i]]->Draw("esame");
  val2width_vs_eta[names[i]]->Draw("esame");

  legs_eta[i]->AddEntry(ref1width_vs_eta[names[i]], legs[0].c_str(), "lep");
  legs_eta[i]->AddEntry(ref2width_vs_eta[names[i]], legs[2].c_str(), "lep");
  legs_eta[i]->AddEntry(val1width_vs_eta[names[i]], legs[1].c_str(), "lep");
  legs_eta[i]->AddEntry(val2width_vs_eta[names[i]], legs[3].c_str(), "lep");
  legs_eta[i]->SetLineStyle(0);
  legs_eta[i]->SetBorderSize(0);
  legs_eta[i]->SetFillStyle(0);
  legs_eta[i]->Draw();

  lines_eta[i]->SetLineWidth(2);
  lines_eta[i]->SetLineStyle(kDashed);
  if(not plotResidual){
  lines_eta[i]->Draw();
  }

  cs_eta[i]->Update();
  if(fit){
  cs_eta[i]->SaveAs(Form("%s%swidth_vs_eta_%s.pdf", path.c_str(),
  detector.c_str(), names[i].c_str())); } else {
  cs_eta[i]->SaveAs(Form("%s%s_RMS_vs_eta_%s.pdf", path.c_str(),
  detector.c_str(), names[i].c_str()));
  }
  }

  // pull width vs pt
  for(int i=0; i<names.size(); i++){
  ref1width_vs_pt[names[i]]->GetXaxis()->SetTitle("Truth p_{T} [GeV]");

  //std::string leg = (names[i]=="phi" or names[i] == "theta")? "#"+names[i]:
  names[i]; std::string leg = plotResidual? resTitles[i]:pullTitles[i]; if(fit){
  ref1width_vs_pt[names[i]]->GetYaxis()->SetTitle(Form("Fitted #sigma of %s",
  leg.c_str())); }else{
  ref1width_vs_pt[names[i]]->GetYaxis()->SetTitle(Form("RMS of %s",
  leg.c_str()));
  }

  if(plotResidual){
  ref1width_vs_pt[names[i]]->GetYaxis()->SetRangeUser(yRange_reswidth_vs_pt[i].first,
  yRange_reswidth_vs_pt[i].second); } else {
  ref1width_vs_pt[names[i]]->GetYaxis()->SetRangeUser(yRange_pullwidth_vs_pt[i].first,
  yRange_pullwidth_vs_pt[i].second);
  }
  //setYRange(ref1width_vs_pt[names[i]], 0.2, 2.5);
  //setYRange(val1width_vs_pt[names[i]], 0.2, 2.5);


  setThisHistStyle(ref1width_vs_pt[names[i]], 64, 24, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
  setThisHistStyle(val1width_vs_pt[names[i]], 95, 26, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
  setThisHistStyle(ref2width_vs_pt[names[i]], 64, 20, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
  setThisHistStyle(val2width_vs_pt[names[i]], 95, 22, xTitleSize, yTitleSize,
  xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);

  cs_pt[i]->cd();
  ref1width_vs_pt[names[i]]->Draw("e");
  val1width_vs_pt[names[i]]->Draw("esame");
  ref2width_vs_pt[names[i]]->Draw("esame");
  val2width_vs_pt[names[i]]->Draw("esame");
  //gPad->Update();
  //gPad->SetLogx();

  legs_pt[i]->AddEntry(ref1width_vs_pt[names[i]], legs[0].c_str(), "lep");
  legs_pt[i]->AddEntry(ref2width_vs_pt[names[i]], legs[2].c_str(), "lep");
  legs_pt[i]->AddEntry(val1width_vs_pt[names[i]], legs[1].c_str(), "lep");
  legs_pt[i]->AddEntry(val2width_vs_pt[names[i]], legs[3].c_str(), "lep");
  legs_pt[i]->SetLineStyle(0);
  legs_pt[i]->SetBorderSize(0);
  legs_pt[i]->SetFillStyle(0);
  legs_pt[i]->Draw();

  lines_pt[i]->SetLineWidth(2);
  lines_pt[i]->SetLineStyle(kDashed);
  if(not plotResidual){
  lines_pt[i]->Draw();
  }

  cs_pt[i]->Update();
  if(fit){
  cs_pt[i]->SaveAs(Form("%s%swidth_vs_pt_%s.pdf", path.c_str(),
  detector.c_str(), names[i].c_str())); } else {
  cs_pt[i]->SaveAs(Form("%s%s_RMS_vs_pt_%s.pdf", path.c_str(), detector.c_str(),
  names[i].c_str()));
  }
  }
  */
}
