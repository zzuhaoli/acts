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
  //hist->SetLineColor(0);
  //hist->SetLineColorAlpha(0, 1); 
  //hist->SetLineStyle(0); 
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

void anaHisto(TH1F* inputHist, int j, TH1F* meanHist, TH1F* widthHist,
              bool fit = false, double scale=1) {
  // evaluate mean and width via the Gauss fit
  assert(inputHist != nullptr);
  if (inputHist->GetEntries() > 0) {
    if (fit) {
      //TFitResultPtr r = inputHist->Fit("gaus", "QS0");
      TFitResultPtr r = inputHist->Fit("gaus", "QS0");
      r = inputHist->Fit("gaus", "QS0");
      if (r.Get() and ((r->Status() % 1000) == 0)) {
      //if (r.Get()) {
        // fill the mean and width into 'j'th bin of the meanHist and widthHist,
        // respectively
        meanHist->SetBinContent(j, r->Parameter(1));
        meanHist->SetBinError(j, r->ParError(1));
        widthHist->SetBinContent(j, r->Parameter(2)*scale);
        widthHist->SetBinError(j, r->ParError(2));
        if(scale!=1) 
	  std::cout<<"hist "<<widthHist->GetName() <<" width = "<< r->Parameter(2)*scale << std::endl; 
        }else {
         std::cout<<"Fitting failed " << std::endl;	
	}
    } else {
      meanHist->SetBinContent(j, inputHist->GetMean());
      meanHist->SetBinError(j, inputHist->GetMeanError());
      // meanHist->SetBinError(j, inputHist->GetRMS());

      widthHist->SetBinContent(j, inputHist->GetRMS()*scale);
      widthHist->SetBinError(j, inputHist->GetRMSError());
    }
  }
}


void comparePerigee_v3(
       std::string inputPath = "/nfs/xiaocong/acts_workspace/scan/CKF.estimated.chi2Cut15.maxPropSteps330.NoStepAjustError.maxSeeds2",
       std::string particle = "mu-",
       std::vector<std::string> degs = {"90","60","30"},
       //std::vector<std::string> degs = {"90"},
       //in GeV
       bool usePt = true, 
       std::vector<double> ps = {0.1, 0.125, 0.15, 0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8},
       //std::vector<double> ps = {0.175, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8},
       bool savePlot  = false, 
    double absEtaMin = 0, double absEtaMax = 1.75, double ptMin = 0.05,
    double ptMax = 1.8, bool saveAs = false, bool showEta = false,
    bool showPt = true, bool fit = true, bool plotResidual = true,
    // plotType 0: mean,   1:width,   2: mean and width
    int plotType = 1, bool plotResidualRatio = false, bool absEta = true,
    bool variablePtBin = true, 
    std::map<std::string, std::string> tags = {
     {"90","|cos#theta|=0"},
     {"60","|cos#theta|=0.5"},
     {"30","|cos#theta|=0.87"},
    },
    std::vector<int> colors = {872, 866, 854, 896}, std::vector<int> markers ={24, 26, 20, 22}) {
 
  gStyle->SetOptFit(0000);
  gStyle->SetOptStat(0000);
  gStyle->SetPadLeftMargin(0.20);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetNdivisions(510, "y");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::string saveTag ="STCF_" + particle + "_res"; 

  if (plotType == 2 or not plotResidual)
    plotResidualRatio = false;

  if (plotType != 0 and plotType != 1 and plotType != 2) {
    throw std::invalid_argument("The plotType must be 0, 1 or 2!");
  }

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


  double xTitleSize = 0.05;
  double yTitleSize = 0.05;
  double xLabelSize = 0.05;
  double yLabelSize = 0.05;
  double xTitleOffset = 1.2;
  double yTitleOffset = plotResidual ? 1.4 : 1.7;

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

  std::pair<double, double> pullRange = {-5, 5};
  std::vector<std::pair<double, double>> resRanges = {
      {-4, 4},     {-4, 4},   {-0.02, 0.02},
      {-0.02, 0.02}, {-0.2, 0.2}, {-3.5, 3.5},
  };
  
  //For each deg, a few ps
  std::map<std::string, std::vector<TChain*>> chains;
  // For each deg, nParams hists 
  std::map<std::string, std::vector<TH1F*>> means;
  std::map<std::string, std::vector<TH1F*>> widths; 
  for(int i=0; i<degs.size(); ++i){
    for(int j=0; j<ps.size(); ++j){
      // Load the tree chain
      TChain* treeChain = new TChain("tracksummary");
      std::string psstring = std::to_string(static_cast<int>(ps[j]*1000));
      std::string ptUnit = (usePt)? "Pt":"";
      std::string inFile = inputPath + "/" + particle + "/absThetaDeg_" + degs[i] + "_" + degs[i] + "_momentum" + ptUnit + "Mev_" + psstring + "_" + psstring + "/tracksummary_ckf.root";
      std::cout<<"Reading file " << inFile << std::endl; 
      treeChain->Add(inFile.c_str()); 
  
      if (treeChain->GetEntries() == 0) {
        std::cout << "[x] No entries found ... " << std::endl;
        return -1;
       }
      
      chains[degs[i]].push_back(treeChain); 
    }
  }
   

  std::vector<std::string> names = {"l0", "l1", "phi", "theta", "qop", "t"};
  std::vector<std::string> stores = {"eLOC0", "eLOC1", "ePHI", "eTHETA", "eQOP", "eT"};

  for(int i=0; i<degs.size(); ++i){
    for(int j=0; j<names.size(); ++j){
      means[degs[i]].push_back(new TH1F(Form("%s_%s_mean", names[j].c_str(), degs[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
      widths[degs[i]].push_back(new TH1F(Form("%s_%s_width", names[j].c_str(), degs[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
    }
    means[degs[i]].push_back(new TH1F(Form("pt_%s_mean", degs[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
    widths[degs[i]].push_back(new TH1F(Form("pt_%s_width", degs[i].c_str()), "", ptbins.size()-1, ptbins.data() )); 
  }


  
  std::vector<std::pair<double, double>> yRange_resmean_vs_eta;
  std::vector<std::pair<double, double>> yRange_resmean_vs_pt;
  std::vector<std::pair<double, double>> yRange_reswidth_vs_eta;
  std::vector<std::pair<double, double>> yRange_reswidth_vs_pt;
  std::vector<std::pair<double, double>> yRange_pullmean_vs_eta;
  std::vector<std::pair<double, double>> yRange_pullmean_vs_pt;
  std::vector<std::pair<double, double>> yRange_pullwidth_vs_eta;
  std::vector<std::pair<double, double>> yRange_pullwidth_vs_pt;

  
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
      {0.0, 0.5},        {0.25, 0.6},    {0.0, 0.015},
      {0, 0.015}, {0, 0.02}, {0.95, 1.1},
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

  auto timeInns = [](double timeInmm) -> double {
    return timeInmm / 299792458000.0 * 1e+9;
  };

  // book hists
  for(const auto& [deg, cs] : chains){
    for (int i = 0;i < cs.size(); ++i) {
       auto& chain = cs[i];
       
       double degVal;
       if(deg=="90"){
         degVal= 90./180*M_PI;
       } else if(deg=="60"){
         degVal= 60./180*M_PI;
       } else if (deg=="30"){
         degVal= 30./180*M_PI;
       }
       double p_;
       if(usePt){
          p_ = ps[i]/sin(degVal);
       } else {
          p_ = ps[i];
       }
       
       std::string prefix = plotResidual?"res_":"pull_";

       for (int j = 0;j < names.size(); ++j) {
         const auto& name = names[j];

         // The hist y range
         std::pair<double, double> yRange = pullRange;
         if (plotResidual) {
           yRange = resRanges[j];
         }

	 if(name=="qop" and ps[i] < 0.13){
           yRange = {-0.5, 0.5};
	 }
	 //if(name=="qop"){
         //  std::cout<<yRange.first  << std::endl;
         //}


         std::string draw = prefix + stores[j] + "_fit";
         
         //TH1F* temp = new TH1F(Form("%s_ptbin%i_%s", deg.c_str(), i, name.c_str()), "", 100, yRange.first, yRange.second);
         //cs[i]->Draw(Form("%s>>temp(200,%f,%f)", draw.c_str(), yRange.first, yRange.second),"nMeasurements>=5&&nHoles<=1");
         cs[i]->Draw(Form("%s>>temp(200,%f,%f)", draw.c_str(), yRange.first, yRange.second),"nMeasurements>=5");
         TH1F* temp = (TH1F*)gDirectory->Get("temp"); 

	 double scale = (name!="qop")? 1 : p_*100; 
         anaHisto(temp, i+1, means[deg][j], widths[deg][j], fit, scale);
      }

      {
        std::string draw_ = "abs(1./eQOP_fit)*sin(eTHETA_fit) - t_p*sin(t_theta)";
        std::pair<double, double> yRange_ = {-0.01, 0.01};
        if(ps[i]>0.13) yRange_ = {-0.03, 0.03};	
         std::cout<<yRange_.first  << std::endl;
        cs[i]->Draw(Form("%s>>temp_(100,%f,%f)", draw_.c_str(), yRange_.first, yRange_.second),"nMeasurements>=5");
        TH1F* temp_ = (TH1F*)gDirectory->Get("temp_");
	double scale_ = (usePt)? 100./ps[i] : 100./(ps[i]*sin(degVal));
	anaHisto(temp_, i+1, means[deg][6], widths[deg][6], fit, scale_);
      }
    }
  }


  std::vector<TLegend*> legs;
  for(int i=0; i<names.size()+1; ++i){
   legs.push_back(new TLegend(0.3, 0.7, 0.6, 0.85));
   legs[i]->SetLineStyle(0);
   legs[i]->SetBorderSize(0);
   //legs[i]->SetFillStyle(0);
  }

  std::string xAxisTitle = (usePt)? "p_{T} [GeV]" :  "p [GeV]";
  std::string plotName = (usePt)? "pt":"p";

  TCanvas *c1 = new TCanvas("c1", "", 600, 500); 
  c1->SetGrid();  
  int i=0; 
  for(const auto& [deg, hists] : widths){
    hists[0]->Draw("EsameX0");
    setThisHistStyle(hists[0], colors[i], markers[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[0]->GetXaxis()->SetTitle(xAxisTitle.c_str()); 
    hists[0]->GetYaxis()->SetTitle("#sigma(d_{0}) [mm]");
    hists[0]->GetYaxis()->SetRangeUser(0, 1.4); 
    legs[0]->AddEntry(hists[0],tags[deg].c_str(), "APL"); 
    i++; 
  }
  legs[0]->Draw();
  if(savePlot) 
  c1->SaveAs(Form("%s_d0_vs_%s.pdf", saveTag.c_str(), plotName.c_str()));

  i=0; 
  TCanvas *c2 = new TCanvas("c2", "", 600, 500); 
  c2->SetGrid();  
  for(const auto& [deg, hists] : widths){
    hists[1]->Draw("EsameX0");
    setThisHistStyle(hists[1], colors[i], markers[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[1]->GetXaxis()->SetTitle(xAxisTitle.c_str()); 
    hists[1]->GetYaxis()->SetTitle("#sigma(z_{0}) [mm]");
    hists[1]->GetYaxis()->SetRangeUser(0.2, 2.0); 
    legs[1]->AddEntry(hists[1],tags[deg].c_str(), "APL"); 
    i++; 
  }
  legs[1]->Draw();
  if(savePlot) 
  c2->SaveAs(Form("%s_z0_vs_%s.pdf", saveTag.c_str(), plotName.c_str()));

  i=0;
  TCanvas *c3 = new TCanvas("c3", "", 600, 500);
  c3->SetGrid();  
  for(const auto& [deg, hists] : widths){
    hists[4]->Draw("EsameX0");
    setThisHistStyle(hists[4], colors[i], markers[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[4]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[4]->GetYaxis()->SetTitle("#sigma(p)/p [%]");
    hists[4]->GetYaxis()->SetRangeUser(0.2, 1.6);
    legs[4]->AddEntry(hists[4],tags[deg].c_str(), "APL");
    i++;
  }
  legs[4]->Draw();
  if(savePlot) 
  c3->SaveAs(Form("%s_p_vs_%s.pdf", saveTag.c_str(), plotName.c_str()));

  i=0;
  TCanvas *c4 = new TCanvas("c4", "", 600, 500);
  c4->SetGrid();
  for(const auto& [deg, hists] : widths){
    hists[6]->Draw("EsameX0");
    setThisHistStyle(hists[6], colors[i], markers[i], xTitleSize, yTitleSize,
                         xLabelSize, yLabelSize, xTitleOffset, yTitleOffset);
    hists[6]->GetXaxis()->SetTitle(xAxisTitle.c_str());
    hists[6]->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T} [%]");
    hists[6]->GetYaxis()->SetRangeUser(0.2, 1.6);
    legs[6]->AddEntry(hists[4],tags[deg].c_str(), "APL");
    i++;
  }
  legs[6]->Draw();
  if(savePlot) 
  c4->SaveAs(Form("%s_pt_vs_%s.pdf", saveTag.c_str(), plotName.c_str()));


}
