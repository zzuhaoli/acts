#pragma once

#include "TLatex.h"
#include "TStyle.h"

void setStyleAtlas(double tsize = 0.05) {
  Int_t icol = 0;
  gStyle->SetFrameBorderMode(icol);
  gStyle->SetFrameFillColor(icol);
  gStyle->SetCanvasBorderMode(icol);
  gStyle->SetCanvasColor(icol);
  gStyle->SetPadBorderMode(icol);
  gStyle->SetPadColor(icol);
  gStyle->SetStatColor(icol);
  // gStyle->SetFillColor(icol); // don't use: white fill color for *all*
  // objects

  // set the paper & margin sizes
  gStyle->SetPaperSize(20, 26);

  // set margin sizes
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  gStyle->SetTitleXOffset(1.4);
  gStyle->SetTitleYOffset(1.4);

  // use large fonts
  // Int_t font=72; // Helvetica italics
  Int_t font = 42;  // Helvetica
  //    Double_t tsize = 0.05;
  gStyle->SetTextFont(font);

  gStyle->SetTextSize(tsize);
  gStyle->SetLabelFont(font, "x");
  gStyle->SetTitleFont(font, "x");
  gStyle->SetLabelFont(font, "y");
  gStyle->SetTitleFont(font, "y");
  gStyle->SetLabelFont(font, "z");
  gStyle->SetTitleFont(font, "z");

  gStyle->SetLabelSize(tsize, "x");
  gStyle->SetTitleSize(tsize, "x");
  gStyle->SetLabelSize(tsize, "y");
  gStyle->SetTitleSize(tsize, "y");
  gStyle->SetLabelSize(tsize, "z");
  gStyle->SetTitleSize(tsize, "z");

  // use bold lines and markers
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1.2);
  gStyle->SetHistLineWidth((Width_t)2.);
  gStyle->SetLineStyleString(2, "[12 12]");  // postscript dashes

  // get rid of X error bars
  // gStyle->SetErrorX(0.001);
  // get rid of error bar caps
  // gStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  gStyle->SetOptTitle(0);
  // gStyle->SetOptStat(1111);
  gStyle->SetOptStat(0);
  // gStyle->SetOptFit(1111);
  gStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
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

void ATLAS_LABEL(Double_t x, Double_t y, Color_t color, Double_t tsize = 0.05) {
  // Double_t tsize=0.05;
  TLatex l;  //
  l.SetTextAlign(12);
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x, y, "ATLAS");
}

void ATLASLabel(Double_t x, Double_t y, Color_t color, Double_t tsize = 0.05,
                std::string text = "Internal", double delx = 0) {
  TLatex l;
  // l.SetTextAlign(12);
  l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextAlign(12);
  l.SetTextFont(72);
  //  l.SetTextFont(73);
  l.SetTextColor(color);
  l.DrawLatex(x, y, "ATLAS");
  if (delx == 0)
    delx = 0.115 * 696 * gPad->GetWh() / (472 * gPad->GetWw());
  // std::cout<<"delx = " << delx << std::endl;
  if (not text.empty()) {
    TLatex p;
    p.SetNDC();
    p.SetTextAlign(12);
    p.SetTextSize(tsize);
    p.SetTextFont(42);
    //    p.SetTextFont(43);
    p.SetTextColor(color);
    p.DrawLatex(x + delx, y, text.c_str());
  }
}

void ATLASLabelOld(Double_t x, Double_t y, bool Preliminary, Color_t color) {
  TLatex l;  // l.SetTextAlign(12); l.SetTextSize(tsize);
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x, y, "ATLAS");
  if (Preliminary) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x + 0.115, y, "Preliminary");
  }
}

void ATLASLabelNew(Double_t x, Double_t y, const char* text, Color_t color,
                   float text_size, float delx) {
  TLatex l;
  l.SetNDC();
  l.SetTextFont(73);
  l.SetTextColor(color);
  l.SetTextSize(text_size);

  l.DrawLatex(x, y, "ATLAS");
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextFont(43);
    p.SetTextColor(color);
    p.SetTextSize(text_size);
    p.DrawLatex(x + delx, y, text);
  }
  return;
}

void ATLASVersion(const char* version, Double_t x, Double_t y, Color_t color) {
  if (version) {
    char versionString[100];
    sprintf(versionString, "Version %s", version);
    TLatex l;
    l.SetTextAlign(22);
    l.SetTextSize(0.04);
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextColor(color);
    l.DrawLatex(x, y, versionString);
  }
}
