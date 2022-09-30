#include "DD4hep/DetFactoryHelper.h"
#include "RecGeoData.h"
#include "TMath.h"
#include "XML/Layering.h"
#include <exception>
#include <vector>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;
using namespace TMath;

double rDefocus(double beta, double gamma, double rmin, double b) {
  double r = (rmin / Sin(beta) * Cos(gamma) - (b - 28) / 2.0) *
             Sin(M_PI / 2 - beta - gamma);
  return r;
}

double RDefocus(double beta, double gamma, double rmin, double b,
                double alpha) {
  double r = rDefocus(beta, gamma, rmin, b);
  double R = r * Cos(alpha);
  return R;
}

double rmax(double r, double h, double alpha) {
  double R = r * Cos(alpha);
  double a = r * Sin(alpha);
  double b = (R + h) / Cos(alpha);
  double max = Sqrt(a * a + b * b + 2 * a * b * Sin(alpha));
  return max;
}

Trap *MyTrap(double theta0, double dtheta, double dphi, double b, double rmin,
             double r, double R, double offsetw, double cthickness,
             bool iswrapping, bool isLeft) {
  double a, c, h1, h2, alpha1, alpha2, h3, h4, mtheta, mphi;
  a = rmin / Sin(theta0) * Sin(dtheta) *
      ((rmin / Sin(theta0) * Cos(dtheta) - (b - 28) / 2.0) /
       (rmin / Sin(theta0) * Cos(dtheta)));
  c = a + b * Tan(dtheta);
  h1 = ((rmin / Sin(theta0) *
         ((rmin / Sin(theta0) * Cos(dtheta) - (b - 28) / 2.0) /
          (rmin / Sin(theta0) * Cos(dtheta)))) *
            Cos(theta0) -
        r + R) *
       Tan(dphi / 2);
  h2 = ((rmin / Sin(theta0) * Cos(dtheta) - (b - 28) / 2.0) *
            Sin(M_PI / 2 - theta0 - dtheta) -
        r + R) *
       Tan(dphi / 2);
  alpha1 = 0;
  alpha2 = 0;
  h3 =
      ((rmin / Tan(theta0) + (b - (b - 28) / 2.0) / Cos(dtheta) * Cos(theta0)) -
       r + R) *
      Tan(dphi / 2);
  h4 = ((rmin / Sin(theta0) * Cos(dtheta) * Sin(M_PI / 2 - theta0 - dtheta) +
         (b - (b - 28) / 2.0) * Sin(M_PI / 2 - theta0 - dtheta)) -
        r + R) *
       Tan(dphi / 2);
  mtheta = ATan((c - a) / b / 2);
  mphi = M_PI / 2;

  if (offsetw > 0.0) {
    double tan1 = (h1 - h2) / a;
    double tan2 = (h3 - h4) / c;
    h2 = h2 + offsetw * tan1 - offsetw;
    h1 = h1 - offsetw * tan1 - offsetw;
    h4 = h4 + offsetw * tan2 - offsetw;
    h3 = h3 - offsetw * tan2 - offsetw;
    a = a - 2 * offsetw;
    c = a + b * Tan(dtheta);
  }

  if (iswrapping) {
    a = a + (offsetw - cthickness) * (c - a) / b;
    h1 = h1 +
         (offsetw - cthickness) / Cos(dtheta) * (h3 - h1) / (b / Cos(dtheta));
    h2 = h2 + (offsetw - cthickness) * (h4 - h2) / b;
    b = b - offsetw + cthickness;
  }

  Trap *MyTrap1 = new Trap(b / 2, mtheta, mphi, a / 2, h2, h1, alpha1, c / 2,
                           h4, h3, alpha2);
  Trap *MyTrap2 = new Trap(b / 2, -mtheta, mphi, c / 2, h4, h3, alpha1, a / 2,
                           h2, h1, alpha2);
  if (isLeft)
    return MyTrap2;
  else
    return MyTrap1;
}

static Ref_t create_element(Detector &lcdd, xml_h e, SensitiveDetector sens) {
  // Set generalized detector infomation based on XML file
  xml_det_t x_det = e;
  Layering layering(e);
  string det_name = x_det.nameStr();
  xml_dim_t dim = x_det.dimensions();
  int det_id = x_det.id();
  double rmin = dim.rmin();
  double zmax = dim.z();
  int nLadders1 = dim.inner_r();
  int nLadders2 = dim.inner_z();
  int nLadders3 = dim.outer_r();
  int nRow = dim.number();
  double dphi1 = 2 * Pi() / nLadders1;
  double dphi2 = 2 * Pi() / nLadders2;
  double dphi3 = 2 * Pi() / nLadders3;

  DetElement ECALEndcap(det_name, x_det.id());
  Material air = lcdd.air();
  Volume motherVol = lcdd.pickMotherVolume(ECALEndcap);
  sens.setType("calorimeter");

  double a, b, c, dtheta, theta0, mthickness;

  int module_id;

  rec::ECALERecGeoData *GeoData = new rec::ECALERecGeoData;

  // Set Angle for Crystal
  double d = 4.5;
  double q1 = (1.0 - d / 2.0 / rmin) / (1.0 + d / 2.0 / rmin);
  double Dtheta[22] = {0};
  double Theta0[22] = {0};
  double r1, r2, h1, h2, max1, max2;

  for (int i = 0; i < 15; i++) {
    if (i == 0)
      Dtheta[i] = M_PI / 2.0 - 2.0 * ATan(Power(q1, i + 1));
    Dtheta[i] = 2.0 * ATan(Power(q1, i)) - 2.0 * ATan(Power(q1, i + 1));
    Theta0[i] = 2.0 * ATan(Power(q1, i + 1));
  }
  r1 = rDefocus(Theta0[12], Dtheta[12], rmin, 28 + 1);
  h1 = rmin / Tan(Theta0[14]) - r1;
  max1 = rmax(r1, h1, dphi1 / 2.0);

  double q2 = Power(
      Sqrt((max1 / Cos(dphi1 / 2.0) / rmin) * (max1 / Cos(dphi1 / 2.0) / rmin) +
           1) -
          (max1 / Cos(dphi1 / 2.0) / rmin),
      1 / 15.0);

  for (int i = 15; i < 18; i++) {
    Dtheta[i] = 2.0 * ATan(Power(q2, i)) - 2.0 * ATan(Power(q2, i + 1));
    Theta0[i] = 2.0 * ATan(Power(q2, i + 1));
  }
  r2 = rDefocus(Theta0[15], Dtheta[15], rmin, 28 + 1);
  h2 = rmin / Tan(Theta0[17]) - r2;
  max2 = rmax(r2, h2, dphi2 / 2.0);

  double q3 = Power(
      Sqrt((max2 / Cos(dphi2 / 2.0) / rmin) * (max2 / Cos(dphi2 / 2.0) / rmin) +
           1) -
          (max2 / Cos(dphi2 / 2.0) / rmin),
      1 / 18.0);

  for (int i = 18; i < 22; i++) {
    Dtheta[i] = 2.0 * ATan(Power(q3, i)) - 2.0 * ATan(Power(q3, i + 1));
    Theta0[i] = 2.0 * ATan(Power(q3, i + 1));
  }

  for (int i = 0; i < 10; i++) {
    GeoData->Theta.push_back(Theta0[i + 12]);
    GeoData->DTheta.push_back(Dtheta[i + 12]);
  }

  // Create two Cone-like envelopes representing the whole detector volume
  Cone envelope1(14.5, (rmin - 0.5) / Tan(Theta0[11]) - 0.05,
                 (rmin - 0.5) / Tan(Theta0[21]) + 0.1,
                 (rmin + 28.5) / Tan(Theta0[11]) - 0.05,
                 (rmin + 28.5) / Tan(Theta0[21]) + 0.1);
  Cone envelope2(14.5, (rmin + 28.5) / Tan(Theta0[11]) - 0.05,
                 (rmin + 28.5) / Tan(Theta0[21]) + 0.1,
                 (rmin - 0.5) / Tan(Theta0[11]) - 0.05,
                 (rmin - 0.5) / Tan(Theta0[21]) + 0.1);
  Volume envelopeVol1(det_name, envelope1, air);
  Volume envelopeVol2(det_name, envelope2, air);
  envelopeVol1.setAttributes(lcdd, x_det.regionStr(), x_det.limitsStr(),
                             x_det.visStr());
  envelopeVol2.setAttributes(lcdd, x_det.regionStr(), x_det.limitsStr(),
                             x_det.visStr());

  // Get each module thickness
  double thickness[4];
  for (xml_coll_t clayer(x_det, _U(layer)); clayer; ++clayer) {
    xml_comp_t x_layer = clayer;
    int i = 0;
    for (xml_coll_t cmodule(x_layer, _U(module)); cmodule; cmodule++) {
      xml_comp_t x_module(cmodule);
      if (!x_module.isSensitive()) {
        mthickness = x_module.thickness();
        thickness[i] = mthickness;
        i++;
      }
    }
  }

  // Loop over all the layers in the detector
  // Right EndCap
  for (xml_coll_t clayer(x_det, _U(layer)); clayer; ++clayer) {
    xml_comp_t x_layer = clayer;
    double lthickness = x_layer.thickness();

    for (int row = 12; row < nRow; row++) {
      module_id = row;
      dtheta = Dtheta[module_id];
      theta0 = Theta0[module_id];

      double dphi, beta, gamma;
      if (module_id >= 12 && module_id < 15) {
        dphi = dphi1;
        beta = Theta0[12];
        gamma = Dtheta[12];
      }
      if (module_id >= 15 && module_id < 18) {
        dphi = dphi2;
        beta = Theta0[15];
        gamma = Dtheta[15];
      }
      if (module_id >= 18 && module_id < 22) {
        dphi = dphi3;
        beta = Theta0[18];
        gamma = Dtheta[18];
      }

      // Definition of geometry parameters
      a = rmin / Sin(theta0) * Sin(dtheta);
      b = lthickness;
      c = a + b * Tan(dtheta);

      double z0 = b / 2;
      double y0 = (c - 3 * a) / 4;
      // double x0=0;
      double y = z0 * Sin(M_PI / 2 - theta0 - dtheta) +
                 y0 * Cos(M_PI / 2 - theta0 - dtheta);
      double z = z0 * Cos(M_PI / 2 - theta0 - dtheta) -
                 y0 * Sin(M_PI / 2 - theta0 - dtheta);
      y = y + rmin / Tan(theta0);
      z = z - z0;

      // Creat a Trap-shaped envelope for layer
      string rowStr = _toString(row, "Rightrow%02d");
      Volume rowVol;
      rowVol = Volume(rowStr,
                      *MyTrap(theta0, dtheta, dphi, b + 1, rmin,
                              rDefocus(beta, gamma, rmin, b + 1),
                              RDefocus(beta, gamma, rmin, b + 1, dphi / 2), 0,
                              0, false, false),
                      air);
      DetElement rowdet(ECALEndcap, rowStr, det_id);
      rowdet.setAttributes(lcdd, rowVol, x_layer.regionStr(),
                           x_layer.limitsStr(), x_layer.visStr());
      // Definition of each module in layer
      double awidth = 0;
      for (xml_coll_t cmodule(x_layer, _U(module)); cmodule; cmodule++) {
        xml_comp_t x_module(cmodule);
        Material moduleMat = lcdd.material(x_module.materialStr());
        mthickness = x_module.thickness();
        string moduleStr = _toString(x_module.id(), "REmodule%02d");
        Volume *moduleVol;

        if (!x_module.isSensitive()) {
          awidth += mthickness;
          if (x_module.id() == 1)
            moduleVol = new Volume(
                moduleStr,
                SubtractionSolid(
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth - mthickness, 0, false, false),
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth, 0, false, false)),
                moduleMat);
          if (x_module.id() == 2)
            moduleVol = new Volume(
                moduleStr,
                SubtractionSolid(
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth - mthickness, 0, false, false),
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth, thickness[0], true, false),
                    Position(0, 0, mthickness / 2)),
                moduleMat);
          if (x_module.id() == 3 || x_module.id() == 4)
            moduleVol = new Volume(
                moduleStr,
                SubtractionSolid(
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth - mthickness, thickness[0], true, false),
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth, thickness[0], true, false),
                    Position(0, 0, mthickness / 2)),
                moduleMat);
        } else {
          if (x_module.id() == 5)
            moduleVol =
                new Volume(moduleStr,
                           *MyTrap(theta0, dtheta, dphi, b, rmin,
                                   rDefocus(beta, gamma, rmin, b + 1),
                                   RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                                   awidth, thickness[0], true, false),
                           moduleMat);
          if (x_module.id() == 6)
            moduleVol =
                new Volume(moduleStr, Box(0.25, 0.5, 0.0025), moduleMat);

          moduleVol->setSensitiveDetector(sens);
        }

        DetElement moduledet(rowdet, moduleStr, det_id);
        moduledet.setAttributes(lcdd, *moduleVol, x_module.regionStr(),
                                x_module.limitsStr(), x_module.visStr());
        // Place the module in the proper position
        PlacedVolume modulePlv;
        if (x_module.id() == 1 || x_module.id() == 2)
          modulePlv = rowVol.placeVolume(*moduleVol, Position(0, 0, 0));
        if (x_module.id() == 3)
          modulePlv = rowVol.placeVolume(*moduleVol,
                                         Position(0, 0, thickness[1] / 2.0));
        if (x_module.id() == 4)
          modulePlv = rowVol.placeVolume(
              *moduleVol, Position(0, 0, (thickness[1] + thickness[2]) / 2.0));
        if (x_module.id() == 5)
          modulePlv = rowVol.placeVolume(
              *moduleVol,
              Position(0, 0,
                       (thickness[1] + thickness[2] + thickness[3]) / 2.0));
        if (x_module.id() == 6)
          modulePlv = rowVol.placeVolume(*moduleVol, Position(0, 0, 14.0025));
        modulePlv.addPhysVolID("module", x_module.id());
        moduledet.setPlacement(modulePlv);
      }

      // Loop placing in phi direction
      int n;
      double offsety = rDefocus(beta, gamma, rmin, b + 1);
      if (module_id >= 12 && module_id < 15)
        n = nLadders1;
      if (module_id >= 15 && module_id < 18)
        n = nLadders2;
      if (module_id >= 18 && module_id < 22)
        n = nLadders3;
      for (int i = 0; i < n; i++) {
        double phi = dphi * i;
        double posx =
            offsety * Sin(-phi) + (y - offsety) * Sin(-phi - dphi / 2);
        double posy = offsety * Cos(phi) + (y - offsety) * Cos(phi + dphi / 2);

        rec::ECALERecGeoData::RightRowPos rightrowpos;
        rightrowpos.x = posx * 10;
        rightrowpos.y = posy * 10;
        rightrowpos.z = z * 10;
        GeoData->RightRows.push_back(rightrowpos);

        Transform3D tr(RotationZ(phi + dphi / 2) *
                           RotationZYX(0, 0, -M_PI / 2 + theta0 + dtheta),
                       Position(posx, posy, z));
        PlacedVolume rowPlv = envelopeVol1.placeVolume(rowVol, tr);
        rowPlv.addPhysVolID("row", i + 1);
        rowdet.setPlacement(rowPlv);
      }
    }
  }
  PlacedVolume envelopePlv1 =
      motherVol.placeVolume(envelopeVol1, Position(0, 0, zmax + 14));
  envelopePlv1.addPhysVolID("system", x_det.id());
  envelopePlv1.addPhysVolID("Endcap", 0);
  ECALEndcap.setPlacement(envelopePlv1);

  // Left EndCap
  for (xml_coll_t clayer(x_det, _U(layer)); clayer; ++clayer) {
    xml_comp_t x_layer = clayer;
    double lthickness = x_layer.thickness();

    for (int row = 12; row < nRow; row++) {
      module_id = row;
      dtheta = Dtheta[module_id];
      theta0 = Theta0[module_id];

      double dphi, beta, gamma;
      if (module_id >= 12 && module_id < 15) {
        dphi = dphi1;
        beta = Theta0[12];
        gamma = Dtheta[12];
      }
      if (module_id >= 15 && module_id < 18) {
        dphi = dphi2;
        beta = Theta0[15];
        gamma = Dtheta[15];
      }
      if (module_id >= 18 && module_id < 22) {
        dphi = dphi3;
        beta = Theta0[18];
        gamma = Dtheta[18];
      }

      // Definition of geometry parameters
      a = rmin / Sin(theta0) * Sin(dtheta);
      b = lthickness;
      c = a + b * Tan(dtheta);

      double z0 = -b / 2;
      double y0 = (c - 3 * a) / 4;
      // double x0=0;
      double y = -z0 * Sin(M_PI / 2 - theta0 - dtheta) +
                 y0 * Cos(M_PI / 2 - theta0 - dtheta);
      double z = z0 * Cos(M_PI / 2 - theta0 - dtheta) +
                 y0 * Sin(M_PI / 2 - theta0 - dtheta);
      y = y + rmin / Tan(theta0);
      z = z - z0;

      // Creat a Trap-shaped envelope for layer
      string rowStr = _toString(row, "Leftrow%02d");
      Volume rowVol;
      rowVol = Volume(rowStr,
                      *MyTrap(theta0, dtheta, dphi, b + 1, rmin,
                              rDefocus(beta, gamma, rmin, b + 1),
                              RDefocus(beta, gamma, rmin, b + 1, dphi / 2), 0,
                              0, false, true),
                      air);
      DetElement rowdet(ECALEndcap, rowStr, det_id);
      rowdet.setAttributes(lcdd, rowVol, x_layer.regionStr(),
                           x_layer.limitsStr(), x_layer.visStr());

      // Definition of each module in layer
      double awidth = 0;
      for (xml_coll_t cmodule(x_layer, _U(module)); cmodule; cmodule++) {
        xml_comp_t x_module(cmodule);
        Material moduleMat = lcdd.material(x_module.materialStr());
        mthickness = x_module.thickness();
        string moduleStr = _toString(x_module.id(), "LEmodule%02d");
        Volume *moduleVol;

        if (!x_module.isSensitive()) {
          awidth += mthickness;
          if (x_module.id() == 1)
            moduleVol = new Volume(
                moduleStr,
                SubtractionSolid(
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth - mthickness, 0, false, true),
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth, 0, false, true)),
                moduleMat);
          if (x_module.id() == 2)
            moduleVol = new Volume(
                moduleStr,
                SubtractionSolid(
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth - mthickness, 0, false, true),
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth, thickness[0], true, true),
                    Position(0, 0, -mthickness / 2)),
                moduleMat);
          if (x_module.id() == 3 || x_module.id() == 4)
            moduleVol = new Volume(
                moduleStr,
                SubtractionSolid(
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth - mthickness, thickness[0], true, true),
                    *MyTrap(theta0, dtheta, dphi, b, rmin,
                            rDefocus(beta, gamma, rmin, b + 1),
                            RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                            awidth, thickness[0], true, true),
                    Position(0, 0, -mthickness / 2)),
                moduleMat);
        } else {
          if (x_module.id() == 5)
            moduleVol =
                new Volume(moduleStr,
                           *MyTrap(theta0, dtheta, dphi, b, rmin,
                                   rDefocus(beta, gamma, rmin, b + 1),
                                   RDefocus(beta, gamma, rmin, b + 1, dphi / 2),
                                   awidth, thickness[0], true, true),
                           moduleMat);
          if (x_module.id() == 6)
            moduleVol =
                new Volume(moduleStr, Box(0.25, 0.5, 0.0025), moduleMat);

          moduleVol->setSensitiveDetector(sens);
        }
        DetElement moduledet(rowdet, moduleStr, det_id);
        moduledet.setAttributes(lcdd, *moduleVol, x_module.regionStr(),
                                x_module.limitsStr(), x_module.visStr());
        // Place the module in the proper position
        PlacedVolume modulePlv;
        if (x_module.id() == 1 || x_module.id() == 2)
          modulePlv = rowVol.placeVolume(*moduleVol, Position(0, 0, 0));
        if (x_module.id() == 3)
          modulePlv = rowVol.placeVolume(*moduleVol,
                                         Position(0, 0, -thickness[1] / 2.0));
        if (x_module.id() == 4)
          modulePlv = rowVol.placeVolume(
              *moduleVol, Position(0, 0, -(thickness[1] + thickness[2]) / 2.0));
        if (x_module.id() == 5)
          modulePlv = rowVol.placeVolume(
              *moduleVol,
              Position(0, 0,
                       -(thickness[1] + thickness[2] + thickness[3]) / 2.0));
        if (x_module.id() == 6)
          modulePlv = rowVol.placeVolume(*moduleVol, Position(0, 0, -14.0025));
        modulePlv.addPhysVolID("module", x_module.id());
        moduledet.setPlacement(modulePlv);
      }

      // Loop placing in phi direction
      int n;
      double offsety = rDefocus(beta, gamma, rmin, b + 1);
      if (module_id >= 12 && module_id < 15)
        n = nLadders1;
      if (module_id >= 15 && module_id < 18)
        n = nLadders2;
      if (module_id >= 18 && module_id < 22)
        n = nLadders3;
      for (int i = 0; i < n; i++) {
        double phi = dphi * i;
        double posx =
            offsety * Sin(-phi) + (y - offsety) * Sin(-phi - dphi / 2);
        double posy = offsety * Cos(phi) + (y - offsety) * Cos(phi + dphi / 2);

        rec::ECALERecGeoData::LeftRowPos leftrowpos;
        leftrowpos.x = posx * 10;
        leftrowpos.y = posy * 10;
        leftrowpos.z = z * 10;
        GeoData->LeftRows.push_back(leftrowpos);

        Transform3D tr(RotationZ(phi + dphi / 2) *
                           RotationZYX(0, 0, M_PI / 2 - theta0 - dtheta),
                       Position(posx, posy, z));
        PlacedVolume rowPlv = envelopeVol2.placeVolume(rowVol, tr);
        rowPlv.addPhysVolID("row", i + 1);
        rowdet.setPlacement(rowPlv);
      }
    }
  }

  PlacedVolume envelopePlv2 =
      motherVol.placeVolume(envelopeVol2, Position(0, 0, -zmax - 14));
  envelopePlv2.addPhysVolID("system", x_det.id());
  envelopePlv2.addPhysVolID("Endcap", 1);
  ECALEndcap.setPlacement(envelopePlv2);

  ECALEndcap.addExtension<rec::ECALERecGeoData>(GeoData);

  return ECALEndcap;
}
DECLARE_DETELEMENT(STCFEEMC_03, create_element)
