#include "DD4hep/DetFactoryHelper.h"
#include "TMath.h"
#include "XML/Layering.h"
#include <exception>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
using namespace TMath;

Trap *MyTrap(double theta0, double dtheta, double dphi, double b, double rmin,
             double offsetw, double cthickness, bool iswrapping, bool isMiddle,
             bool isLeft) {
  double a, c, h1, h2, alpha1, alpha2, h3, h4, mtheta, mphi;
  if (isMiddle) {
    a = 5 - 2 * offsetw;
    c = a;
    h1 = (rmin - (b - 28) / 2.0) * Tan(dphi / 2) - offsetw;
    h2 = h1;
    alpha1 = 0;
    alpha2 = 0;
    h3 = (b + rmin - (b - 28) / 2.0) * Tan(dphi / 2) - offsetw;
    h4 = h3;
    mtheta = 0;
    mphi = 0;
    if (iswrapping) {
      b = b - offsetw + cthickness;
      h1 = h1 + (offsetw - cthickness) * Tan(dphi / 2);
      h2 = h1;
    }
  } else {
    a = (rmin - (b - 28) / 2.0 * Cos(theta0)) * Tan(dtheta) / Cos(theta0);
    c = a + b * Tan(dtheta);
    h1 = (rmin - (b - 28) / 2.0 * Cos(theta0)) * Tan(dphi / 2);
    h2 = ((rmin - (b - 28) / 2.0 * Cos(theta0)) + a * Sin(theta0)) *
         Tan(dphi / 2);
    alpha1 = 0;
    alpha2 = 0;
    h3 = (b * Cos(theta0) + (rmin - (b - 28) / 2.0 * Cos(theta0))) *
         Tan(dphi / 2);
    h4 = (c * Sin(theta0) + b * Cos(theta0) +
          (rmin - (b - 28) / 2.0 * Cos(theta0))) *
         Tan(dphi / 2);
    mtheta = ATan((c - a) / b / 2);
    mphi = M_PI / 2;

    if (offsetw > 0.0) {
      double tan1 = (h2 - h1) / a;
      double tan2 = (h4 - h3) / c;
      h2 = h2 - offsetw * tan1 - offsetw;
      h1 = h1 + offsetw * tan1 - offsetw;
      h4 = h4 - offsetw * tan2 - offsetw;
      h3 = h3 + offsetw * tan2 - offsetw;
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
  }
  Trap *MyTrap1 = new Trap(b / 2, mtheta, mphi, a / 2, h1, h2, alpha1, c / 2,
                           h3, h4, alpha2);
  Trap *MyTrap2 = new Trap(b / 2, mtheta, -mphi, a / 2, h2, h1, alpha1, c / 2,
                           h4, h3, alpha2);
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
  xml_comp_t x_stave = x_det.staves();
  int det_id = x_det.id();
  double rmin = dim.rmin();
  int nLadders = dim.nsides();
  int nRow = dim.number();
  double dphi = 2 * Pi() / nLadders;
  double zmax = dim.z();
  double sthickness = x_stave.thickness();
  Material air = lcdd.air();
  DetElement ECALBarrel(det_name, x_det.id());
  Volume motherVol = lcdd.pickMotherVolume(ECALBarrel);
  sens.setType("calorimeter");
  double a, b, c, dtheta, theta0, mthickness;
  int module_id;
  double offsetr = 0;

  // Create a Tube-like envelope representing the whole detector volume
  Tube envelope(rmin - 0.5 - offsetr, rmin + 30 + offsetr, zmax, 0, 2 * M_PI);
  Volume envelopeVol(det_name, envelope, air);
  envelopeVol.setAttributes(lcdd, x_det.regionStr(), x_det.limitsStr(),
                            x_det.visStr());

  // Creat a Trapezoid-shaped stave
  Trapezoid stave((rmin - 0.5 - offsetr) * Tan(dphi / 2),
                  (rmin + 0.5 + sthickness + offsetr) * Tan(dphi / 2), zmax,
                  zmax, (sthickness + 1) / 2.0 + offsetr);
  Volume staveVol("stave", stave, air);
  DetElement stavedet(ECALBarrel, "stave", x_det.id());
  stavedet.setAttributes(lcdd, staveVol, x_stave.regionStr(),
                         x_stave.limitsStr(), x_stave.visStr());

  double Dtheta[25] = {2.72, 2.71, 2.70, 2.68, 2.65, 2.62, 2.58, 2.54, 2.50,
                       2.45, 2.40, 2.34, 2.28, 2.22, 2.16, 2.10, 2.04, 1.97,
                       1.91, 1.85, 1.78, 1.72, 1.66, 1.60, 1.54};
  double Theta0[25] = {2.72,  5.43,  8.13,  10.81, 13.46, 16.08, 18.66,
                       21.20, 23.70, 26.15, 28.55, 30.89, 33.17, 35.39,
                       37.55, 39.65, 41.69, 43.66, 45.57, 47.42, 49.20,
                       50.92, 52.58, 54.18, 55.72};
  for (int i = 0; i < 25; i++) {
    Dtheta[i] = Dtheta[i] * M_PI / 180.0;
    Theta0[i] = Theta0[i] * M_PI / 180.0;
  }

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
  for (xml_coll_t clayer(x_det, _U(layer)); clayer; ++clayer) {

    xml_comp_t x_layer = clayer;
    double lthickness = x_layer.thickness();
    // double dtheta  = x_layer.angle();
    // dtheta = dtheta*M_PI/180;

    for (int row = 1; row <= nRow; row++) {
      module_id = row;
      if (module_id > 26) {
        dtheta = Dtheta[module_id - 27];
        theta0 = Theta0[module_id - 27];
      }
      if (module_id < 26) {
        dtheta = Dtheta[25 - module_id];
        theta0 = Theta0[25 - module_id];
      }
      if (module_id == 26) {
        dtheta = 0;
        theta0 = 0;
      }

      // Definition of geometry parameters
      a = rmin * Tan(dtheta) / Cos(theta0);
      b = lthickness;
      c = a + b * Tan(dtheta);
      double y0 = (a + c) / 4;
      // double x0=0;
      double z0 = b / 2;
      double x = 0;
      double y = -z0 * Sin(theta0) + y0 * Cos(theta0);
      double z = z0 * Cos(theta0) + y0 * Sin(theta0);
      y = y - rmin * Tan(theta0) - 2.5;
      z = z - sthickness / 2;
      double theta1;
      if (module_id > 26) {
        theta1 = theta0;
      } else {
        theta1 = -theta0;
      }
      if (module_id < 26)
        y = -y;
      if (module_id == 26) {
        y = 0;
        z = z0 - sthickness / 2;
      }

      // Creat a Trap-shaped envelope for layer
      string rowStr = _toString(row, "row%02d");
      Volume rowVol;
      //	Volume rowVol(rowStr, *MyTrap(theta0, dtheta, dphi, b, rmin,
      //0,0,false, true,false), air);
      if (module_id == 26)
        rowVol = Volume(rowStr,
                        *MyTrap(theta0, dtheta, dphi, b + 1, rmin, 0, 0, false,
                                true, false),
                        air);
      if (module_id > 26)
        rowVol = Volume(rowStr,
                        *MyTrap(theta0, dtheta, dphi, b + 1, rmin, 0, 0, false,
                                false, false),
                        air);
      if (module_id < 26)
        rowVol = Volume(rowStr,
                        *MyTrap(theta0, dtheta, dphi, b + 1, rmin, 0, 0, false,
                                false, true),
                        air);
      DetElement rowdet(stavedet, rowStr, det_id);
      rowdet.setAttributes(lcdd, rowVol, x_layer.regionStr(),
                           x_layer.limitsStr(), x_layer.visStr());

      // Definition of each module in layer
      if (module_id == 26) {
        double awidth = 0;
        for (xml_coll_t cmodule(x_layer, _U(module)); cmodule; cmodule++) {
          xml_comp_t x_module(cmodule);
          Material moduleMat = lcdd.material(x_module.materialStr());
          mthickness = x_module.thickness();
          string moduleStr = _toString(x_module.id(), "Bmodule%02d");
          Volume *moduleVol;
          if (!x_module.isSensitive()) {
            awidth += mthickness;
            if (x_module.id() == 1)
              moduleVol = new Volume(
                  moduleStr,
                  SubtractionSolid(*MyTrap(theta0, dtheta, dphi, b, rmin,
                                           awidth - mthickness, 0, false, true,
                                           false),
                                   *MyTrap(theta0, dtheta, dphi, b, rmin,
                                           awidth, 0, false, true, false)),
                  moduleMat);
            if (x_module.id() == 2)
              moduleVol = new Volume(
                  moduleStr,
                  SubtractionSolid(
                      *MyTrap(theta0, dtheta, dphi, b, rmin,
                              awidth - mthickness, 0, false, true, false),
                      *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                              thickness[0], true, true, false),
                      Position(0, 0, mthickness / 2)),
                  moduleMat);
            if (x_module.id() == 3 || x_module.id() == 4)
              moduleVol =
                  new Volume(moduleStr,
                             SubtractionSolid(
                                 *MyTrap(theta0, dtheta, dphi, b, rmin,
                                         awidth - mthickness, thickness[0],
                                         true, true, false),
                                 *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                                         thickness[0], true, true, false),
                                 Position(0, 0, mthickness / 2)),
                             moduleMat);
          } else {
            if (x_module.id() == 5)
              moduleVol =
                  new Volume(moduleStr,
                             *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                                     thickness[0], true, true, false),
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
                *moduleVol,
                Position(0, 0, (thickness[1] + thickness[2]) / 2.0));
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
      } else {
        double awidth = 0;
        for (xml_coll_t cmodule(x_layer, _U(module)); cmodule; cmodule++) {
          xml_comp_t x_module(cmodule);
          Material moduleMat = lcdd.material(x_module.materialStr());
          mthickness = x_module.thickness();
          string moduleStr = _toString(x_module.id(), "Bmodule%02d");
          Volume *moduleVol;
          if (!x_module.isSensitive()) {
            awidth += mthickness;
            if (module_id > 26) {
              if (x_module.id() == 1)
                moduleVol = new Volume(
                    moduleStr,
                    SubtractionSolid(*MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth - mthickness, 0, false,
                                             false, false),
                                     *MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth, 0, false, false, false)),
                    moduleMat);
              if (x_module.id() == 2)
                moduleVol = new Volume(
                    moduleStr,
                    SubtractionSolid(
                        *MyTrap(theta0, dtheta, dphi, b, rmin,
                                awidth - mthickness, 0, false, false, false),
                        *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                                thickness[0], true, false, false),
                        Position(0, 0, mthickness / 2)),
                    moduleMat);
              if (x_module.id() == 3 || x_module.id() == 4)
                moduleVol = new Volume(
                    moduleStr,
                    SubtractionSolid(*MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth - mthickness, thickness[0],
                                             true, false, false),
                                     *MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth, thickness[0], true, false,
                                             false),
                                     Position(0, 0, mthickness / 2)),
                    moduleMat);
            } else {
              if (x_module.id() == 1)
                moduleVol = new Volume(
                    moduleStr,
                    SubtractionSolid(*MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth - mthickness, 0, false,
                                             false, true),
                                     *MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth, 0, false, false, true)),
                    moduleMat);
              if (x_module.id() == 2)
                moduleVol = new Volume(
                    moduleStr,
                    SubtractionSolid(
                        *MyTrap(theta0, dtheta, dphi, b, rmin,
                                awidth - mthickness, 0, false, false, true),
                        *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                                thickness[0], true, false, true),
                        Position(0, 0, mthickness / 2)),
                    moduleMat);
              if (x_module.id() == 3 || x_module.id() == 4)
                moduleVol = new Volume(
                    moduleStr,
                    SubtractionSolid(*MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth - mthickness, thickness[0],
                                             true, false, true),
                                     *MyTrap(theta0, dtheta, dphi, b, rmin,
                                             awidth, thickness[0], true, false,
                                             true),
                                     Position(0, 0, mthickness / 2)),
                    moduleMat);
            }
          } else {
            if (x_module.id() == 5) {
              if (module_id > 26)
                moduleVol =
                    new Volume(moduleStr,
                               *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                                       thickness[0], true, false, false),
                               moduleMat);
              else
                moduleVol =
                    new Volume(moduleStr,
                               *MyTrap(theta0, dtheta, dphi, b, rmin, awidth,
                                       thickness[0], true, false, true),
                               moduleMat);
            }
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
                *moduleVol,
                Position(0, 0, (thickness[1] + thickness[2]) / 2.0));
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
      }

      // After putting all the module into the layer,palce the row into the
      // stave
      PlacedVolume rowPlv = staveVol.placeVolume(
          rowVol, Transform3D(RotationZYX(0, 0, theta1), Position(x, y, z)));
      rowPlv.addPhysVolID("row", row);
      rowdet.setPlacement(rowPlv);
    }
  }

  // Loop placing in phi direction
  double offsetx = 0;
  double offsety = rmin + sthickness / 2;
  for (int i = 0; i < nLadders; i++) {
    double phi1 = dphi * i;
    double posx = offsetx * Cos(phi1) - offsety * Sin(phi1);
    double posy = offsetx * Sin(phi1) + offsety * Cos(phi1);
    Transform3D tr(RotationZYX(0, -phi1, -M_PI / 2), Position(posx, posy, 0));
    PlacedVolume stavePlv = envelopeVol.placeVolume(staveVol, tr);
    stavePlv.addPhysVolID("stave", i + 1);
    stavedet.setPlacement(stavePlv);
  }
  // put the envelope into world
  PlacedVolume envelopePlv =
      motherVol.placeVolume(envelopeVol, Position(0, 0, 0));
  envelopePlv.addPhysVolID("system", x_det.id());
  envelopePlv.addPhysVolID("barrel", 0);
  ECALBarrel.setPlacement(envelopePlv);
  return ECALBarrel;
}
DECLARE_DETELEMENT(STCFBEMC_02, create_element)
