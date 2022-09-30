//====================================================================
//  Simple tracking detector made from planar sensors that are parallel
//  to the z-axis. There are two materials per ladder: one sensitive
//  and one support.
//--------------------------------------------------------------------
//
//  Author     : F.Gaede
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"

#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hep/DetType.h"
#include <exception>

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {

  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  xml_dim_t dim = x_det.dimensions();
  std::string gasvis = dim.visStr();

  DetElement innerTracker(name, x_det.id());
  DetType type( DetType::TRACKER);
  innerTracker.setTypeFlag(type.to_ulong());
 
  PlacedVolume pv;

  double rmin = dim.inner_r();
  double rmax = dim.outer_r();
  double zmax = dim.length();
  Tube TubeInner(rmin, rmax, zmax);

  Volume SubSolidInner_vol(name + "_Gas", TubeInner,
                           description.material("Air"));
  SubSolidInner_vol.setVisAttributes(description.visAttributes("gasVis"));
  sens.setType("innerTracker");

  double r1, r2;
  int i = 1;

  for (xml_coll_t c(e, _U(layer)); c; ++c, i++) {
    xml_comp_t x_layer(c);

    double startr = x_layer.attr<double>(_Unicode(startr));
    double halfz = x_layer.attr<double>(_Unicode(halfz));

    xml_comp_t x_strip = x_layer.child(_Unicode(strip));
    xml_comp_t x_uRwell = x_layer.child(_Unicode(uRwell));

    double Cu1_th = x_uRwell.attr<double>(_Unicode(Cu1_th));
    std::string Cu1_ma = x_uRwell.attr<std::string>(_Unicode(Cu1_ma));
    double PI1_th = x_uRwell.attr<double>(_Unicode(PI1_th));
    std::string PI1_ma = x_uRwell.attr<std::string>(_Unicode(PI1_ma));
    double sup_th = x_uRwell.attr<double>(_Unicode(sup_th));
    std::string sup_ma = x_uRwell.attr<std::string>(_Unicode(sup_ma));
    double PI2_th = x_uRwell.attr<double>(_Unicode(PI2_th));
    std::string PI2_ma = x_uRwell.attr<std::string>(_Unicode(PI2_ma));
    double Cu2_th = x_uRwell.attr<double>(_Unicode(Cu2_th));
    std::string Cu2_ma = x_uRwell.attr<std::string>(_Unicode(Cu2_ma));
    double Drift_Gas_th = x_uRwell.attr<double>(_Unicode(Drift_Gas_th));
    std::string Drift_Gas_ma =
        x_uRwell.attr<std::string>(_Unicode(Drift_Gas_ma));
    double Copper_th = x_uRwell.attr<double>(_Unicode(Copper_th));
    std::string Copper_ma = x_uRwell.attr<std::string>(_Unicode(Copper_ma));
    double well_th = x_uRwell.attr<double>(_Unicode(well_th));
    std::string well_ma = x_uRwell.attr<std::string>(_Unicode(well_ma));
    double DLC_th = x_uRwell.attr<double>(_Unicode(DLC_th));
    std::string DLC_ma = x_uRwell.attr<std::string>(_Unicode(DLC_ma));
    double base_th = x_uRwell.attr<double>(_Unicode(base_th));
    std::string base_ma = x_uRwell.attr<std::string>(_Unicode(base_ma));
    double endRing_h = x_uRwell.attr<double>(_Unicode(endRing_h));
    std::string endRing_ma = x_uRwell.attr<std::string>(_Unicode(endRing_ma));

    std::string strip_ma = x_strip.attr<std::string>(_Unicode(material));
    double Strip_th = x_strip.attr<double>(_Unicode(Strip_th));
    double topStrip_wid = x_strip.attr<double>(_Unicode(topStrip_wid));
    double mediumStrip_th = x_strip.attr<double>(_Unicode(mediumStrip_th));
    //double mediumStrip_wid = x_strip.attr<double>(_Unicode(mediumStrip_wid));
    double quan_top = x_strip.attr<double>(_Unicode(quan_top));
    double dr1 = x_strip.attr<double>(_Unicode(Strip1_dr));
    double dr2 = x_strip.attr<double>(_Unicode(Strip2_dr));

    r1 = startr;
    r2 = r1 + Cu1_th;
    Tube Tube1(r1, r2, halfz);
    Volume Cu1("Cu1_" + _toString(i), Tube1, description.material(Cu1_ma));
    Cu1.setVisAttributes(description.visAttributes("CuVis"));
    pv = SubSolidInner_vol.placeVolume(Cu1);

    r1 = r2;
    r2 = r1 + PI1_th;
    Tube Tube2(r1, r2, halfz);
    Volume PI1("PI1_" + _toString(i), Tube2, description.material(PI1_ma));
    PI1.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(PI1);

    r1 = r2;
    r2 = r1 + sup_th;
    Tube Tube3(r1, r2, halfz);
    Volume support("support_" + _toString(i), Tube3,
                   description.material(sup_ma));
    support.setVisAttributes(description.visAttributes("SupportVis"));
    pv = SubSolidInner_vol.placeVolume(support);

    r1 = r2;
    r2 = r1 + PI2_th;
    Tube Tube4(r1, r2, halfz);
    Volume PI2("PI2_" + _toString(i), Tube4, description.material(PI2_ma));
    PI2.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(PI2);

    r1 = r2;
    r2 = r1 + Cu2_th;
    Tube Tube5(r1, r2, halfz);
    Volume Cu2("Cu2_" + _toString(i), Tube5, description.material(Cu2_ma));
    Cu2.setVisAttributes(description.visAttributes("CuVis"));
    pv = SubSolidInner_vol.placeVolume(Cu2);

    r1 = r2;
    r2 = r1 + Drift_Gas_th;
    Tube Tube6(r1, r2, halfz);
    Volume Drift_Gas("Drift_Gas_" + _toString(i), Tube6,
                     description.material(Drift_Gas_ma));
    Drift_Gas.setVisAttributes(description.visAttributes("gasVis"));
    pv = SubSolidInner_vol.placeVolume(Drift_Gas);

    r1 = r2;
    r2 = r1 + Copper_th;
    Tube Tube7(r1, r2, halfz);
    Volume Copper("Copper_" + _toString(i), Tube7,
                  description.material(Copper_ma));
    Copper.setVisAttributes(description.visAttributes("CuVis"));
    pv = SubSolidInner_vol.placeVolume(Copper);

    r1 = r2;
    r2 = r1 + well_th;
    Tube Tube8(r1, r2, halfz);
    Volume well("well_" + _toString(i), Tube8, description.material(well_ma));
    well.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(well);

    r1 = r2;
    r2 = r1 + DLC_th;
    Tube Tube9(r1, r2, halfz);
    Volume DLC("DLC_" + _toString(i), Tube9, description.material(DLC_ma));
    DLC.setVisAttributes(description.visAttributes("DLCVis"));
    pv = SubSolidInner_vol.placeVolume(DLC);

    r1 = r2;
    r2 = r1 + base_th;
    Tube Tube10(r1, r2, halfz);
    Volume base("base_" + _toString(i), Tube10, description.material(base_ma));
    base.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(base);

    Tube topStripTube(r1 + dr1, r1 + dr1 + Strip_th, halfz,
                      -0.5 * topStrip_wid / (r1 + dr1) * rad,
                      0.5 * topStrip_wid / (r1 + dr1) * rad);
    Volume topStrip("topStrip_", topStripTube, description.material(strip_ma));
    topStrip.setVisAttributes(description.visAttributes("StripVis"));
    for (int j = 0; j < quan_top; j++) {
      RotationZYX rotate(-2. * M_PI * double(j) / quan_top * rad, 0, 0);
      pv = base.placeVolume(topStrip, rotate);
    }

    Tube mediumStripTube(r1 + dr2, r1 + dr2 + mediumStrip_th, halfz);
    Volume mediumStrip("mediumStrip_", mediumStripTube,
                       description.material(strip_ma));
    mediumStrip.setVisAttributes(description.visAttributes("StripVis"));
    pv = base.placeVolume(mediumStrip);

    Tube endRingTube(startr - 2. * mm, r2 + 2. * mm, endRing_h / 2.);
    Volume endRing("endRing_", endRingTube, description.material(endRing_ma));
    endRing.setVisAttributes(description.visAttributes("EpoxyVis"));
    Position pos1(0., 0., halfz + endRing_h / 2.);
    Position pos2(0., 0., -halfz - endRing_h / 2.);
    pv = SubSolidInner_vol.placeVolume(endRing, pos1);
    pv = SubSolidInner_vol.placeVolume(endRing, pos2);
  }

  Volume mother = description.pickMotherVolume(innerTracker);
  pv = mother.placeVolume(SubSolidInner_vol);
  pv.addPhysVolID("system", x_det.id());
  innerTracker.setPlacement(pv);
  return innerTracker;
}

DECLARE_DETELEMENT(CuRwell02, create_element)
