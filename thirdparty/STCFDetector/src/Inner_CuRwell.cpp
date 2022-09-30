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

  xml_comp_t x_endcap = x_det.child(_Unicode(endcap));

  xml_comp_t x_innerWall = x_det.child(_Unicode(inWall));

  double iw_th = x_innerWall.attr<double>(_Unicode(thickness));
  double ir = x_innerWall.attr<double>(_Unicode(rmin));
  double iw_l = x_innerWall.attr<double>(_Unicode(halfz));
  std::string iw_ma = x_innerWall.attr<std::string>(_Unicode(material));

  Tube Wall(ir, ir + iw_th, iw_l);
  Volume innerWall("innerWall_", Wall, description.material(iw_ma));
  innerWall.setVisAttributes(description.visAttributes("InnerVis"));
  pv = SubSolidInner_vol.placeVolume(innerWall);

  double end_rmin = x_endcap.attr<double>(_Unicode(rmin));
  double end_rmax = x_endcap.attr<double>(_Unicode(rmax));
  double end_th = x_endcap.attr<double>(_Unicode(thickness));
  double end_l = x_endcap.attr<double>(_Unicode(halfz));
  double end_pos = x_endcap.attr<double>(_Unicode(zpos));
  std::string end_ma = x_endcap.attr<std::string>(_Unicode(material));

  Cone Endcapcone("endcap", end_l, end_rmin, end_rmin + end_th, end_rmax,
                  end_rmax + end_th);
  Volume endcap("endcap_", Endcapcone, description.material(end_ma));
  endcap.setVisAttributes(description.visAttributes("InnerVis"));
  pv = SubSolidInner_vol.placeVolume(endcap, Position(0., 0., end_pos));
  Position pos(0, 0, -end_pos);
  RotationZYX rot(0, 0, 180 * deg);
  pv = SubSolidInner_vol.placeVolume(endcap, Transform3D(rot, pos));

  for (xml_coll_t c(e, _U(layer)); c; ++c, i++) {
    xml_comp_t x_layer(c);

    double startr = x_layer.attr<double>(_Unicode(startr));
    double halfz = x_layer.attr<double>(_Unicode(halfz));

    xml_comp_t x_strip = x_layer.child(_Unicode(strip));
    xml_comp_t x_uRwell = x_layer.child(_Unicode(uRwell));

    double Cathode_Kapton_th =
        x_uRwell.attr<double>(_Unicode(Cathode_Kapton_th));
    std::string Cathode_Kapton_ma =
        x_uRwell.attr<std::string>(_Unicode(Cathode_Kapton_ma));
    double Cathode_DLC_th = x_uRwell.attr<double>(_Unicode(Cathode_DLC_th));
    std::string Cathode_DLC_ma =
        x_uRwell.attr<std::string>(_Unicode(Cathode_DLC_ma));
    double Drift_Gas_th = x_uRwell.attr<double>(_Unicode(Drift_Gas_th));
    std::string Drift_Gas_ma =
        x_uRwell.attr<std::string>(_Unicode(Drift_Gas_ma));
    double PCB_DLC1_th = x_uRwell.attr<double>(_Unicode(PCB_DLC1_th));
    std::string PCB_DLC1_ma = x_uRwell.attr<std::string>(_Unicode(PCB_DLC1_ma));
    double PCB_Well_Kapton_th =
        x_uRwell.attr<double>(_Unicode(PCB_Well_Kapton_th));
    std::string PCB_Well_Kapton_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_Well_Kapton_ma));
    double PCB_DLC2_th = x_uRwell.attr<double>(_Unicode(PCB_DLC2_th));
    std::string PCB_DLC2_ma = x_uRwell.attr<std::string>(_Unicode(PCB_DLC2_ma));
    double PCB_Epoxy1_th = x_uRwell.attr<double>(_Unicode(PCB_Epoxy1_th));
    std::string PCB_Epoxy1_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_Epoxy1_ma));
    double PCB_Prepreg_th = x_uRwell.attr<double>(_Unicode(PCB_Prepreg_th));
    std::string PCB_Prepreg_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_Prepreg_ma));
    double PCB_Kapton1_th = x_uRwell.attr<double>(_Unicode(PCB_Kapton1_th));
    std::string PCB_Kapton1_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_Kapton1_ma));
    double PCB_topEpoxy_th = x_uRwell.attr<double>(_Unicode(PCB_topEpoxy_th));
    std::string PCB_topEpoxy_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_topEpoxy_ma));
    double PCB_Kapton2_th = x_uRwell.attr<double>(_Unicode(PCB_Kapton2_th));
    std::string PCB_Kapton2_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_Kapton2_ma));
    double PCB_mediumEpoxy_th =
        x_uRwell.attr<double>(_Unicode(PCB_mediumEpoxy_th));
    std::string PCB_mediumEpoxy_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_mediumEpoxy_ma));
    double PCB_Kapton3_th = x_uRwell.attr<double>(_Unicode(PCB_Kapton3_th));
    std::string PCB_Kapton3_ma =
        x_uRwell.attr<std::string>(_Unicode(PCB_Kapton3_ma));

    double glue_th = x_strip.attr<double>(_Unicode(glue_th));
    std::string strip_ma = x_strip.attr<std::string>(_Unicode(material));
    double topStrip_wid = x_strip.attr<double>(_Unicode(topStrip_wid));
    double mediumStrip_wid = x_strip.attr<double>(_Unicode(mediumStrip_wid));
    double quan_top = x_strip.attr<double>(_Unicode(quan_top));
    double quan_medium = x_strip.attr<double>(_Unicode(quan_medium));

    r1 = startr;
    r2 = r1 + Cathode_Kapton_th;
    Tube Tube1(r1, r2, halfz);
    Volume Cathode_Kapton("Cathode_Kapton_" + _toString(i), Tube1,
                          description.material(Cathode_Kapton_ma));
    Cathode_Kapton.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(Cathode_Kapton);

    r1 = r2;
    r2 = r1 + Cathode_DLC_th;
    Tube Tube2(r1, r2, halfz);
    Volume Cathode_DLC("Cathode_DLC_" + _toString(i), Tube2,
                       description.material(Cathode_DLC_ma));
    Cathode_DLC.setVisAttributes(description.visAttributes("DLCVis"));
    pv = SubSolidInner_vol.placeVolume(Cathode_DLC);

    r1 = r2;
    r2 = r1 + Drift_Gas_th;
    Tube Tube3(r1, r2, halfz);
    Volume Drift_Gas("Drift_Gas_" + _toString(i), Tube3,
                     description.material(Drift_Gas_ma));
    Drift_Gas.setVisAttributes(description.visAttributes("gasVis"));
    pv = SubSolidInner_vol.placeVolume(Drift_Gas);

    r1 = r2;
    r2 = r1 + PCB_DLC1_th;
    Tube Tube4(r1, r2, halfz);
    Volume PCB_DLC1("PCB_DLC1_" + _toString(i), Tube4,
                    description.material(PCB_DLC1_ma));
    PCB_DLC1.setVisAttributes(description.visAttributes("DLCVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_DLC1);

    r1 = r2;
    r2 = r1 + PCB_Well_Kapton_th;
    Tube Tube5(r1, r2, halfz);
    Volume PCB_Well_Kapton("PCB_Well_Kapton_" + _toString(i), Tube5,
                           description.material(PCB_Well_Kapton_ma));
    PCB_Well_Kapton.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_Well_Kapton);

    r1 = r2;
    r2 = r1 + PCB_DLC2_th;
    Tube Tube6(r1, r2, halfz);
    Volume PCB_DLC2("PCB_DLC2_" + _toString(i), Tube6,
                    description.material(PCB_DLC2_ma));
    PCB_DLC2.setVisAttributes(description.visAttributes("DLCVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_DLC2);

    r1 = r2;
    r2 = r1 + PCB_Epoxy1_th;
    Tube Tube7(r1, r2, halfz);
    Volume PCB_Epoxy1("PCB_Epoxy1_" + _toString(i), Tube7,
                      description.material(PCB_Epoxy1_ma));
    PCB_Epoxy1.setVisAttributes(description.visAttributes("EpoxyVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_Epoxy1);

    r1 = r2;
    r2 = r1 + PCB_Prepreg_th;
    Tube Tube8(r1, r2, halfz);
    Volume PCB_Prepreg("PCB_Prepreg_" + _toString(i), Tube8,
                       description.material(PCB_Prepreg_ma));
    PCB_Prepreg.setVisAttributes(description.visAttributes("PrepregVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_Prepreg);

    r1 = r2;
    r2 = r1 + PCB_Kapton1_th;
    Tube Tube9(r1, r2, halfz);
    Volume PCB_Kapton1("PCB_Kapton1_" + _toString(i), Tube9,
                       description.material(PCB_Kapton1_ma));
    PCB_Kapton1.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_Kapton1);

    r1 = r2;
    r2 = r1 + PCB_topEpoxy_th;
    Tube Tube10(r1, r2, halfz);
    Volume PCB_topEpoxy("PCB_topEpoxy_" + _toString(i), Tube10,
                        description.material(PCB_topEpoxy_ma));
    PCB_topEpoxy.setVisAttributes(description.visAttributes("EpoxyVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_topEpoxy);

    Tube topStripTube(r1 + glue_th, r2 - glue_th, halfz,
                      -0.5 * topStrip_wid / (r1 + glue_th) * rad,
                      0.5 * topStrip_wid / (r1 + glue_th) * rad);
    Volume topStrip("topStrip_", topStripTube, description.material(strip_ma));
    topStrip.setVisAttributes(description.visAttributes("StripVis"));
    for (int j = 0; j < quan_top; j++) {
      RotationZYX rotate(-2. * M_PI * double(i) / quan_top * rad, 0, 0);
      pv = PCB_topEpoxy.placeVolume(topStrip, rotate);
    }

    r1 = r2;
    r2 = r1 + PCB_Kapton2_th;
    Tube Tube11(r1, r2, halfz);
    Volume PCB_Kapton2("PCB_Kapton2_" + _toString(i), Tube11,
                       description.material(PCB_Kapton2_ma));
    PCB_Kapton2.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_Kapton2);

    r1 = r2;
    r2 = r1 + PCB_mediumEpoxy_th;
    Tube Tube12(r1, r2, halfz);
    Volume PCB_mediumEpoxy("PCB_mediumEpoxy_" + _toString(i), Tube12,
                           description.material(PCB_mediumEpoxy_ma));
    PCB_mediumEpoxy.setVisAttributes(description.visAttributes("EpoxyVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_mediumEpoxy);

    Tube mediumStripTube(r1 + glue_th, r2 - glue_th, mediumStrip_wid / 2.);
    Volume mediumStrip("mediumStrip_", mediumStripTube,
                       description.material(strip_ma));
    mediumStrip.setVisAttributes(description.visAttributes("StripVis"));
    for (int j = 0; j < quan_medium; j++) {
      double z_in = halfz * 2. / quan_medium;
      double z_pos = -z_in * (quan_medium / 2. - 0.5);
      Position pos(0, 0, (z_pos + z_in * double(j)) * mm);
      pv = PCB_mediumEpoxy.placeVolume(mediumStrip, pos);
    }

    r1 = r2;
    r2 = r1 + PCB_Kapton3_th;
    Tube Tube13(r1, r2, halfz);
    Volume PCB_Kapton3("PCB_Kapton3_" + _toString(i), Tube13,
                       description.material(PCB_Kapton3_ma));
    PCB_Kapton3.setVisAttributes(description.visAttributes("KaptonVis"));
    pv = SubSolidInner_vol.placeVolume(PCB_Kapton3);
  }

  Volume mother = description.pickMotherVolume(innerTracker);
  pv = mother.placeVolume(SubSolidInner_vol);
  pv.addPhysVolID("system", x_det.id());
  innerTracker.setPlacement(pv);
  return innerTracker;
}

DECLARE_DETELEMENT(CuRwell, create_element)
