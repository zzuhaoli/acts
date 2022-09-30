#include "DD4hep/DetFactoryHelper.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include "DD4hep/DetType.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

// --------------------------------------------------------------

// Coding conventions:
// Global variables start with upper character "G".

static Ref_t create_detector(Detector &lcdd, xml_h e, SensitiveDetector sens) {
  // #################### init ##################################
  xml_det_t G_x_det = e;
  string G_det_name = G_x_det.nameStr();
  Material G_air = lcdd.air();
  DetElement G_det(G_det_name, G_x_det.id());
  DetType type( DetType::TRACKER);
  G_det.setTypeFlag(type.to_ulong());
 
  Volume G_motherVol = lcdd.pickMotherVolume(G_det);

  xml_coll_t G_constant(G_x_det, _U(constant));
  xml_comp_t G_x_constant = G_constant;
  xml_coll_t G_barrel(G_x_det, _U(barrel));
  xml_comp_t G_x_barrel = G_barrel;
  xml_coll_t G_endcap(G_x_det, _U(endcap));
  xml_comp_t G_x_endcap = G_endcap;

  // #################### read constants ########################
  double G_barrel_rmin = G_x_constant.rmin();
  double G_zhalf = G_x_constant.zhalf();
  double G_endcap_angle = G_x_constant.angle();
  int G_numsides = G_x_constant.numsides();

  // #################### declare global variables ##############
  double G_barrel_total_thickness = 0;
  double G_endcap_inner_thickness = 0;
  double G_endcap_total_thickness = 0;

  // caculate global variables
  for (xml_coll_t i(G_x_barrel, _U(slice)); i; i++) {
    xml_comp_t x_slice = i;
    G_barrel_total_thickness += x_slice.thickness();
  }
  // ---------------
  xml_coll_t aa(G_x_endcap, _U(slice));
  xml_comp_t aa_slice = aa;
  G_endcap_inner_thickness = aa_slice.thickness();
  // ---------------
  for (xml_coll_t i(G_x_endcap, _U(slice)); i; i++) {
    xml_comp_t x_slice = i;
    G_endcap_total_thickness += x_slice.thickness();
  }

  // ************************************************************
  // #################### define envelopeVol ####################
  // ---------- define barrel geometry ----------------
  PolyhedraRegular G_barrel_geo(G_numsides, G_barrel_rmin,
                                G_barrel_rmin + G_barrel_total_thickness,
                                2 * G_zhalf);

  // ---------- define endcap geometry ----------------
  PolyhedraRegular G_endcap_geo(G_numsides, G_zhalf * tan(G_endcap_angle),
                                G_barrel_rmin + G_barrel_total_thickness,
                                G_endcap_total_thickness);

  // ---------- place endcap to barrel ----------------
  Transform3D env_trans_zp(
      RotationZ(0),
      Translation3D(0, 0, G_zhalf + G_endcap_total_thickness / 2));
  Transform3D env_trans_zn(
      RotationZ(0),
      Translation3D(0, 0, -G_zhalf - G_endcap_total_thickness / 2));

  UnionSolid G_env_geo(G_barrel_geo, G_endcap_geo, env_trans_zp);
  G_env_geo = UnionSolid(G_env_geo, G_endcap_geo, env_trans_zn);

  Volume G_envelope_vol(G_det_name, G_env_geo, G_air);
  G_envelope_vol.setAttributes(lcdd, G_x_det.regionStr(), G_x_det.limitsStr(),
                               G_x_det.visStr());

  // ************************************************************
  // #################### define endtrap geometry ###############
  double endcap_thickness = 0;
  UnionSolid G_endtrap_geo_zp;
  UnionSolid G_endtrap_geo_zn;
  Trap slice_tmp_geo;

  for (xml_coll_t i(G_x_endcap, _U(slice)); i; ++i) {
    xml_comp_t x_slice = i;

    // slice's size
    double slice_inner_width = 2 * (G_zhalf + endcap_thickness) *
                               tan(G_endcap_angle) * tan(M_PI / G_numsides);
    double slice_outer_width =
        2 * (G_barrel_rmin + G_barrel_total_thickness) * tan(M_PI / G_numsides);
    double slice_length = G_barrel_rmin + G_barrel_total_thickness -
                          (G_zhalf + endcap_thickness) * tan(G_endcap_angle);
    double slice_thickness = x_slice.thickness();

    Trap slice_endtrap_geo(slice_thickness / 2, 0, 0, slice_length / 2,
                           slice_inner_width / 2, slice_outer_width / 2, 0,
                           slice_length / 2, slice_inner_width / 2,
                           slice_outer_width / 2, 0);

    if (!G_endtrap_geo_zp) {
      if (!slice_tmp_geo) {
        slice_tmp_geo = slice_endtrap_geo;
      } else {
        Transform3D slice_endtrap_trans_zp(
            Translation3D(0, endcap_thickness * tan(G_endcap_angle) / 2,
                          endcap_thickness + slice_thickness / 2 -
                              G_endcap_inner_thickness / 2));
        Transform3D slice_endtrap_trans_zn(
            Translation3D(0, endcap_thickness * tan(G_endcap_angle) / 2,
                          -endcap_thickness - slice_thickness / 2 +
                              G_endcap_inner_thickness / 2));

        G_endtrap_geo_zp = UnionSolid(slice_tmp_geo, slice_endtrap_geo,
                                      slice_endtrap_trans_zp);
        G_endtrap_geo_zn = UnionSolid(slice_tmp_geo, slice_endtrap_geo,
                                      slice_endtrap_trans_zn);
      }
    } else {
      Transform3D slice_endtrap_trans_zp(
          Translation3D(0, endcap_thickness * tan(G_endcap_angle) / 2,
                        endcap_thickness + slice_thickness / 2 -
                            G_endcap_inner_thickness / 2));
      Transform3D slice_endtrap_trans_zn(
          Translation3D(0, endcap_thickness * tan(G_endcap_angle) / 2,
                        -endcap_thickness - slice_thickness / 2 +
                            G_endcap_inner_thickness / 2));

      G_endtrap_geo_zp = UnionSolid(G_endtrap_geo_zp, slice_endtrap_geo,
                                    slice_endtrap_trans_zp);
      G_endtrap_geo_zn = UnionSolid(G_endtrap_geo_zn, slice_endtrap_geo,
                                    slice_endtrap_trans_zn);
    }

    endcap_thickness += slice_thickness;
  }

  // ************************************************************
  // #################### read barrel construction ##############
  for (int n = 1; n <= G_numsides; n++) { // construct subdetector
    double barrel_thickness = 0;

    string barrel_name = G_x_barrel.nameStr() + "_" + to_string(n);
    DetElement barrel_det(G_det, barrel_name, n);
    Trapezoid barrel_geo(G_barrel_rmin * tan(M_PI / G_numsides),
                         (G_barrel_rmin + G_barrel_total_thickness) *
                             tan(M_PI / G_numsides),
                         G_zhalf, G_zhalf, G_barrel_total_thickness / 2);
    Volume barrel_vol(barrel_name, barrel_geo, G_air);
    barrel_vol.setAttributes(lcdd, G_x_barrel.regionStr(),
                             G_x_barrel.limitsStr(), G_x_barrel.visStr());

    for (xml_coll_t i(G_x_barrel, _U(slice)); i; ++i) {
      xml_comp_t x_slice = i;

      double slice_length = 2 * G_zhalf;
      double slice_width =
          2 * (G_barrel_rmin + barrel_thickness) * tan(M_PI / G_numsides);
      double slice_thickness = x_slice.thickness();
      Material slice_material = lcdd.material(x_slice.materialStr());

      string slice_name = G_x_barrel.nameStr() + "_" + to_string(n) + "_" +
                          x_slice.nameStr() + "_" + to_string(x_slice.number());
      DetElement slice_det(barrel_det, slice_name, x_slice.number());
      Box slice_geo(slice_width / 2, slice_length / 2, slice_thickness / 2);
      Volume slice_vol;

      // cout<< " >> " << slice_name << " : " << slice_width << " : " <<
      // slice_length << endl;

      xml_coll_t c_strip(x_slice, _Unicode(strip));
      xml_comp_t x_strip(c_strip);
      if (x_strip) {
        int strip_num = x_strip.nmodules();
        double strip_gap = x_strip.gap();
        string strip_type = x_strip.typeStr();
        double strip_pos = 0;

        slice_vol = Volume(slice_name, slice_geo, G_air);
        slice_vol.setAttributes(lcdd, x_slice.regionStr(), x_slice.limitsStr(),
                                x_slice.visStr());

        for (int k = 1; k <= strip_num; k++) {
          double strip_width = 0;
          double strip_length = 0;

          string strip_name;
          if (strip_type == "x")
            strip_name = G_x_barrel.nameStr() + "_" + to_string(n) + "_" +
                         x_slice.nameStr() + "_" + to_string(x_slice.number()) +
                         "_stripX_" + to_string(k);
          else if (strip_type == "y")
            strip_name = G_x_barrel.nameStr() + "_" + to_string(n) + "_" +
                         x_slice.nameStr() + "_" + to_string(x_slice.number()) +
                         "_stripY_" + to_string(k);
          else
            cout << "WARNING: strips are not placed correctly!!!!!!!!!!"
                 << endl;
          DetElement strip_det(slice_det, strip_name, k);

          Box strip_tmp_geo;
          Transform3D strip_tmp_trans;

          if (strip_type == "x") {
            strip_width = (slice_width + strip_gap) / strip_num - strip_gap;
            strip_length = slice_length;
            strip_tmp_geo =
                Box(strip_width / 2, strip_length / 2, slice_thickness);
            strip_tmp_trans = Transform3D(Translation3D(
                -slice_width / 2 + strip_width / 2 + strip_pos, 0, 0));
            strip_pos += strip_width + strip_gap;
          } else if (strip_type == "y") {
            strip_width = slice_width;
            strip_length = (slice_length + strip_gap) / strip_num - strip_gap;
            strip_tmp_geo =
                Box(strip_width / 2, strip_length / 2, slice_thickness);
            strip_tmp_trans = Transform3D(Translation3D(
                0, -slice_length / 2 + strip_length / 2 + strip_pos, 0));
            strip_pos += strip_length + strip_gap;
          }

          IntersectionSolid strip_geo(slice_geo, strip_tmp_geo,
                                      strip_tmp_trans);
          Volume strip_vol(strip_name, strip_geo, slice_material);
          strip_vol.setAttributes(lcdd, x_slice.regionStr(),
                                  x_slice.limitsStr(), x_slice.visStr());
          if (x_strip.isSensitive()) {
            sens.setType("PS");
            strip_vol.setSensitiveDetector(sens);
          }

          Transform3D strip_trans(Translation3D(0, 0, 0));
          PlacedVolume strip_phv =
              slice_vol.placeVolume(strip_vol, strip_trans);
          if (x_slice.nameStr() == "OuterReadout")
            strip_phv.addPhysVolID("sublayer", 2).addPhysVolID("strip", k);
          else
            strip_phv.addPhysVolID("sublayer", 1).addPhysVolID("strip", k);
          if (x_slice.nameStr() == "InnerReadout")
            strip_phv.addPhysVolID("direction", 1);
          else
            strip_phv.addPhysVolID("direction", 2);
          strip_det.setPlacement(strip_phv);
        }
      } else {
        slice_vol = Volume(slice_name, slice_geo, slice_material);
        slice_vol.setAttributes(lcdd, x_slice.regionStr(), x_slice.limitsStr(),
                                x_slice.visStr());
        if (x_slice.isSensitive()) {
          sens.setType("RPC");
          slice_vol.setSensitiveDetector(sens);
        }
      }

      Transform3D slice_trans(Translation3D(0, 0,
                                            -G_barrel_total_thickness / 2 +
                                                barrel_thickness +
                                                slice_thickness / 2));
      PlacedVolume slice_phv = barrel_vol.placeVolume(slice_vol, slice_trans);
      slice_phv.addPhysVolID("layer", x_slice.number());

      slice_det.setPlacement(slice_phv);

      barrel_thickness += slice_thickness;
    }

    double barrel_radius = G_barrel_rmin + G_barrel_total_thickness / 2;
    double barrel_thetax = -M_PI / 2;
    double barrel_thetaz = 2 * (n - 1) * M_PI / G_numsides + M_PI / G_numsides;
    Transform3D barrel_trans(
        Rotation3D(cos(barrel_thetaz), -sin(barrel_thetaz) * cos(barrel_thetax),
                   sin(barrel_thetaz) * sin(barrel_thetax), sin(barrel_thetaz),
                   cos(barrel_thetaz) * cos(barrel_thetax),
                   -cos(barrel_thetaz) * sin(barrel_thetax), 0,
                   sin(barrel_thetax), cos(barrel_thetax)),
        Translation3D(-barrel_radius * sin(barrel_thetaz),
                      barrel_radius * cos(barrel_thetaz), 0));
    PlacedVolume barrel_phv =
        G_envelope_vol.placeVolume(barrel_vol, barrel_trans);
    barrel_phv.addPhysVolID("part", 0).addPhysVolID("sector", n);

    barrel_det.setPlacement(barrel_phv);
  }

  // ************************************************************
  // #################### read endcap construction ##############
  for (int n = 1; n <= G_numsides; n++) { // construct subdetector
    double endcap_thickness = 0;

    string endcap_name_zp = G_x_endcap.nameStr() + "_" + to_string(n);
    string endcap_name_zn =
        G_x_endcap.nameStr() + "_" + to_string(n + G_numsides);

    DetElement endcap_det_zp(G_det, endcap_name_zp, n);
    DetElement endcap_det_zn(G_det, endcap_name_zn, n + G_numsides);
    Volume endcap_vol_zp(endcap_name_zp, G_endtrap_geo_zp, G_air);
    Volume endcap_vol_zn(endcap_name_zn, G_endtrap_geo_zn, G_air);
    endcap_vol_zp.setAttributes(lcdd, G_x_endcap.regionStr(),
                                G_x_endcap.limitsStr(), G_x_endcap.visStr());
    endcap_vol_zn.setAttributes(lcdd, G_x_endcap.regionStr(),
                                G_x_endcap.limitsStr(), G_x_endcap.visStr());

    for (xml_coll_t i(G_x_endcap, _U(slice)); i; ++i) {
      xml_comp_t x_slice = i;

      double slice_inner_width = 2 * (G_zhalf + endcap_thickness) *
                                 tan(G_endcap_angle) * tan(M_PI / G_numsides);
      double slice_outer_width = 2 *
                                 (G_barrel_rmin + G_barrel_total_thickness) *
                                 tan(M_PI / G_numsides);
      double slice_length = G_barrel_rmin + G_barrel_total_thickness -
                            (G_zhalf + endcap_thickness) * tan(G_endcap_angle);
      double slice_thickness = x_slice.thickness();
      Material slice_material = lcdd.material(x_slice.materialStr());

      string slice_name_zp = G_x_endcap.nameStr() + "_" + to_string(n) + "_" +
                             x_slice.nameStr() + "_" +
                             to_string(x_slice.number());
      string slice_name_zn =
          G_x_endcap.nameStr() + "_" + to_string(n + G_numsides) + "_" +
          x_slice.nameStr() + "_" + to_string(x_slice.number());

      DetElement slice_det_zp(endcap_det_zp, slice_name_zp, x_slice.number());
      DetElement slice_det_zn(endcap_det_zn, slice_name_zn, x_slice.number());
      Trap slice_geo(slice_thickness / 2, 0, 0, slice_length / 2,
                     slice_inner_width / 2, slice_outer_width / 2, 0,
                     slice_length / 2, slice_inner_width / 2,
                     slice_outer_width / 2, 0);

      Volume slice_vol_zp;
      Volume slice_vol_zn;

      // cout<< " >> " << slice_name_zp << " : " << slice_outer_width << " : "
      // << slice_length << endl;

      xml_coll_t c_strip(x_slice, _Unicode(strip));
      xml_comp_t x_strip(c_strip);
      if (x_strip) {
        int strip_num = x_strip.nmodules();
        double strip_gap = x_strip.gap();
        string strip_type = x_strip.typeStr();
        double strip_pos = 0;

        slice_vol_zp = Volume(slice_name_zp, slice_geo, G_air);
        slice_vol_zn = Volume(slice_name_zn, slice_geo, G_air);
        slice_vol_zp.setAttributes(lcdd, x_slice.regionStr(),
                                   x_slice.limitsStr(), x_slice.visStr());
        slice_vol_zn.setAttributes(lcdd, x_slice.regionStr(),
                                   x_slice.limitsStr(), x_slice.visStr());

        for (int k = 1; k <= strip_num; k++) {
          double strip_width = 0;
          double strip_length = 0;

          string strip_name_zp;
          string strip_name_zn;
          if (strip_type == "x") {
            strip_name_zp = G_x_endcap.nameStr() + "_" + to_string(n) + "_" +
                            x_slice.nameStr() + "_" +
                            to_string(x_slice.number()) + "_stripX_" +
                            to_string(k);
            strip_name_zn =
                G_x_endcap.nameStr() + "_" + to_string(n + G_numsides) + "_" +
                x_slice.nameStr() + "_" + to_string(x_slice.number()) +
                "_stripX_" + to_string(k);
          } else if (strip_type == "y") {
            strip_name_zp = G_x_endcap.nameStr() + "_" + to_string(n) + "_" +
                            x_slice.nameStr() + "_" +
                            to_string(x_slice.number()) + "_stripY_" +
                            to_string(k);
            strip_name_zn =
                G_x_endcap.nameStr() + "_" + to_string(n + G_numsides) + "_" +
                x_slice.nameStr() + "_" + to_string(x_slice.number()) +
                "_stripY_" + to_string(k);
          } else
            cout << "WARNING: strips are not placed correctly!!!!!!!!!!"
                 << endl;
          DetElement strip_det_zp(slice_det_zp, strip_name_zp, k);
          DetElement strip_det_zn(slice_det_zn, strip_name_zn, k);

          Box strip_tmp_geo;
          Transform3D strip_tmp_trans;

          if (strip_type == "x") {
            strip_width =
                (slice_outer_width + strip_gap) / strip_num - strip_gap;
            strip_length = slice_length;
            strip_tmp_geo =
                Box(strip_width / 2, strip_length / 2, slice_thickness);
            strip_tmp_trans = Transform3D(Translation3D(
                -slice_outer_width / 2 + strip_width / 2 + strip_pos, 0, 0));
            strip_pos += strip_width + strip_gap;
          } else if (strip_type == "y") {
            strip_width = slice_outer_width;
            strip_length = (slice_length + strip_gap) / strip_num - strip_gap;
            strip_tmp_geo =
                Box(strip_width / 2, strip_length / 2, slice_thickness);
            strip_tmp_trans = Transform3D(Translation3D(
                0, -slice_length / 2 + strip_length / 2 + strip_pos, 0));
            strip_pos += strip_length + strip_gap;
          }

          IntersectionSolid strip_geo(slice_geo, strip_tmp_geo,
                                      strip_tmp_trans);
          Volume strip_vol_zp(strip_name_zp, strip_geo, slice_material);
          Volume strip_vol_zn(strip_name_zn, strip_geo, slice_material);
          strip_vol_zp.setAttributes(lcdd, x_strip.regionStr(),
                                     x_strip.limitsStr(), x_slice.visStr());
          strip_vol_zn.setAttributes(lcdd, x_strip.regionStr(),
                                     x_strip.limitsStr(), x_slice.visStr());
          if (x_strip.isSensitive()) {
            sens.setType("PS");
            strip_vol_zp.setSensitiveDetector(sens);
            strip_vol_zn.setSensitiveDetector(sens);
          }

          Transform3D strip_trans(Translation3D(0, 0, 0));
          PlacedVolume strip_phv_zp =
              slice_vol_zp.placeVolume(strip_vol_zp, strip_trans);
          PlacedVolume strip_phv_zn =
              slice_vol_zn.placeVolume(strip_vol_zn, strip_trans);
          if (x_slice.nameStr() == "OuterReadout" ||
              x_slice.nameStr() == "PSY") {
            strip_phv_zp.addPhysVolID("sublayer", 2).addPhysVolID("strip", k);
            strip_phv_zn.addPhysVolID("sublayer", 2).addPhysVolID("strip", k);
          } else {
            strip_phv_zp.addPhysVolID("sublayer", 1).addPhysVolID("strip", k);
            strip_phv_zn.addPhysVolID("sublayer", 1).addPhysVolID("strip", k);
          }
          if (x_slice.nameStr() == "InnerReadout" ||
              x_slice.nameStr() == "PSX") {
            strip_phv_zp.addPhysVolID("direction", 1);
            strip_phv_zn.addPhysVolID("direction", 1);
          } else {
            strip_phv_zp.addPhysVolID("direction", 2);
            strip_phv_zn.addPhysVolID("direction", 2);
          }
          strip_det_zp.setPlacement(strip_phv_zp);
          strip_det_zn.setPlacement(strip_phv_zn);
        }
      } else {
        slice_vol_zp = Volume(slice_name_zp, slice_geo, slice_material);
        slice_vol_zn = Volume(slice_name_zn, slice_geo, slice_material);
        slice_vol_zp.setAttributes(lcdd, x_slice.regionStr(),
                                   x_slice.limitsStr(), x_slice.visStr());
        slice_vol_zn.setAttributes(lcdd, x_slice.regionStr(),
                                   x_slice.limitsStr(), x_slice.visStr());
        if (x_slice.isSensitive()) {
          sens.setType("RPC");
          slice_vol_zp.setSensitiveDetector(sens);
          slice_vol_zn.setSensitiveDetector(sens);
        }
      }

      Transform3D slice_trans_zp(
          Translation3D(0, endcap_thickness * tan(G_endcap_angle) / 2,
                        endcap_thickness + slice_thickness / 2 -
                            G_endcap_inner_thickness / 2));
      Transform3D slice_trans_zn(
          Translation3D(0, endcap_thickness * tan(G_endcap_angle) / 2,
                        -endcap_thickness - slice_thickness / 2 +
                            G_endcap_inner_thickness / 2));
      PlacedVolume slice_phv_zp =
          endcap_vol_zp.placeVolume(slice_vol_zp, slice_trans_zp);
      PlacedVolume slice_phv_zn =
          endcap_vol_zn.placeVolume(slice_vol_zn, slice_trans_zn);
      slice_phv_zp.addPhysVolID("layer", x_slice.number());
      slice_phv_zn.addPhysVolID("layer", x_slice.number());

      slice_det_zp.setPlacement(slice_phv_zp);
      slice_det_zn.setPlacement(slice_phv_zn);

      endcap_thickness += slice_thickness;
    }

    double endcap_radius = (G_barrel_rmin + G_barrel_total_thickness +
                            G_zhalf * tan(G_endcap_angle)) /
                           2;
    double endcap_thetaz = 2 * (n - 1) * M_PI / G_numsides + M_PI / G_numsides;
    Transform3D endcap_trans_zp(
        RotationZ(endcap_thetaz),
        Translation3D(-endcap_radius * sin(endcap_thetaz),
                      endcap_radius * cos(endcap_thetaz),
                      G_zhalf + G_endcap_inner_thickness / 2));
    Transform3D endcap_trans_zn(
        RotationZ(endcap_thetaz),
        Translation3D(-endcap_radius * sin(endcap_thetaz),
                      endcap_radius * cos(endcap_thetaz),
                      -G_zhalf - G_endcap_inner_thickness / 2));
    PlacedVolume endcap_phv_zp =
        G_envelope_vol.placeVolume(endcap_vol_zp, endcap_trans_zp);
    PlacedVolume endcap_phv_zn =
        G_envelope_vol.placeVolume(endcap_vol_zn, endcap_trans_zn);
    endcap_phv_zp.addPhysVolID("part", 1).addPhysVolID("sector", n);
    endcap_phv_zn.addPhysVolID("part", 2).addPhysVolID("sector", n);

    endcap_det_zp.setPlacement(endcap_phv_zp);
    endcap_det_zn.setPlacement(endcap_phv_zn);
  }

  // ************************************************************
  double env_thetaz = -M_PI / G_numsides;
  Transform3D env_trans(RotationZ(env_thetaz), Translation3D(0, 0, 0));
  PlacedVolume env_phv = G_motherVol.placeVolume(G_envelope_vol, env_trans);
  env_phv.addPhysVolID("system", G_x_det.id());
  G_det.setPlacement(env_phv);

  return G_det;
}

DECLARE_DETELEMENT(DetectorMUC_02, create_detector)
