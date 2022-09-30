/**********************************************************
 * Author        : Hang Zhou
 * Email         : zh25@mail.ustc.edu.cn
 * Filename      : Tracker_geo.cpp
 * Description   : for Main drift chamber construction
 * *******************************************************/
#include "DD4hep/DetFactoryHelper.h"

using namespace dd4hep;
using namespace std;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  DetElement DriftChamber(name, x_det.id());
  PlacedVolume pv;

  xml_dim_t dim = x_det.dimensions();
  double rmin = dim.inner_r();
  double rmax = dim.outer_r();
  double zmax = dim.length();
  Tube TubeMDC(rmin, rmax, zmax / 2.0);

  Volume SubSolidMDC_vol(name + "_Gas", TubeMDC,
                         description.material("He/C4H10(90/10)"));
  SubSolidMDC_vol.setVisAttributes(description.visAttributes("SupportVis"));
  sens.setType("tracker");
  // SubSolidMDC_vol.setSensitiveDetector(sens);

  double D, r1, r2, alpha1, alpha2;
  double l, alpha, r_cen;
  int copy_no = 1;
  int i = 0;
  for (xml_coll_t c(e, _U(layer)); c; ++c, i++) {
    xml_comp_t x_layer(c);
    // read tags
    int ncell = x_layer.attr<int>(_Unicode(ncell));
    int layer_id = x_layer.id();
    double radius = x_layer.attr<double>(_Unicode(radius));
    double dr = x_layer.attr<double>(_Unicode(dr));
    double length = x_layer.attr<double>(_Unicode(length));
    double twist_angle = x_layer.attr<double>(_Unicode(twisted));
    int edge = x_layer.attr<int>(_Unicode(edge));
    int uv = x_layer.attr<int>(_Unicode(uvflag));
    std::string vislayer = x_layer.attr<std::string>(_Unicode(vis));
    const double segment_angle = 2 * M_PI / double(ncell);

    std::string layername = x_layer.nameStr() + _toString(layer_id, "_layer%d");

    D = (radius)*sin(twist_angle / 2) * 2;
    r1 = (radius)*cos(twist_angle / 2);
    alpha1 = atan(D / length);

    D = (radius + dr) * sin(twist_angle / 2) * 2;
    r2 = (radius + dr) * cos(twist_angle / 2);
    alpha2 = atan(D / length);
    /*
                    Assembly layer_vol( layername );
                    DetElement layerDE(DriftChamber,layername,layer_id);
                    // put the layer to the whole assembly
                    // pv = assembly.placeVolume(  layer_assembly );
                    pv = SubSolidMDC_vol.placeVolume(  layer_vol );
                    pv.addPhysVolID("layer", layer_id);
                    layerDE.setPlacement(pv);
    */

    Hyperboloid layer(r1, alpha1, r2, alpha2, length / 2.);
    Volume layer_vol(layername, layer, description.material("He/C4H10(90/10)"));
    layer_vol.setVisAttributes(description.visAttributes(vislayer));

    // DetElement layerDE(DriftChamber,layername,layer_id);
    pv = SubSolidMDC_vol.placeVolume(layer_vol);
    // pv.addPhysVolID("layer", layer_id);
    // layerDE.setPlacement(pv);

    xml_comp_t x_sense = x_layer.child(_Unicode(sense));
    xml_comp_t x_field = x_layer.child(_Unicode(field));

    // for field wires
    double fthick = x_field.attr<double>(_Unicode(cladthick));
    double radius_f = x_field.attr<double>(_Unicode(radius));
    std::string fmaterial_f = x_field.attr<std::string>(_Unicode(material));
    std::string fvisattr_f = x_field.attr<std::string>(_Unicode(vis));
    std::string fmaterial_clad = x_field.attr<std::string>(_Unicode(clad));
    Tube fieldwire_clad(0., fthick + radius_f, length / 2.);
    Volume wire_f("field_wire", fieldwire_clad,
                  description.material(fmaterial_clad));
    Tube fieldwire(0., radius_f, length / 2.);
    Volume wiref("fieldwire", fieldwire, description.material(fmaterial_f));
    pv = wire_f.placeVolume(wiref);

    Tube fieldwire_clad1(0., fthick + radius_f, length / 2.,
                         -M_PI + segment_angle / 2., 0 + segment_angle / 2.);
    Volume wire_f1("field_wire", fieldwire_clad1,
                   description.material(fmaterial_clad));
    Tube fieldwire1(0., radius_f, length / 2., -M_PI + segment_angle / 2.,
                    0 + segment_angle / 2.);
    Volume wiref1("fieldwire", fieldwire1, description.material(fmaterial_f));
    pv = wire_f1.placeVolume(wiref1);

    Tube fieldwire_clad2(0., fthick + radius_f, length / 2.,
                         0 - segment_angle / 2., M_PI - segment_angle / 2.);
    Volume wire_f2("field_wire", fieldwire_clad2,
                   description.material(fmaterial_clad));
    Tube fieldwire2(0., radius_f, length / 2., 0 - segment_angle / 2.,
                    M_PI - segment_angle / 2.);
    Volume wiref2("fieldwire", fieldwire2, description.material(fmaterial_f));
    pv = wire_f2.placeVolume(wiref2);

    double sthick = x_sense.attr<double>(_Unicode(cladthick));
    double radius_s = x_sense.attr<double>(_Unicode(radius));
    std::string smaterial_s = x_sense.attr<std::string>(_Unicode(material));
    std::string svisattr_s = x_sense.attr<std::string>(_Unicode(vis));
    std::string smaterial_clad = x_sense.attr<std::string>(_Unicode(clad));
    Tube sensewire_clad(0., sthick + radius_s, length / 2.);
    Volume wire_s("sense_wire", sensewire_clad,
                  description.material(smaterial_clad));
    Tube sensewire(0., radius_s, length / 2.);
    Volume wires("sensewire", sensewire, description.material(smaterial_s));
    pv = wire_s.placeVolume(wires);

    if (twist_angle == 0) {
      Tube cell(radius, radius + dr, length / 2., -segment_angle / 2.,
                segment_angle / 2.);
      string cellname = "layer_" + to_string(i) + "axialcell";
      // Volume
      // cell_vol("axialcell",cell,description.material("He/C4H10(90/10)"));
      // /
      Volume cell_vol(cellname, cell, description.material("He/C4H10(90/10)"));
      cell_vol.setSensitiveDetector(sens);

      Position pos(radius + dr / 2., 0, 0);
      pv = cell_vol.placeVolume(wire_s, pos);

      Position pos1((radius + 0.375 * dr) * cos(segment_angle / 2.),
                    (radius + 0.375 * dr) * sin(segment_angle / 2.), 0);
      Position pos2((radius + 0.625 * dr) * cos(segment_angle / 2.),
                    (radius + 0.625 * dr) * sin(segment_angle / 2.), 0);
      Position pos3((radius + 0.375 * dr) * cos(-segment_angle / 2.),
                    (radius + 0.375 * dr) * sin(-segment_angle / 2.), 0);
      Position pos4((radius + 0.625 * dr) * cos(-segment_angle / 2.),
                    (radius + 0.625 * dr) * sin(-segment_angle / 2.), 0);

      pv = cell_vol.placeVolume(wire_f1, pos1);
      pv = cell_vol.placeVolume(wire_f1, pos2);
      pv = cell_vol.placeVolume(wire_f2, pos3);
      pv = cell_vol.placeVolume(wire_f2, pos4);

      Position pos5((radius + 1.2 * radius_f) *
                        cos(segment_angle / 2. - 1.2 * radius_f / radius),
                    (radius + 1.2 * radius_f) *
                        sin(segment_angle / 2. - 1.2 * radius_f / radius),
                    0);
      Position pos6(
          (radius + dr - 1.2 * radius_f) *
              cos(segment_angle / 2. - 1.2 * radius_f / (dr + radius)),
          (radius + dr - 1.2 * radius_f) *
              sin(segment_angle / 2. - 1.2 * radius_f / (dr + radius)),
          0);

      pv = cell_vol.placeVolume(wire_f, pos5);
      pv = cell_vol.placeVolume(wire_f, pos6);

      Position pos7((radius + 1.2 * radius_f), 0, 0);
      Position pos8((radius + dr - 1.2 * radius_f), 0, 0);

      pv = cell_vol.placeVolume(wire_f, pos7);
      pv = cell_vol.placeVolume(wire_f, pos8);

      if (uv == 0) {
        for (int j = 0; j < ncell; j++) {
          RotationZYX rot(segment_angle * double(j) * rad, 0, 0);
          //	pv= SubSolidMDC_vol.placeVolume( cell_vol,copy_no++,rot );
          pv = layer_vol.placeVolume(cell_vol, copy_no++, rot);
          cell_vol.setVisAttributes(description.visAttributes("CellVis"));
        }
      } else {
        for (int j = 0; j < ncell; j++) {
          RotationZYX rot(segment_angle * (double(j) + 0.25) * rad, 0, 0);
          // pv= SubSolidMDC_vol.placeVolume( cell_vol,rot );
          pv = layer_vol.placeVolume(cell_vol, copy_no++, rot);
          cell_vol.setVisAttributes(description.visAttributes("CellVis"));
        }
      }

    } else {

      if (uv == 1) {
        // TwistedTube
        // cell(-twist_angle*180/M_PI,radius,radius+dr,length/2.,segment_angle);
        TwistedTube cell(
            -twist_angle * 180. / M_PI, radius * cos(twist_angle / 2.),
            (radius + dr) * cos(twist_angle / 2.), length / 2., segment_angle);

        string cellname = "layer_" + to_string(i) + "stereocellv";
        Volume cell_vol(cellname, cell,
                        description.material("He/C4H10(90/10)"));

        cell_vol.setSensitiveDetector(sens);

        D = (radius + 0.375 * dr) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 0.375 * dr) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot1(
            0,
            asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                 length) *
                rad);
        Position pos1(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        // Position pos1(100*mm,100*mm,0);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot1, pos1));

        D = (radius + 0.625 * dr) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 0.625 * dr) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot2(
            0,
            asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                 length) *
                rad);
        Position pos2(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot2, pos2));

        D = (radius + 0.5 * dr) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 0.5 * dr) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot3(0, 0, atan(D * cos(0) / length) * rad);
        Position pos3(r_cen, 0., 0.);

        pv = cell_vol.placeVolume(wire_s, Transform3D(rot3, pos3));

        D = (radius + 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot4(0, 0, atan(D * cos(0) / length) * rad);
        Position pos4((r_cen), 0, 0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot4, pos4));

        D = ((radius + dr) - 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = ((radius + dr) - 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot5(0, 0, atan(D * cos(0) / length) * rad);
        Position pos5((r_cen), 0, 0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot5, pos5));

        D = (radius + 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot6(
            0,
            asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                 length) *
                rad);
        Position pos6(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot6, pos6));

        D = (radius + dr - 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = (radius + dr - 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot7(
            0,
            asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                 length) *
                rad);
        Position pos7(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot7, pos7));

        for (int j = 0; j < ncell; j++) {

          RotationZYX rot(segment_angle * (0.25 + double(j)) * rad, 0, 0);
          // pv= SubSolidMDC_vol.placeVolume( cell_vol,rot );
          pv = layer_vol.placeVolume(cell_vol, copy_no++, rot);
          cell_vol.setVisAttributes(description.visAttributes("CellVis"));
        }

        //	pv.addPhysVolID("cell", j);
      }

      else {
        TwistedTube cell(
            twist_angle * 180. / M_PI, radius * cos(twist_angle / 2.),
            (radius + dr) * cos(twist_angle / 2.), length / 2., segment_angle);

        D = (radius + 0.375 * dr) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 0.375 * dr) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        string cellname = "layer_" + to_string(i) + "stereocellu";
        Volume cell_vol(cellname, cell,
                        description.material("He/C4H10(90/10)"));

        cell_vol.setSensitiveDetector(sens);

        RotationZYX rot1(
            0,
            -asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            -atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                  length) *
                rad);
        Position pos1(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        // Position pos1(100*mm,100*mm,0);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot1, pos1));

        D = (radius + 0.625 * dr) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 0.625 * dr) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot2(
            0,
            -asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            -atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                  length) *
                rad);
        Position pos2(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot2, pos2));

        D = (radius + 0.5 * dr) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 0.5 * dr) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot3(0, 0, -atan(D * cos(0) / length) * rad);
        Position pos3(r_cen, 0, 0.);

        pv = cell_vol.placeVolume(wire_s, Transform3D(rot3, pos3));

        // Position pos1(100*mm,100*mm,0);

        D = (radius + 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot4(0, 0, -atan(D * cos(0) / length) * rad);
        Position pos4((r_cen), 0, 0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot4, pos4));

        D = ((radius + dr) - 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = ((radius + dr) - 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot5(0, 0, -atan(D * cos(0) / length) * rad);
        Position pos5((r_cen), 0, 0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot5, pos5));

        D = (radius + 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = (radius + 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot6(
            0,
            -asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            -atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                  length) *
                rad);
        Position pos6(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot6, pos6));

        D = (radius + dr - 1.5 * radius_f) * sin(twist_angle / 2) * 2;
        r_cen = (radius + dr - 1.5 * radius_f) * cos(twist_angle / 2);
        alpha = atan(D / length);
        l = sqrt(D * D + length * length);

        RotationZYX rot7(
            0,
            -asin(D * sin(segment_angle / 2. - radius_f / r_cen * 1.5) / l) *
                rad,
            -atan(D * cos(segment_angle / 2. - radius_f / r_cen * 1.5) /
                  length) *
                rad);
        Position pos7(r_cen * cos(segment_angle / 2. - radius_f / r_cen * 1.5),
                      r_cen * sin(segment_angle / 2. - radius_f / r_cen * 1.5),
                      0.);

        pv = cell_vol.placeVolume(wire_f, Transform3D(rot7, pos7));

        for (int j = 0; j < ncell; j++) {
          RotationZYX rot(segment_angle * double(j) * rad, 0, 0);
          //	pv= SubSolidMDC_vol.placeVolume( cell_vol,copy_no++,rot );
          pv = layer_vol.placeVolume(cell_vol, copy_no++, rot);
          cell_vol.setVisAttributes(description.visAttributes("CellVis"));
        }
      }
    }
  }

  Volume mother = description.pickMotherVolume(DriftChamber);
  pv = mother.placeVolume(SubSolidMDC_vol);
  pv.addPhysVolID("system", x_det.id());
  DriftChamber.setPlacement(pv);

  return DriftChamber;
}
DECLARE_DETELEMENT(DriftChamber_02, create_element)
