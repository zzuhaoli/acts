/**********************************************************
 * Author        :  Zhou Hang
 * Email         : zh25@mail.ustc.edu.cn
 * Filename      : Tracker_geo.cpp
 * Description   : for Main drift chamber construction
 * *******************************************************/
#include "DD4hep/DetFactoryHelper.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hep/DetType.h"
#include "RecGeoData.h"
#include <exception>

using namespace dd4hep;
using namespace std;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {
  std::cout<<"Calling detectorMDC_03.cpp" << std::endl; 
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  DetElement DriftChamber(name, x_det.id());
  DetType type( DetType::TRACKER);
  DriftChamber.setTypeFlag(type.to_ulong());
 
  PlacedVolume pv;

  xml_dim_t dim = x_det.dimensions();
  double rmin = dim.inner_r();
  double rmax = dim.outer_r();
  double zmax = dim.length();
  std::string gas = dim.attr<std::string>(_Unicode(gas));
  Tube TubeMDC(rmin, rmax, zmax / 2.0);
  Cone ConeMDC1(412.5 * mm, 0, 180 * mm, 0, 450 * mm);
  Cone ConeMDC2(412.5 * mm, 0, 450 * mm, 0, 180 * mm);
  SubtractionSolid StMDC1(TubeMDC, ConeMDC1, Position(0, 0, 987.5 * mm));
  SubtractionSolid StMDC(StMDC1, ConeMDC2, Position(0, 0, -987.5 * mm));

  Volume SubSolidMDC_vol(name + "_Gas", StMDC, description.material(gas));
  SubSolidMDC_vol.setVisAttributes(description.visAttributes("GasVis"));
  sens.setType("tracker");
  // SubSolidMDC_vol.setSensitiveDetector(sens);

  int copy_no = 0;

  const double normal_l = 2476. * mm;

  int i = 0;

  rec::MDCRecGeoData *MdcGeoData = new rec::MDCRecGeoData;
  MdcGeoData->normal_l = 2476.;

  xml_comp_t ce = x_det.child(_Unicode(cell));
  std::string viscell = ce.attr<std::string>(_Unicode(vis));
  for (xml_coll_t c(ce, _U(layer)); c; ++c, i++) {
    xml_comp_t x_layer(c);
    // read tags
    int ncell = x_layer.attr<int>(_Unicode(ncell));
    int layer_id = x_layer.id();
    double sr = x_layer.attr<double>(_Unicode(radius));
    double dr = x_layer.attr<double>(_Unicode(dr));
    double length = x_layer.attr<double>(_Unicode(length));
    double twist = x_layer.attr<double>(_Unicode(twisted));
    int edge = x_layer.attr<int>(_Unicode(edge));
    int uv = x_layer.attr<int>(_Unicode(uvflag));
    std::string vislayer = x_layer.attr<std::string>(_Unicode(vis));
    double lr1, lr2, angle1, angle2;
    double cr1, cr2, twangle; // cell(twistedtube or tube) parameters
    const double segment_angle = 2 * M_PI / double(ncell);

    std::string layername = x_layer.nameStr() + _toString(layer_id, "_layer%d");

    lr1 = sr * cos(twist * M_PI / double(ncell));
    lr2 = (sr + dr) * cos(twist * M_PI / double(ncell));
    angle1 = atan(2. * sr * sin(twist * M_PI / double(ncell)) / normal_l);
    angle2 =
        atan(2. * (sr + dr) * sin(twist * M_PI / double(ncell)) / normal_l);

    rec::MDCRecGeoData::LayerLayout thisLayer;

    thisLayer.SuperL_id = layer_id / 6;
    thisLayer.Layer_id = layer_id;
    thisLayer.CellPerLayer = ncell;
    thisLayer.inner_Radius = lr1;
    thisLayer.outer_Radius = lr2;
    thisLayer.Zlength = length;
    thisLayer.SwireR_z0 = (lr1 + lr2) / 2.;
    thisLayer.Twist = twist;
    thisLayer.Stereo = 0;
    thisLayer.Type = twist == 0 ? 0 : (twist > 0 ? 1 : -1);
    thisLayer.StartShift = uv * 0.5;
    thisLayer.LayerInSup = 6;

    Hyperboloid layer(lr1, angle1, lr2, angle2, length / 2.);
    Volume layer_vol(layername, layer, description.material(gas));
    layer_vol.setVisAttributes(description.visAttributes(vislayer));
    layer_vol.setRegion(description, x_det.regionStr());

    // cout<<"region name  "<<x_det.regionStr()<<endl;
    // cout<<"layer id "<<layer_id<<endl;

    pv = SubSolidMDC_vol.placeVolume(layer_vol);

    xml_comp_t x_sense = x_layer.child(_Unicode(sense));
    xml_comp_t x_field = x_layer.child(_Unicode(field));

    double sense_angle;
    double tran_angle;
    double field_angle;

    // for field wires
    double fthick = x_field.attr<double>(_Unicode(cladthick));
    double r_field = x_field.attr<double>(_Unicode(radius));
    std::string fmaterial_f = x_field.attr<std::string>(_Unicode(material));
    std::string fvisattr_f = x_field.attr<std::string>(_Unicode(vis));
    std::string fmaterial_clad = x_field.attr<std::string>(_Unicode(clad));
    std::string visfield = x_field.attr<std::string>(_Unicode(vis));

    Tube fieldwire_clad(0., fthick + r_field, length / 2.);
    Volume wire_f("field_wire", fieldwire_clad,
                  description.material(fmaterial_clad));
    Tube fieldwire(0., r_field, length / 2.);
    Volume wiref("fieldwire", fieldwire, description.material(fmaterial_f));
    pv = wire_f.placeVolume(wiref);
    wire_f.setVisAttributes(description.visAttributes(visfield));

    Tube fieldwire_clad1(0., fthick + r_field, length / 2., 0., M_PI);
    Volume wire_f1("field_wire_half", fieldwire_clad1,
                   description.material(fmaterial_clad));
    Tube fieldwire1(0., r_field, length / 2., 0., M_PI);
    Volume wiref1("fieldwire_half", fieldwire1,
                  description.material(fmaterial_f));
    pv = wire_f1.placeVolume(wiref1); // field wire half
    wire_f1.setVisAttributes(description.visAttributes(visfield));

    // sense_wire
    double sthick = x_sense.attr<double>(_Unicode(cladthick));
    double r_sense = x_sense.attr<double>(_Unicode(radius));
    std::string smaterial_s = x_sense.attr<std::string>(_Unicode(material));
    std::string svisattr_s = x_sense.attr<std::string>(_Unicode(vis));
    std::string smaterial_clad = x_sense.attr<std::string>(_Unicode(clad));
    std::string vissense = x_sense.attr<std::string>(_Unicode(vis));

    Tube sensewire_clad(0., sthick + r_sense, length / 2.);
    Volume wire_s("sense_wire", sensewire_clad,
                  description.material(smaterial_clad));
    Tube sensewire(0., r_sense, length / 2.);
    Volume wires("sensewire", sensewire, description.material(smaterial_s));
    pv = wire_s.placeVolume(wires);
    wire_f.setVisAttributes(description.visAttributes(vissense));

    if (twist == 0) {
      Tube cell(lr1, lr2, length / 2., -segment_angle / 2., segment_angle / 2.);
      string cellname = "layer_" + to_string(i) + "axialcell";
      // Volume cell_vol("axialcell",cell,description.material(gas));	/
      Volume cell_vol(cellname, cell, description.material(gas));
      cell_vol.setSensitiveDetector(sens);
      cell_vol.setVisAttributes(description.visAttributes(viscell));
      Position pos((lr1 + lr2) / 2., 0, 0);
      pv = cell_vol.placeVolume(wire_s, pos);

      RotationZYX rot1(-segment_angle / 2. * rad, 0, 0);
      Position pos1((lr1 + lr2) * cos(-segment_angle / 2. * rad) / 2.,
                    (lr1 + lr2) * sin(-segment_angle / 2. * rad) / 2., 0.);
      pv = cell_vol.placeVolume(wire_f1, Transform3D(rot1, pos1));

      RotationZYX rot2((M_PI + segment_angle / 2.) * rad, 0, 0);
      Position pos2((lr1 + lr2) * cos(segment_angle / 2. * rad) / 2.,
                    (lr1 + lr2) * sin(segment_angle / 2. * rad) / 2., 0.);
      pv = cell_vol.placeVolume(wire_f1, Transform3D(rot2, pos2));

      tran_angle = 1.1 * r_field / lr1;
      Position posa((lr1 + r_field + fthick) *
                        cos((segment_angle / 2. - tran_angle) * rad),
                    (lr1 + r_field + fthick) *
                        sin((segment_angle / 2. - tran_angle) * rad),
                    0.);
      Position posb(lr1 + r_field + fthick, 0, 0);
      pv = cell_vol.placeVolume(wire_f, posa);
      pv = cell_vol.placeVolume(wire_f, posb);

      if (edge == -1) {
        tran_angle = 1.1 * r_field / lr2;
        Position posc((lr2 - r_field - fthick) *
                          cos((segment_angle / 2. - tran_angle) * rad),
                      (lr2 - r_field - fthick) *
                          sin((segment_angle / 2. - tran_angle) * rad),
                      0.);
        Position posd(lr2 - r_field - fthick, 0, 0);
        pv = cell_vol.placeVolume(wire_f, posc);
        pv = cell_vol.placeVolume(wire_f, posd);
      }
      for (int j = 0; j < ncell; j++) {
        RotationZYX rot9;
        if (i % 2 == 0) {
          rot9 = RotationZYX(segment_angle * double(j) * rad, 0, 0);
        } else {
          rot9 = RotationZYX(segment_angle * (double(j) + 0.5) * rad, 0, 0);
        }
        pv = layer_vol.placeVolume(cell_vol, copy_no++, rot9);
        cell_vol.setVisAttributes(description.visAttributes("CellVis"));
      }
    }

    else {
      cr1 = lr1;
      cr2 = lr2;
      twangle = 2 * atan(tan(twist * M_PI / double(ncell)) * length / normal_l);
      sense_angle = -atan(2. * (sr + dr / 2.) *
                          sin(twist * M_PI / double(ncell)) / normal_l);

      thisLayer.Stereo = -sense_angle;
      /*
                      TwistedTube test(60*deg,500*mm,520*mm,800*mm,5*deg);
                      Volume test_vol("idi",test,description.material(gas));
                      test_vol.setVisAttributes(description.visAttributes("CellVis"));
                      SubSolidMDC_vol.placeVolume(test_vol);
*/

      TwistedTube cells(twangle * 180 / M_PI, cr1, cr2, length / 2.,
                        segment_angle);

      string cellname = "layer_" + to_string(i) + "stereocell";
      Volume cell_vol(cellname, cells, description.material(gas));

      cell_vol.setSensitiveDetector(sens);

      RotationZYX rot3(0, 0, sense_angle * rad);
      Position pos3((lr1 + lr2) / 2., 0., 0.);
      pv = cell_vol.placeVolume(wire_s, Transform3D(rot3, pos3));
      tran_angle = 1.1 * r_field / ((lr1 + lr2) / 2.);

      RotationZYX rot4(
          0, asin(sin(sense_angle) * sin(segment_angle / 2. - tran_angle)),
          atan(tan(sense_angle) * cos(segment_angle / 2. - tran_angle)));
      Position pos4(
          (lr1 + lr2) * cos((segment_angle / 2. - tran_angle) * rad) / 2.,
          (lr1 + lr2) * sin((segment_angle / 2. - tran_angle) * rad) / 2., 0.);
      pv = cell_vol.placeVolume(wire_f, Transform3D(rot4, pos4));

      field_angle = -atan(2. * (sr + r_field + fthick) *
                          sin(twist * M_PI / double(ncell)) / normal_l);
      tran_angle = 1.1 * r_field / (lr1 + r_field + fthick);
      RotationZYX rot5(
          0, asin(sin(field_angle) * sin(segment_angle / 2. - tran_angle)),
          atan(tan(field_angle) * cos(segment_angle / 2. - tran_angle)));
      Position pos5(
          (lr1 + 1.1 * r_field) * cos((segment_angle / 2. - tran_angle) * rad),
          (lr1 + 1.1 * r_field) * sin((segment_angle / 2. - tran_angle) * rad),
          0.);
      pv = cell_vol.placeVolume(wire_f, Transform3D(rot5, pos5));

      RotationZYX rot6(0, 0., atan(tan(field_angle)));
      Position pos6((lr1 + r_field + fthick), 0., 0.);
      pv = cell_vol.placeVolume(wire_f, Transform3D(rot6, pos6));

      if (edge == -1) {
        field_angle = -atan(2. * (sr + dr - r_field - fthick) *
                            sin(twist * M_PI / double(ncell)) / normal_l);
        tran_angle = 1.1 * r_field / (lr2 - r_field - fthick);
        RotationZYX rot7(
            0, asin(sin(field_angle) * sin(segment_angle / 2. - tran_angle)),
            atan(tan(field_angle) * cos(segment_angle / 2. - tran_angle)));
        Position pos7((lr2 - 1.1 * r_field) *
                          cos((segment_angle / 2. - tran_angle) * rad),
                      (lr2 - 1.1 * r_field) *
                          sin((segment_angle / 2. - tran_angle) * rad),
                      0.);
        pv = cell_vol.placeVolume(wire_f, Transform3D(rot7, pos7));

        RotationZYX rot8(0, 0., atan(tan(field_angle)));
        Position pos8((lr2 - r_field - fthick), 0., 0.);
        pv = cell_vol.placeVolume(wire_f, Transform3D(rot8, pos8));
      }

      for (int j = 0; j < ncell; j++) {
        RotationZYX rot9;
        if (i % 2 == 0) {
          rot9 = RotationZYX(segment_angle * double(j) * rad, 0, 0);
        } else {
          rot9 = RotationZYX(segment_angle * (double(j) + 0.5) * rad, 0, 0);
        }
        pv = layer_vol.placeVolume(cell_vol, copy_no++, rot9);
        cell_vol.setVisAttributes(description.visAttributes("CellVis"));
      }
    }

    MdcGeoData->layers.push_back(thisLayer);
  }

  // end cap support, divide to two parts
  xml_comp_t endcap = x_det.child(_Unicode(endcap));
  // inner region, a Cone shape
  xml_comp_t cone = endcap.child(_Unicode(cone1));
  xml_comp_t cap = endcap.child(_Unicode(end3));

  std::string cone_ma = cone.attr<std::string>(_Unicode(material));
  std::string viscone = cone.attr<std::string>(_Unicode(vis));
  for (xml_coll_t c(cone, _U(layer)); c; ++c) {
    xml_comp_t x_segtube(c);
    double cone_rmin = x_segtube.attr<double>(_Unicode(rmin));
    double cone_rmax = x_segtube.attr<double>(_Unicode(rmax));
    double z = x_segtube.attr<double>(_Unicode(th));
    double zpos = x_segtube.attr<double>(_Unicode(zpos));
    Tube coneTube(cone_rmin, cone_rmax, z / 2.0);
    Volume coneVol(name + "Tube_support", coneTube,
                   description.material(cone_ma));
    coneVol.setVisAttributes(description.visAttributes(viscone));
    SubSolidMDC_vol.placeVolume(
        coneVol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., zpos)));
    SubSolidMDC_vol.placeVolume(
        coneVol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., -zpos)));
  }

  std::string end_ma = cap.attr<std::string>(_Unicode(material));
  std::string viscap = cap.attr<std::string>(_Unicode(vis));
  for (xml_coll_t c(cap, _U(layer)); c; ++c) {
    xml_comp_t x_segtube(c);
    double cap_rmin = x_segtube.attr<double>(_Unicode(rmin));
    double cap_rmax = x_segtube.attr<double>(_Unicode(rmax));
    double z = x_segtube.attr<double>(_Unicode(th));
    double zpos = x_segtube.attr<double>(_Unicode(zpos));
    Tube endcapTube(cap_rmin, cap_rmax, z / 2.0);
    Volume endcapVol(name + "Tube_support", endcapTube,
                     description.material(end_ma));
    endcapVol.setVisAttributes(description.visAttributes(viscap));
    SubSolidMDC_vol.placeVolume(
        endcapVol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., zpos)));
    SubSolidMDC_vol.placeVolume(
        endcapVol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., -zpos)));
  }

  xml_comp_t barrel = x_det.child(_Unicode(barrel));
  xml_comp_t innertube = barrel.child(_Unicode(innertube));
  std::string visattr = barrel.attr<std::string>(_Unicode(vis));
  double BarrelTubeLength = innertube.attr<double>(_Unicode(length));
  double BarrelTubeRmin = innertube.attr<double>(_Unicode(rmin));
  double BarrelTubeRmax = innertube.attr<double>(_Unicode(rmax));
  double film0_th = innertube.attr<double>(_Unicode(film0_th));
  double film1_th = innertube.attr<double>(_Unicode(film1_th));
  std::string iwall_ma = innertube.attr<std::string>(_Unicode(wall_material));
  std::string film_ma = innertube.attr<std::string>(_Unicode(film_ma));
  Tube TubeSolidInnerBarrel_Al1(BarrelTubeRmin - film0_th, BarrelTubeRmin,
                                BarrelTubeLength / 2.0);
  Volume TubeSolidInnerBarrel_Al1Vol(name + "barrel_support",
                                     TubeSolidInnerBarrel_Al1,
                                     description.material(film_ma));
  SubSolidMDC_vol.placeVolume(
      TubeSolidInnerBarrel_Al1Vol,
      Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));
  Tube TubeSolidInnerBarrel_Al2(BarrelTubeRmax, BarrelTubeRmax + film1_th,
                                BarrelTubeLength / 2.0);
  Volume TubeSolidInnerBarrel_Al2Vol(name + "barrel_support",
                                     TubeSolidInnerBarrel_Al2,
                                     description.material(film_ma));
  SubSolidMDC_vol.placeVolume(
      TubeSolidInnerBarrel_Al2Vol,
      Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));
  Tube TubeSolidInnerBarrel(BarrelTubeRmin, BarrelTubeRmax,
                            BarrelTubeLength / 2.0);
  Volume TubeSolidInnerBarrelVol(name + "barrel_support", TubeSolidInnerBarrel,
                                 description.material(iwall_ma));
  TubeSolidInnerBarrelVol.setVisAttributes(description.visAttributes(visattr));
  SubSolidMDC_vol.placeVolume(
      TubeSolidInnerBarrelVol,
      Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));

  // barrel part, inner and outer
  xml_comp_t outertube = barrel.child(_Unicode(outertube));
  BarrelTubeLength = outertube.attr<double>(_Unicode(length));
  BarrelTubeRmin = outertube.attr<double>(_Unicode(rmin));
  BarrelTubeRmax = outertube.attr<double>(_Unicode(rmax));
  double film2_th = outertube.attr<double>(_Unicode(film2_th));
  double film3_th = outertube.attr<double>(_Unicode(film3_th));
  std::string owall_ma = outertube.attr<std::string>(_Unicode(wall_material));
  film_ma = outertube.attr<std::string>(_Unicode(film_ma));
  Tube TubeSolidouterBarrel_Al1(BarrelTubeRmin - film2_th, BarrelTubeRmin,
                                BarrelTubeLength / 2.0);
  Volume TubeSolidouterBarrel_Al1Vol(name + "barrel_support",
                                     TubeSolidouterBarrel_Al1,
                                     description.material(film_ma));
  SubSolidMDC_vol.placeVolume(
      TubeSolidouterBarrel_Al1Vol,
      Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));
  Tube TubeSolidouterBarrel_Al2(BarrelTubeRmax, BarrelTubeRmax + film3_th,
                                BarrelTubeLength / 2.0);
  Volume TubeSolidouterBarrel_Al2Vol(name + "barrel_support",
                                     TubeSolidouterBarrel_Al2,
                                     description.material(film_ma));
  SubSolidMDC_vol.placeVolume(
      TubeSolidouterBarrel_Al2Vol,
      Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));
  Tube TubeSolidouterBarrel(BarrelTubeRmin, BarrelTubeRmax,
                            BarrelTubeLength / 2.0);
  Volume TubeSolidouterBarrelVol(name + "barrel_support", TubeSolidouterBarrel,
                                 description.material(iwall_ma));
  TubeSolidouterBarrelVol.setVisAttributes(description.visAttributes(visattr));
  SubSolidMDC_vol.placeVolume(
      TubeSolidouterBarrelVol,
      Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));

  xml_comp_t board = x_det.child(_Unicode(board));
  double board_x = board.attr<double>(_Unicode(x));
  double board_y = board.attr<double>(_Unicode(y));
  double board_z = board.attr<double>(_Unicode(z));
  double Cu_x = board.attr<double>(_Unicode(Cu_x));
  double Cu_y = board.attr<double>(_Unicode(Cu_y));
  double Cu_z = board.attr<double>(_Unicode(Cu_z));
  std::string board_ma = board.attr<std::string>(_Unicode(material));
  std::string copper_ma = board.attr<std::string>(_Unicode(Cu_ma));
  std::string visboard = board.attr<std::string>(_Unicode(vis));
  Box board_box((board_x + Cu_x) / 2., board_y / 2., board_z / 2.);
  Volume boardVol(name + "board_box", board_box, description.material("Air"));
  boardVol.setVisAttributes(description.visAttributes(visboard));
  Box base_box(board_x / 2., board_y / 2., board_z / 2.);
  Volume baseVol(name + "base_box", base_box, description.material(board_ma));
  Box copper_box(Cu_x / 2., Cu_y / 2., Cu_z / 2.);
  Volume copperVol(name + "copper_box", copper_box,
                   description.material(copper_ma));
  boardVol.placeVolume(
      baseVol, Transform3D(RotationZYX(0, 0., 0), Position(-Cu_x / 2., 0., 0)));
  boardVol.placeVolume(
      copperVol, Transform3D(RotationZYX(0, 0., 0),
                             Position(board_x / 2., (board_y - Cu_y) / 2., 0)));

  for (xml_coll_t c(board, _U(layer)); c; ++c) {
    xml_comp_t x_slayer(c);
    int n_board = x_slayer.attr<int>(_Unicode(nboard));
    double r_place = x_slayer.attr<double>(_Unicode(rpos));
    double z_place = x_slayer.attr<double>(_Unicode(zpos));
    for (int k = 0; k < n_board; k++) {
      SubSolidMDC_vol.placeVolume(
          boardVol,
          Transform3D(
              RotationZYX(2 * M_PI * double(k) / double(n_board), 0., 0),
              Position(r_place * cos(2 * M_PI * double(k) / double(n_board)),
                       r_place * sin(2 * M_PI * double(k) / double(n_board)),
                       z_place + board_z / 2.)));
      SubSolidMDC_vol.placeVolume(
          boardVol,
          Transform3D(
              RotationZYX(2 * M_PI * double(k) / double(n_board), 0., 0),
              Position(r_place * cos(2 * M_PI * double(k) / double(n_board)),
                       r_place * sin(2 * M_PI * double(k) / double(n_board)),
                       -z_place - board_z / 2.)));
    }
  }

  DriftChamber.addExtension<rec::MDCRecGeoData>(MdcGeoData);
  std::cout << "mytest " << MdcGeoData->layers.size() << endl;

  Volume mother = description.pickMotherVolume(DriftChamber);
  pv = mother.placeVolume(SubSolidMDC_vol);
  pv.addPhysVolID("system", x_det.id());
  DriftChamber.setPlacement(pv);

  return DriftChamber;
}
DECLARE_DETELEMENT(DriftChamber02, create_element)
