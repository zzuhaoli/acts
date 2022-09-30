/**********************************************************
 * Author        : Dong Liu
 * Email         : dliu13@ustc.edu.cn
 * Last modified : 2018-11-29 16:06
 * Filename      : Tracker_geo.cpp
 * Description   : for Main drift chamber construction
 * *******************************************************/
#include "DD4hep/DetFactoryHelper.h"

using namespace dd4hep;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  // Assembly assembly( name+"_assembly" );
  DetElement DriftLayer(name, x_det.id());
  PlacedVolume pv;

  xml_dim_t dim = x_det.dimensions();
  double rmin = dim.inner_r();
  double rmax = dim.outer_r();
  double mdcLength = dim.length();
  double zmax = mdcLength / 2;
  double z1 = dim.z1();
  double z2 = dim.z2();
  std::string gasvis = dim.visStr();
  Tube TubeMDC(rmin, rmax, zmax);

  double rcnt = dim.attr<double>(_Unicode(rcnt));
  Cone ConeInnerF((z2 - z1) / 2, rmin, rmin, rmin, rcnt);
  Cone ConeOuterF((zmax - z2) / 2, rmin, rcnt, rmin, rmax);
  Cone ConeInnerB((z2 - z1) / 2, rmin, rcnt, rmin, rmin);
  Cone ConeOuterB((zmax - z2) / 2, rmin, rmax, rmin, rcnt);

  SubtractionSolid TubeSub1(TubeMDC, ConeOuterF,
                            Position(0, 0, (zmax + z2) / 2));
  SubtractionSolid TubeSub2(TubeSub1, ConeInnerF,
                            Position(0, 0, (z1 + z2) / 2));
  SubtractionSolid TubeSub3(TubeSub2, ConeInnerB,
                            Position(0, 0, -(z1 + z2) / 2));
  SubtractionSolid TubeSub4(TubeSub3, ConeOuterB,
                            Position(0, 0, -(z2 + zmax) / 2));

  Volume SubSolidMDC_vol(name + "_Gas", TubeSub4,
                         description.material("G4_He"));
  sens.setType("tracker");
  SubSolidMDC_vol.setSensitiveDetector(sens);
  SubSolidMDC_vol.setVisAttributes(description.visAttributes(gasvis));

  // loop over layers
  double positioncnt = 0;
  double lengthcnt = 0;
  for (xml_coll_t c(e, _U(layer)); c; ++c) {
    // std::cout<<"ccccc "<< c <<std::endl;
    //  get component of a layer
    xml_comp_t x_layer(c);
    // read tags
    int nwire = x_layer.attr<int>(_Unicode(nwire));
    int layer_id = x_layer.id();
    double position = x_layer.attr<double>(_Unicode(position));
    double dpos = x_layer.attr<double>(_Unicode(dpos));
    double length = x_layer.attr<double>(_Unicode(length));
    double dlen = x_layer.attr<double>(_Unicode(dlen));
    // double rotation = x_layer.attr<double>(  _Unicode(rotation) ) ;

    // create a collection for a layer, all subdetector will be placed in it
    // std::string layername =  name +_toString(layer_id,"_layer%d");
    std::string layername = x_layer.nameStr() + _toString(layer_id, "_layer%d");
    Assembly layer_assembly(layername);
    DetElement layerDE(DriftLayer, layername, layer_id);
    // put the layer to the whole assembly
    // pv = assembly.placeVolume(  layer_assembly );
    pv = SubSolidMDC_vol.placeVolume(layer_assembly);
    pv.addPhysVolID("layer", layer_id);
    layerDE.setPlacement(pv);

    double dphi = 2. * M_PI / double(nwire);

    xml_comp_t x_signal = x_layer.child(_Unicode(signal));
    xml_comp_t x_field = x_layer.child(_Unicode(field));

    // for field wires
    double phi0_f = x_field.phi0();
    double radius_f = x_field.attr<double>(_Unicode(radius));
    double rot_f = x_field.attr<double>(_Unicode(rotation));
    std::string material_f = x_field.attr<std::string>(_Unicode(material));
    std::string visattr_f = x_field.attr<std::string>(_Unicode(vis));

    // create logical volume, not a realistic one.
    if (length > 0.1)
      lengthcnt = length;
    else
      lengthcnt += dlen;
    Tube tubeSolid_f(0, radius_f, lengthcnt / 2.0);
    // Tube   tubeSolid_f (0, radius_f, lengthcnt/2.0-0.5 ) ;
    Volume tube_vol_f(name + "WiresField", tubeSolid_f,
                      description.material(material_f.c_str()));
    // tube_vol_f.setSensitiveDetector(sens);
    tube_vol_f.setVisAttributes(description.visAttributes(visattr_f));

    if (position > 0.1)
      positioncnt = position;
    else
      positioncnt += dpos;

    rot_f = rot_f / 180 * M_PI; // deg to rad
    double dlphi = 0.5 * lengthcnt * tan(rot_f);
    double centerpos = sqrt(pow(positioncnt, 2) - pow(dlphi, 2));
    for (int j = 0; j < nwire; ++j) {
      std::string laddername = layername + _toString(j, "Wire%dField");
      Assembly wire_assembly(laddername.c_str());

      double phi = phi0_f + j * dphi;
      // RotationZYX rot(0.0, 0.0, rot_f);
      // pv = wire_assembly.placeVolume(tube_vol_f,Transform3D(rot,
      // Position(centerpos, 0, 0. ) )); RotationZYX rot2(phi, 0.0, 0); pv =
      // layer_assembly.placeVolume(wire_assembly, Transform3D(rot2, Position(0,
      // 0, 0. ) ));
      RotationZYX rot2(0, rot_f * sin(phi), rot_f * cos(phi));
      pv = layer_assembly.placeVolume(
          tube_vol_f, Transform3D(rot2, Position(centerpos * cos(phi),
                                                 centerpos * sin(phi), 0.)));
      pv.addPhysVolID("sub", 0);
      pv.addPhysVolID("wire", j);
      DetElement wireDE(layerDE, laddername, j);
      wireDE.setPlacement(pv);
    }

    // signal wires
    double phi0_s = x_signal.phi0();
    double radius_s = x_signal.attr<double>(_Unicode(radius));
    double rot_s = x_signal.attr<double>(_Unicode(rotation));
    std::string material_s = x_signal.attr<std::string>(_Unicode(material));
    std::string visattr_s = x_signal.attr<std::string>(_Unicode(vis));

    // create logical volume, not a realistic one.
    // if (length > 0.1 ) lengthcnt = length;
    // else lengthcnt += dlen;
    lengthcnt += dlen;
    Tube tubeSolid_s(0, radius_s, lengthcnt / 2.0);
    Volume tube_vol_s(name + "WiresSignal", tubeSolid_s,
                      description.material(material_s.c_str()));
    tube_vol_s.setVisAttributes(description.visAttributes(visattr_s));

    // tube_vol_s.setSensitiveDetector(sens);

    ////if (position > 0.1) positioncnt = position;
    ////else positioncnt += dpos;
    positioncnt += dpos;
    rot_s = rot_s / 180 * M_PI; // deg to rad
    dlphi = 0.5 * lengthcnt * tan(rot_s);
    centerpos = sqrt(pow(positioncnt, 2) - pow(dlphi, 2));
    for (int j = 0; j < nwire; j += 2) {
      std::string laddername = layername + _toString(j, "Wire%dSignal");

      Assembly wire_assembly(laddername.c_str());
      // pv = assembly.placeVolume(  layer_assembly );

      double phi = phi0_s + j * dphi;
      // RotationZYX rot(0.0, 0.0, rot_s);
      // pv = wire_assembly.placeVolume(tube_vol_s,Transform3D(rot,
      // Position(centerpos, 0, 0. ) )); RotationZYX rot2(phi, 0.0, 0); pv =
      // layer_assembly.placeVolume(wire_assembly, Transform3D(rot2, Position(0,
      // 0, 0. ) ));
      RotationZYX rot2(0, rot_s * sin(phi), rot_s * cos(phi));
      pv = layer_assembly.placeVolume(
          tube_vol_s, Transform3D(rot2, Position(centerpos * cos(phi),
                                                 centerpos * sin(phi), 0.)));
      pv.addPhysVolID("sub", 1);
      pv.addPhysVolID("wire", j);
      DetElement wireDE(layerDE, laddername, j);
      wireDE.setPlacement(pv);
    }

    // place field wire in the signal layer
    Tube tubeSolid_f2(0, radius_f, lengthcnt / 2.0);
    Volume tube_vol_f2(name + "WiresFieldInSignalLayer", tubeSolid_f2,
                       description.material(material_f.c_str()));
    tube_vol_f2.setVisAttributes(description.visAttributes(visattr_f));
    for (int j = 1; j < nwire; j += 2) {
      std::string laddername =
          layername + _toString(j, "Wire%dFieldInSignalLayer");

      // Assembly wire_assembly( laddername.c_str() );
      // pv = assembly.placeVolume(  layer_assembly );
      double phi = phi0_s + j * dphi;
      // RotationZYX rot(0.0, 0.0, rot_s);
      // pv = wire_assembly.placeVolume(tube_vol_f2,Transform3D(rot,
      // Position(centerpos, 0, 0. ) )); RotationZYX rot2(phi, 0.0, 0); pv =
      // layer_assembly.placeVolume(wire_assembly, Transform3D(rot2, Position(0,
      // 0, 0. ) ));
      RotationZYX rot2(0, rot_s * sin(phi), rot_s * cos(phi));
      pv = layer_assembly.placeVolume(
          tube_vol_f2, Transform3D(rot2, Position(centerpos * cos(phi),
                                                  centerpos * sin(phi), 0.)));
      pv.addPhysVolID("sub", 2);
      pv.addPhysVolID("wire", j);
      DetElement wireDE(layerDE, laddername, j);
      wireDE.setPlacement(pv);
    }

    bool comp2 = x_layer.hasChild(_Unicode(field2));
    if (comp2) {
      // if (false) {
      xml_comp_t x_field2 = x_layer.child(_Unicode(field2));
      // for field wires
      double phi0_f2 = x_field2.phi0();
      double len_f2 = x_field2.attr<double>(_Unicode(length));
      double pos_f2 = x_field2.attr<double>(_Unicode(position));
      double radius_f2 = x_field2.attr<double>(_Unicode(radius));
      double rot_f2 = x_field2.attr<double>(_Unicode(rotation));
      std::string material_f2 = x_field2.attr<std::string>(_Unicode(material));
      std::string visattr_f2 = x_field2.attr<std::string>(_Unicode(vis));

      // create logical volume, not a realistic one.
      lengthcnt = len_f2;
      Tube tubeSolid_f2(0, radius_f2, lengthcnt / 2.0);
      // Tube   tubeSolid_f2 (0, radius_f2, lengthcnt/2.0-1 ) ;
      Volume tube_vol_f2(name + "OuterWiresField", tubeSolid_f2,
                         description.material(material_f2.c_str()));
      tube_vol_f2.setVisAttributes(description.visAttributes(visattr_f2));

      positioncnt = pos_f2;
      rot_f2 = rot_f2 / 180 * M_PI; // deg to rad
      double dlphi = 0.5 * lengthcnt * tan(rot_f2);
      double centerpos = sqrt(pow(positioncnt, 2) - pow(dlphi, 2));
      for (int j = 0; j < nwire; ++j) {
        std::string laddername = layername + _toString(j, "Wire%dField2");
        // Assembly wire_assembly( laddername.c_str() );

        double phi = phi0_f + j * dphi;
        ////RotationZYX rot(0.0, 0.0, rot_f2);
        ////pv = wire_assembly.placeVolume(tube_vol_f,Transform3D(rot,
        ///Position(centerpos, 0, 0. ) )); /RotationZYX rot2(phi, 0.0, 0); /pv =
        ///layer_assembly.placeVolume(wire_assembly, Transform3D(rot2,
        ///Position(0, 0, 0. ) ));
        RotationZYX rot2(0, rot_f2 * sin(phi), rot_f2 * cos(phi));
        pv = layer_assembly.placeVolume(
            tube_vol_f2, Transform3D(rot2, Position(centerpos * cos(phi),
                                                    centerpos * sin(phi), 0.)));
        pv.addPhysVolID("sub", 3);
        pv.addPhysVolID("wire", j);
        DetElement wireDE(layerDE, laddername, j);
        wireDE.setPlacement(pv);
      }
    }
  }

  Volume mother = description.pickMotherVolume(DriftLayer);
  // pv = mother.placeVolume( assembly ) ;
  pv = mother.placeVolume(SubSolidMDC_vol);
  pv.addPhysVolID("system", x_det.id()); //.addPhysVolID("side",0 )  ;
  // pv.addPhysVolID( "system", x_det.id() ).addPhysVolID("side",0 )  ;
  DriftLayer.setPlacement(pv);

  return DriftLayer;
}
DECLARE_DETELEMENT(DriftChamber, create_element)
