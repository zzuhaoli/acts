//====================================================================
//  // Simple tube filled with air
//  // used for tracking purposes ...
//
//--------------------------------------------------------------------
//
//  Author     : F.Gaede
//
//====================================================================
#include "DD4hep/DetFactoryHelper.h"

#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hep/DetType.h"

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

Tube construct_tube(xml_comp_t comp);
Polycone construct_polycone(xml_comp_t comp);
Box construct_box(xml_comp_t comp);
Trd2 construct_trd(xml_comp_t comp);

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector /* sens */) {

  xml_det_t x_det = e;
  std::string name = x_det.nameStr();

  DetElement beamPipe(name, x_det.id());
  DetType type( DetType::TRACKER | DetType::BEAMPIPE);
  beamPipe.setTypeFlag(type.to_ulong());

  Assembly assembly(name + "_assembly");

  PlacedVolume pv;

  // ----- read xml ----------------------
  int id = 0;
  double frontDis = 0.0;
  // double pi = 3.1415926;
  for (xml_coll_t i(x_det, _Unicode(sector)); i; i++) {
    id++;
    xml_comp_t sect = i;
    std::string sectName = sect.nameStr();
    std::string typeN = sect.attr<std::string>(_U(type));
    std::string material = sect.attr<std::string>(_U(material));
    std::string visStr = sect.visStr();
    if (typeN == "Tube") {
      // double rmin = sect.rmin();
      // double rmax = sect.rmax();
      double zmin = sect.zmin();
      double zmax = sect.zmax();
      // double zhalf = (zmax - zmin)/2;
      double zpos = (zmax - zmin) / 2.0 + zmin;

      // Tube tubeSolid(rmin, rmax, zhalf);
      Tube tubeSolid = construct_tube(sect);
      Volume tube_vol(sectName, tubeSolid, description.material(material));
      tube_vol.setVisAttributes(description.visAttributes(visStr));
      // Volume mother =  description.pickMotherVolume( beamPipe ) ;
      pv = assembly.placeVolume(
          tube_vol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., zpos)));
      pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
      id++;
      pv = assembly.placeVolume(
          tube_vol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., -zpos)));
      pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
    } else if (typeN == "CrossTube") {
      double rmin = sect.rmin();
      double rmax = sect.rmax();
      double zmin = sect.zmin();
      double zmax = sect.zmax();
      double cdis = sect.attr<double>(_Unicode(centerDistance));
      double angle = sect.angle();
      double zhalf = (zmax - zmin) / 2;
      double zpos = (zmax - zmin) / 2.0 + zmin;
      if (cdis < 1e-6) {
        cdis = frontDis;
      }
      double centerDisx = cdis * cos(angle) + zhalf * sin(angle * 2);
      double centerDisz = cdis * sin(angle) + zhalf * (1 - cos(angle * 2));
      double centerDis = cdis + zhalf * sin(angle) * 2;
      frontDis = cdis + 2 * zhalf * sin(angle) * 2;

      Tube tube1(rmin, rmax, zhalf);
      Tube tube2(rmin, rmax, zhalf);

      UnionSolid crossTube(tube1, tube2,
                           Transform3D(RotationZYX(0, angle * 2, 0),
                                       Position(centerDisx, 0, centerDisz)));
      Volume tube_vol(sectName, crossTube, description.material(material));

      tube_vol.setVisAttributes(description.visAttributes(visStr));
      pv = assembly.placeVolume(
          tube_vol, Transform3D(RotationZYX(0, -angle, 0),
                                Position(0. - centerDis / 2, 0, zpos)));
      pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
      id++;
      pv = assembly.placeVolume(
          tube_vol, Transform3D(RotationZYX(0, -angle + pi, 0),
                                Position(0. + centerDis / 2, 0, -zpos)));
      pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);

    } else if (typeN == "Polycone") {
      // std::vector<double> rmin;
      // std::vector<double> rmax;
      // std::vector<double> z;
      // for (xml_coll_t no(sect, _Unicode(node)); no; no++) {
      //     xml_comp_t node = no;
      //     double r1 = node.rmin();
      //     double r2 = node.rmax();
      //     double z1 = node.z();
      //     rmin.push_back(r1);
      //     rmax.push_back(r2);
      //     z.push_back(z1);
      // }
      double zpos = 0;

      // Polycone pcone(0, 360*deg, rmin, rmax, z); // std::vector<double> rmin,
      // rmax, z;
      Polycone pcone = construct_polycone(sect);
      Volume pc_vol(sectName, pcone, description.material(material));
      pc_vol.setVisAttributes(description.visAttributes(visStr));
      pv = assembly.placeVolume(
          pc_vol, Transform3D(RotationZYX(0, 0., 0), Position(0, 0., zpos)));
      pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
      id++;
      pv = assembly.placeVolume(
          pc_vol, Transform3D(RotationZYX(0, pi, 0), Position(0, 0., -zpos)));
      pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
    } else if (typeN == "UnionSolid") {
      UnionSolid *pUniSolid = 0;
      Solid *firstSolid = 0;
      for (xml_coll_t unicol(sect, _Unicode(block)); unicol; unicol++) {
        xml_comp_t blk = unicol;

        // std::string blkName  = blk.nameStr();
        std::string typeN2 = blk.attr<std::string>(_U(type));
        Solid *blkSolid = 0;
        double angle = blk.angle();
        double transx = blk.attr<double>(_Unicode(transx));
        double transy = blk.attr<double>(_Unicode(transy));
        double transz = blk.attr<double>(_Unicode(transz));
        Transform3D transform3DUni(RotationZYX(0, angle, 0),
                                   Position(transx, transy, transz));
        if (typeN2 == "SubtractionSolid") {
          xml_comp_t mS = blk.child(_Unicode(mother));
          std::string mType = mS.attr<std::string>(_U(type));
          Solid *pMother = 0;
          // std::cout<<"subtractionSolid mother type "<<mType<<std::endl;
          if (mType == "Trd") {
            pMother = new Solid(construct_trd(mS));
          }
          if (mType == "Box") {
            pMother = new Solid(construct_box(mS));
          }
          if (mType == "Tube") {
            pMother = new Solid(construct_tube(mS));
          }
          if (mType == "Polycone") {
            pMother = new Solid(construct_polycone(mS));
          }

          xml_comp_t dS = blk.child(_Unicode(daughter));
          std::string dType = dS.attr<std::string>(_U(type));
          Solid *pDaughter = 0;
          if (dType == "Trd") {
            pDaughter = new Solid(construct_trd(dS));
          }
          if (dType == "Box") {
            pDaughter = new Solid(construct_box(dS));
          }
          if (dType == "Tube") {
            pDaughter = new Solid(construct_tube(dS));
          }
          if (dType == "Polycone") {
            pDaughter = new Solid(construct_polycone(dS));
          }

          if (pMother != 0 && pDaughter != 0) {
            blkSolid = new Solid(SubtractionSolid(*pMother, *pDaughter));
          } else {
            std::cout << "No mother or No daughter in the SubtractionSolid"
                      << std::endl;
          }
        } else if (typeN2 == "Tube") {
          blkSolid = new Solid(construct_tube(blk));
        } else if (typeN2 == "Polycone") {
          blkSolid = new Solid(construct_polycone(blk));
        } else if (typeN2 == "Box") {
          blkSolid = new Solid(construct_box(blk));
        } else if (typeN2 == "Trd") {
          blkSolid = new Solid(construct_trd(blk));
        }

        if (pUniSolid != 0 && blkSolid != 0) {
          pUniSolid = new UnionSolid(*pUniSolid, *blkSolid, transform3DUni);
        } else if (pUniSolid == 0 && firstSolid != 0 && blkSolid != 0) {
          pUniSolid = new UnionSolid(*firstSolid, *blkSolid, transform3DUni);
        } else if (firstSolid == 0 && blkSolid != 0) {
          firstSolid = blkSolid;
        }
      }
      if (pUniSolid != 0) {
        Volume union_vol(sectName, *pUniSolid, description.material(material));
        union_vol.setVisAttributes(description.visAttributes(visStr));
        double transx = sect.attr<double>(_Unicode(transx));
        double transy = sect.attr<double>(_Unicode(transy));
        double transz = sect.attr<double>(_Unicode(transz));
        // Volume mother =  description.pickMotherVolume( beamPipe ) ;
        pv = assembly.placeVolume(
            union_vol, Transform3D(RotationZYX(0, 0., 0),
                                   Position(transx, transy, transz)));
        pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
        id++;
        pv = assembly.placeVolume(
            union_vol, Transform3D(RotationZYX(0, pi, 0),
                                   Position(transx, transy, -transz)));
        pv.addPhysVolID("system", x_det.id()).addPhysVolID(sectName, id);
      }

    } else {
      std::cout << "Warning: DD4hep does not know how to construct shape "
                << typeN << std::endl;
    }
  }

  Volume mother = description.pickMotherVolume(beamPipe);
  pv = mother.placeVolume(assembly);
  pv.addPhysVolID("system", x_det.id());
  beamPipe.setPlacement(pv);

  return beamPipe;
}

Tube construct_tube(xml_comp_t comp) {
  double rmin = comp.rmin();
  double rmax = comp.rmax();
  double zmin = comp.zmin();
  double zmax = comp.zmax();
  double zhalf = (zmax - zmin) / 2;
  // double zpos  = (zmax-zmin)/2.0+zmin;

  Tube tubeSolid(rmin, rmax, zhalf);
  return tubeSolid;
}
Polycone construct_polycone(xml_comp_t comp) {
  std::vector<double> rmin;
  std::vector<double> rmax;
  std::vector<double> z;
  for (xml_coll_t no(comp, _Unicode(node)); no; no++) {
    xml_comp_t node = no;
    double r1 = node.rmin();
    double r2 = node.rmax();
    double z1 = node.z();
    rmin.push_back(r1);
    rmax.push_back(r2);
    z.push_back(z1);
  }
  // double zpos = 0;

  Polycone pcone(0, 360 * deg, rmin, rmax,
                 z); // std::vector<double> rmin, rmax, z;

  return pcone;
}
Box construct_box(xml_comp_t comp) {
  double xhalf = comp.x();
  double yhalf = comp.y();
  double zhalf = comp.z();

  Box box(xhalf, yhalf, zhalf);
  return box;
}
Trd2 construct_trd(xml_comp_t comp) {
  double x1 = comp.x1();
  double x2 = comp.x2();
  double y1 = comp.y1();
  double y2 = comp.y2();
  double z = comp.z();

  Trd2 trd(x1, x2, y1, y2, z);
  return trd;
}

DECLARE_DETELEMENT(supportMDI01, create_element)
