#include "DD4hep/DetFactoryHelper.h"

#include "DD4hep/Detector.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"

#include <exception>
#include <iostream>
#include <map>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {
  xml_det_t x_det = e;

  DetElement Dirc(x_det.nameStr(), x_det.id());
  Volume mother = description.pickMotherVolume(Dirc);

  //---------------------Dirc and Sector-----------------------
  xml_comp_t x_sector(x_det.child(_Unicode(Sector)));
  const int SectorNu = x_sector.attr<int>(_Unicode(nSectors));
  double SectorR1 = x_sector.attr<double>(_Unicode(radius1));
  double SectorR2 = x_sector.attr<double>(_Unicode(radius2));
  double SectorZ = x_sector.attr<double>(_Unicode(positionZ));
  double SectorT = x_sector.thickness();
  double SectorH = SectorR2 * cos(360. * deg / 2. / SectorNu) - SectorR1;
  double SectorW1 = SectorR1 * tan(360. * deg / 2. / SectorNu) * 2.;
  double SectorW2 = SectorR2 * sin(360. * deg / 2. / SectorNu) * 2.;

  Tube SolidDircP1(SectorR1, SectorR2, SectorZ + SectorT / 2., 0., 360. * deg);
  Tube SolidDircP2(SectorR1, SectorR2, SectorZ - SectorT / 2., 0., 360. * deg);
  SubtractionSolid SolidDirc(SolidDircP1, SolidDircP2);
  Volume DircVol(x_det.nameStr(), SolidDirc,
                 description.material(x_det.attr<string>(_U(material))));
  DircVol.setVisAttributes(description, x_det.visStr());
  PlacedVolume pvDirc = mother.placeVolume(DircVol, x_det.id());
  pvDirc.addPhysVolID(x_det.nameStr(), x_det.id());
  Dirc.setPlacement(pvDirc);

  Trapezoid SolidSector(SectorW1 / 2., SectorW2 / 2., SectorT / 2.,
                        SectorT / 2., SectorH / 2.);
  Volume SectorVol(x_sector.nameStr(), SolidSector,
                   description.material(x_sector.materialStr()));
  SectorVol.setVisAttributes(description, x_sector.visStr());
  PlacedVolume pvSector[SectorNu * 2];
  for (int i = 0; i < SectorNu * 2; i++) {
    int tag = i / SectorNu == 0 ? 1 : -1;
    int j = i % SectorNu;
    double Sector_posZ = -SectorH / 2. - SectorR1;
    double theta = 360. * deg / SectorNu * tag;
    double theta0 = 90. * deg * (1 - tag);
    pvSector[i] = DircVol.placeVolume(
        SectorVol, i,
        Transform3D(RotationZYX(0, theta * j + theta0, -90 * deg * tag),
                    Position(-Sector_posZ * sin(theta * j + theta0),
                             -Sector_posZ * cos(theta * j + theta0) * tag,
                             SectorZ * tag)));
    pvSector[i].addPhysVolID(x_sector.nameStr(), i);
  }
  //-----------------------------------------------------------

  //---------------------Shield--------------------------------
  xml_comp_t x_shield(x_det.child(_Unicode(Shield)));
  double Shieldthick = x_shield.thickness();
  double ShieldW1 =
      SectorW1 + (tan(360. * deg / 2. / SectorNu) - 1.) * Shieldthick;
  double ShieldW2 =
      SectorW2 - (tan(360. * deg / 2. / SectorNu) + 1.) * Shieldthick;
  double ShieldT = SectorT - Shieldthick;
  double ShieldH = SectorH - Shieldthick;

  Trapezoid SolidShieldP1(ShieldW1 / 2., ShieldW2 / 2., ShieldT / 2.,
                          ShieldT / 2., ShieldH / 2.);
  SubtractionSolid SolidShield(SolidSector, SolidShieldP1);
  Volume ShieldVol(x_shield.nameStr(), SolidShield,
                   description.material(x_shield.materialStr()));
  ShieldVol.setVisAttributes(description, x_shield.visStr());
  PlacedVolume pvShield = SectorVol.placeVolume(ShieldVol, 0);
  pvShield.addPhysVolID(x_shield.nameStr(), 0);
  //-----------------------------------------------------------

  //---------------------Quartz and Focus----------------------
  xml_comp_t x_quartz(x_det.child(_Unicode(Quartz)));
  double QuartzR = x_quartz.radius();
  double QuartzT = x_quartz.thickness();
  double QuartzW1 = x_quartz.attr<double>(_Unicode(width1));
  double QuartzW2 = x_quartz.attr<double>(_Unicode(width2));
  double QuartzH = x_quartz.height();

  Trapezoid SolidQuartz(QuartzW1 / 2., QuartzW2 / 2., QuartzT / 2.,
                        QuartzT / 2., QuartzH / 2.);
  Volume QuartzVol(x_quartz.nameStr(), SolidQuartz,
                   description.material(x_quartz.materialStr()));
  QuartzVol.setVisAttributes(description, x_quartz.visStr());
  PlacedVolume pvQuartz =
      SectorVol.placeVolume(QuartzVol, 0, Position(0, 0, -QuartzR / 2.));
  pvQuartz.addPhysVolID(x_quartz.nameStr(), 0);

  Tube SolidFocus(0., QuartzR, QuartzW2 / 2., 90. * deg, 180. * deg);
  Volume FocusVol(x_quartz.nameStr() + "Focus", SolidFocus,
                  description.material(x_quartz.materialStr()));
  FocusVol.setVisAttributes(description, x_quartz.visStr());
  PlacedVolume pvFocus =
      SectorVol.placeVolume(FocusVol, 0,
                            Transform3D(RotationZYX(0, 90. * deg, 0),
                                        Position(0, QuartzT / 2. - QuartzR,
                                                 QuartzH / 2. - QuartzR / 2.)));
  pvFocus.addPhysVolID(x_quartz.nameStr() + "Focus", 0);
  //-----------------------------------------------------------

  //---------------------Coating-------------------------------
  xml_comp_t x_coating(x_det.child(_Unicode(Coating)));
  double CoatingT = x_coating.thickness();

  Box SolidCoatingTop(QuartzW1 / 2., QuartzT / 2., CoatingT / 2.);
  Volume CoatingTopVol(x_coating.nameStr() + "Top", SolidCoatingTop,
                       description.material(x_coating.materialStr()));
  CoatingTopVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingTop = SectorVol.placeVolume(
      CoatingTopVol, 0,
      Position(0, 0., -QuartzH / 2. - QuartzR / 2. - CoatingT / 2.));
  pvCoatingTop.addPhysVolID(x_coating.nameStr() + "Top", 0);

  Box SolidCoatingBottom(QuartzW2 / 2., QuartzR / 2. - QuartzT / 2.,
                         CoatingT / 2.);
  Volume CoatingBottomVol(x_coating.nameStr() + "Bottom", SolidCoatingBottom,
                          description.material(x_coating.materialStr()));
  CoatingBottomVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingBottom = SectorVol.placeVolume(
      CoatingBottomVol, 0,
      Position(0, -QuartzR / 2., QuartzH / 2. - QuartzR / 2. - CoatingT / 2.));
  pvCoatingBottom.addPhysVolID(x_coating.nameStr() + "Bottom", 0);

  Tube SolidCoatingFocusP1(0., QuartzR + CoatingT, QuartzW2 / 2., 90. * deg,
                           180. * deg);
  SubtractionSolid SolidCoatingFocus(SolidCoatingFocusP1, SolidFocus);
  Volume CoatingFocusVol(x_coating.nameStr() + "Focus", SolidCoatingFocus,
                         description.material(x_coating.materialStr()));
  CoatingFocusVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingFocus =
      SectorVol.placeVolume(CoatingFocusVol, 0,
                            Transform3D(RotationZYX(0., 90. * deg, 0.),
                                        Position(0., QuartzT / 2. - QuartzR,
                                                 QuartzH / 2. - QuartzR / 2.)));
  pvCoatingFocus.addPhysVolID(x_coating.nameStr() + "Focus", 0);

  Tube SolidCoatingSideFP1(0., QuartzR, QuartzW2 / 2. + CoatingT, 90. * deg,
                           180. * deg);
  SubtractionSolid SolidCoatingSideF(SolidCoatingSideFP1, SolidFocus);
  Volume CoatingSideFVol(x_coating.nameStr() + "SideF", SolidCoatingSideF,
                         description.material(x_coating.materialStr()));
  CoatingSideFVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingSideF =
      SectorVol.placeVolume(CoatingSideFVol, 0,
                            Transform3D(RotationZYX(0., 90. * deg, 0.),
                                        Position(0., QuartzT / 2. - QuartzR,
                                                 QuartzH / 2. - QuartzR / 2.)));
  pvCoatingSideF.addPhysVolID(x_coating.nameStr() + "SideF", 0);

  Trapezoid SolidCoatingSideQP1(QuartzW1 / 2. + CoatingT,
                                QuartzW2 / 2. + CoatingT, QuartzT / 2.,
                                QuartzT / 2., QuartzH / 2.);
  SubtractionSolid SolidCoatingSideQ(SolidCoatingSideQP1, SolidQuartz);
  Volume CoatingSideQVol(x_coating.nameStr() + "sideQ", SolidCoatingSideQ,
                         description.material(x_coating.materialStr()));
  CoatingSideQVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingSideQ =
      SectorVol.placeVolume(CoatingSideQVol, 0, Position(0, 0, -QuartzR / 2.));
  pvCoatingSideQ.addPhysVolID(x_coating.nameStr(), 0);
  //-----------------------------------------------------------

  //---------------------PhotonDet and Cathode-----------------
  xml_comp_t x_PhotonDet(x_det.child(_Unicode(PhotonDet)));
  xml_comp_t x_cathode(x_det.child(_Unicode(Cathode)));
  const int PhotonDetNu = x_PhotonDet.attr<int>(_Unicode(nPhotonDets));
  double PhotonDetW = x_PhotonDet.width();
  double PhotonDetT = x_PhotonDet.thickness();
  double PhotonDetH = x_PhotonDet.height();
  double PhotonDetGap = x_PhotonDet.attr<double>(_Unicode(gap));
  const int CathodeNu = x_cathode.attr<int>(_Unicode(nCathodes));
  double CathodeW = x_cathode.width();
  double CathodeT = x_cathode.thickness();
  double CathodeH = x_cathode.height();

  Box SolidPhotonDet(PhotonDetW / 2., PhotonDetT / 2., PhotonDetH / 2.);
  Volume PhotonDetVol(x_PhotonDet.nameStr(), SolidPhotonDet,
                      description.material(x_PhotonDet.materialStr()));
  PhotonDetVol.setVisAttributes(description, x_PhotonDet.visStr());
  PlacedVolume pvPhotonDet[PhotonDetNu];
  for (int i = 0; i < PhotonDetNu; i++) {
    double PhotonDet_posX =
        -(PhotonDetW + PhotonDetGap) * (PhotonDetNu - 1) / 2. +
        (PhotonDetW + PhotonDetGap) * i;
    double PhotonDet_posY = QuartzT / 2. - PhotonDetT / 2. - QuartzR;
    double PhotonDet_posZ = QuartzH / 2. - CathodeH / 2. + QuartzR / 2.;
    pvPhotonDet[i] = SectorVol.placeVolume(
        PhotonDetVol, i,
        Position(PhotonDet_posX, PhotonDet_posY, PhotonDet_posZ));
    pvPhotonDet[i].addPhysVolID(x_PhotonDet.nameStr(), i);
  }

  Box SolidCathode(CathodeW / 2., CathodeT / 2., CathodeH / 2.);
  Volume CathodeVol(x_cathode.nameStr(), SolidCathode,
                    description.material(x_cathode.materialStr()));
  CathodeVol.setVisAttributes(description, x_cathode.visStr());
  PlacedVolume pvCathode[CathodeNu];
  for (int i = 0; i < CathodeNu; i++) {
    double Cathode_posX = -CathodeW * (CathodeNu - 1) / 2. + CathodeW * i;
    double Cathode_posY = PhotonDetT / 2. - CathodeT / 2.;
    pvCathode[i] = PhotonDetVol.placeVolume(
        CathodeVol, i, Position(Cathode_posX, Cathode_posY, 0));
    pvCathode[i].addPhysVolID(x_cathode.nameStr(), i);
  }
  //-----------------------------------------------------------

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 17, 0)
  //---------------------OpticalSurface------------------------
  OpticalSurfaceManager surfMgr = description.surfaceManager();

  OpticalSurface QuartzAirSurf = surfMgr.opticalSurface("QuartzAirSurface");
  cout << "11111" << endl;
  for (int i = 0; i < SectorNu * 2; i++) {
    BorderSurface(description, Dirc, _toString(i, "QuartzAirSurf_%d"),
                  QuartzAirSurf, pvQuartz, pvSector[i]);
    BorderSurface(description, Dirc, _toString(i, "FocusAirSurf_%d"),
                  QuartzAirSurf, pvFocus, pvSector[i]);
  }

  OpticalSurface QuartzCoatingFocusSurf =
      surfMgr.opticalSurface("QuartzCoatingFocusSurface");
  BorderSurface(description, Dirc, "QuartzCoatingFocusSurface",
                QuartzCoatingFocusSurf, pvFocus, pvCoatingFocus);

  OpticalSurface QuartzCoatingTopSurf =
      surfMgr.opticalSurface("QuartzCoatingTopSurface");
  BorderSurface(description, Dirc, "QuartzCoatingTopSurface",
                QuartzCoatingTopSurf, pvQuartz, pvCoatingTop);

  OpticalSurface QuartzCoatingBottomSurf =
      surfMgr.opticalSurface("QuartzCoatingBottomSurface");
  BorderSurface(description, Dirc, "QuartzCoatingBottomSurface",
                QuartzCoatingBottomSurf, pvFocus, pvCoatingBottom);

  OpticalSurface QuartzCoatingSideFSurf =
      surfMgr.opticalSurface("QuartzCoatingSideFSurface");
  BorderSurface(description, Dirc, "QuartzCoatingSideFSurface",
                QuartzCoatingSideFSurf, pvFocus, pvCoatingSideF);

  OpticalSurface QuartzCoatingSideQSurf =
      surfMgr.opticalSurface("QuartzCoatingSideQSurface");
  BorderSurface(description, Dirc, "QuartzCoatingSideQSurface",
                QuartzCoatingSideQSurf, pvQuartz, pvCoatingSideQ);

  OpticalSurface QuartzPhotonDetSurf =
      surfMgr.opticalSurface("QuartzPhotonDetSurface");
  for (int i = 0; i < PhotonDetNu; i++)
    BorderSurface(description, Dirc, _toString(i, "QuartzPhotonDetSurface_%d"),
                  QuartzPhotonDetSurf, pvFocus, pvPhotonDet[i]);

  OpticalSurface QuartzCathodeSurf =
      surfMgr.opticalSurface("QuartzCathodeSurface");
  for (int i = 0; i < CathodeNu; i++)
    BorderSurface(description, Dirc, _toString(i, "QuartzCathodeSurface_%d"),
                  QuartzCathodeSurf, pvFocus, pvCathode[i]);

  OpticalSurface PhotonDetSurf = surfMgr.opticalSurface("PhotonDetSurface");
  SkinSurface(description, Dirc, "PhotonDetSurface", PhotonDetSurf,
              PhotonDetVol);

  OpticalSurface CoatingSurf = surfMgr.opticalSurface("CoatingSurface");
  SkinSurface(description, Dirc, "CoatingTopSurface", CoatingSurf,
              CoatingTopVol);
  SkinSurface(description, Dirc, "CoatingBottomSurface", CoatingSurf,
              CoatingBottomVol);
  SkinSurface(description, Dirc, "CoatingSideFSurface", CoatingSurf,
              CoatingSideFVol);
  SkinSurface(description, Dirc, "CoatingSideQSurface", CoatingSurf,
              CoatingSideQVol);
  SkinSurface(description, Dirc, "CoatingFocusSurface", CoatingSurf,
              CoatingFocusVol);
  //-----------------------------------------------------------
#endif

  return Dirc;
}
DECLARE_DETELEMENT(StcfDircGeoV1, create_element)
