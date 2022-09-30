#include "DD4hep/DetFactoryHelper.h"

#include "DD4hep/Detector.h"
#include "DD4hep/OpticalSurfaces.h"
#include "DD4hep/Printout.h"
#include "DDRec/DetectorData.h"
#include "DDRec/Surface.h"
#include "DD4hep/DetType.h"

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

  DetElement Ftof(x_det.nameStr(), x_det.id());
  Volume mother = description.pickMotherVolume(Ftof);
  DetType type( DetType::TRACKER);
  Ftof.setTypeFlag(type.to_ulong());



  //---------------------Ftof and Sector-----------------------
  xml_comp_t x_sector(x_det.child(_Unicode(Sector)));
  const int SectorNu = x_sector.attr<int>(_Unicode(nSectors));
  double SectorR1 = x_sector.attr<double>(_Unicode(radius1)) / sin(75 * deg);
  double SectorR2 = x_sector.attr<double>(_Unicode(radius2));
  double SectorZ = x_sector.attr<double>(_Unicode(positionZ));
  double SectorT = x_sector.thickness();
  const vector<double> SectorTArray = {-0.5 * SectorT, 0.5 * SectorT};
  const vector<double> SectorR1Array = {SectorR1 * sin(75 * deg),
                                        SectorR1 * sin(75 * deg)};
  const vector<double> SectorR2Array = {SectorR2 * sin(75 * deg),
                                        SectorR2 * sin(75 * deg)};

  Tube SolidFtofP1(SectorR1 * sin(75 * deg), SectorR2, SectorZ + SectorT / 2.,
                   0., 360. * deg);
  Tube SolidFtofP2(SectorR1 * sin(75 * deg), SectorR2, SectorZ - SectorT / 2.,
                   0., 360. * deg);
  SubtractionSolid SolidFtof(SolidFtofP1, SolidFtofP2);
  Volume FtofVol(x_det.nameStr(), SolidFtof,
                 description.material(x_det.attr<string>(_U(material))));
  FtofVol.setVisAttributes(description, x_det.visStr());
  PlacedVolume pvFtof = mother.placeVolume(FtofVol, x_det.id());
  pvFtof.addPhysVolID(x_det.nameStr(), x_det.id());
  Ftof.setPlacement(pvFtof);

  Polyhedra SolidSector(3, 0, 90 * deg, SectorTArray, SectorR1Array,
                        SectorR2Array);
  Volume SectorVol(x_sector.nameStr(), SolidSector,
                   description.material(x_sector.materialStr()));
  SectorVol.setVisAttributes(description, x_sector.visStr());
  PlacedVolume pvSector[SectorNu * 2];
  for (int i = 0; i < SectorNu; i++) // at +Z
  {
    double theta = 360. * deg / SectorNu;
    double theta0 = 0. * deg;
    pvSector[i] =
        FtofVol.placeVolume(SectorVol, i,
                            Transform3D(RotationZYX(theta * i + theta0, 0, 0),
                                        Position(0, 0, SectorZ)));
    pvSector[i].addPhysVolID(x_sector.nameStr(), i);
  }
  for (int i = 0; i < SectorNu; i++) // at -Z
  {
    double theta = 360. * deg / SectorNu;
    double theta0 = 0. * deg;
    pvSector[i + SectorNu] = FtofVol.placeVolume(
        SectorVol, i + SectorNu,
        Transform3D(RotationZYX(theta * i + theta0, 0, 180 * deg),
                    Position(0, 0, -SectorZ)));
    pvSector[i + SectorNu].addPhysVolID(x_sector.nameStr(), i);
  }
  //-----------------------------------------------------------

  //---------------------Shield--------------------------------
  xml_comp_t x_shield(x_det.child(_Unicode(Shield)));
  double ShieldT = x_shield.thickness();
  double ShieldR1 = SectorR1;
  double ShieldR2 = SectorR2 - 2 * ShieldT;
  const vector<double> ShieldTArray = {-0.5 * SectorT + ShieldT,
                                       0.5 * SectorT - ShieldT};
  const vector<double> ShieldR1Array = {ShieldR1 * sin(75 * deg),
                                        ShieldR1 * sin(75 * deg)};
  const vector<double> ShieldR2Array = {ShieldR2 * sin(75 * deg),
                                        ShieldR2 * sin(75 * deg)};

  Polyhedra SolidShieldOut(3, 0, 90 * deg, SectorTArray, SectorR1Array,
                           SectorR2Array);
  Polyhedra SolidShieldIn(3, 0, 90 * deg, ShieldTArray, ShieldR1Array,
                          ShieldR2Array);
  SubtractionSolid SolidShield(SolidShieldOut, SolidShieldIn,
                               Position(ShieldT, ShieldT, 0));
  Volume ShieldVol(x_shield.nameStr(), SolidShield,
                   description.material(x_shield.materialStr()));
  ShieldVol.setVisAttributes(description, x_shield.visStr());
  PlacedVolume pvShield = SectorVol.placeVolume(ShieldVol, 0);
  pvShield.addPhysVolID(x_shield.nameStr(), 0);
  //-----------------------------------------------------------

  //---------------------Quartz -----------------------------
  xml_comp_t x_quartz(x_det.child(_Unicode(Quartz)));
  double QuartzT = x_quartz.thickness();
  double interval = x_quartz.attr<double>(_Unicode(interval));
  double QuartzR1 = SectorR1;
  double QuartzR2 = SectorR2 - 2 * interval;
  double QuartzZ = -SectorT / 2 + ShieldT + interval + QuartzT / 2;
  const vector<double> QuartzTArray = {-0.5 * QuartzT, 0.5 * QuartzT};
  const vector<double> QuartzR1Array = {QuartzR1 * sin(75 * deg),
                                        QuartzR1 * sin(75 * deg)};
  const vector<double> QuartzR2Array = {QuartzR2 * sin(75 * deg),
                                        QuartzR2 * sin(75 * deg)};

  Polyhedra SolidQuartz(3, 0, 90 * deg, QuartzTArray, QuartzR1Array,
                        QuartzR2Array);
  Volume QuartzVol(x_quartz.nameStr(), SolidQuartz,
                   description.material(x_quartz.materialStr()));
  QuartzVol.setVisAttributes(description, x_quartz.visStr());
  PlacedVolume pvQuartz = SectorVol.placeVolume(
      QuartzVol, 0, Position(interval, interval, QuartzZ));
  pvQuartz.addPhysVolID(x_quartz.nameStr(), 0);

  //-----------------------------------------------------------

  ////---------------------Coating-------------------------------
  xml_comp_t x_coating(x_det.child(_Unicode(Coating)));
  double CoatingT = x_coating.thickness();
  const vector<double> UpCoatingTArray = {-0.5 * QuartzT, 0.5 * QuartzT};
  const vector<double> UpCoatingR1Array = {QuartzR2 * sin(75 * deg),
                                           QuartzR2 * sin(75 * deg)};
  const vector<double> UpCoatingR2Array = {QuartzR2 * sin(75 * deg) + CoatingT,
                                           QuartzR2 * sin(75 * deg) + CoatingT};
  const vector<double> DownCoatingTArray = {-0.5 * QuartzT, 0.5 * QuartzT};
  const vector<double> DownCoatingR1Array = {
      QuartzR1 * sin(75 * deg) - CoatingT, QuartzR1 * sin(75 * deg) - CoatingT};
  const vector<double> DownCoatingR2Array = {QuartzR1 * sin(75 * deg),
                                             QuartzR1 * sin(75 * deg)};

  Polyhedra SolidCoatingUp(3, 0, 90 * deg, UpCoatingTArray, UpCoatingR1Array,
                           UpCoatingR2Array);
  Volume CoatingUpVol(x_coating.nameStr() + "Up", SolidCoatingUp,
                      description.material(x_coating.materialStr()));
  CoatingUpVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingUp = SectorVol.placeVolume(
      CoatingUpVol, 0, Position(interval, interval, QuartzZ));
  pvCoatingUp.addPhysVolID(x_coating.nameStr() + "Up", 0);

  Polyhedra SolidCoatingDown(3, 0, 90 * deg, DownCoatingTArray,
                             DownCoatingR1Array, DownCoatingR2Array);
  Volume CoatingDownVol(x_coating.nameStr() + "Down", SolidCoatingDown,
                        description.material(x_coating.materialStr()));
  CoatingDownVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingDown = SectorVol.placeVolume(
      CoatingDownVol, 0, Position(interval, interval, QuartzZ));
  pvCoatingDown.addPhysVolID(x_coating.nameStr() + "Down", 0);

  Box SolidCoatingSideBox((QuartzR2 - QuartzR1) / 2, CoatingT / 2, QuartzT / 2);
  UnionSolid SolidCoatingSide(
      SolidCoatingSideBox, SolidCoatingSideBox,
      Transform3D(RotationZYX(90 * deg, 0, 0),
                  Position(-QuartzR1 / 2 - QuartzR2 / 2 - CoatingT / 2,
                           QuartzR1 / 2 + QuartzR2 / 2 + CoatingT / 2, 0)));
  Volume CoatingSideVol(x_coating.nameStr() + "side", SolidCoatingSide,
                        description.material(x_coating.materialStr()));
  CoatingSideVol.setVisAttributes(description, x_coating.visStr());
  PlacedVolume pvCoatingSide =
      SectorVol.placeVolume(CoatingSideVol, 0,
                            Position(interval + QuartzR1 / 2 + QuartzR2 / 2,
                                     interval - CoatingT / 2, QuartzZ));
  pvCoatingSide.addPhysVolID(x_coating.nameStr(), 0);
  //-----------------------------------------------------------

  ////---------------------PhotonDet and Cathode-----------------
  xml_comp_t x_PhotonDet(x_det.child(_Unicode(PhotonDet)));
  xml_comp_t x_cathode(x_det.child(_Unicode(Cathode)));
  const int PhotonDetNu = x_PhotonDet.attr<int>(_Unicode(nPhotonDets));
  double PhotonDetW = x_PhotonDet.width();
  double PhotonDetT = x_PhotonDet.thickness();
  double PhotonDetH = x_PhotonDet.height();
  double PhotonDetR = QuartzR2 * sin(75 * deg) - PhotonDetH / 2.;
  double QuartzW2 = PhotonDetR * tan(360 * deg / SectorNu / 3. / 2.) * 2.;
  double PhotonDetMargin = PhotonDetH / 2 * tan(360 * deg / SectorNu / 3. / 2.);
  double PhotonDetGap =
      (QuartzW2 - PhotonDetMargin * 2. - PhotonDetW * PhotonDetNu) /
      PhotonDetNu;
  PhotonDetMargin += PhotonDetGap / 2.;
  const int CathodeNuR = x_cathode.attr<int>(_Unicode(Nr));
  const int CathodeNuC = x_cathode.attr<int>(_Unicode(Nc));
  double CathodeW = x_cathode.width();
  double CathodeT = x_cathode.thickness();
  double CathodeH = x_cathode.height();

  Box SolidPhotonDet(PhotonDetW / 2., PhotonDetT / 2., PhotonDetH / 2.);
  Volume PhotonDetVol(x_PhotonDet.nameStr(), SolidPhotonDet,
                      description.material(x_PhotonDet.materialStr()));
  PhotonDetVol.setVisAttributes(description, x_PhotonDet.visStr());
  double PDW, theta0, posX0, posY0;
  PlacedVolume pvPhotonDet[PhotonDetNu * 3];
  for (int i = 0; i < PhotonDetNu * 3; i++) {
    if (i < PhotonDetNu * 1) {
      theta0 = 75 * deg;
      posY0 = interval + PhotonDetR / sin(75 * deg);
      posX0 = interval;
    } else if (i < PhotonDetNu * 2) {
      theta0 = 45 * deg;
      posY0 = interval + PhotonDetR / sin(75 * deg) - QuartzW2 * cos(75 * deg);
      posX0 = interval + QuartzW2 * sin(75 * deg);
    } else {
      theta0 = 15 * deg;
      posY0 = interval + PhotonDetR / sin(75 * deg) - QuartzW2 * cos(75 * deg) -
              QuartzW2 * cos(45 * deg);
      posX0 = interval + QuartzW2 * sin(75 * deg) + QuartzW2 * sin(45 * deg);
    }
    PDW = PhotonDetMargin + PhotonDetW / 2. +
          (PhotonDetW + PhotonDetGap) * (i % PhotonDetNu);

    double PhotonDet_posX = posX0 + PDW * sin(theta0);
    double PhotonDet_posY = posY0 - PDW * cos(theta0);
    double PhotonDet_posZ = QuartzZ + PhotonDetT / 2. + QuartzT / 2.;

    pvPhotonDet[i] = SectorVol.placeVolume(
        PhotonDetVol, i,
        Transform3D(RotationZYX(0, 90 * deg - theta0, -90 * deg),
                    Position(PhotonDet_posX, PhotonDet_posY, PhotonDet_posZ)));
    pvPhotonDet[i].addPhysVolID(x_PhotonDet.nameStr(), i);
  }

  Box SolidCathode(CathodeW / 2., CathodeT / 2., CathodeH / 2.);
  Volume CathodeVol(x_cathode.nameStr(), SolidCathode,
                    description.material(x_cathode.materialStr()));
  CathodeVol.setVisAttributes(description, x_cathode.visStr());
  PlacedVolume pvCathode[CathodeNuR * CathodeNuC];
  for (int i = 0; i < CathodeNuR; i++)
    for (int j = 0; j < CathodeNuC; j++) {
      double Cathode_posX = -CathodeW * (CathodeNuR - 1) / 2. + CathodeW * i;
      double Cathode_posY = PhotonDetT / 2. - CathodeT / 2.;
      double Cathode_posZ = CathodeH * (CathodeNuC - 1) / 2. - CathodeH * j;
      pvCathode[CathodeNuC * i + j] = PhotonDetVol.placeVolume(
          CathodeVol, CathodeNuC * i + j,
          Position(Cathode_posX, Cathode_posY, Cathode_posZ));
      pvCathode[CathodeNuC * i + j].addPhysVolID(x_cathode.nameStr(),
                                                 CathodeNuC * i + j);
    }
    //-----------------------------------------------------------

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 17, 0)
  //---------------------OpticalSurface------------------------
  OpticalSurfaceManager surfMgr = description.surfaceManager();

  OpticalSurface QuartzAirSurf = surfMgr.opticalSurface("QuartzAirSurface");
  for (int i = 0; i < SectorNu * 2; i++) {
    BorderSurface(description, Ftof, _toString(i, "QuartzAirSurf_%d"),
                  QuartzAirSurf, pvQuartz, pvSector[i]);
  }

  OpticalSurface QuartzCoatingUpSurf =
      surfMgr.opticalSurface("QuartzCoatingUpSurface");
  BorderSurface(description, Ftof, "QuartzCoatingUpSurface",
                QuartzCoatingUpSurf, pvQuartz, pvCoatingUp);

  OpticalSurface QuartzCoatingDownSurf =
      surfMgr.opticalSurface("QuartzCoatingDownSurface");
  BorderSurface(description, Ftof, "QuartzCoatingDownSurface",
                QuartzCoatingDownSurf, pvQuartz, pvCoatingDown);

  OpticalSurface QuartzCoatingSideSurf =
      surfMgr.opticalSurface("QuartzCoatingSideSurface");
  BorderSurface(description, Ftof, "QuartzCoatingSideSurface",
                QuartzCoatingSideSurf, pvQuartz, pvCoatingSide);

  OpticalSurface QuartzCathodeSurf =
      surfMgr.opticalSurface("QuartzCathodeSurface");
  for (int i = 0; i < CathodeNuR * CathodeNuC; i++)
    BorderSurface(description, Ftof, _toString(i, "QuartzCathodeSurface_%d"),
                  QuartzCathodeSurf, pvQuartz, pvCathode[i]);
#endif

  return Ftof;
}
DECLARE_DETELEMENT(StcfFtofGeoV1, create_element)
