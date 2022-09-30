#include "DD4hep/DetFactoryHelper.h"

using namespace dd4hep;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {
  xml_det_t x_det = e;
  std::string name = x_det.nameStr();
  Assembly assembly(name + "_assembly");
  DetElement DriftLayer(name, x_det.id());
  PlacedVolume pv;

  std::string visattr = x_det.attr<std::string>(_Unicode(vis));
  std::string supmat = x_det.attr<std::string>(_Unicode(material));
  // end cap support, divide to two parts
  xml_comp_t endcap = x_det.child(_Unicode(endcap));
  // inner region, a Cone shape
  xml_comp_t InnerCone = endcap.child(_Unicode(innercone));

  double Conez = InnerCone.attr<double>(_Unicode(z));
  double Coner1min = InnerCone.attr<double>(_Unicode(r1_min));
  double Coner1max = InnerCone.attr<double>(_Unicode(r1_max));
  double Coner2min = InnerCone.attr<double>(_Unicode(r2_min));
  double Coner2max = InnerCone.attr<double>(_Unicode(r2_max));
  double Positionmin = InnerCone.attr<double>(_Unicode(positionmin));
  Cone ConeInner(Conez / 2, Coner1min, Coner1max, Coner2min, Coner2max);
  Volume ConeInner_vol(name + "endcap_support", ConeInner,
                       description.material(supmat));
  ConeInner_vol.setVisAttributes(description.visAttributes(visattr));
  assembly.placeVolume(ConeInner_vol,
                       Transform3D(RotationZYX(0, 0., 0),
                                   Position(0, 0., Positionmin + Conez / 2)));
  assembly.placeVolume(
      ConeInner_vol, Transform3D(RotationZYX(0.0, M_PI, 0),
                                 Position(0, 0., -(Positionmin + Conez / 2))));

  // outer region, a Cone shape
  if (endcap.hasChild(_Unicode(outercone))) {
    xml_comp_t OuterCone = endcap.child(_Unicode(outercone));

    Conez = OuterCone.attr<double>(_Unicode(z));
    Coner1min = OuterCone.attr<double>(_Unicode(r1_min));
    Coner1max = OuterCone.attr<double>(_Unicode(r1_max));
    Coner2min = OuterCone.attr<double>(_Unicode(r2_min));
    Coner2max = OuterCone.attr<double>(_Unicode(r2_max));
    Positionmin = OuterCone.attr<double>(_Unicode(positionmin));
    Cone ConeOuter(Conez / 2, Coner1min, Coner1max, Coner2min, Coner2max);
    Volume ConeOuter_vol(name + "endcap_support", ConeOuter,
                         description.material(supmat));
    ConeOuter_vol.setVisAttributes(description.visAttributes(visattr));
    assembly.placeVolume(ConeOuter_vol,
                         Transform3D(RotationZYX(0, 0., 0),
                                     Position(0, 0., Positionmin + Conez / 2)));
    assembly.placeVolume(
        ConeOuter_vol,
        Transform3D(RotationZYX(0.0, M_PI, 0),
                    Position(0, 0., -(Positionmin + Conez / 2))));
  }

  // outer region, a tube shape
  if (endcap.hasChild(_Unicode(outertube))) {
    xml_comp_t OuterTube = endcap.child(_Unicode(outertube));
    double EndcapTubePosition = OuterTube.attr<double>(_Unicode(position));
    double EndcapTubeLength = OuterTube.attr<double>(_Unicode(length));
    double EndcapTubeRmin = OuterTube.attr<double>(_Unicode(rmin));
    double EndcapTubeRmax = OuterTube.attr<double>(_Unicode(rmax));
    Tube tubeSolid_endcap(EndcapTubeRmin, EndcapTubeRmax,
                          EndcapTubeLength / 2.0);
    Volume tubeSolid_endcap_vol(name + "endcap_support", tubeSolid_endcap,
                                description.material(supmat));
    tubeSolid_endcap_vol.setVisAttributes(description.visAttributes(visattr));
    assembly.placeVolume(
        tubeSolid_endcap_vol,
        Transform3D(
            RotationZYX(0, 0., 0),
            Position(0, 0., EndcapTubePosition + EndcapTubeLength / 2)));
    assembly.placeVolume(
        tubeSolid_endcap_vol,
        Transform3D(
            RotationZYX(0, 0., 0),
            Position(0, 0., -(EndcapTubePosition + EndcapTubeLength / 2))));
  }

  // barrel part, inner and outer
  xml_comp_t barrel = x_det.child(_Unicode(barrel));
  xml_comp_t innertube = barrel.child(_Unicode(innertube));
  double BarrelTubeLength = innertube.attr<double>(_Unicode(length));
  double BarrelTubeRmin = innertube.attr<double>(_Unicode(rmin));
  double BarrelTubeRmax = innertube.attr<double>(_Unicode(rmax));
  Tube TubeSolidInnerBarrel(BarrelTubeRmin, BarrelTubeRmax,
                            BarrelTubeLength / 2.0);
  Volume TubeSolidInnerBarrelVol(name + "barrel_support", TubeSolidInnerBarrel,
                                 description.material(supmat));
  TubeSolidInnerBarrelVol.setVisAttributes(description.visAttributes(visattr));
  assembly.placeVolume(TubeSolidInnerBarrelVol,
                       Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));

  // barrel part, inner and outer
  xml_comp_t OuterTubeBarrel = barrel.child(_Unicode(outertube));
  BarrelTubeLength = OuterTubeBarrel.attr<double>(_Unicode(length));
  BarrelTubeRmin = OuterTubeBarrel.attr<double>(_Unicode(rmin));
  BarrelTubeRmax = OuterTubeBarrel.attr<double>(_Unicode(rmax));
  Tube TubeSolidOuterBarrel(BarrelTubeRmin, BarrelTubeRmax,
                            BarrelTubeLength / 2.0);
  Volume TubeSolidOuterBarrelVol(name + "barrel_support", TubeSolidOuterBarrel,
                                 description.material(supmat));
  TubeSolidOuterBarrelVol.setVisAttributes(description.visAttributes(visattr));
  assembly.placeVolume(TubeSolidOuterBarrelVol,
                       Transform3D(RotationZYX(0, 0., 0), Position(0, 0., 0)));

  Volume mother = description.pickMotherVolume(DriftLayer);
  pv = mother.placeVolume(assembly);
  pv.addPhysVolID("system", x_det.id()).addPhysVolID("side", 0);
  DriftLayer.setPlacement(pv);

  return DriftLayer;
}
DECLARE_DETELEMENT(MDCSupport, create_element)
