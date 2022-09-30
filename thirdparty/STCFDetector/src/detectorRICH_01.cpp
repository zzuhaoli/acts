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
#include "DD4hep/DetType.h"
#include "DDRec/Surface.h"
#include <exception>

using namespace dd4hep;
using namespace dd4hep::detail;
using namespace dd4hep::rec;

static Ref_t create_element(Detector &description, xml_h e,
                            SensitiveDetector sens) {

  xml_det_t x_det = e;
  std::string name = x_det.nameStr();

  // put the whole detector into an assembly
  //  - should be replaced by an envelope volume ...

  Assembly assembly(name + "_assembly");

  DetElement tracker(name, x_det.id());
  DetType type( DetType::TRACKER);
  tracker.setTypeFlag(type.to_ulong());


  PlacedVolume pv;

  // rec::ZPlanarData*  zPlanarData = new rec::ZPlanarData ;

  double minRadius = 1e99;
  double minZhalf = 1e99;

  bool isStripDetector = false;
  try {
    isStripDetector = x_det.attr<bool>(_Unicode(isStripDetector));

  } catch (std::runtime_error) {
  }

  //=========  loop over layer elements in xml
  //======================================

  for (xml_coll_t c(e, _U(layer)); c; ++c) {

    xml_comp_t x_layer(c);

    int layer_id = x_layer.id();
    int nLadders = x_layer.attr<double>(_Unicode(nLadders));
    int interval = x_det.attr<double>(_Unicode(interval));
    int initial = x_det.attr<double>(_Unicode(initial));

    double dphi = 2. * M_PI / double(nLadders);

    std::string layername = name + _toString(layer_id, "_layer%d");

    // --- create an assembly and DetElement for the layer

    Assembly layer_assembly("layer_assembly" + _toString(layer_id, "_%d"));

    DetElement layerDE(tracker, _toString(layer_id, "layer_%d"), x_det.id());

    pv = assembly.placeVolume(layer_assembly);

    pv.addPhysVolID("layer", layer_id);

    layerDE.setPlacement(pv);

    double phi0 = x_layer.phi0();
    // -------- create a measurement plane for the tracking surface attched to
    // the sensitive volume -----

    //    Vector3D o( 0. , 0. , 0. ) ;
    Vector3D u(0., 1., 0.);
    Vector3D v(0., 0., 1.);
    Vector3D n(1., 0., 0.);
    // child elements: ladder and sensitive
    for (xml_coll_t i(c, _U(sensitive)); i; ++i) {
      xml_comp_t x_sensitive(i);
      int sensitive_id = x_sensitive.id();
      std::string sensitivename =
          layername + _toString(sensitive_id, "_sensitive%d");

      double sens_zhalf = x_sensitive.length();
      double sens_offset = x_sensitive.offset();
      double sens_distance = x_sensitive.distance();
      double sens_thickness = x_sensitive.thickness();
      double sens_width = x_sensitive.width();
      // double sens_initial   = x_sensitive.initial();
      double sens_initial = x_sensitive.attr<double>(_Unicode(initial));

      std::string sens_vis = x_sensitive.visStr();
      std::string sens_matS = x_sensitive.materialStr();

      if (sens_distance < minRadius)
        minRadius = sens_distance;
      if (sens_zhalf < minZhalf)
        minZhalf = sens_zhalf;

      Material sens_mat = description.material(sens_matS);
      Box sens_box(sens_thickness / 2., sens_width / 2., sens_zhalf / 2.);

      Volume sens_vol(layername + "_sens", sens_box, sens_mat);
      sens.setType("tracker");
      sens_vol.setSensitiveDetector(sens);

      sens_vol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(),
                             sens_vis);
      for (int j = 0; j < nLadders; ++j) {

        double phi = phi0 + j * dphi;
        RotationZYX rot(phi, 0, 0);
        std::string nsensitivename =
            sensitivename + _toString(j, "_nsensitive%d");
        // --- place sensitive -----
        double lthick = sens_thickness;
        double radius = sens_distance;
        double offset = sens_offset;

        pv = layer_assembly.placeVolume(
            sens_vol,
            Transform3D(
                rot,
                Position((radius + lthick / 2.0) * cos(phi) - offset * sin(phi),
                         (radius + lthick / 2.0) * sin(phi) + offset * cos(phi),
                         (sens_initial))));
        pv.addPhysVolID("module", j).addPhysVolID("sensor", 0);
        DetElement sensitiveDE(layerDE, nsensitivename, x_det.id());
        sensitiveDE.setPlacement(pv);
      }
    }

    for (xml_coll_t k(c, _U(ladder)); k; ++k) {

      xml_comp_t x_ladder(k);
      int ladder_id = x_ladder.id();

      std::string laddername = layername + _toString(ladder_id, "_ladder%d");

      //--------------------------------

      double supp_zhalf = x_ladder.attr<double>(_Unicode(length));
      double supp_offset = x_ladder.attr<double>(_Unicode(offset));
      double supp_distance = x_ladder.attr<double>(_Unicode(distance));
      double supp_thickness = x_ladder.attr<double>(_Unicode(thickness));
      double supp_width = x_ladder.attr<double>(_Unicode(width));
      double supp_initial = x_ladder.attr<double>(_Unicode(initial));

      std::string supp_vis = x_ladder.visStr();
      std::string supp_matS = x_ladder.materialStr();

      if (supp_distance < minRadius)
        minRadius = supp_distance;

      if (supp_zhalf < minZhalf)
        minZhalf = supp_zhalf;

      //-----------------------------------
      //  store the data in an extension to be used for reconstruction
      // rec::ZPlanarData::LayerLayout thisLayer ;

      // thisLayer.sensorsPerLadder =  1 ; // for now only one planar sensor
      // thisLayer.lengthSensor     =  2. * sens_zhalf ;
      //
      // thisLayer.distanceSupport  = supp_distance ;
      // thisLayer.offsetSupport    = supp_offset ;
      // thisLayer.thicknessSupport = supp_thickness ;
      // thisLayer.zHalfSupport     = supp_zhalf ;
      // thisLayer.widthSupport     = supp_width ;
      //
      // thisLayer.distanceSensitive  = sens_distance ;
      // thisLayer.offsetSensitive    = sens_offset ;
      // thisLayer.thicknessSensitive = sens_thickness ;
      // thisLayer.zHalfSensitive     = sens_zhalf ;
      // thisLayer.widthSensitive     = sens_width ;
      //
      // thisLayer.ladderNumber =  nLadders ;
      // thisLayer.phi0         =  phi0 ;
      //
      // zPlanarData->layers.push_back( thisLayer ) ;
      //-----------------------------------

      Material supp_mat = description.material(supp_matS);

      //-------

      Box supp_box(supp_thickness / 2., supp_width / 2., supp_zhalf / 2.);

      Volume supp_vol(layername + "_supp", supp_box, supp_mat);

      // compute the inner and outer thicknesses that need to be assigned to the
      // tracking surface depending on wether the support is above or below the
      // sensor
      // double inner_thickness = ( sens_distance > supp_distance ?  (
      // sens_distance - supp_distance ) + sens_thickness/2  : sens_thickness/2
      // ) ; double outer_thickness = ( sens_distance > supp_distance ?
      // sens_thickness/2  :  ( supp_distance - sens_distance ) + supp_thickness
      // - sens_thickness/2   ) ;

      // SurfaceType type( SurfaceType::Sensitive ) ;

      // if( isStripDetector )
      //   type.setProperty( SurfaceType::Measurement1D , true ) ;

      // VolPlane surf( sens_vol , type , inner_thickness , outer_thickness ,
      // u,v,n ) ; //,o ) ;

      //--------------------------------------------

      supp_vol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(),
                             supp_vis);

      //--------- loop over ladders ---------------------------

      for (int j = 0; j < nLadders; ++j) {

        double phi = phi0 + j * dphi;

        std::string nladdername = laddername + _toString(j, "_nladder%d");

        RotationZYX rot(phi, 0, 0);

        // --- place support -----
        double lthick = supp_thickness;
        double radius = supp_distance;
        double offset = supp_offset;

        /* pv = */ layer_assembly.placeVolume(
            supp_vol,
            Transform3D(
                rot,
                Position((radius + lthick / 2.) * cos(phi) - offset * sin(phi),
                         (radius + lthick / 2.) * sin(phi) + offset * cos(phi),
                         supp_initial)));

        //      pv.addPhysVolID("layer", layer_id ).addPhysVolID( "module" , j
        //      ).addPhysVolID("sensor", 0 )   ;

        pv.addPhysVolID("module", j).addPhysVolID("sensor", 0);
        DetElement ladderDE(layerDE, nladdername, x_det.id());
        ladderDE.setPlacement(pv);
      }
    }

    //  volSurfaceList( ladderDE )->push_back( surf ) ;

    //    tracker.setVisAttributes(description, x_det.visStr(),laddervol);

    // is this needed ??
    layer_assembly->GetShape()->ComputeBBox();
  }

#if 0 //-------- add an inscribing cylinder of air for tracking purposes
      //-----------------
      //  this screw up the geometry and the material scan does not work anymore
      //  !!!!!?????

	double tube_thick =  1.0 * dd4hep::mm ;
	double inner_r    =  minRadius - 1.1 * tube_thick ;
	double outer_r    =  inner_r + tube_thick ;
	double z_half     =  minZhalf ; 

	Tube   tubeSolid (inner_r, outer_r, z_half ) ;
	Volume tube_vol( name+"_inner_cylinder_air", tubeSolid ,  description.material("Air") ) ;

	assembly.placeVolume( tube_vol , Transform3D() ) ;

	Vector3D ocyl(  inner_r + 0.5*tube_thick , 0. , 0. ) ;

	VolCylinder cylSurf( tube_vol , SurfaceType( SurfaceType::Helper ) , 0.5*tube_thick  , 0.5*tube_thick , ocyl ) ;

	volSurfaceList( tracker )->push_back( cylSurf ) ;

#endif //----------------------------------------------------------------------------------

  // tracker.addExtension< rec::ZPlanarData >( zPlanarData ) ;

  Volume mother = description.pickMotherVolume(tracker);

  pv = mother.placeVolume(assembly);

  pv.addPhysVolID("system", x_det.id()).addPhysVolID("side", 0);

  tracker.setPlacement(pv);

  assembly->GetShape()->ComputeBBox();

  return tracker;
}

DECLARE_DETELEMENT(BarrelDIRC, create_element)
