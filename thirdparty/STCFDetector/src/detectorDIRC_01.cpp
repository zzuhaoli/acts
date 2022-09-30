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

  Assembly assembly(name + "_assembly");
  DetElement DIRC(name, x_det.id());
  PlacedVolume pv;

  for (xml_coll_t c(e, _U(layer)); c; ++c) {
    xml_comp_t x_layer(c);

    xml_comp_t x_quartz(x_layer.child(_Unicode(sensitive)));
    xml_comp_t x_PMT(x_layer.child(_Unicode(sensPMT)));
    xml_comp_t x_MCP(x_layer.child(_Unicode(sensMCP)));

    int layer_id = x_layer.id();
    int nLadders = x_layer.attr<int>(_Unicode(nLadders));
    double layer_Z = x_layer.attr<double>(_Unicode(positionZ));
    std::string layername = name + _toString(layer_id, "_layer%d");

    Assembly layer_assembly("layer_assembly" + _toString(layer_id, "_%d"));
    DetElement layerDE(DIRC, _toString(layer_id, "layer_%d"), x_det.id());
    pv = assembly.placeVolume(layer_assembly);
    pv.addPhysVolID("layer", layer_id);
    layerDE.setPlacement(pv);

    double quartz_distance = x_quartz.distance();
    double quartz_thickness = x_quartz.thickness();
    double quartz_width1 = x_quartz.attr<double>(_Unicode(width1));
    double quartz_width2 = x_quartz.attr<double>(_Unicode(width2));
    double quartz_length = x_quartz.length();
    std::string quartz_vis = x_quartz.visStr();
    std::string quartz_mats = x_quartz.materialStr();
    Material quartz_mat = description.material(quartz_mats);

    int nPMT = x_PMT.attr<int>(_Unicode(nPMT));
    double PMT_thickness = x_PMT.thickness();
    double PMT_length = x_PMT.length();
    double PMT_width = x_PMT.width();
    double PMT_gap = x_PMT.gap();
    std::string PMT_mats = x_PMT.materialStr();
    std::string PMT_vis = x_PMT.visStr();
    Material PMT_mat = description.material(PMT_mats);

    int nMCP = x_MCP.attr<int>(_Unicode(nMCP));
    double MCP_thickness = x_MCP.thickness();
    double MCP_length = x_MCP.length();
    double MCP_width = x_MCP.width();
    std::string MCP_mats = x_MCP.materialStr();
    std::string MCP_vis = x_MCP.visStr();
    Material MCP_mat = description.material(MCP_mats);

    for (int j = 0; j < nLadders; ++j) {
      double phi0 = x_layer.phi0();
      double dphi = 2. * M_PI / double(nLadders);
      double phi = phi0 + j * dphi;
      double y = layer_Z;
      double z = -quartz_distance - quartz_length / 2.0;
      std::string sectorname = layername + _toString(j, "_sector%d");

      Assembly sector_assembly("sector_assembly" + _toString(j, "_%d"));
      pv = layer_assembly.placeVolume(
          sector_assembly,
          Transform3D(RotationZYX(0, -phi, 0),
                      Position(-z * sin(phi), y, z * cos(phi))));
      pv.addPhysVolID("sector", j);
      DetElement sectorDE(layerDE, sectorname, x_det.id());
      sectorDE.setPlacement(pv);

      Trapezoid quartzTrd(quartz_width1 / 2., quartz_width2 / 2.,
                          quartz_thickness / 2., quartz_thickness / 2.,
                          quartz_length / 2.);
      Volume quartzVol(layername + "_quartz", quartzTrd, quartz_mat);
      quartzVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(),
                              quartz_vis);
      pv = sector_assembly.placeVolume(quartzVol);
      pv.addPhysVolID("quartz", j);
      DetElement quartzDE(sectorDE, sectorname, x_det.id());
      quartzDE.setPlacement(pv);

      for (int k = 0; k < nPMT; k++) {
        double x0 = -(PMT_gap + PMT_width) * (nPMT - 1) / 2;
        double dx = PMT_gap + PMT_width;
        double y = layer_Z > 0 ? (quartz_thickness + PMT_thickness) / 2
                               : -(quartz_thickness + PMT_thickness) / 2;
        double z = -quartz_length / 2 + PMT_length / 2;
        double x = x0 + dx * k;

        Box PMTBox(PMT_width / 2.0, PMT_thickness / 2.0, PMT_length / 2.0);
        Volume PMTVol(layername + "_PMT", PMTBox, PMT_mat);
        PMTVol.setAttributes(description, x_det.regionStr(), x_det.limitsStr(),
                             PMT_vis);

        pv = sector_assembly.placeVolume(PMTVol, Position(x, y, z));
        pv.addPhysVolID("module", k);
        DetElement PMTDE(sectorDE, sectorname + _toString(k, "_module%d"),
                         x_det.id());
        PMTDE.setPlacement(pv);

        for (int l = 0; l < nMCP; ++l) {
          double x0 = -MCP_width * (nMCP - 1) / 2;
          double dx = MCP_width;
          double y = layer_Z > 0 ? (-PMT_thickness + MCP_thickness) / 2
                                 : (PMT_thickness - MCP_thickness) / 2;
          double x = x0 + dx * l;

          Box MCPBox(MCP_width / 2.0, MCP_thickness / 2.0, MCP_length / 2.0);
          Volume MCPVol(layername + "MCP", MCPBox, MCP_mat);
          MCPVol.setAttributes(description, x_det.regionStr(),
                               x_det.limitsStr(), MCP_vis);

          pv = PMTVol.placeVolume(MCPVol, Position(x, y, 0));
          pv.addPhysVolID("senser", l);
          DetElement MCPDE(PMTDE,
                           sectorname + _toString(k, "_module%d") +
                               _toString(l, "_senser%d"),
                           x_det.id());
          MCPDE.setPlacement(pv);
        }
      }
    }
  }

  Volume mother = description.pickMotherVolume(DIRC);
  // pv=mother.placeVolume(assembly,RotationZYX(M_PI/2.0,0,0));
  pv = mother.placeVolume(assembly, RotationZYX(0, M_PI / 2.0, 0));
  pv.addPhysVolID("system", x_det.id());
  DIRC.setPlacement(pv);
  assembly->GetShape()->ComputeBBox();
  return DIRC;
}
DECLARE_DETELEMENT(DIRCPID, create_element)
