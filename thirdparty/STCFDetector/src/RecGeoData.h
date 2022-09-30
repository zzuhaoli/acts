#ifndef REC_DATA_H_
#define REC_DATA_H_

#include <map>
#include <bitset>
#include <ostream>

#include "DD4hep/DetElement.h"
#include "DDRec/DetectorData.h"

namespace dd4hep{
 namespace rec{
   struct MDCRecGeoStruct{

   int nLayer ; // The number of Layer
   int nSuperLayer ; // The number of SuperLayer
   int nWire ; // The number of Wire
   double normal_l;//normalized length in z direction

   struct LayerLayout{
   //  default c'tor with zero initialization
    LayerLayout() :
    SuperL_id(0),
    Layer_id(0),
    CellPerLayer(0),
    inner_Radius(0),
    outer_Radius(0),
    Zlength(0),
    SwireR_z0(0),
    Twist(0),
    Stereo(0),
    Type(0),
    StartShift(0),
    LayerInSup(0){
    }
                    
    int SuperL_id ;  
    int Layer_id ;
    int CellPerLayer ;
    double inner_Radius ;
    double outer_Radius ;
    double Zlength ;
    double SwireR_z0 ;
    int Twist ;
    double Stereo ;
    int Type ;
    double StartShift ;
    int LayerInSup ;
    };
                                                  
    std::vector<LayerLayout> layers ;
    
    };
    typedef StructExtension<MDCRecGeoStruct> MDCRecGeoData ;

    struct ECALBRecGeoStruct{

      struct RowPos{
        RowPos() :
          x(0),
          y(0),
          z(0){

          }

        double x;
        double y;
        double z;
      };

      struct StavePos{
        StavePos() :
          x(0),
          y(0),
          z(0){

          }

          double x;
          double y;
          double z;
      };

      std::vector<double> Theta;
      std::vector<RowPos> Rows;
      std::vector<StavePos> Staves;
    };

    typedef StructExtension<ECALBRecGeoStruct> ECALBRecGeoData ;

    struct ECALERecGeoStruct{

      struct LeftRowPos{
        LeftRowPos() :
          x(0),
          y(0),
          z(0){

          }

        double x;
        double y;
        double z;
      };
      
      struct RightRowPos{
        RightRowPos() :
          x(0),
          y(0),
          z(0){

          }

        double x;
        double y;
        double z;
      };

      std::vector<double> Theta;
      std::vector<double> DTheta;
      std::vector<LeftRowPos> LeftRows;
      std::vector<RightRowPos> RightRows;
    };

    typedef StructExtension<ECALERecGeoStruct> ECALERecGeoData ;
                                                              
    std::ostream& operator<<( std::ostream& io , const MDCRecGeoData& d ) ;
  }
}
#endif
