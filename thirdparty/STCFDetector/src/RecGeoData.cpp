//#include "RecGeoData/RecGeoData.h"
#include "RecGeoData.h"

#include <boost/io/ios_state.hpp>

namespace dd4hep{
 namespace rec{  
  std::ostream& operator<<( std::ostream& io , const MDCRecGeoData& d ) {
   boost::io::ios_base_all_saver ifs(io);

   io <<  " -- MDCData: "  << std::scientific << std::endl ;
   io <<  " nLayer"   << d.nLayer << "   " ;
   io <<  " nSup"     << d.nSuperLayer << "   ";
   io <<  " nWire"    << d.nWire << std::endl ;
   
   std::vector<MDCRecGeoData::LayerLayout> layers = d.layers ;

   io <<  " superL   Layer   CellNo     Rin     Rout     Zlength     SwireR     Twist     Stereo     Type    StartShift     LayerInSup" << std::endl ;

   for(unsigned i=0,N=layers.size() ; i<N ; ++i){

    MDCRecGeoData::LayerLayout l = layers[i] ;

    io << " " 
       << l.SuperL_id << "       "
       << l.Layer_id  << "       "
       << l.CellPerLayer << "      "
       << l.inner_Radius << "      "
       << l.outer_Radius << "      "
       << l.Zlength      << "      "
       << l.SwireR_z0    << "      "
       << l.Twist        << "      "
       << l.Stereo       << "      "
       << l.Type         << "      "
       << l.StartShift   << "      "
       << l.LayerInSup   << "      "
       << std::endl ;
   }

   return io ;
  }
 }
}
