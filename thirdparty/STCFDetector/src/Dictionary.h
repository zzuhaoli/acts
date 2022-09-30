#ifndef DICTIONARY_H
#define DICTIONARY_H

#include "RecGeoData.h"

namespace {
  class Dictionary {};
}

#if defined(__CINT__) || defined(__MAKECINT__) || defined(__CLING__) || defined(__ROOTCLING__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace dd4hep;
#pragma link C++ namespace dd4hep::detail;
#pragma link C++ namespace dd4hep::rec;

using namespace dd4hep;
using namespace dd4hep::rec;

#pragma link C++ class MDCRecGeoStruct+;
#pragma link C++ class MDCRecGeoStruct::LayerLayout+;
#pragma link C++ class StructExtension<MDCRecGeoStruct>+;

#pragma link C++ class ECALBRecGeoStruct+;
#pragma link C++ class ECALBRecGeoStruct::RowPos+;
#pragma link C++ class ECALBRecGeoStruct::StavePos+;
#pragma link C++ class StructExtension<ECALBRecGeoStruct>+;
#pragma link C++ class ECALERecGeoStruct+;
#pragma link C++ class ECALERecGeoStruct::LeftRowPos+;
#pragma link C++ class ECALERecGeoStruct::RightRowPos+;
#pragma link C++ class StructExtension<ECALERecGeoStruct>+;

#endif

#endif /*DICTIONARY_H  */
