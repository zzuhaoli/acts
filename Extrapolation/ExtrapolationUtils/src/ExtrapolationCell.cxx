///////////////////////////////////////////////////////////////////
// ExtrapolationCell.cxx 
///////////////////////////////////////////////////////////////////

#include "ExtrapolationUtils/ExtrapolationCell.h"

std::vector<std::string> Acts::ExtrapolationCode::s_ecodeNames = { "Unset",                  
                                                                  "InProgress",
                                                                  "SuccessDestination",     
                                                                  "SuccessBoundaryReached", 
                                                                  "SuccessPathLimit", 
                                                                  "SuccessMaterialLimit",  
                                                                  "Recovered",             
                                                                  "FailureDestination",     
                                                                  "FailureLoop",     
                                                                  "FailureNavigation",      
                                                                  "FailureUpdateKill",      
                                                                  "FailureConfiguration",
                                                                  "LeftKnownWorld" };
