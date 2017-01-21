#include <pulsar/modulemanager/ModuleCreationFuncs.hpp>
#include "FManIIPulsarAPI.hpp"

using pulsar::ModuleCreationFuncs;

extern "C" {

ModuleCreationFuncs insert_supermodule(void){
    ModuleCreationFuncs cf;
    cf.add_cpp_creator<FManII::FFTermPulsar>("FFTerm");
    return cf;
}

}
