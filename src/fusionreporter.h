#ifndef FUSION_REPORTER_H
#define FUSION_REPORTER_H

#include "fusionopt.h"
#include "software.h"
#include "CLI.hpp"

struct FusionReporter{
    FusionOptions* fuseOpt;///< fusion report options
    Software* softEnv;///< software environments

    /** FusionReporter constructor */
    FusionReporter();

    /** FusionReporter destructor */
    ~FusionReporter();

    /** update options after command arguments parsed
    * @param argc number of arguments feed to main
    * @param argv array of arguments feed to main
    */
    void update(int argc, char** argv);
    
    /** do fusion report job */
    void report();
};

#endif
