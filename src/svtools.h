#ifndef SVTOOLS_H
#define SVTOOLS_H

#include <cmath>
#include <unordered_map>
#include "fuseflags.h"
#include "svbam.h"
#include "svmerge.h"
#include "svrnadb.h"
#include "svdnadb.h"
#include "svdebug.h"
#include "svfilter.h"
#include "fusewlist.h"
#include "fusionreporter.h"

struct SVToolOpts{
    Software* softEnv;

    SVToolOpts(){
        softEnv = new Software();
    }

    ~SVToolOpts(){
        if(softEnv){
            delete softEnv;
        }
    }

    void update(int argc, char** argv){
        softEnv->cwd = util::cwd();
        for(int i = 0; i < argc; ++i){
            softEnv->cmd.append(argv[i]);
            softEnv->cmd.append(" ");
        }
    }
};

#endif
