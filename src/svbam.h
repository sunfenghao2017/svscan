#ifndef SVBAM_H
#define SVBAM_H

#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <htslib/sam.h>
#include "software.h"
#include <util.h>

struct SVBAMOpt{
    std::string ibam;
    int32_t svid;
    std::string obam;
    Software* softEnv;
    
    SVBAMOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->ver + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        obam = "./svi.bam";
    }

    ~SVBAMOpt(){
        if(softEnv) delete softEnv;
    }

    void update(int argc, char** argv){
        softEnv->cwd = util::cwd();
        for(int i = 0; i < argc; ++i){
            softEnv->cmd.append(argv[i]);
            softEnv->cmd.append(" ");
        }
    }

    void getBam();
};

#endif
