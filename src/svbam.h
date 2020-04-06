#ifndef SVBAM_H
#define SVBAM_H

#include <iostream>
#include <string>
#include <cassert>
#include <vector>
#include <htslib/sam.h>
#include "software.h"
#include "onlinebwa.h"
#include <util.h>
#include <bamutil.h>

struct SVBAMOpt{
    std::string ibam;
    int32_t svid;
    std::string obam;
    std::string obam2;
    std::string ref;
    bool outalt;
    Software* softEnv;
    OnlineBWA* obwa;
    
    SVBAMOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->ver + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        obam = "./svi.bam";
        obam2 = "";
        ref = "";
        outalt = false;
        obwa = NULL;
    }

    ~SVBAMOpt(){
        if(softEnv) delete softEnv;
        if(obwa) delete obwa;
    }

    void update(int argc, char** argv){
        softEnv->cwd = util::cwd();
        for(int i = 0; i < argc; ++i){
            softEnv->cmd.append(argv[i]);
            softEnv->cmd.append(" ");
        }
        if(obam2.size() && ref.size()){
            outalt = true;
            obwa = new OnlineBWA();
            obwa->loadIndex(ref);
        }
    }

    void getBam();
};

#endif
