#ifndef SVANNO_H
#define SVANNO_H

#include <map>
#include <set>
#include "util.h"
#include <stdio.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include "software.h"

struct SVAnnoDBOpt{
    std::string refGeneDB;      ///< refgene database with transcript version added
    std::string primaryTrsList; ///< canonical transcripts list
    std::string svAnnoDB;       ///< sv annotation database output file
    Software* softEnv;          ///< software env
    
    SVAnnoDBOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->version + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        svAnnoDB = "./svandb.gz";
    }

    ~SVAnnoDBOpt(){
    }

    void update(int argc, char** argv){
        softEnv->cwd = util::cwd();
        for(int i = 0; i < argc; ++i){
            softEnv->cmd.append(argv[i]);
            softEnv->cmd.append(" ");
        }
    }

    void prepDB();
};

#endif
