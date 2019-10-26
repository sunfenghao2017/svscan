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

struct SVDNADBOpt{
    std::string refSeqDB;   ///< refseq database with transcript version added
    std::string gene2cnc; ///< gene to canonical transcripts list
    std::string svAnnoDB;   ///< sv annotation database output file
    Software* softEnv;      ///< software env
    
    SVDNADBOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->ver + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        svAnnoDB = "./dnaAnno.gz";
    }

    ~SVDNADBOpt(){
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
