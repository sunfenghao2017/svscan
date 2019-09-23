#ifndef FUSEWLIST_H
#define FUSEWLIST_H

#include <fstream>
#include <iostream>
#include <set>
#include <map>
#include <string>
#include "util.h"
#include "software.h"

struct FuseWOpt{
    std::string genelist;
    std::string fusedb;
    std::string whitelist;
    Software* softEnv;

    FuseWOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->version + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        whitelist = "./gwhite.tsv";
    }

    ~FuseWOpt(){}

    void update(int argc, char** argv){
        softEnv->cwd = util::cwd();
        for(int i = 0; i < argc; ++i){
            softEnv->cmd.append(argv[i]);
            softEnv->cmd.append(" ");
        }
    }

    void prepWlist();
};

#endif
