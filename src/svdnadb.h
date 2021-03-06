#ifndef SVANNO_H
#define SVANNO_H

#include <map>
#include <set>
#include <util.h>
#include <stdio.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include "software.h"

struct SVDNADBOpt{
    std::string refSeqDB;   ///< refseq database with transcript version added
    std::string gene2cnc;   ///< gene to canonical transcripts list
    std::string svAnnoDB;   ///< sv annotation database output file
    std::string ref2gene;   ///< refseq accession number to gene name
    std::string bedCDS;     ///< cds bed region file of each transcript
    std::string bedUnits;   ///< exon/utr region file of each transcript
    std::string geneRng;    ///< gene canonical transcript coordinate range
    Software* softEnv;      ///< software env
    
    SVDNADBOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->ver + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        svAnnoDB = "anno.gz";
        ref2gene = "ref2gene.tsv";
        bedCDS = "ref.cds.bed";
        bedUnits = "ref.unit.bed";
        geneRng = "gene.coord.tsv";
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
