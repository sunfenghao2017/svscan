#ifndef SVRNA_H
#define SVRNA_H

#include <map>
#include <set>
#include "util.h"
#include <stdio.h>
#include <cassert>
#include <fstream>
#include <sstream>
#include <htslib/bgzf.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include "software.h"

struct SVRNADBOpt{
    std::string genome;     ///< genome fasta file
    std::string refSeqDB;   ///< refseq database with transcript version added
    std::string cncTrsList; ///< canonical transcripts list
    std::string svAnnoDB;   ///< sv annotation database output file
    std::string gene2trs;   ///< genename to used transcript output file
    std::string refMrna;    ///< refmrna fasta file
    std::string bedCDS;     ///< cds bed region file of each gene
    std::string bedUnits;   ///< exon/utr region file of each gene
    std::string utr3len;    ///< 3'utr length list of each gene
    bool keepNCRna;         ///< keep nc rna to final result if true
    Software* softEnv;      ///< software env
    
    SVRNADBOpt(){
        softEnv = new Software();
        softEnv->cmp += "version: " + softEnv->ver + "\n";
        softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
        svAnnoDB = "anno.gz";
        gene2trs = "gene2trs.tsv";
        refMrna = "ref.fa.gz";
        bedCDS = "ref.cds.bed";
        bedUnits = "ref.unit.bed";
        utr3len = "ref.3utr.len";
        keepNCRna = false;
    }

    ~SVRNADBOpt(){
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
