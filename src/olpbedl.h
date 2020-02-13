#ifndef OVERLAP_BEDS_H
#define OVERLAP_BEDS_H

#include <CLI.hpp>
#include <map>
#include <string>
#include <vector>
#include <zlib.h>
#include <util.h>
#include <kseq.h>
#include "cgranges.h"

KSEQ_INIT(gzFile, gzread)

typedef struct {
    int32_t l;
    char *s;
}bed_rest1_t;

typedef struct {
    int64_t n, m;
    bed_rest1_t *a;
}bed_rest_t;

struct OlpStore{
    cgranges_t* olpret;
    std::vector<int32_t> idx;
};

class BedRegs{
    public:
        std::vector<cgranges_t*> mRegs;     ///< cgranges_t
        std::vector<std::string> mNames;    ///< bed names
        std::vector<std::string> mBedPaths; ///< bed paths
        std::string outdir;                 ///< output directory
        std::string bedlist;                ///< bed list file
        
        /** constructor of BedRegs */
        BedRegs(){
            outdir = "anaout";
        }

        /** destructor of BedRegs */
        ~BedRegs(){}

    public:
        cgranges_t* loadOneBed(const std::string& bedFile);
        void loadBeds();
        void olpAna();

    private:
        static char *parse_bed3b(char *s, int32_t *st_, int32_t *en_, char **r);
        static char *parse_bed3(char *s, int32_t *st_, int32_t *en_);
        static cgranges_t *read_bed3b(const char *fn, bed_rest_t *r);
        static cgranges_t *read_bed3(const char *fn);

};
#endif
