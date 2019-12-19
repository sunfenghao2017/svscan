#ifndef SVFILTER_H
#define SVFILTER_H

#include <map>
#include <set>
#include <vector>
#include <cstdint>
#include <string>
#include <htslib/vcf.h>
#include "software.h"
#include "svinfo.h"
#include "svutil.h"
#include "CLI.hpp"
#include "util.h"

typedef std::vector<std::vector<SVInfo>> SVList;

/** class to store SV filter options */
struct SVFilter{
    int32_t minsr = 3;       ///< min SR support for an valid SV
    int32_t mindp = 3;       ///< min DP support for an valid SV
    int32_t minsu = 3;       ///< min SU support for an valid SV
    float minaf = 0.05;            ///< min VAF for an valid SV
    int32_t maxBpOffset = 10;       ///< max breakpoint offset of an SV against background SV to be excluded
    std::string mBgBCF;              ///< background BCF file
    std::string infile;             ///< input file of svscan sv tsv format result file
    std::string outfile = "svf.tsv"; ///< output file of reported SV
    SVList bgSVs;                   ///< to store background SVs
    bool initialized = false;       ///< SVFilter is initialized if true
    Software* softEnv;               ///< software environment

    /** SVFilter constructor */
    SVFilter();

    /** SVFilter destructor */
    ~SVFilter();

    /** update options after command arguments parsed
    * @param argc number of arguments feed to main
    * @param argv array of arguments feed to main
    */
    void update(int argc, char** argv);

    /** initialize filter options */
    void init();

    /** test whether an SV breakpoint is not in background
     * @param svt SV type
     * @param chr1 big chr of SV
     * @param chr2 lite chr of SV
     * @param start starting position
     * @param end ending position
     * @return true if this SV breakpoint it not in background
     */
    bool validSV(int32_t svt, const std::string& chr1, const std::string& chr2, int32_t start, int32_t end);

    /** get background SV events */
    void getBgSVs();

    /** filter and get result */
    void filter();
};

#endif
