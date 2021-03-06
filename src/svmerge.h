#ifndef MERGEVCF_H
#define MERGEVCF_H

#include <map>
#include <set>
#include <string>
#include <vector>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <htslib/vcf.h>
#include "software.h"
#include "interval.h"
#include "svutil.h"
#include <CLI.hpp>
#include <util.h>

/** type to store Interval of each contigs */
typedef std::vector<std::vector<Interval>> ContigIntervals;

/** type to store contig and id value pairs */
typedef std::map<std::string, int32_t> ContigMap;

/** class to store merge options */
struct VCFMerger{
    bool filterForPass = false;        ///< keep PASS record if true
    bool filterForPrecise = false;     ///< keep PRECISE record if true
    uint32_t chunksize = 100;          ///< max input files to merge at one time
    int32_t svcounter = 0;             ///< output SV record counter
    int32_t bpoffset = 300;            ///< max difference allowed for two duplicated SV pos
    int32_t minsize = 0;               ///< min SV size to output
    int32_t maxsize = (1 << 20);       ///< max SV size to output
    int32_t coverage = 0;              ///< min coverage for SV to output
    float recoverlap = 0.8;            ///< min overlap ratio for SV to output
    float vaf = 0.0;                   ///< min VAF for SV to output
    std::string tmpdir = "./.tmpdir";  ///< temp directory used during merging
    std::string outfile = "merge.bcf"; ///< output bcf file path
    std::string infilelist;            ///< input bcf file list
    std::vector<std::string> infiles;  ///< input bcf file list parsed
    Software* softEnv = NULL;          ///< pointer to Software options

    /** VCFMerger constructor */
    VCFMerger();

    /** VCFMerger destructor */
    ~VCFMerger();

    /** update options after command arguments parsed
    * @param argc number of arguments feed to main
    * @param argv array of arguments feed to main
    */
    void update(int argc, char** argv);

    /** get unique tid for each contigs across samples
     * @param cmap map to store contig name ~ contig id value pairs
     */
    void getContigMap(ContigMap& cmap);

    /** get Interval representation of each SV of some type
     * @param ci reference of ContigIntervals to store SV on each contig
     * @param cmap contig name ~ contig id value pairs
     * @param svt sv type
     */   
    void getIntervals(ContigIntervals& ci, ContigMap& cmap, int32_t svt);

    /** merge Interval supporting same SV
     * @param ici referene of ContigIntervals to store input Intervals of SV
     * @param oci reference of ContigIntervals to store output Interval of SV
     * @param svt sv type
     */
    void mergeIntervals(ContigIntervals& ici, ContigIntervals& oci, int32_t svt);

    /** write specific type SV in all input bcf files to output file(dup excluded)
     * @param ci valid contig intervals covering SV
     * @param cmap contig name ~ contig id value pairs
     * @param svt sv type
     */
    void writeIntervals(ContigIntervals& ci, ContigMap& cmap, int32_t svt);

    /** merge a list of BCF files into one
     * @param ifiles list of BCF files
     */
    void mergeBCFs(std::vector<std::string>& ifiles);

    /** merge one type of SV from all bcf files to output file(dup excluded)
     * @param svt SV type
     */
    void mergeOneType(int32_t svt);

    /** merge all types of SV from all bcf files to outut file(dup excluded) */
    void mergeAllType();

};

#endif
