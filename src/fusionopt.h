#ifndef FUSION_OPT_H
#define FUSION_OPT_H

#include <map>
#include <set>
#include <vector>
#include <cstdint>
#include <string>
#include <htslib/vcf.h>
#include "svinfo.h"
#include "svutil.h"
#include "util.h"

/** class to store detect range of one gene */
struct DetectRange{
    std::string mGene;           /// < gene name
    std::set<int32_t> mExonList; /// < exon list to detect
    std::set<int32_t> mSVT;      /// < sv type to detect
};

/** type to store predefined genes to detect sv events in same gene */
typedef std::map<std::string, DetectRange> DetectRangeMaps;

/** type to store each kind of structural variants */
typedef std::vector<std::vector<SVInfo>> SVList;

/** type to stroe fusion gene partner information hgene ~ tgene list */
typedef std::map<std::string, std::set<std::string>> FusePairs;

/** filter options for fusion event */
struct FilterOptions{
    int32_t mMinIntraGeneSVSize = 2000; ///< min intra-gene SV size for a valid fusion
    int32_t mMinDepth = 3;              ///< min depth covering breakpoint needed
    int32_t mMinSRSupport = 3;          ///< min SR support for a valid fusion
    int32_t mMinDPSupport = 3;          ///< min DP support for a valid fusion
    int32_t mMinSUSupport = 3;          ///< min SU support for a valid fusion
    int32_t mMinSupport = 3;            ///< min supporting reads needed for a valid fusion
    float mMinVAF = 0;                  ///< min VAF needed for a valid fusion

    /** FilterOptions constructor */
    FilterOptions(){};

    /** FilterOptions destructor */
    ~FilterOptions(){};
};

/** class to store fusion report options */
struct FusionOptions{
    FilterOptions mWhiteFilter;        ///< filter options for fusion events in whitelist
    FilterOptions mUsualFilter;        ///< filter options for fusion event not in whitelist
    uint32_t mFsMaskInclude;           ///< result must match this fusion mask;
    uint32_t mFsMaskExclude;           ///< result must not match this fusion mask;
    int32_t mMaxBpOffset = 10;         ///< max breakpoint offset of an SV against background SV to be excluded
    std::string mBgBCF;                ///< background BCF file
    std::string mWhiteList;            ///< fusion event which will keep always if found
    std::string mBlackList;            ///< fusion event which will drop always if found
    std::string mSameGeneSVList;       ///< fusion event in same gene to be reported
    std::string mInfile;               ///< input file of sver sv tsv format result file
    std::string mOutFile = "fs.tsv";   ///< output file of reported fusion
    std::string mSupFile = "ss.tsv";   ///< output file of supplementary fusions
    std::string mSVModFile;            ///< output sv tsv file with fsmask updated
    SVList mBgSVs;                     ///< to store background SVs
    FusePairs mWhiteFusions;           ///< to store fusion events in fusion whitelist
    FusePairs mBlackFusions;           ///< to store fusion events in fusion blacklist
    std::set<std::string> m5Partners;  ///< genes commonly acts as 5' partner
    std::set<std::string> m3Partners;  ///< genes commonly acts as 3' partner
    std::set<std::string> mWhiteGenes; ///< to store hot gene in whitelist
    std::set<std::string> mBlackGenes; ///< to store black gene which should excluded by all fusion events
    DetectRangeMaps mDetectRngMap;     ///< detect range of sv event in the same gene
    bool mInitialized = false;         ///< FusionOptions is initialized if true

    /** FusionOptions constructor */
    FusionOptions(){
        mWhiteFilter.mMinVAF = 0;
        mWhiteFilter.mMinSupport = 3;
        mWhiteFilter.mMinDepth = 3;
        mUsualFilter.mMinVAF = 0.01;
        mUsualFilter.mMinSupport = 5;
        mUsualFilter.mMinDepth = 5;
    }
        
    /** FusionOptions destructor */
    ~FusionOptions(){};

    /** initialize filter options */
    void init();

    /** test whether an fusion event match hot gene common fusion direction
     * @param hgene head gene
     * @param tgene tail gene
     * @return true if fusion event match hot gene common fusion direction
     */
    bool matchHotDirec(const std::string& hgene, const std::string& tgene);

    /** test whether an fusion event contains hot gene partner
     * @param hgene head gene
     * @param tgene tail gene
     * @return true if fusion event contains hot gene partner
     */
    bool hasWhiteGene(const std::string& hgene, const std::string& tgene);

    /** test whether an fusion event contains black gene partner
     * @param hgene head gene
     * @param tgene tail gene
     * @return true if fusion event contains black gene partner
     */
    bool hasBlackGene(const std::string& hgene, const std::string& tgene);

    /** test whether an fusion event is in whitelist
     * @param hgene head gene
     * @param tgene tail gene
     * @return true if fusion is in whitelist
     */
    bool inWhiteList(const std::string& hgene, const std::string& tgene);

    /** test whether an fusion event is in whitelist
     * @param hgene head gene
     * @param tgene tail gene
     * @return true if fusion is in whitelist
     */
    bool inBlackList(const std::string& hgene, const std::string& tgene);
   
    /** test whether an fusion event in same gene is in same sv detection range
     * @param gene gene name
     * @param exon exon list
     * @param svt sv type
     * @return true if an fusion event is in same sv detection range
     */ 
    bool inSameSVRngMap(const std::string& gene, const std::vector<int32_t>& exon, int32_t svt);

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

    /** get fusion events in fusion whitelist */
    void parseWhiteList();

    /** get fusion events in fusion blacklist */
    void parseBlackList();

    /** get fusion events to report in same event */
    void parseSameGeneEventList();
};

#endif
