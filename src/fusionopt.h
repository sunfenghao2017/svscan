#ifndef FUSION_OPT_H
#define FUSION_OPT_H

#include <map>
#include <set>
#include <vector>
#include <cstdint>
#include <string>
#include <htslib/vcf.h>
#include "extanno.h"
#include "svinfo.h"
#include "svutil.h"
#include "util.h"

struct GeneRange{
    std::string mGene;
    std::string mChr;
    int32_t mStart;
    int32_t mEnd;

    bool operator<(const GeneRange& other) const {
        return (mChr < other.mChr) || 
               (mChr == other.mChr && mStart < other.mStart) ||
               (mChr == other.mChr && mStart == other.mStart && mEnd < other.mEnd);
    }
    
};

typedef std::vector<GeneRange> GeneRangeVector;

/** class to store a fusion gene fusion range */
struct FusionRange{
    std::string mHgene;                               ///< hgene
    std::string mTgene;                               ///< tgene
    std::vector<std::pair<int32_t, int32_t>> mExPair; ///< hgene exon and tgene exon pairs
    std::vector<std::pair<int32_t, int32_t>> mInPair; ///< hgene intron and tgene intron pairs
    std::vector<std::pair<int32_t, int32_t>> mExInts; ///< hgene exon and tgene intron pairs
    std::vector<std::pair<int32_t, int32_t>> mIntExs; ///< hgene intron and tgene exon pairs

    static void showPairs(std::ostream& os, const std::vector<std::pair<int32_t, int32_t>>& p){
        for(uint32_t i = 0; i < p.size(); ++i){
            os << "(" << p[i].first << ", " << p[i].second << ")\n";
        }
    }

    inline friend std::ostream& operator<<(std::ostream& os, const FusionRange& fr){
        os << "Hgene: " << fr.mHgene << "\n";
        os << "Tgene: " << fr.mTgene << "\n";
        os << "ExonPair: " << "\n";
        showPairs(os, fr.mExPair);
        os << "IntronPair: " << "\n";
        showPairs(os, fr.mInPair);
        os << "ExInPair: " << "\n";
        showPairs(os, fr.mExInts);
        os << "InExPair: " << "\n";
        showPairs(os, fr.mIntExs);
        return os;
        return os;
    }
    /** add pair */
    inline void addp(int32_t h, int32_t t, std::string uu){
        if(uu == "ee"){
            mExPair.push_back({h, t});
        }else if(uu == "ei"){
            mExInts.push_back({h, t});
        }else if(uu == "ii"){
            mInPair.push_back({h, t});
        }else if(uu == "ie"){
            mIntExs.push_back({h, t});
        }
    }
    
    /** return true if exon exon fusion got */
    inline bool eegot(int32_t h, int32_t t){
        return uugot(mExPair, h, t);
    }
    
    /** return true if exon intron fusion got */
    inline bool eigot(int32_t h, int32_t t){
        return uugot(mExInts, h, t);
    }

    /** return true if intron intron fusion got */
    inline bool iigot(int32_t h, int32_t t){
        return uugot(mInPair, h, t);
    }

    /** return true if intron exon fusion got */
    inline bool iegot(int32_t h, int32_t t){
        return uugot(mIntExs, h, t);
    } 

    /** magic func */
    inline bool uugot(const std::vector<std::pair<int32_t, int32_t>>& plist, int32_t h, int32_t t){
        for(uint32_t i = 0; i < plist.size(); ++i){
            if(plist[i].first == h && plist[i].second == t) return true;
        }
        return false;
    }
};

/** type to store a series of fusion pair ranges */
typedef std::map<std::string, FusionRange> FusionReportRange;

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
    int32_t mMinSRSeed = 3;             ///< min SR seed for a valid fusion
    int32_t mMinDPSeed = 3;             ///< min DP seed for a valid fusion
    int32_t mMinSUSeed = 3;             ///< min SU seed for a valid fusion
    int32_t mMaxRepHit = 4;             ///< max hit of two part of fusion seq
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
    ExtraAnno mExtraAnnotator;         ///< extra gene annotator 
    uint32_t mFsMaskInclude;           ///< result must match this fusion mask;
    uint32_t mFsMaskExclude;           ///< result must not match this fusion mask;
    int32_t mMaxBpOffset = 10;         ///< max breakpoint offset of an SV against background SV to be excluded
    std::string mRef;                  ///< reference file used in alignment of bam
    std::string mBgBCF;                ///< background BCF file
    std::string mWhiteList;            ///< fusion event which will keep always if found
    std::string mBlackList;            ///< fusion event which will drop always if found
    std::string mFsRptList;            ///< if provided, report fusions in this list only
    std::string mGeneCrdList;          ///< gene coordinate range list
    std::string mSameGeneSVList;       ///< fusion event in same gene to be reported
    std::string mExtraAnnoList;        ///< fusion gene which should be annotated seperately
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
    FusionReportRange mFusionRptMap;   ///< report fusion event only in this range
    GeneRangeVector mGeneRangeVec;     ///< gene coordinates vector
    bool mInitialized = false;         ///< FusionOptions is initialized if true

    /** FusionOptions constructor */
    FusionOptions(){
        mWhiteFilter.mMinVAF = 0;
        mWhiteFilter.mMinSupport = 3;
        mWhiteFilter.mMinDepth = 30;
        mWhiteFilter.mMinSRSeed = 3;
        mWhiteFilter.mMinDPSeed = 3;
        mWhiteFilter.mMaxRepHit = 4;
        mUsualFilter.mMinVAF = 0.01;
        mUsualFilter.mMinSupport = 5;
        mUsualFilter.mMinDepth = 50;
        mUsualFilter.mMinSRSeed = 5;
        mUsualFilter.mMinDPSeed = 5;
        mUsualFilter.mMaxRepHit = 4;
    }
        
    /** FusionOptions destructor */
    ~FusionOptions(){};

    /** initialize filter options */
    void init();

    /** initialize gene coordinate range list */
    void initGeneCrdRange();

    /** test whether two gene is near each other
     * @param g1 gene1 name
     * @param chr1 gene1 chr
     * @param pos1 one pos of gene1
     * @param g2 geng2 name
     * @param chr2 gene2 chr
     * @return positive distance of two near gene, -1 if two gene not near
     */
    int32_t geneNear(const std::string& g1, const std::string& chr1, int32_t pos1, const std::string& g2, const std::string& chr2);

    /** initialize fusion report range */
    void initFusionRptRange();

    /** test whether an fusion event is in report range
     * @param hg hgene name
     * @param tg tgene name
     * @param hu hgene unit number
     * @param tu tgene unit number
     * @param uu unit type 
     */
    bool inFsRptRange(std::string hgene, std::string tgene,  int32_t hu, int32_t tu, std::string uu = "ee");

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
