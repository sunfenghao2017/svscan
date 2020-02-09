#ifndef SVDEBUG_H
#define SVDEBUG_H

#include <map>
#include <util.h>
#include <string>
#include <sstream>
#include <bamutil.h>
#include <ThreadPool.h>
#include <htslib/tbx.h>
#include <htslib/sam.h>


/** class to represent gene region on chrosome */
struct GeneRegion{
    std::string beg;    ///< beg pos
    std::string end;    ///< end pos
    std::string chr;    ///< chr name
    std::string strand; ///< strand

    /** operator to output an GeneRegion to ostream
     * @param os reference of ostream
     * @param reg reference of GeneRegion
     * @return ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const GeneRegion& reg){
        os << reg.chr << ":" << reg.beg << "-" << reg.end;
        return os;
    }

    /** get string representation of GeneRegion
     * @return string representation of GeneRegion
     */
    std::string toString(){
        std::stringstream ss;
        ss << chr << ":" << beg << "-" << end;
        return ss.str();
    }
};

/** class to store one fusion gene evidences */
struct FusionDetail{
    std::string hgene;
    std::string tgene;
    int32_t srcnt = 0;
    int32_t dpcnt = 0;
    int32_t svid = -1;

    inline friend std::ostream& operator <<(std::ostream& os, const FusionDetail& fd){
        os << fd.svid << "\t" << fd.hgene << "\t" << fd.tgene << "\t" << fd.srcnt << "\t" << fd.dpcnt << "\n";
        return os;
    }

    static void outheader(std::ostream& os){
        os << "svid\tgene1\tgene2\tsrcnt\tdpcnt\n";
    }
};

/** class to do sv calling debug */
struct SVDebug{
    std::string inbam;                         ///< input bam to do debug work on
    std::string annodb;                        ///< annotation database generated by svtools rnadb
    std::string gene2trans;                    ///< gene to transcript table generated by svtools rnadb
    std::vector<std::string> hgl;              ///< hgl names
    std::vector<std::string> tgl;              ///< tgl names
    std::map<std::string, std::string> g2tmap; ///< gene to transcript map
    std::string fglist;                        ///< hgene->tgene list
    std::string svbam;                         ///< output supporting bam file
    samFile* svbfp;                            ///< file handler of out bam
    std::string table;                         ///< output tatistical information of all evidences
    ThreadPool::ThreadPool* tp;                ///< thread pool
    size_t nthread = 8;                        ///< thread number
    std::mutex logmtx, wmtx;                   ///< log mutex, writting mutex
    bool rnamode;                              ///< input bam is from rnaseq if true

    /** SVDebug constructor */
    SVDebug(){
        svbam = "sv.bam";
        table = "stat.tsv";
        tp = NULL;
        svbfp = NULL;
        rnamode = false;
    }

    /** SVDebug destructor */
    ~SVDebug(){
        if(tp){
            delete tp;
            tp = NULL;
        }
        if(svbfp){
            sam_close(svbfp);
        }
    }

    /** init */
    void init(){
        // parse h/t pairs
        if(fglist.size()){
            std::ifstream fr(fglist);
            std::string line;
            std::vector<std::string> vstr;
            while(std::getline(fr, line)){
                util::split(line, vstr, "\t");
                if(vstr[0] == vstr[1]){
                    util::loginfo("Skip scan same h/t gene:" + vstr[0]);
                    continue;
                }
                hgl.push_back(vstr[0]);
                tgl.push_back(vstr[1]);
            }
            fr.close();
        }
        // construct thread pool
        tp = new ThreadPool::ThreadPool(std::min(hgl.size(), nthread));
        // make g2t
        util::makeMapPairFromFileByLine(gene2trans, g2tmap);
        // open out bam for write
        samFile* ifp = sam_open(inbam.c_str(), "r");
        bam_hdr_t* hdr = sam_hdr_read(ifp);
        svbfp = sam_open(svbam.c_str(), "wb");
        assert(sam_hdr_write(svbfp, hdr) >= 0);
        sam_close(ifp);
        sam_hdr_destroy(hdr);
    }

    /** debug dna sv */
    void debugDNA();
    
    /** debug rna sv */
    void debugRNA();

    /** debug one pair of fusion gene in dna */
    void debugOnePairDNA(FusionDetail& ft);

    /** debug one pair of fusion gene in rna */
    void debugOnePairRNA(FusionDetail& ft);

    /** output FusionDetails to ofstream */
    void outStat(const std::vector<FusionDetail>& fdtv);

    /** get gene region on genome 
     * @param trs gene transcript name
     * @param tbx tabix of annotabion db
     * @param reg region object to store result
     */
    static void getReg(const std::string& trs, tbx_t* tbx, htsFile* tfp, GeneRegion& reg);
};
#endif
