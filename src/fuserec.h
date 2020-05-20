#ifndef FUSEREC_H
#define FUSEREC_H

#include <string>
#include <cstdint>
#include <iostream>
#include "svutil.h"
#include "fusionopt.h"

/** web URI prefix */
const std::string URI_PRE = "http://10.100.2.10:9999/igv/?";

/** class to store fusion record schema */
struct FusionRecSchema{
    int fusegene = 0;
    int fusionmols = 1;
    int totalmols = 2;
    int fuserate = 3;
    int indb = 4;
    int fusepattern = 5;
    int status = 6;
    int distance = 7;
    int gene1 = 8;
    int chr1 = 9;
    int jctpos1 = 10;
    int strand1 = 11;
    int transcript1 = 12;
    int gene2 = 13;
    int chr2 = 14;
    int jctpos2 = 15;
    int strand2 = 16;
    int transcript2 = 17;
    int fusionsequence = 18;
    int fseqbp = 19;
    int svt = 20;
    int svsize = 21;
    int srcount = 22;
    int dpcount = 23;
    int srrescued = 24;
    int dprescued = 25;
    int srrefcount = 26;
    int dprefcount = 27;
    int srsrescued = 28;
    int srsmalncnt = 29;
    int srsmrate = 30;
    int insbp = 31;
    int insseq = 32;
    int svid = 33;
    int svint = 34;
    int fsmask = 35;
    int fsHits = 36;
    int exon1 = 37;
    int exon2 = 38;
    int ts1name = 37;
    int ts1pos = 38;
    int ts2name = 39;
    int ts2pos = 40;
    int cigar = 41;
    int urld = 39;
    int urlr = 42;
};

/** class to store an fusion record */
struct FusionRecord{
    std::string fusegene;       ///< fusion gene hgene->tgene
    std::string fusepattern;    ///< fusion gene strand e.g., +->-, -->-, +->+
    int32_t fusionmols;         ///< # of molecules supporting fusion event
    int32_t totalmols;          ///< # of molecules supporting fusion event and refrence type
    float fuserate;             ///< rate of fusionmols in totalmols
    std::string status;         ///< status str
    std::string gene1;          ///< hgene
    std::string chr1;           ///< chr on which hgene situated
    int32_t jctpos1;            ///< junction position of hgene on chr1
    std::string strand1;        ///< strand hgene situated on chr1
    std::string transcript1;    ///< transcript of hgene
    std::string gene2;          ///< tgene
    std::string chr2;           ///< chr on which tgene situated
    int32_t jctpos2;            ///< junction position of tgene on chr2
    std::string strand2;        ///< strand tgene situated on chr2
    std::string transcript2;    ///< transcript of tgene
    std::string fusionsequence; ///< fusion concensus sequence
    int32_t fseqbp;             ///< fusion sequence length of one partner in fusiongene
    std::string indb;           ///< this fusion event is in public database if true
    std::string svt;            ///< svtype, raw str representation
    int32_t svsize;             ///< sv size of this fusion event
    int32_t srcount;            ///< seed split reads supporting this fusion event
    int32_t dpcount;            ///< seed discordant pair of reads supporting this fusion event
    int32_t srrescued;          ///< all split reads supporting this fusion event
    int32_t dprescued;          ///< all siscordant pair of reads supporting this fusion event
    int32_t srrefcount;         ///< all reads supporting ref type
    int32_t dprefcount;         ///< all paired reads supporting ref type
    int32_t srsrescued;         ///< number of sr seed rescued
    int32_t srsmalncnt;         ///< number of sr rescued seed with partner multiple realn
    double srsmrate;            ///< rate of sr seeds with partner multiple realn
    int32_t insbp;              ///< inserted sequence length around breakpoint
    std::string insseq;         ///< inserted sequence
    int32_t svid;               ///< sv id
    int32_t svint;              ///< svtype, precise integer representation
    TFUSION_FLAG fsmask;        ///< fusion mask
    int32_t exon1;              ///< approximate gene1 exon(DNA only)
    int32_t exon2;              ///< approximate gene2 exon(DNA only)
    int32_t ie1;                ///< intron/exon count 
    int32_t ie2;                ///< intron/exon count
    std::string ts1name;        ///< transcript1 name(RNA only)
    int32_t ts1pos;             ///< transcript1 breakpoint(RNA only)
    std::string ts2name;        ///< transcript2 name(RNA only)
    int32_t ts2pos;             ///< transcript2 breakpoint(RNA only)
    std::string cigar;          ///< cigar string of bp(RNA only)
    int32_t distance;           ///< distance of two gene
    int32_t fsHits;             ///< fusion seq hits int mask
    bool report;                ///< this fusion will be reported if true
    std::string url;            ///< fusion event online igv address

    /** construct an FusionRecord */
    FusionRecord(){
        fsmask = 0;
        url = "-";
        report = false;
    }

    /** destroy an FusionRecord */
    ~FusionRecord(){}

    /** update url */
    void get_url(const std::string& bam){
        std::stringstream ss;
        ss << URI_PRE;
        ss << "b=" << bam << "&r=" << chr1 << ":" << jctpos1 + 1 << "-" << jctpos1 + 1 << ",";
        ss << chr2 << ":" << jctpos2 + 1 << "-" << jctpos2 + 2;
        ss << "&fsid=" << svid;
        url = ss.str();
    }

    /** operator to compare two FusionRecord */
    inline bool operator<(const FusionRecord& other) const {
        return (svint < other.svint) ||
               (svint == other.svint && chr1 < other.chr1) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 < other.chr2) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 < other.gene1) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 < other.gene2) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 == other.gene2 && exon1 < other.exon1) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 == other.gene2 && exon1 == other.exon1 && exon2 < other.exon2) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 == other.gene2 && exon1 == other.exon1 && exon2 == other.exon2 && srrescued < other.srrescued) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 == other.gene2 && exon1 == other.exon1 && exon2 == other.exon2 && srrescued == other.srrescued && dprescued < other.dprescued) ||
               (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 == other.gene2 && exon1 == other.exon1 && exon2 == other.exon2 && srrescued == other.srrescued && dprescued == other.dprescued && fuserate < other.fuserate);

    }

    /** test whether two FusionRecord is the same fusion */
    inline bool samefs(const FusionRecord& other){
        return (svint == other.svint && chr1 == other.chr1 && chr2 == other.chr2 && gene1 == other.gene1 && gene2 == other.gene2 && exon1 == other.exon1 && exon2 == other.exon2);
    } 
    /** output an FusionRecord to ostream
     * @param os reference of ostream
     * @param fsr reference of FusionRecord
     * @return reference of ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const FusionRecord& fsr){
        os << fsr.fusegene << "\t" << fsr.fusionmols << "\t" << fsr.totalmols << "\t" << fsr.fuserate << "\t";
        os << fsr.indb << "\t" << fsr.fusepattern << "\t" << fsr.getStatus() << "\t" << fsr.distance << "\t";
        os << fsr.gene1 << "\t" << fsr.chr1 << "\t" << fsr.jctpos1 + 1 << "\t" << fsr.strand1 << "\t" << fsr.transcript1 << "\t";
        os << fsr.gene2 << "\t" << fsr.chr2 << "\t" << fsr.jctpos2 + 1<< "\t" << fsr.strand2 << "\t" << fsr.transcript2 << "\t";
        os << fsr.fusionsequence << "\t" << fsr.fseqbp + 1 << "\t" << fsr.svt << "\t" << fsr.svsize << "\t";
        os << fsr.srcount << "\t" << fsr.dpcount << "\t" << fsr.srrescued << "\t" << fsr.dprescued << "\t";
        os << fsr.srrefcount << "\t" << fsr.dprefcount << "\t";
        os << fsr.srsrescued << "\t" << fsr.srsmalncnt << "\t" << fsr.srsmrate << "\t";
        os << fsr.insbp << "\t" << fsr.insseq << "\t" << fsr.svid << "\t" << fsr.svint << "\t" << fsr.fsmask << "\t" << fsr.fsHits;
        if(fsr.fsmask & FUSION_FCALLFROMRNASEQ){
            os << "\t" << fsr.ts1name << "\t" << fsr.ts1pos + 1 << "\t" << fsr.ts2name << "\t" << fsr.ts2pos + 1 << "\t" << fsr.cigar;
        }else{
            os << "\t" << fsr.exon1 << "\t" << fsr.exon2;
        }
        os << fsr.url << "\n";
       return os;
    }

    /** get headline of final output
     * @param rnamode rnamode or not
     * @return headline
     */
    static std::string gethead(bool rnamode = false){
        std::string header = "FusionGene\tFusionMols\tTotalMols\tFusionRate\t"; //[0-3]
        header.append("inDB\tfsPattern\tstatus\tfpDist\t"); //[4-7]
        header.append("Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t"); //[8-12]
        header.append("Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t"); //[13-17]
        header.append("FusionSequence\tfseqBp\tsvType\tsvSize\t"); //[18-21]
        header.append("srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"); //[22-27]
        header.append("srSRescued\tsrSResMaln\tsrSResMalnRate\t"); //[28,30]
        header.append("insBp\tinsSeq\tsvID\tsvtInt\tfsMask\tfsHits"); //[31-36]
        if(rnamode) header.append("\tts1Name\tts1Pos\tts2Name\tts2Pos\tfsCigar\t");
        else header.append("\texon1\texon2\t");
        header.append("URL\n");
        return header;
    }

    /** get status of fusion event */
    std::string getStatus() const {
        std::string sts;
        if(fsmask & FUSION_FMINDB) sts.append("M");
        if(!(fsmask & FUSION_FNORMALCATDIRECT)) sts.append("D");
        if(!(fsmask & FUSION_FCOMMONHOTDIRECT)) sts.append("C");
        if(gene1 == gene2) sts.append("S");
        if(!(fsmask & FUSION_FALLGENE)) sts.append("G");
        if(sts.empty()) sts.append("Y");
        return sts;
    }

    /** update seed/rescue/depth/.. determined fusion mask */
    inline void maskFusion(FusionOptions* fsopt){
        fsmask &= (~(FUSION_FLOWSUPPORT | FUSION_FLOWAF | FUSION_FLOWDEPTH | FUSION_FINREPORTRNG | FUSION_FSRSEEDERROR)); //clear some mask
        if(srcount > 0 && srsrescued == srsmalncnt) fsmask |= FUSION_FSRSEEDERROR; // fusion srseed are err-prone
        if(fsmask & (FUSION_FINDB | FUSION_FMINDB)){ // fusion/mirror in public db
            if(fusionmols < fsopt->mWhiteFilter.mMinSupport){ // total molecule support
                fsmask |= FUSION_FLOWSUPPORT;
            }
            if(srcount < fsopt->mWhiteFilter.mMinSRSeed && dpcount < fsopt->mWhiteFilter.mMinDPSeed){ // seed 
                fsmask |= FUSION_FLOWSUPPORT;
            }
            if(srrescued < fsopt->mWhiteFilter.mMinSRSupport && dprescued < fsopt->mWhiteFilter.mMinDPSupport){ // rescue
                fsmask |= FUSION_FLOWSUPPORT;
            }
            if(fuserate < fsopt->mWhiteFilter.mMinVAF){ // vaf
                fsmask |= FUSION_FLOWAF;
            }
            if(totalmols < fsopt->mWhiteFilter.mMinDepth){ // depth
                fsmask |= FUSION_FLOWDEPTH;
            }
            if(svint != 4 && (fsmask & FUSION_FINSAMEGENE) && (!(fsmask & FUSION_FCALLFROMRNASEQ))){ // size
                if(svsize < fsopt->mWhiteFilter.mMinIntraGeneSVSize){
                    fsmask |= FUSION_FTOOSMALLSIZE;
                }
            }
            if(svint < 4 && (!(fsmask & FUSION_FALLGENE))){// size of gene->nongene
                if(svsize < fsopt->mWhiteFilter.mMinIntraGeneSVSize){
                    fsmask |= FUSION_FTOOSMALLSIZE;
                }
            }
        }else if(fsmask & FUSION_FHOTGENE){ // fusion only with gene in target region
            if(fsmask & FUSION_FPRECISE){
                if(fusionmols < fsopt->mUsualFilter.mMinSupport){
                    fsmask |= FUSION_FLOWSUPPORT;
                }
                if(srcount < fsopt->mUsualFilter.mMinSRSeed && dpcount < fsopt->mUsualFilter.mMinDPSeed){
                    fsmask |= FUSION_FLOWSUPPORT;
                }
                if(srrescued < fsopt->mUsualFilter.mMinSRSupport && dprescued < fsopt->mUsualFilter.mMinDPSupport){
                    fsmask |= FUSION_FLOWSUPPORT;
                }
                if(fuserate < fsopt->mUsualFilter.mMinVAF){
                    fsmask |= FUSION_FLOWAF;
                }
                if(totalmols < fsopt->mUsualFilter.mMinDepth){
                    fsmask |= FUSION_FLOWDEPTH;
                }
            }else{
                fsmask |= FUSION_FLOWSUPPORT;
            }
            if(svint != 4 && (fsmask & FUSION_FINSAMEGENE) && (!(fsmask & FUSION_FCALLFROMRNASEQ))){ // size
                if(svsize < fsopt->mUsualFilter.mMinIntraGeneSVSize){
                    fsmask |= FUSION_FTOOSMALLSIZE;
                }
            }
            if(svint < 4 && (!(fsmask & FUSION_FALLGENE))){// size of gene->nongene
                if(svsize < fsopt->mUsualFilter.mMinIntraGeneSVSize){
                    fsmask |= FUSION_FTOOSMALLSIZE;
                }
            }
        }
        if(fsopt->mFsRptList.empty()){
            fsmask |= FUSION_FINREPORTRNG;
        }else{
            std::vector<std::string> vstr1, vstr2;
            util::split(transcript1, vstr1, ",");
            util::split(transcript2, vstr2, ",");
            int32_t unum1 = std::atoi(vstr1[3].c_str());
            int32_t unum2 = std::atoi(vstr2[3].c_str());
            std::string uumk = "";
            if(vstr1[2] == "intron") uumk.append("i");
            else if(vstr1[2] == "exon") uumk.append("e");
            else uumk.append("-");
            if(vstr2[2] == "intron") uumk.append("i");
            else if(vstr2[2] == "exon") uumk.append("e");
            else uumk.append("-");
            if(fsopt->inFsRptRange(gene1, gene2, exon1, exon2, "ee")){
                fsmask |= FUSION_FINREPORTRNG;
            }else if(fsopt->inFsRptRange(gene1, gene2, unum1, unum2, uumk)){
                fsmask |= FUSION_FINREPORTRNG;
            }
        }
        if(distance < 0) fsmask |= FUSION_FFBG; // overlap gene leading fusion treat as background
        if((fsmask & FUSION_FPRECISE) && srrescued == 0) fsmask |= FUSION_FLOWSUPPORT; // rescue failed fusion drop
        if(util::startsWith(transcript1, "NR") || util::startsWith(transcript2, "NR")) fsmask |= FUSION_FWITHNCRNA; // ncrna participated
    }

    /** get inter gene nomenclature */
    inline void getInterGeneInfo(FusionOptions* fsopt){
        if(fsmask & FUSION_FALLGENE) return;
        InterGeneInfo ig1, ig2;
        ig1.fsgene = gene1; ig1.fsinfo = gene1;
        ig2.fsgene = gene2; ig2.fsinfo = gene2;
        if(gene1 == "-"){
            fsopt->getInterGeneInfo(ig1, chr1, jctpos1);
            ig1.getFuseInfo();
            ig1.getFuseGeneName();
        }
        if(gene2 == "-"){
            fsopt->getInterGeneInfo(ig2, chr2, jctpos2);
            ig2.getFuseInfo();
            ig2.getFuseGeneName();
        }
        fusegene = ig1.fsgene + "->" + ig2.fsgene;
        gene1 = ig1.fsinfo;
        gene2 = ig2.fsinfo;
    }

    /** line2rec */
    void line2rec(const std::string& line, const FusionRecSchema& s){
        std::vector<std::string> vstr;
        util::split(line, vstr, "\t");
        fusegene = vstr[s.fusegene];
        fusionmols = std::atoi(vstr[s.fusionmols].c_str());
        totalmols = std::atoi(vstr[s.totalmols].c_str());
        fuserate = std::atof(vstr[s.fuserate].c_str());
        indb = vstr[s.indb];
        fusepattern = vstr[s.fusepattern];
        status = vstr[s.status];
        distance = std::atoi(vstr[s.distance].c_str());
        gene1 = vstr[s.gene1];
        chr1 = vstr[s.chr1];
        jctpos1 = std::atoi(vstr[s.jctpos1].c_str());
        strand1 = vstr[s.strand1];
        transcript1 = vstr[s.transcript1];
        gene2 = vstr[s.gene2];
        chr2 = vstr[s.chr2];
        jctpos2 = std::atoi(vstr[s.jctpos2].c_str());
        strand2 = vstr[s.strand2];
        transcript2 = vstr[s.transcript2];
        fusionsequence = vstr[s.fusionsequence];
        fseqbp = std::atoi(vstr[s.fseqbp].c_str());
        svt = vstr[s.svt];
        svsize = std::atoi(vstr[s.svsize].c_str());
        srcount = std::atoi(vstr[s.srcount].c_str());
        dpcount = std::atoi(vstr[s.dpcount].c_str());
        srrescued = std::atoi(vstr[s.srrescued].c_str());
        dprescued = std::atoi(vstr[s.dprescued].c_str());
        srrefcount = std::atoi(vstr[s.srrefcount].c_str());
        dprefcount = std::atoi(vstr[s.dprefcount].c_str());
        srsrescued = std::atoi(vstr[s.srsrescued].c_str());
        srsmalncnt = std::atoi(vstr[s.srsmalncnt].c_str());
        srsmrate = std::atof(vstr[s.srsmrate].c_str());
        insbp = std::atoi(vstr[s.insbp].c_str());
        insseq = vstr[s.insseq];
        svid = std::atoi(vstr[s.svid].c_str());
        svint = std::atoi(vstr[s.svint].c_str());
        fsmask = std::atoi(vstr[s.fsmask].c_str());
        fsHits = std::atoi(vstr[s.fsHits].c_str());
        if(fsmask & FUSION_FCALLFROMRNASEQ){
            ts1name = vstr[s.ts1name];
            ts1pos = std::atoi(vstr[s.ts1pos].c_str());
            ts2name = vstr[s.ts2name];
            ts2pos = std::atoi(vstr[s.ts2pos].c_str());
            cigar = vstr[s.cigar];
            url = vstr[s.urlr];
        }else{
            exon1 = std::atoi(vstr[s.exon1].c_str());
            exon2 = std::atoi(vstr[s.exon2].c_str());
            url = vstr[s.urld];
        }
    }
};

/** class to store a list of FusionRecords */
typedef std::vector<FusionRecord> FusionRecordList;

/** mark mirror fusion event which came from mirror structural evnet to report the desired one */
inline void markFusionMirrorFromMirrorSVEvent(FusionRecordList& frl){
    if(frl.empty()) return;
    // construct index
    std::vector<int32_t> idxs(9, 0), idxe(9, 0);
    for(uint32_t i = 1; i < frl.size(); ++i){
        if(frl[i].svint != frl[i - 1].svint){
            idxs[frl[i].svint] = i;
            idxe[frl[i - 1].svint] = i - 1;
        }
    }
    idxe[frl[frl.size() - 1].svint] = frl.size() - 1;
    // construct mark pair
    std::map<int32_t, int32_t> mkp = {{0, 1}, {5, 6}, {7, 8}};
    // begin mark
    for(auto iter = mkp.begin(); iter != mkp.end(); ++iter){
        for(int32_t firidx = idxs[iter->first]; firidx <= idxe[iter->first]; ++firidx){
            if(frl[firidx].report && (frl[firidx].fsmask & FUSION_FWITHMIRROR) && (frl[firidx].fsmask & FUSION_FALLGENE)){
                for(int32_t secidx = idxs[iter->second]; secidx <= idxe[iter->second]; ++secidx){
                    if(frl[secidx].report && (frl[secidx].fsmask & FUSION_FWITHMIRROR) && (frl[secidx].fsmask & FUSION_FALLGENE)){
                        if(frl[firidx].gene1 == frl[secidx].gene2 && 
                           frl[firidx].gene2 == frl[secidx].gene1 &&
                           frl[firidx].exon1 == frl[secidx].exon2 &&
                           frl[firidx].exon2 == frl[secidx].exon1){
                            bool kf = false, ks = false;
                            if(frl[firidx].fsmask & FUSION_FCOMMONHOTDIRECT) kf = true;
                            if(frl[secidx].fsmask & FUSION_FCOMMONHOTDIRECT) ks = true;
                            if(kf ^ ks){
                                if(kf) frl[secidx].report = false;
                                else frl[firidx].report = false;
                            }else{
                                if(frl[firidx].fusionmols > frl[secidx].fusionmols){
                                    frl[secidx].report = false;
                                }else{
                                    frl[firidx].report = false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

/** mark mirror fusion event to keep the PRIMARY only if mirror SECONDARY exists */
inline void markMirrorFusionEvent(FusionRecordList& frl){
    std::set<std::string> fsoutset;
    for(auto& e: frl){
        if(e.report && (e.fsmask & FUSION_FPRIMARY)) fsoutset.insert(e.fusegene);
    }
    for(auto& e: frl){
        if(e.report && (e.fsmask & FUSION_FSUPPLEMENTARY)){
            std::string revfg = e.gene2 + "->" + e.gene1;
            if(fsoutset.find(e.fusegene) != fsoutset.end() ||
               fsoutset.find(revfg) != fsoutset.end()){
                e.report = false;
            }
        }
    }
}

inline std::string adjustPattern(const int& mask, const std::string& pat, const std::string& ss){
    if(pat == ss) return pat; // same pattern
    if(pat[0] == ss[1] && pat[1] == ss[0]) return ss; // swapped pattern
    bool needSwapt = false;
    if(mask & 1){//h is hot
        if(pat[0] != ss[0]){
            needSwapt = true;
        }
    }else if(mask & 2){//t is hot
        if(pat[1] != ss[1]){
            needSwapt = true;
        }
    }
    if(needSwapt){
        if(pat == "++") return "--";
        if(pat == "+-") return "-+";
        if(pat == "-+") return "+-";
        if(pat == "--") return "++";
    }
    return pat;
}

#endif
