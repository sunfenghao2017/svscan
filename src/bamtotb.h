#ifndef BAMTOETB_H
#define BAMTOETB_H

#include <map>
#include <regex>
#include <string>
#include <fstream>
#include <htslib/sam.h>
#include <xlsxwriter.h>
#include <bamutil.h>
#include <lxwutil.h>
#include <util.h>
#include "bamsorter.h"

/** hit got */
struct HitPat{
    bool isread1;
    int32_t index;
};

/** fusion information corresponding to a svid */
struct FRExtraInfo{
    int32_t svid;
    std::string fsgene;
};

/** bam record items to output */
struct BamRec{
    std::string chr;     ///< chr
    int32_t pos;         ///< pos
    char strand;         ///< strand
    std::string mchr;    ///< mate chr
    int32_t mpos;        ///< mate pos
    char mstrand;        ///< mate strand
    std::string cigar;   ///< cigar
    std::string mcigar;  ///< mate cigar
    std::string sa;      ///< sa
    std::string seq;     ///< seq
    std::string lseq;    ///< leading seq
    std::string tseq;    ///< tailing seq
    std::string barcode; ///< barcode
    std::string qname;   ///< read name
    std::string sname;   ///< sample name
    int32_t rbp;         ///< read breakpoint position
    int32_t sbp;         ///< supplementary read breakpoint position
    int32_t svid;        ///< sv id
    bool read1;          ///< read1 if true
    int32_t svrt;        ///< sv read type
    std::string fsgene;  ///< fusion gene
    int32_t lhit;        ///< lseq hit times on genome
    int32_t thit;        ///< tseq hit times on genome
    
    /** constructor */
    BamRec(){
        sa = '-';
        mcigar = "-";
        barcode = "-";
        lseq = "";
        tseq = "";
        fsgene = "-";
        lhit = 0;
        thit = 0;
    }

    /** destructor */
    ~BamRec(){}

    /** convert BamRec to string */
    inline std::string toStr(){
        if(lseq.empty()) lseq = "-";
        if(tseq.empty()) tseq = "-";
        std::ostringstream oss;
        oss << svid << "\t";
        oss << fsgene << "\t";
        oss << chr  << "," << pos  << "," << strand  << "," << cigar  << "\t";
        oss << mchr << "," << mpos << "," << mstrand << "," << mcigar << "\t"; 
        oss << sa << "\t" << rbp << "\t" << sbp << "\t" << lhit << "\t" << thit << "\t" << lseq << "\t" << tseq << "\t";
        // aseq output beg
        if(svrt) oss << "-";
        else oss << lseq << tseq;
        oss << "\t";
        // aseq output end
        oss << barcode << "\t" << qname << "\t" << std::boolalpha << read1 << "\t" << svrt << "\n";
        return oss.str();
    }

    /** get header items of str rec */
    static std::string getHeader(){
        return "svid\tfsgene\trmap\tmmap\tsa\trbp\tsbp\tlhit\tthit\tlseq\ttseq\taseq\tbarcode\tqname\tread1\tsvrt\n";
    }

    /** compare two BamRec */
    inline bool operator<(const BamRec& other) const {
        return (sname < other.sname) ||
               (sname == other.sname && svid < other.svid) ||
               (sname == other.sname && svid == other.svid && chr < other.chr) ||
               (sname == other.sname && svid == other.svid && chr == other.chr && sa.length() > other.sa.length()) || 
               (sname == other.sname && svid == other.svid && chr == other.chr && sa.length() == other.sa.length() && pos < other.pos) || 
               (sname == other.sname && svid == other.svid && chr == other.chr && sa.length() == other.sa.length() && pos == other.pos && mchr < other.mchr) ||
               (sname == other.sname && svid == other.svid && chr == other.chr && sa.length() == other.sa.length() && pos == other.pos && mchr == other.mchr && 
                mpos < other.mpos);
    }
};

/** type to store a list of BamRec */
typedef std::vector<BamRec> BamRecVector;

/** class to extract sv supporting bam to table */
class BamToTable{
    public:
        std::string svbam;                 ///< sv supporting bam
        std::string newbam;                ///< newly updated sv supporting bam
        std::string fstsv;                 ///< fusion result tsv
        std::vector<int32_t> usrid;        ///< fusion id list
        std::string bamtb;                 ///< bam output table(excel format)
        std::string bamtt;                 ///< bam output txt(tsv format)
        std::vector<std::string> fsgene;   ///< fusion gene of each svid
        std::map<std::string, HitPat> peout; ///< reads output pe support records
        int32_t svidf = 33;                ///< svid column index in tsv
        int32_t fsidf = 0;                 ///< fusion column index in tsv
        bool refinedp = false;             ///< refine dp if true

    /** BamToTable constructor */
    BamToTable(){
        newbam = "fs.bam";
        bamtb = "bt.xlsx";
        bamtt = "bt.tsv";
        refinedp = false;
    }

    /** BamToTable destructor */
    ~BamToTable(){
    }

    /** extract sv supporting sv bam records into excel */
    void b2t();

    /** convert an bam record to BamRec*/
    void b2r(const bam1_t* b, bam_hdr_t* h, BamRec& br, int32_t id);

    /** get sv info from fs tsvs */
    void getFRExtraInfo(std::map<int32_t, FRExtraInfo>& fim);
};

#endif
