#ifndef BAMTOETB_H
#define BAMTOETB_H

#include <map>
#include <string>
#include <fstream>
#include <htslib/sam.h>
#include <xlsxwriter.h>
#include "bamutil.h"
#include "util.h"

/** bam record items to output */
struct BamRec{
    std::string chr;     ///< chr
    int32_t pos;         ///< pos
    std::string mchr;    ///< mate chr
    int32_t mpos;        ///< mate pos
    std::string cigar;   ///< cigar
    std::string mcigar;  ///< mate cigar
    std::string sa;      ///< sa
    std::string seq;     ///< seq
    std::string barcode; ///< barcode
    std::string qname;   ///< read name
    std::string sname;   ///< sample name
    int32_t svid;        ///< sv id
    
    /** constructor */
    BamRec(){}

    /** destructor */
    ~BamRec(){}

    /** convert BamRec to string */
    inline std::string toStr(){
        std::ostringstream oss;
        oss << chr << "\t" << pos << "\t" << mchr << "\t" << mpos << "\t";
        oss << cigar << "\t" << mcigar << "\t" << sa << "\t" << seq << "\t";
        oss << barcode << "\t" << qname;
        return oss.str();
    }

    /** get header items of str rec */
    static std::string getHeader(){
        return "chr\tpos\tmchr\tmpos\tcigar\tmcigar\tsa\tseq\tbarcode\tqname";
    }

    /** compare two BamRec */
    inline bool operator<(const BamRec& other) const {
        return (sname < other.sname) ||
               (sname == other.sname && svid < other.svid) ||
               (sname == other.sname && svid == other.svid && chr < other.chr) ||
               (sname == other.sname && svid == other.svid && chr == other.chr && sa.length() > other.sa.length()) || 
               (sname == other.sname && svid == other.svid && chr == other.chr && sa.length() == other.sa.length() && pos < other.pos);
    }
};

/** type to store a list of BamRec */
typedef std::vector<BamRec> BamRecVector;

/** class to extract sv supporting bam to table */
class BamToTable{
    public:
        std::string svbam;  ///< sv supporting bam
        std::string fstsv;  ///< primary fusion result tsv
        std::string sstsv;  ///< secondary fusion result tsv
        std::string bamtb;  ///< bam output table
        int32_t svidf = 28; ///< svid column index in tsv

    /** BamToTable constructor */
    BamToTable(){
        bamtb = "b2t.xlsx";
    }

    /** BamToTable destructor */
    ~BamToTable(){}

    /** extract sv supporting sv bam records into excel */
    void b2t();

    /** convert an bam record to BamRec*/
    void b2r(bam1_t* b, bam_hdr_t* h, BamRec& br, int32_t id);

    /** get svids from fs/ss tsvs */
    void getsvid(std::set<int32_t>& svids);

    /** output line buffers to a sheet */
    static void lines2sheet(lxw_worksheet* sheet, const std::string& buf, lxw_format* fmt = NULL);
};
#endif