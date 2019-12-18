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
    int32_t svid;        ///< sv id
    bool read1;          ///< read1 if true
    int32_t svrt;        ///< sv read type
    
    /** constructor */
    BamRec(){
        sa = '-';
        mcigar = "-";
        barcode = "-";
        lseq = "-";
        tseq = "-";
    }

    /** destructor */
    ~BamRec(){}

    /** convert BamRec to string */
    inline std::string toStr(){
        std::ostringstream oss;
        oss << svid << "\t";
        oss << chr << "\t" << pos << "\t" << strand << "\t";
        oss << mchr << "\t" << mpos << "\t" << mstrand << "\t";
        oss << cigar << "\t" << mcigar << "\t" << sa << "\t" << lseq << "\t" << tseq << "\t";
        oss << barcode << "\t" << qname << "\t" << std::boolalpha << read1 << "\t" << svrt << "\n";
        return oss.str();
    }

    /** get header items of str rec */
    static std::string getHeader(){
        return "svid\tchr\tpos\tstrand\tmchr\tmpos\tmstrand\tcigar\tmcigar\tsa\tlseq\ttseq\tbarcode\tqname\tread1\tsvrt\n";
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
        std::string svbam;          ///< sv supporting bam
        std::string fstsv;          ///< primary fusion result tsv
        std::string sstsv;          ///< secondary fusion result tsv
        std::vector<int32_t> usrid; ///< fusion id list
        std::string bamtb;          ///< bam output table
        int32_t svidf = 29;         ///< svid column index in tsv

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
    static int lines2sheet(lxw_worksheet* sheet, const std::string& buf, lxw_format* fmt = NULL);
};
#endif
