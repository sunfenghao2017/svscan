#ifndef BAMTOETB_H
#define BAMTOETB_H

#include <map>
#include <string>
#include <fstream>
#include <htslib/sam.h>
#include <xlsxwriter.h>
#include "bamutil.h"
#include "util.h"

/** class to extract sv supporting bam to table */
class BamToTable{
    public:
        std::string svbam;  ///< sv supporting bam
        std::string fstsv;  ///< primary fusion result tsv
        std::string sstsv;  ///< secondary fusion result tsv
        std::string bamtb;  ///< bam output table
        int32_t svidf = 18; ///< svid column index in tsv

    /** BamToTable constructor */
    BamToTable(){
        bamtb = "b2t.xlsx";
    }

    /** BamToTable destructor */
    ~BamToTable(){}

    /** extract sv supporting sv bam records into excel */
    void b2t();

    /** convert an bam record to string */
    std::string b2s(bam1_t* b, bam_hdr_t* h);

    /** get svids from fs/ss tsvs */
    void getsvid(std::set<int32_t>& svids);

    /** output line buffers to a sheet */
    static void lines2sheet(lxw_worksheet* sheet, const std::string& buf, lxw_format* fmt = NULL);
};
#endif
