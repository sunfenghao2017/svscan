#ifndef BEDTOOL_H
#define BEDTOOL_H

#include <string>
#include <zlib.h>
#include <bwa/kseq.h>
#include "cgranges.h"

typedef struct {
    int32_t l;
    char *s;
}bed_rest1_t;

typedef struct {
    int64_t n, m;
    bed_rest1_t *a;
}bed_rest_t;

class BedRegs{
    public:
        cgranges_t* mCR; ///< pointer to cgranges_t

        /** constructor of BedRegs */
        BedRegs();

        /** destructor of BedRegs */
        ~BedRegs();

    public:
        void loadBed(const std::string& bedFile);
        bool overlap(const std::string& chr, const int32_t& beg, const int32_t& end);

    private:
        static char *parse_bed3b(char *s, int32_t *st_, int32_t *en_, char **r);
        static char *parse_bed3(char *s, int32_t *st_, int32_t *en_);
        static cgranges_t *read_bed3b(const char *fn, bed_rest_t *r);
        static cgranges_t *read_bed3(const char *fn);

};
#endif
