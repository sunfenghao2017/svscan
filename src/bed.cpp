#include "bed.h"

BedRegs::BedRegs(){
    mCR = NULL;
}

BedRegs::~BedRegs(){
    if(mCR){
        cr_destroy(mCR);
        mCR = NULL;
    }
}

char* BedRegs::parse_bed3b(char *s, int32_t *st_, int32_t *en_, char **r){
    char *p, *q, *ctg = 0;
    int32_t i, st = -1, en = -1;
    if (r) *r = 0;
    for (i = 0, p = q = s;; ++q) {
        if (*q == '\t' || *q == '\0') {
            int c = *q;
            *q = 0;
            if (i == 0) ctg = p;
            else if (i == 1) st = atol(p);
            else if (i == 2) {
                en = atol(p);
                if (r && c != 0) *r = q, *q = c;
            }
            ++i, p = q + 1;
            if (i == 3 || c == '\0') break;
        }
    }
    *st_ = st, *en_ = en;
    return i >= 3? ctg : 0;
}

char* BedRegs::parse_bed3(char *s, int32_t *st_, int32_t *en_){
    return parse_bed3b(s, st_, en_, 0);
}

cgranges_t* BedRegs::read_bed3b(const char *fn, bed_rest_t *r){
    gzFile fp;
    cgranges_t *cr;
    kstream_t *ks;
    kstring_t str = {0,0,0};
    int64_t k = 0;
    fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
    if (fp == 0) {
        fprintf(stderr, "ERROR: failed to open the input file\n");
        return 0;
    }
    ks = ks_init(fp);
    cr = cr_init();
    if (r) r->m = r->n = 0, r->a = 0;
    while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
        char *ctg, *rest;
        int32_t st, en;
        ctg = parse_bed3b(str.s, &st, &en, &rest);
        if (ctg) {
            cr_add(cr, ctg, st, en, k);
            if (r) {
                bed_rest1_t *p;
                if (r->n == r->m) {
                    int64_t old_m = r->m;
                    r->m = r->m? r->m<<1 : 16;
                    r->a = (bed_rest1_t*)realloc(r->a, r->m * sizeof(bed_rest1_t));
                    memset(&r->a[old_m], 0, (r->m - old_m) * sizeof(bed_rest1_t));
                }
                p = &r->a[r->n++];
                p->l = rest? str.l - (rest - str.s) : 0;
                if (rest) {
                    p->s = (char*)malloc(p->l + 1);
                    memcpy(p->s, rest, p->l);
                    p->s[p->l] = 0;
                }
            }
            ++k;
        }
    }
    if (k > INT32_MAX)
        fprintf(stderr, "WARNING: more than %d records; some functionality may not work!\n", INT32_MAX);
    free(str.s);
    ks_destroy(ks);
    gzclose(fp);
    return cr;
}

cgranges_t* BedRegs::read_bed3(const char* fn){
    return read_bed3b(fn, 0);
}

void BedRegs::loadBed(const std::string& bedFile){
    mCR = read_bed3(bedFile.c_str());
    cr_index2(mCR, 1);
}

bool BedRegs::overlap(const std::string &chr, const int32_t &beg, const int32_t &end){
    int64_t n, *b = 0, max_b = 0;
    n = cr_overlap(mCR, chr.c_str(), beg, end, &b, &max_b);
    return n > 0;
}
