#include "olpbedl.h"

BedRegs::BedRegs(){
}

BedRegs::~BedRegs(){
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

cgranges_t* BedRegs::loadOneBed(const std::string& bedFile){
    cgranges_t* cgr =read_bed3(bedFile.c_str());
    cr_index2(cgr, 1);
    return cgr;
}

void BedRegs::loadBeds(){
    mRegs.resize(mBedPaths.size());
    mNames.resize(mBedPaths.size());
    std::string bname;
    for(uint32_t i = 0; i < mBedPaths.size(); ++i){
        bname = util::basename(mBedPaths[i]);
        bname = util::replace(bname, ".bed", "");
        mNames[i] = bname;
        mRegs[i] = loadOneBed(mBedPaths[i]);
    }
}

void BedRegs::olpAna(){
    std::map<std::string, OlpStore*> olpret;
    std::vector<cgranges_t*> twoins;
    std::string key, okey, lk;
    std::vector<std::string> lks(mRegs.size());
    // store each one regs into olpret
    for(uint32_t i = 0; i < mRegs.size(); ++i){
        OlpStore* osr = new OlpStore();
        osr->idx.push_back(i);
        osr->olpret = mRegs[i];
        lks[i] = std::to_string(i);
        olpret[lks[i]] = osr;
    }
    // process overlap of all combinations
    for(uint32_t i = 2; i <= mRegs.size(); ++i){
        std::vector<std::vector<int>> cmb;
        util::combSet(mRegs.size(), i, cmb);
        for(auto& e: cmb){
            std::stringstream ss;
            ss << e[0];
            for(uint32_t k = 1; k < i - 1; ++k){
                ss << "-" << e[k];
            }
            okey = ss.str();
            ss << "-" << e[i - 1];
            key = ss.str();
            lk = lks[e[i - 1]];
            OlpStore* osr = new OlpStore();
            osr->idx = e;
            osr->olpret = cr_overlap2(olpret[okey]->olpret, olpret[lk]->olpret);
            olpret[key] = osr;
            if(i == 2) twoins.push_back(osr->olpret);
        }
    }
    // output results of all overlaps
    for(auto iter = olpret.begin(); iter != olpret.end(); ++iter){
        if(iter->second->idx.size() > 1){
            std::string outname = outdir + "/" + mNames[iter->second->idx[0]];
            for(uint32_t i = 1; i < iter->second->idx.size(); ++i){
                outname.append("_" + mNames[iter->second->idx[i]]);
            }
            outname.append(".bed");
            FILE* fp = fopen(outname.c_str(), "w");
            cr_iter_indexed(iter->second->olpret, fp, true);
            fclose(fp);
        }
    }
    // get uniq to each one
    cgranges_t* cra = cr_init();
    for(auto& e: twoins){
        cr_copyone(cra, e);
    }
    cr_index2(cra, 1);
    cgranges_t* qcr = NULL;
    for(uint32_t i = 0; i < mRegs.size(); ++i){
        qcr = mRegs[i];
        std::string outname = outdir + "/uniq_" + mNames[i] + ".bed";
        FILE* fp = fopen(outname.c_str(), "w");
        for(int32_t ctg_id = 0; ctg_id < qcr->n_ctg; ++ctg_id){
            int64_t i, *b = 0, max_b = 0, n = 0;
            n = cr_overlap_int(qcr, ctg_id, 0, INT_MAX, &b, &max_b);
            for(i = 0; i < n; ++i){
                int64_t j, *bb = 0, max_bb = 0, nn =0;
                char* ctg = qcr->ctg[ctg_id].name;
                int32_t st1 = cr_start(qcr, b[i]);
                int32_t en1 = cr_end(qcr, b[i]);
                int32_t x = 0;
                nn = cr_overlap(cra, ctg, st1, en1, &bb, &max_bb);
                for(j = 0, x = st1; j < nn; ++j){
                    cr_intv_t *r = &cra->r[bb[j]];
                    int32_t st0 = cr_st(r), en0 = cr_en(r);
                    if(st0 < st1) st0 = st1;
                    if(en0 > en1) en0 = en1;
                    if(st0 > x) fprintf(fp, "%s\t%d\t%d\n", ctg, x, st0);
                    x = en0;
                }
                if(x < en1) fprintf(fp, "%s\t%d\t%d\n", ctg, x, en1);
                free(bb);
            }
            free(b);
        }
        fclose(fp);
    }
    // release resources
    cr_destroy(cra);
    for(auto& e: olpret){
        cr_destroy(e.second->olpret);
        e.second->olpret = NULL;
        delete e.second;
        e.second = NULL;
    }
}

int main(int argc, char** argv){
    if(argc == 1){
        std::string helpCMD = std::string(argv[0]) + " -h";
        std::system(helpCMD.c_str());
        return 0;
    }
    BedRegs* pbr = new BedRegs();
    CLI::App app("program: " + std::string(argv[0]));
    app.get_formatter()->column_width(50);
    CLI::Option* pbed = app.add_option("-b,--beds", pbr->mBedPaths, "bed files to do ana")->required(false);
    CLI::Option* plst = app.add_option("-l,--list", pbr->bedlist, "bed list file")->excludes(pbed)->required(false)->check(CLI::ExistingFile);
    app.add_option("-o,--outdir", pbr->outdir, "output directory")->required(true);
    CLI11_PARSE(app, argc, argv);
    util::makedir(pbr->outdir);
    if(pbed->count() == 0 && plst->count() == 0){
        util::errorExit("-b and -l must be provides one");
    }else{
        if(plst->count()) util::makeListFromFileByLine(pbr->bedlist, pbr->mBedPaths);
    }
    pbr->loadBeds();
    pbr->olpAna();
}
