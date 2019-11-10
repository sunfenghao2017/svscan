#include "svcreg.h"

void SVCreg::getCreg(){
    // fetch gene names
    std::set<std::string> gnameSet;
    util::makeSetFromFileByField(glist, gnameSet, 0);
    // fetch all valid regions and merge
    cgranges_t *cr = cr_init();
    std::set<std::string> geneGot;
    std::ifstream fr(ibed);
    std::string line;
    std::vector<std::string> vstr;
    int64_t i = 0;
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        std::string gname = vstr[5];
        auto iter = gnameSet.find(gname);
        if(iter != gnameSet.end()){
            cr_add(cr, vstr[0].c_str(), std::atoi(vstr[1].c_str()), std::atoi(vstr[2].c_str()), i++);
            geneGot.insert(gname);
        }
    }
    fr.close();
    if(!cr_is_sorted(cr)) cr_sort(cr);
    cr_merge_pre_index(cr);
    // write output
    std::ofstream fw(obed);
    for(i = 0; i < cr->n_r; ++i){
        const cr_intv_t *r = &cr->r[i];
        fw << cr->ctg[r->x>>32].name << "\t" << (int32_t)r->x << "\t" << (int32_t)r->y << "\n";
    }
    fw.close();
    cr_destroy(cr);
}
