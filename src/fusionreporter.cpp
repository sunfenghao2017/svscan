#include "fusionreporter.h"

FusionReporter::FusionReporter(){
    fuseOpt = new FusionOptions();
    softEnv = new Software();
    softEnv->cmp += "version: " + softEnv->version + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
}

FusionReporter::~FusionReporter(){
    if(softEnv) delete softEnv;
    if(fuseOpt) delete fuseOpt;
}

void FusionReporter::update(int argc, char** argv){
    // update software environment records
    softEnv->cwd = util::cwd();
    for(int i = 0; i < argc; ++i){
        softEnv->cmd.append(argv[i]);
        softEnv->cmd.append(" ");
    }
}

void FusionReporter::report(){
    sv2fs();
    // mask fusion pair which has not proper directions
    std::map<std::string, std::set<std::string>> fpairs;
    for(uint32_t i = 0; i < fuseList.size(); ++i){
        if(fuseList[i].gene1 != fuseList[i].gene2){
            fpairs[fuseList[i].gene1].insert(fuseList[i].gene2);
        }
    }
    for(uint32_t i = 0; i < fuseList.size(); ++i){
        if(fuseList[i].gene1 == fuseList[i].gene2) continue;
        std::string hg = fuseList[i].gene1;
        std::string tg = fuseList[i].gene2;
        auto titer = fpairs.find(tg);
        if(titer == fpairs.end() || titer->second.find(hg) == titer->second.end()) continue; // no mirror fusion
        if((fuseOpt->m5Partners.find(hg) == fuseOpt->m5Partners.end()) &&
           (fuseOpt->m3Partners.find(tg) == fuseOpt->m3Partners.end())){
            fuseList[i].report = false;
        }
    }
    std::ofstream fw(fuseOpt->mOutFile);
    fw << "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t"; //[0-4]
    fw << "Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t"; //[5-9]
    fw << "Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t"; //[10-14]
    fw << "FusionSequence\tinDB\tsvType\tsvSize\t"; //[15-18]
    fw << "srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"; //[19-24]
    fw << "insBp\tinsSeq\tsvID\tsvtInt\n";//[25-28]
    for(auto& e: fuseList){
        if(e.report){
            fw << e;
        }
    }
    fw.close();
}

void FusionReporter::sv2fs(){
    fuseOpt->init();
    std::ifstream fr(fuseOpt->mInfile);
    std::string tmpstr;
    std::vector<std::string> vstr;
    std::getline(fr, tmpstr);
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        int32_t svt = std::atoi(vstr[29].c_str());
        int32_t start = std::stoi(vstr[11].c_str());
        int32_t end = std::atoi(vstr[14].c_str());
        std::string chr1 = vstr[10];
        std::string chr2 = vstr[13];
        std::string fusegene = vstr[3];
        std::string hgene = vstr[4];
        std::string tgene = vstr[7];
        std::string gene1 = vstr[12];
        std::string hend = vstr[5];
        std::string tend = vstr[8];
        std::string hstrand = vstr[6];
        std::string tstrand = vstr[9];
        // keep only(hgene+5'->tgene+3') fusion
        if(hstrand[0] != '+' || tstrand[0] != '+' || hend != "5" || tend != "3") continue;
        // skip fusion gene which has no partner in white gene list
        if(!fuseOpt->hasWhiteGene(hgene, tgene)) continue;
        // skip fusion gene which has any partner in black gene list
        if(fuseOpt->hasBlackGene(hgene, tgene)) continue;
        int32_t sr = std::atoi(vstr[18].c_str());
        int32_t dp = std::atoi(vstr[19].c_str());
        int32_t srr = std::atoi(vstr[20].c_str());
        int32_t dpr = std::atoi(vstr[21].c_str());
        float af = std::atof(vstr[22].c_str());
        if(fuseOpt->inBlackList(hgene, tgene)) continue; // skip fusion in blacklist
        bool inWhiteList = fuseOpt->inWhiteList(hgene, tgene);
        bool keep = false;
        int32_t svsize = std::atoi(vstr[1].c_str());
        if(inWhiteList){// fusion in whitelist
           if((sr >= fuseOpt->mWhiteFilter.mMinSupport || dp >= fuseOpt->mWhiteFilter.mMinSupport) &&
              (af >= fuseOpt->mWhiteFilter.mMinVAF)) keep = true;
           if((sr + srr) < fuseOpt->mWhiteFilter.mMinDepth && ((dpr + dp) < fuseOpt->mWhiteFilter.mMinDepth)) keep = false;
           if(hgene == tgene && svsize < fuseOpt->mWhiteFilter.mMinIntraGeneSVSize) keep = false;
        }else{// fusion not in whitelist
            if((sr >= fuseOpt->mUsualFilter.mMinSupport) && (af > fuseOpt->mUsualFilter.mMinVAF)) keep = true;
            if((sr + srr) < fuseOpt->mUsualFilter.mMinDepth) keep = false;
            if(hgene == tgene && svsize < fuseOpt->mUsualFilter.mMinIntraGeneSVSize) keep = false;
        }
        if(!keep) continue;
        // skip fusion in background
        if(!fuseOpt->validSV(svt, chr1, chr2, start, end)) continue;
        FusionRecord fsr;
        fsr.fusegene = vstr[3]; // FusionGene
        int cnt[2] = {0, 0};
        for(auto& e: vstr[25]){
            if(e == '+') cnt[0] += 1;
            if(e == '-') cnt[1] += 1;
        }
        std::string strand1 = (cnt[0] > cnt[1] ? "+" : "-");
        cnt[0] = 0;
        cnt[1] = 0;
        for(auto& e: vstr[26]){
            if(e == '+') cnt[0] += 1;
            if(e == '-') cnt[1] += 1;
        }
        std::string strand2 = (cnt[0] > cnt[1] ? "+" : "-");
        if(hgene == gene1){//FusionPattern
            fsr.fusepattern.append(strand1);
            fsr.fusepattern.append(strand2);
        }else{
            fsr.fusepattern.append(strand2);
            fsr.fusepattern.append(strand1);
        }
        if(sr){//FusionReads TotalReads
            fsr.fusionreads = sr;
            fsr.totalreads = std::atoi(vstr[20].c_str()) + sr;
        }else{
            fsr.fusionreads = dp;
            fsr.totalreads = std::atoi(vstr[21].c_str()) + dp;
        }
        fsr.fuserate = af;                    // FusionRate
        fsr.gene1 = hgene;
        fsr.gene2 = tgene;
        if(gene1 == hgene){
            fsr.chr1 = chr1;                  // Chr1
            fsr.junctionposition1 = start;    // JunctionPosition1
            fsr.strand1 = strand1;            // Strand1;
            fsr.transcript1 = vstr[25];       // Transcript1
            fsr.chr2 = chr2;                  // Chr2
            fsr.junctionposition2 = end;      // JunctionPosition2
            fsr.strand2 = strand2;            // Strand2
            fsr.transcript2 = vstr[26];       // Transcript2
        }else{
            fsr.chr1 = chr2;                  // Chr1
            fsr.junctionposition1 = end;      // JunctionPosition1
            fsr.strand1 = strand2;            // Strand1;
            fsr.transcript1 = vstr[26];       // Transcript1
            fsr.chr2 = chr1;                  // Chr2
            fsr.junctionposition2 = start;    // JunctionPosition2
            fsr.strand2 = strand1;            // Strand2
            fsr.transcript2 = vstr[25];       // Transcript2
        }
        fsr.fusionsequence = vstr[27];        // FusionSequence
        fsr.indb = (inWhiteList ? "Y" : "N"); // inDB
        fsr.svt = vstr[0];                    // svType
        fsr.svsize = vstr[1];                 // svSize
        fsr.srcount = vstr[16];               // srCount
        fsr.dpcount = vstr[17];               // dpCount
        fsr.srrescued = vstr[18];             // srRescued
        fsr.dprescued = vstr[19];             // dpRescued
        fsr.srrefcount  = vstr[20];           // srRefCount
        fsr.dprefcount = vstr[21];            // dpRefCount
        fsr.insbp = vstr[23];                 // insBp
        fsr.insseq = vstr[24];                // insSeq
        fsr.svid = vstr[28];                  // svID
        fsr.svint = vstr[29];                 // svInt
        fsr.report = true;
        fuseList.push_back(fsr);
    }
    fr.close();
}

int main(int argc, char** argv){
    if(argc == 1){
        std::string helpCMD = std::string(argv[0]) + " -h";
        std::system(helpCMD.c_str());
        return 0;
    }
    // parse commandline arguments
    FusionReporter* f = new FusionReporter();
    CLI::App app("program: " + std::string(argv[0]) + "\n" + f->softEnv->cmp);
    app.add_option("-i,--in", f->fuseOpt->mInfile, "input tsv sv result of sver")->required(true)->check(CLI::ExistingFile);
    app.add_option("-o,--out", f->fuseOpt->mOutFile, "output tsv fusion result", true);
    app.add_option("--whitemindep", f->fuseOpt->mWhiteFilter.mMinDepth, "min depth for an valid fusion break point in whitelist", true);
    app.add_option("--usualmindep", f->fuseOpt->mUsualFilter.mMinDepth, "min depth for an valid fusion break point ont in whitelist", true);
    app.add_option("--whiteminr", f->fuseOpt->mWhiteFilter.mMinSupport, "min reads support for an valid fusion in whitelist", true);
    app.add_option("--usualminr", f->fuseOpt->mUsualFilter.mMinSupport, "min reads support for an valid fusion not in whitelist", true);
    app.add_option("--whiteminaf", f->fuseOpt->mWhiteFilter.mMinVAF, "min VAF for an valid fusion in whitelist", true);
    app.add_option("--usualminaf", f->fuseOpt->mUsualFilter.mMinVAF, "min VAF for an valid fusion not in whitelist", true); 
    app.add_option("--whiteminigs", f->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion in whitelist", true);
    app.add_option("--usualminigs", f->fuseOpt->mUsualFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion not in whitelist", true);
    app.add_option("--maxbpoffset", f->fuseOpt->mMaxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true);
    app.add_option("--bgbcf", f->fuseOpt->mBgBCF, "background events BCF file");
    app.add_option("--whitelist", f->fuseOpt->mWhiteList, "white list of fusion events")->check(CLI::ExistingFile);
    app.add_option("--blacklist", f->fuseOpt->mBlackList, "black list of fusion events")->check(CLI::ExistingFile);
    // parse arguments
    CLI_PARSE(app, argc, argv);
    f->update(argc, argv);
    util::loginfo("CMD: " + f->softEnv->cmd);
    util::loginfo("Beg report fusions.");
    f->report();
    util::loginfo("End report fusions.");
    delete f;
}
