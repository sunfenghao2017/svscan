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
    fuseOpt->init();
    std::ifstream fr(fuseOpt->mInfile);
    std::string tmpstr;
    std::vector<std::string> vstr;
    std::getline(fr, tmpstr);
    std::ofstream fw(fuseOpt->mOutFile);
    fw << "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t"; //[0-4]
    fw << "Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t"; //[5-9]
    fw << "Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t"; //[10-14]
    fw << "FusionSequence\tHotFusion\tsvType\tsvSize\t"; //[15-18]
    fw << "srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"; //[19-24]
    fw << "insBp\tinsSeq\tsvID\tsvtInt\n";//[25-28]
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
        int32_t sr = std::atoi(vstr[18].c_str());
        int32_t dp = std::atoi(vstr[19].c_str());
        int32_t srr = std::atoi(vstr[20].c_str());
        int32_t dpr = std::atoi(vstr[21].c_str());
        float af = std::atof(vstr[22].c_str());
        if(fuseOpt->inBlackList(hgene, tgene)) continue; // skip fusion in blacklist
        bool inWhiteList = fuseOpt->inWhiteList(hgene, tgene);
        bool keep = false;
        if(inWhiteList){// fusion in whitelist
           if((sr >= fuseOpt->mWhiteFilter.mMinSupport || dp >= fuseOpt->mWhiteFilter.mMinSupport) &&
              (af >= fuseOpt->mWhiteFilter.mMinVAF)) keep = true;
           if((sr + srr) < fuseOpt->mWhiteFilter.mMinDepth && ((dpr + dp) < fuseOpt->mWhiteFilter.mMinDepth)){
               keep = false;
           }
        }else{// fusion not in whitelist
            if((sr >= fuseOpt->mUsualFilter.mMinSupport || dp >= fuseOpt->mUsualFilter.mMinSupport) && 
               (af > fuseOpt->mUsualFilter.mMinVAF)) keep = true;
            if((sr + srr) < fuseOpt->mUsualFilter.mMinDepth && ((dpr + dp) < fuseOpt->mUsualFilter.mMinDepth)){
                keep = false;
            }
        }
        if(!keep) continue;
        // skip fusion in background
        if(!fuseOpt->validSV(svt, chr1, chr2, start, end)) continue;
        fw << vstr[3] << "\t"; //FusionGene
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
            fw << strand1 << strand2 << "\t";
        }else{
            fw << strand2 << strand1 << "\t";
        }
        if(sr){//FusionReads TotalReads
            fw << sr << "\t";
            fw << (sr + std::atoi(vstr[20].c_str())) << "\t";
        }else{
            fw << dp << "\t";
            fw << (dp + std::atoi(vstr[21].c_str())) << "\t";
        }
        fw << af << "\t"; //FusionRate
        if(gene1 == hgene){
            fw << hgene << "\t";    // Gene1
            fw << chr1 << "\t";     // Chr1
            fw << start << "\t";    // JunctionPosition1
            fw << strand1 << "\t";  // Strand1
            fw << vstr[25] << "\t"; // Transcript1
            fw << tgene << "\t";    // Gene2
            fw << chr2 << "\t";     // Chr2
            fw << end << "\t";      // JunctionPosition2
            fw << strand2 << "\t";  // Strand2
            fw << vstr[26] << "\t"; // Transcript2
        }else{
            fw << hgene << "\t";    // Gene2
            fw << chr2 << "\t";     // Chr2
            fw << end << "\t";      // JunctionPosition2
            fw << strand2 << "\t";  // Strand2
            fw << vstr[26] << "\t"; // Transcript2
            fw << tgene << "\t";    // Gene1
            fw << chr1 << "\t";     // Chr1
            fw << start << "\t";    // JunctionPosition1
            fw << strand1 << "\t";  // Strand1
            fw << vstr[25] << "\t"; // Transcript1
        }
        fw << vstr[27] << "\t";  // FusionSequence
        if(inWhiteList) fw << "Y\t"; // HotFusion
        else fw << "N\t";
        fw << vstr[0] << "\t";   // svType
        fw << vstr[1] << "\t";   // svSize
        fw << vstr[16]  << "\t"; // srCount
        fw << vstr[17] << "\t";  // dpCount
        fw << vstr[18] << "\t";  // srRescued
        fw << vstr[19] << "\t";  // dpRescued
        fw << vstr[20] << "\t";  // srRefCount
        fw << vstr[21] << "\t";  // dpRefCount
        fw << vstr[23] << "\t";  // insBp
        fw << vstr[24] << "\t";  // insSeq
        fw << vstr[28] << "\t";  // svID
        fw << vstr[29] << "\n";  // svtInt
    }
    fr.close();
    fw.close();
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
