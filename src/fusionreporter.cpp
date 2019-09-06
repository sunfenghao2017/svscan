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
    fw << "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t";
    fw << "Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t";
    fw << "Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t";
    fw << "FusionSequence\tSVID\n";
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
        std::string hend = vstr[5];
        std::string tend = vstr[8];
        std::string hstrand = vstr[6];
        std::string tstrand = vstr[9];
        if(hstrand[0] != '+' || tstrand[0] != '+' || hend != "5" || tend != "3") continue;
        int32_t sr = std::atoi(vstr[18].c_str());
        int32_t dp = std::atoi(vstr[19].c_str());
        float af = std::atof(vstr[22].c_str());
        if((!fuseOpt->validFusion(hgene, tgene)) && (sr == 0 || af < fuseOpt->mMinVAF)) continue;
        if(!fuseOpt->validSV(svt, chr1, chr2, start, end)) continue;
        fw << vstr[3] << "\t"; //FusionGene
        std::string strand1 = "-";
        std::string strand2 = "-";
        if(vstr[25].find("+") != std::string::npos) strand1 = "+";
        if(vstr[26].find("+") != std::string::npos) strand2 = "+";
        fw << strand1 << strand2 << "\t"; //FusionPattern
        if(sr){//FusionReads TotalReads
            fw << sr << "\t";
            fw << (sr + std::atoi(vstr[18].c_str())) << "\t";
        }else{
            fw << dp << "\t";
            fw << (dp + std::atoi(vstr[19].c_str())) << "\t";
        }
        fw << af << "\t"; //FusionRate
        fw << hgene << "\t"; //Gene1
        fw << chr1 << "\t"; //Chr1
        fw << start << "\t"; //JunctionPosition1
        fw << strand1 << "\t"; //Strand1
        fw << vstr[25] << "\t"; //Transcript1
        fw << tgene << "\t"; //Gene2
        fw << chr2 << "\t"; //Chr2
        fw << end << "\t"; //JunctionPosition2
        fw << strand2 << "\t"; //Strand2
        fw << vstr[26] << "\t";//Transcript2
        fw << vstr[27] << "\t";//FusionSequence
        fw << vstr[29] << "\n";//SVID
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
    app.add_option("--minsr", f->fuseOpt->mMinSRSupport, "min SR support for an valid fusion", true);
    app.add_option("--mindp", f->fuseOpt->mMinDPSupport, "min DP support for an valid fusion", true);
    app.add_option("--minsu", f->fuseOpt->mMinSUSupport, "min SU support for an valid fusion", true);
    app.add_option("--minaf", f->fuseOpt->mMinVAF, "min VAF for an valid fusion", true);
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
