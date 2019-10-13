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
    // output valid fusions
    std::string header = "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t"; //[0-4]
    header.append("Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t");//[5-9]
    header.append("Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t");//[10-14]
    header.append("FusionSequence\tfseqBp\tinDB\tsvType\tsvSize\t"); //[15-19]
    header.append("srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"); //[20-25]
    header.append("insBp\tinsSeq\tsvID\tsvtInt\tfsMask\n"); //[26-30]
    std::ofstream fw(fuseOpt->mOutFile);
    std::ofstream fs(fuseOpt->mSupFile);
    fw << header;
    fs << header;
    for(auto& e: fuseList){
        if(e.fsmask & FUSION_FPRIMARY){
            fw << e;
        }
        if(e.fsmask & FUSION_FSUPPLEMENTARY){
            fs << e;
        }
    }
    fw.close();
    fs.close();
}

void FusionReporter::sv2fs(){
    fuseOpt->init();
    uint32_t REPORT_REQUEST = (FUSION_FINDB | FUSION_FHOTGENE);
    uint32_t ALL_DROP_MASK = (FUSION_FBLACKGENE | FUSION_FBLACKPAIR | FUSION_FFBG);
    uint32_t PRIMARY_DROP_MASK = (FUSION_FLOWAF | FUSION_FLOWSUPPORT | FUSION_FLOWDEPTH | FUSION_FLOWCOMPLEX | FUSION_FTOOSMALLSIZE);
    std::ifstream fr(fuseOpt->mInfile);
    std::string tmpstr;
    std::vector<std::string> vstr;
    std::getline(fr, tmpstr);
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        uint32_t fsmask = std::atoi(vstr[31].c_str());
        int32_t svt = std::atoi(vstr[30].c_str());
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
        if(!(fsmask & FUSION_FALLGENE)) continue;
        if(!fuseOpt->validSV(svt, chr1, chr2, start, end)) fsmask |= FUSION_FFBG;
        if(fuseOpt->hasWhiteGene(hgene, tgene)) fsmask |= FUSION_FHOTGENE;
        if(fuseOpt->inWhiteList(hgene, tgene)) fsmask |= FUSION_FINDB;
        if(fuseOpt->hasBlackGene(hgene, tgene)) fsmask |= FUSION_FBLACKGENE;
        if(fuseOpt->inBlackList(hgene, tgene)) fsmask |= FUSION_FBLACKPAIR;
        if(fuseOpt->matchHotDirec(hgene, tgene)) fsmask |= FUSION_FCOMMONHOTDIRECT;
        int32_t srv = std::atoi(vstr[18].c_str());
        int32_t dpv = std::atoi(vstr[19].c_str());
        int32_t srr = std::atoi(vstr[20].c_str());
        int32_t dpr = std::atoi(vstr[21].c_str());
        float af = std::atof(vstr[22].c_str());
        if(fsmask & FUSION_FINDB){// fusion in whitelist
           if((srv < fuseOpt->mWhiteFilter.mMinSupport) && (dpv < fuseOpt->mWhiteFilter.mMinSupport)){
               fsmask |= FUSION_FLOWSUPPORT;
           }
           if(af < fuseOpt->mWhiteFilter.mMinVAF){
               fsmask |= FUSION_FLOWAF;
           }
           if(((srv + srr) < fuseOpt->mWhiteFilter.mMinDepth) && ((dpr + dpv) < fuseOpt->mWhiteFilter.mMinDepth)){
               fsmask |= FUSION_FLOWDEPTH;
           }
        }else if(fsmask & FUSION_FHOTGENE){// fusion not in whitelist
            if((srv < fuseOpt->mUsualFilter.mMinSupport) && (dpv < fuseOpt->mUsualFilter.mMinSupport)){
                fsmask |= FUSION_FLOWSUPPORT;
            }
            if(af < fuseOpt->mUsualFilter.mMinVAF){
                fsmask |= FUSION_FLOWAF;
            }
            if(((srv + srr) < fuseOpt->mUsualFilter.mMinDepth) && ((dpr + dpv) < fuseOpt->mUsualFilter.mMinDepth)){
                fsmask |= FUSION_FLOWDEPTH;
            }
        }
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
        if(srv){//FusionReads TotalReads
            fsr.fusionreads = srv;
            fsr.totalreads = std::atoi(vstr[20].c_str()) + srv;
        }else{
            fsr.fusionreads = dpv;
            fsr.totalreads = std::atoi(vstr[21].c_str()) + dpv;
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
        fsr.fseqbp = vstr[28];                // fseqBp
        fsr.indb = ((fsmask & FUSION_FINDB) ? "Y" : "N"); // inDB
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
        fsr.svid = vstr[29];                  // svID
        fsr.svint = vstr[30];                 // svInt
        if(fsmask & ALL_DROP_MASK){
            fsmask &= (!(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }
        if(fsmask & REPORT_REQUEST){
            if((fsmask & FUSION_FNORMALCATDIRECT) &&
               (fsmask & FUSION_FCOMMONHOTDIRECT)){
                if(!(fsmask & PRIMARY_DROP_MASK)){
                    fsmask |= FUSION_FPRIMARY;
                }else{
                    fsmask &= (!FUSION_FPRIMARY);
                }
            }else{
                if(!(fsmask & PRIMARY_DROP_MASK)){
                    fsmask |= FUSION_FSUPPLEMENTARY;
                }else{
                    fsmask &= (!FUSION_FSUPPLEMENTARY);
                }
            }
        }else{
            fsmask &= (!(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }
        if((fsmask & fuseOpt->mFsMaskInclude) != fuseOpt->mFsMaskInclude){
            fsmask &= (!(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }
        if(fsmask & fuseOpt->mFsMaskExclude){
            fsmask &= (!(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }
        fsr.fsmask = fsmask;                  // fsMask
        fuseList.push_back(fsr);
    }
    fr.close();
}
