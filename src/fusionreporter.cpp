#include "fusionreporter.h"

FusionReporter::FusionReporter(){
    fuseOpt = new FusionOptions();
    softEnv = new Software();
    softEnv->cmp += "version: " + softEnv->ver + "\n";
    softEnv->cmp += "updated: " + std::string(__TIME__) + " " + std::string(__DATE__);
    rnamode = false;
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
    header.append("insBp\tinsSeq\tsvID\tsvtInt\tfsMask"); //[26-30]
    if(rnamode) header.append("\tts1Name\tts1Pos\tts2Name\tts2Pos\n"); //[31-34]
    else header.append("\n");
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
    std::ofstream fsv;
    if(!fuseOpt->mSVModFile.empty()){ // update sv tsv file if needed
        fsv.open(fuseOpt->mSVModFile.c_str());
    }
    // drop bits mask of all fusion events, if an fusion match any bit in FUSION_DROP_MASK, it will not be reported
    TFUSION_FLAG FUSION_DROP_MASK = (FUSION_FBLACKGENE | FUSION_FBLACKPAIR  | FUSION_FFBG | FUSION_FLOWCOMPLEX |
                                     FUSION_FTOOSMALLSIZE | FUSION_FLOWAF | FUSION_FLOWDEPTH);
    // primary keep bits mask, fusion reported as primary must match all the bits in PRIMARY_KEEP_MASK
    TFUSION_FLAG PRIMARY_KEEP_MASK = (FUSION_FNORMALCATDIRECT | FUSION_FCOMMONHOTDIRECT | FUSION_FINDB);
    // keep bits mask, an fusion to be reported must match all bits in FUSION_KEEP_MASK
    TFUSION_FLAG FUSION_KEEP_MASK = (FUSION_FALLGENE | FUSION_FHOTGENE);
    std::ifstream fr(fuseOpt->mInfile);
    std::string tmpstr;
    std::vector<std::string> vstr;
    std::getline(fr, tmpstr);
    if(!fuseOpt->mSVModFile.empty()) fsv << tmpstr << "\n";
    while(std::getline(fr, tmpstr)){
        util::split(tmpstr, vstr, "\t");
        TFUSION_FLAG fsmask = std::atoi(vstr[31].c_str());
        if(!rnamode && (fsmask & FUSION_FCALLFROMRNASEQ)){
            rnamode = true;
            FUSION_DROP_MASK |= FUSION_FINSAMEGENE;
        }
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
        if(!fuseOpt->validSV(svt, chr1, chr2, start, end)) fsmask |= FUSION_FFBG;
        else fsmask &= (~FUSION_FFBG);
        if(fuseOpt->hasWhiteGene(hgene, tgene)) fsmask |= FUSION_FHOTGENE;
        else fsmask &= (~FUSION_FHOTGENE);
        if(fuseOpt->inWhiteList(hgene, tgene)) fsmask |= FUSION_FINDB;
        else fsmask &= (~FUSION_FINDB);
        if(fuseOpt->inWhiteList(tgene, hgene)) fsmask |= FUSION_FMIRRORINDB;
        else fsmask &= (~FUSION_FMIRRORINDB);
        if(fuseOpt->hasBlackGene(hgene, tgene)) fsmask |= FUSION_FBLACKGENE;
        else fsmask &= (~FUSION_FBLACKGENE);
        if(fuseOpt->inBlackList(hgene, tgene)) fsmask |= FUSION_FBLACKPAIR;
        else fsmask &= (~FUSION_FBLACKPAIR);
        if(fuseOpt->matchHotDirec(hgene, tgene)) fsmask |= FUSION_FCOMMONHOTDIRECT;
        else fsmask &= (~FUSION_FCOMMONHOTDIRECT);
        int32_t srv = std::atoi(vstr[18].c_str());
        int32_t dpv = std::atoi(vstr[19].c_str());
        int32_t srr = std::atoi(vstr[20].c_str());
        int32_t dpr = std::atoi(vstr[21].c_str());
        float af = std::atof(vstr[22].c_str());
        if(fsmask & (FUSION_FINDB | FUSION_FMIRRORINDB)){// fusion in public database
           if((srv < fuseOpt->mWhiteFilter.mMinSupport) && (dpv < fuseOpt->mWhiteFilter.mMinSupport)){
               fsmask |= FUSION_FLOWSUPPORT;
           }else{
               fsmask &= (~FUSION_FLOWSUPPORT);
           }
           if(af < fuseOpt->mWhiteFilter.mMinVAF){
               fsmask |= FUSION_FLOWAF;
           }else{
               fsmask &= (~FUSION_FLOWAF);
           }
           if(((srv + srr) < fuseOpt->mWhiteFilter.mMinDepth) && ((dpr + dpv) < fuseOpt->mWhiteFilter.mMinDepth)){
               fsmask |= FUSION_FLOWDEPTH;
           }else{
               fsmask &= (~FUSION_FLOWDEPTH);
           }
        }else if(fsmask & FUSION_FHOTGENE){// fusion not in public database
            if((srv < fuseOpt->mUsualFilter.mMinSupport) && (dpv < fuseOpt->mUsualFilter.mMinSupport)){
                fsmask |= FUSION_FLOWSUPPORT;
            }else{
                fsmask &= (~FUSION_FLOWSUPPORT);
            }
            if(af < fuseOpt->mUsualFilter.mMinVAF){
                fsmask |= FUSION_FLOWAF;
            }else{
                fsmask &= (~FUSION_FLOWAF);
            }
            if(((srv + srr) < fuseOpt->mUsualFilter.mMinDepth) && ((dpr + dpv) < fuseOpt->mUsualFilter.mMinDepth)){
                fsmask |= FUSION_FLOWDEPTH;
            }else{
                fsmask &= (~FUSION_FLOWDEPTH);
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
        if(fsmask & FUSION_DROP_MASK){
            fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }else{
            if((fsmask & FUSION_KEEP_MASK) != FUSION_KEEP_MASK){
                fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
            }else{
                if((fsmask & PRIMARY_KEEP_MASK) == PRIMARY_KEEP_MASK){
                    fsmask |= FUSION_FPRIMARY;
                }else{
                    fsmask |= FUSION_FSUPPLEMENTARY;
                }
            }
        }
        if((fsmask & fuseOpt->mFsMaskInclude) != fuseOpt->mFsMaskInclude){
            fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }
        if(fsmask & fuseOpt->mFsMaskExclude){
            fsmask &= (~(FUSION_FPRIMARY | FUSION_FSUPPLEMENTARY));
        }
        fsr.fsmask = fsmask;                 // fsMask
        if(fsmask & FUSION_FCALLFROMRNASEQ){
            fsr.ts1name = vstr[32];          // ts1Name
            fsr.ts1pos = vstr[33];           // ts1Pos
            fsr.ts2name = vstr[34];          // ts2Name
            fsr.ts2pos = vstr[35];           // ts2Pos
        }
        fuseList.push_back(fsr);
        if(!fuseOpt->mSVModFile.empty()){
            vstr[31] = std::to_string(fsmask);
            std::string newrec;
            util::join(vstr, newrec, "\t");
            fsv << newrec << "\n";
        }
    }
    if(!fuseOpt->mSVModFile.empty()) fsv.close();
    fr.close();
}
