#include "options.h"

Options::Options(){
    madCutoff = 9;
    nthread = 8;
    batchsvn = 1000;
    tsvOut = "sv.tsv";
    pool = NULL;
    fbamout = NULL;
    overlapRegs = NULL;
    bamheader = NULL;
    filterOpt = new SVFilter();
    softEnv = new Software();
    msaOpt = new MSAOpt();
    passOpt = new PassOptions();
    fuseOpt = new FusionOptions();
    libInfo = NULL;
    realnf = new RealnFilter();
    contigNum = 0;
    rnamode = false;
    debug = false;
}

Options::~Options(){
    if(filterOpt) delete filterOpt;
    if(softEnv) delete softEnv;
    if(libInfo) delete libInfo;
    if(realnf) delete realnf;
    if(pool) delete pool;
    if(overlapRegs) delete overlapRegs;
    if(bamheader) bam_hdr_destroy(bamheader);
}

void Options::validate(){
    // check genome indexed
    if(!util::exists(alnref + ".fai")){
        util::errorExit("Genome must be indexed by `samtools faidx`");
    }
}

void Options::update(int argc, char** argv){
    // update software environment records
    softEnv->update(argc, argv);
    // get library information
    libInfo = getLibInfo(bamfile);
    // update SV types to discover
    std::vector<std::string> svt = {"INV", "DEL", "DUP", "INS", "BND"};
    std::string allSVT;
    if(svtypes.size() == 0){
        for(int i = 0; i < 9; ++i) SVTSet.insert(i);
        allSVT.append("INV DEL DUP INS BND");
    }else{
        for(auto& e: svtypes){
            allSVT.append(svt[e] + " ");
            if(e == 0){// Inversion
                SVTSet.insert(0);
                SVTSet.insert(1);
            }else if(e == 1){// Deletion
                SVTSet.insert(2);
            }else if(e == 2){// Duplication
                SVTSet.insert(3);
            }else if(e == 3){// Insertion
                SVTSet.insert(4);
            }else if(e == 4){// translocation
                SVTSet.insert({5, 6, 7, 8});
            }
        }
    }
    // adjust sr score
    if(rnamode){
        filterOpt->mMinSRResScore = 0.99;
        batchsvn = 500;
    }
    // update contigNum
    contigNum = libInfo->mContigNum;
    // construct threadpool
    pool = new ThreadPool::ThreadPool(std::min(contigNum, nthread));
    // show SV types to discover
    util::loginfo("SV types to discover: " + allSVT);
    // update scan regs
    getScanRegs();
    // update cregs
    getCregs();
    // write bcf or not
    if(bcfOut.empty()) writebcf = false;
    else writebcf = true;
    // parse bam header
    samFile* fp = sam_open(bamfile.c_str(), "r");
    bamheader = sam_hdr_read(fp);
    sam_close(fp);
}

LibraryInfo* Options::getLibInfo(const std::string& bam){
    LibraryInfo* libInfo = new LibraryInfo();
    int32_t nread = 0;
    int32_t ttread = 0;
    samFile* fp = sam_open(bam.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    bam1_t* b = bam_init1();
    std::vector<int32_t> vecISize, readLen;
    const uint16_t BAM_SKIP_RECORD_MASK = (BAM_FREAD2 | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP);
    while(sam_read1(fp, h, b) >= 0){
        if(b->core.flag & BAM_SKIP_RECORD_MASK) continue;
        ++ttread;
        readLen.push_back(b->core.l_qseq);
        if(b->core.flag & BAM_FPAIRED && b->core.tid == b->core.mtid){
            vecISize.push_back(std::abs(b->core.isize));
            if(++nread > libInfo->mMaxSample) break;
        }
    }
    if(!ttread){
        sam_close(fp);
        bam_destroy1(b);
        bam_hdr_destroy(h);
        delete libInfo;
        writeEmptFile();
        util::loginfo("Bam is empty, empty result file written, program will quit now!!!");
        exit(EXIT_SUCCESS);
    }
    libInfo->mReadLen = statutil::median(readLen);
    if(!vecISize.empty()) libInfo->mMedian = statutil::median(vecISize);
    if(!vecISize.empty()) libInfo->mMad = statutil::mad(vecISize, libInfo->mMedian);
    libInfo->mMaxNormalISize = libInfo->mMedian + (5 * libInfo->mMad);
    libInfo->mMinNormalISize = std::max(0, libInfo->mMedian - (5 * libInfo->mMad));
    libInfo->mMinISizeCutoff = std::max(0, libInfo->mMedian - (madCutoff * libInfo->mMad));
    libInfo->mMaxISizeCutoff = std::max(libInfo->mMedian + (madCutoff * libInfo->mMad), 2 * libInfo->mReadLen);
    libInfo->mMaxISizeCutoff = std::max(libInfo->mMaxISizeCutoff, 600);
    libInfo->mContigNum = h->n_targets;
    if(libInfo->mMaxNormalISize == 0) libInfo->mMaxNormalISize = 3 * libInfo->mReadLen;
    libInfo->mVarisize = std::max(libInfo->mReadLen, libInfo->mMaxNormalISize);
    sam_close(fp);
    bam_destroy1(b);
    bam_hdr_destroy(h);
    return libInfo;
}

void Options::getScanRegs(){
    samFile* fp = sam_open(bamfile.c_str(), "r");
    bam_hdr_t* h = sam_hdr_read(fp);
    hts_idx_t* idx = sam_index_load(fp, bamfile.c_str());
    scanRegs.resize(h->n_targets);
    // Parse valid region if exists
    if(!reg.empty()){
        std::ifstream fr(reg);
        std::vector<std::string> vstr;
        std::string tmpStr;
        while(std::getline(fr, tmpStr)){
            util::split(tmpStr, vstr, "\t");
            int32_t tid = bam_name2id(h, vstr[0].c_str());
            scanRegs[tid].insert({std::atoi(vstr[1].c_str()), std::atoi(vstr[2].c_str())});
        }
    }else{// Set valid region to whole genome if valid region list does not exists
        for(int32_t i = 0; i < h->n_targets; ++i){
            uint64_t mapped = 0, unmapped = 0;
            if(hts_idx_get_stat(idx, i, &mapped, &unmapped) >= 0){
                if(mapped > 0){
                    scanRegs[i].insert({0, (int32_t)h->target_len[i] - 1});
                }
            }
        }
    }
    // Clean-up
    sam_close(fp);
    bam_hdr_destroy(h);
}

void Options::getCregs(){
    if(!creg.empty()){
        overlapRegs = new BedRegs();
        overlapRegs->loadBed(creg);
    }
}

void Options::writeEmptFile(){
    // tsv sv
    std::ofstream fw(tsvOut);
    fw << "svType\tsvSize\tbpMark\t";//[0, 2]
    fw << "fuseGene\t"; //[3];
    fw << "hGene\thTrsEnd\thTrsStrand\t";//[4,6]
    fw << "tGene\ttTrsEnd\ttTrsStrand\t";//[7,9]
    fw << "bp1Chr\tbp1Pos\tbp1Gene\t";//[10,12]
    fw << "bp2Chr\tbp2Pos\tbp2Gene\t";//[13,15]
    fw << "srCount\tdpCount\tsrRescued\tdpRescued\t";//[16,19]
    fw << "srRefCount\tdpRefCount\tAF\tinsBp\tinsSeq\t";//[20,24]
    fw << "bp1Trs\tbp2Trs\tsvSeq\tseqBp\tID\tsvtInt\tfsMask";//[25,31]
    if(rnamode){
        fw << "\tts1Name\tts1Pos\tts2Name\tts2Pos\n"; //[32,35]
    }else{
        fw << "\n";
    }
    fw.close();
    // fusion
    std::string header = "FusionGene\tFusionPattern\tFusionReads\tTotalReads\tFusionRate\t"; //[0-4]
    header.append("Gene1\tChr1\tJunctionPosition1\tStrand1\tTranscript1\t");//[5-9]
    header.append("Gene2\tChr2\tJunctionPosition2\tStrand2\tTranscript2\t");//[10-14]
    header.append("FusionSequence\tfseqBp\tinDB\tsvType\tsvSize\t"); //[15-19]
    header.append("srCount\tdpCount\tsrRescued\tdpRescued\tsrRefCount\tdpRefCount\t"); //[20-25]
    header.append("insBp\tinsSeq\tsvID\tsvtInt\tfsMask"); //[26-29]
    if(rnamode){
        header.append("\tts1Name\tts1Pos\tts2Name\tts2Pos\n"); //[30-33]
    }else{
        header.append("\n");
    }
    fw.open(fuseOpt->mOutFile.c_str());
    fw << header;
    fw.close();
    fw.open(fuseOpt->mSupFile.c_str());
    fw << header;
    fw.close();
    // sv bcf
    if(bcfOut.empty()) return;
    samFile* samfp = sam_open(bamfile.c_str(), "r");
    hts_set_fai_filename(samfp, alnref.c_str());
    bam_hdr_t* bamhdr = sam_hdr_read(samfp);
    htsFile* fp = bcf_open(bcfOut.c_str(), "wb");
    bcf_hdr_t* hdr = bcf_hdr_init("w");
    // Output bcf header
    bcf_hdr_append(hdr, "##ALT=<ID=DEL,Description=\"Deletion\">");
    bcf_hdr_append(hdr, "##ALT=<ID=DUP,Description=\"Duplication\">");
    bcf_hdr_append(hdr, "##ALT=<ID=INV,Description=\"Inversion\">");
    bcf_hdr_append(hdr, "##ALT=<ID=BND,Description=\"Translocation\">");
    bcf_hdr_append(hdr, "##ALT=<ID=INS,Description=\"Insertion\">");
    bcf_hdr_append(hdr, "##FILTER=<ID=LowQual,Description=\"Poor quality and insufficient number of PEs and SRs.\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PEMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRALNQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">");
    bcf_hdr_append(hdr, "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=INSLEN,Number=1,Type=Integer,Description=\"Predicted length of the insertion\">");
    bcf_hdr_append(hdr, "##INFO=<ID=HOMLEN,Number=1,Type=Integer,Description=\"Predicted microhomology length using a max. edit distance of 2\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the SV\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RCL,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the left control region\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RCR,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the right control region\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Read-depth based copy-number estimate for autosomal sites\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">");
    // Add reference
    std::string refloc = "##reference=" + alnref;
    bcf_hdr_append(hdr, refloc.c_str());
    for(int i = 0; i < bamhdr->n_targets; ++i){
        std::string ctginfo("##contig=<ID=");
        ctginfo.append(std::string(bamhdr->target_name[i]) + ",length=" + std::to_string(bamhdr->target_len[i]) + ">");
        bcf_hdr_append(hdr, ctginfo.c_str());
    }
    // Add runtime
    std::string dateStr = "##fileDate=" + util::currentTime();
    bcf_hdr_append(hdr, dateStr.c_str());
    // Add software information
    std::string sinf = "##softwoare=sver v" + softEnv->ver;
    std::string scmd = "##command=" + softEnv->cmd;
    std::string scwd = "##cwd=" + softEnv->cwd;
     bcf_hdr_append(hdr, sinf.c_str());
    bcf_hdr_append(hdr, scmd.c_str());
    bcf_hdr_append(hdr, scwd.c_str());
    // Add Samples
    bcf_hdr_add_sample(hdr, "SAMPLE");
    assert(bcf_hdr_write(fp, hdr) >= 0);
    // Add Records
    sam_close(samfp);
    bam_hdr_destroy(bamhdr);
    bcf_hdr_destroy(hdr);
    bcf_close(fp);
    bcf_index_build(bcfOut.c_str(), 14);
}
