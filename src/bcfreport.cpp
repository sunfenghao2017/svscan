#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "stats.h"
#include "util.h"

void Stats::reportSVBCF(const SVSet& svs){
    // Open file handler
    samFile* samfp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_set_fai_filename(samfp, mOpt->genome.c_str());
    bam_hdr_t* bamhdr = sam_hdr_read(samfp);
    htsFile* fp = bcf_open(mOpt->bcfOut.c_str(), "wb");
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
    bcf_hdr_append(hdr, "##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">");
    if(mOpt->libInfo->mIsHaploTagged){
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP1DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs on haplotype 1\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP2DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs on haplotype 2\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP1DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs on haplotype 1\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP2DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs on haplotype 2\">");
    }
    bcf_hdr_append(hdr, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">");
    if(mOpt->libInfo->mIsHaploTagged){
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP1RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads on haplotype 1\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP2RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads on haplotype 2\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP1RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads on haplotype 1\">");
        bcf_hdr_append(hdr, "##FORMAT=<ID=HP2RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads on haplotype 2\">");
    }
    // Add reference
    std::string refloc = "##reference=" + mOpt->genome;
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
    std::string sinf = "##softwoare=sver v" + mOpt->softEnv->ver;
    std::string scmd = "##command=" + mOpt->softEnv->cmd;
    std::string scwd = "##cwd=" + mOpt->softEnv->cwd;
    bcf_hdr_append(hdr, sinf.c_str());
    bcf_hdr_append(hdr, scmd.c_str());
    bcf_hdr_append(hdr, scwd.c_str());
    // Add Samples
    bcf_hdr_add_sample(hdr, "SAMPLE");
    bcf_hdr_write(fp, hdr);
    // Add Records
    if(svs.empty()){// if empty set, just return
        sam_close(samfp);
        bam_hdr_destroy(bamhdr);
        bcf_hdr_destroy(hdr);
        bcf_close(fp);
        bcf_index_build(mOpt->bcfOut.c_str(), 14);
        return;
    }
    // Genotype arrays
    int* gts = (int32_t*)std::calloc(bcf_hdr_nsamples(hdr) * 2,  sizeof(int));
    float* gls = (float*)std::calloc(bcf_hdr_nsamples(hdr) * 3, sizeof(float));
    int* drcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* dvcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp1drcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp2drcount = (int*) std::calloc(bcf_hdr_nsamples(hdr),sizeof(int));
    int* hp1dvcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp2dvcount = (int*) std::calloc(bcf_hdr_nsamples(hdr),sizeof(int));
    int* rrcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* rvcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp1rrcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp2rrcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp1rvcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* hp2rvcount = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    int* gqval = (int*) std::calloc(bcf_hdr_nsamples(hdr), sizeof(int));
    const char* ftarr = {NULL};
    bcf1_t* rec = bcf_init1();
    int32_t tmpi = 0;
    for(auto itsv = svs.begin(); itsv != svs.end(); ++itsv){
        // Prepare Filter field
        int filter = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
        if(itsv->mChr1 == itsv->mChr2){// Intra-chromosomal sv
            if((itsv->mPESupport < mOpt->passOpt->mIntraChrSVMinDPCnt || itsv->mPEMapQuality < mOpt->passOpt->mIntraChrSVMinDPQual) ||
               (itsv->mSRSupport < mOpt->passOpt->mIntraChrSVMinSRCnt || itsv->mSRMapQuality < mOpt->passOpt->mIntraChrSVMinSRQual)){
                filter = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
            }
        }else{// inter-chromosomal sv
            if((itsv->mPESupport < mOpt->passOpt->mInterChrSVMinDPCnt || itsv->mPEMapQuality < mOpt->passOpt->mInterChrSVMinDPQual) ||
               (itsv->mSRSupport < mOpt->passOpt->mInterChrSVMinSRCnt || itsv->mSRMapQuality < mOpt->passOpt->mInterChrSVMinSRQual)){
                filter = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
            }
        }
        rec->rid = bcf_hdr_name2id(hdr, bamhdr->target_name[itsv->mChr1]); // CHROM
        int32_t svStartPos = std::max(1, itsv->mSVStart - 1);
        int32_t svEndPos = std::max(1, itsv-> mSVEnd);
        if(svEndPos >= (int32_t) bamhdr->target_len[itsv->mChr2]) svEndPos = bamhdr->target_len[itsv->mChr2] - 1;
        rec->pos = svStartPos; // POS
        std::string id = svutil::addID(itsv->mSVT) + std::to_string(itsv->mID);
        bcf_update_id(hdr, rec, id.c_str()); // ID
        bcf_update_alleles_str(hdr, rec, itsv->mAlleles.c_str()); // REF and ALT
        bcf_update_filter(hdr, rec, &filter, 1); // FILTER
        // Add INFO fields
        if(itsv->mPrecise) bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
        else bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);
        bcf_update_info_string(hdr, rec, "SVTYPE", svutil::addID(itsv->mSVT).c_str());
        bcf_update_info_string(hdr, rec, "CHR2", bamhdr->target_name[itsv->mChr2]);
        bcf_update_info_int32(hdr, rec, "END", &itsv->mSVEnd, 1);
        bcf_update_info_int32(hdr, rec, "PE", &itsv->mPESupport, 1);
        tmpi = itsv->mPEMapQuality;
        bcf_update_info_int32(hdr, rec, "PEMAPQ", &tmpi, 1);
        bcf_update_info_string(hdr, rec, "CT", svutil::addOrientation(itsv->mSVT).c_str());
        int32_t tmpai[2];
        tmpai[0] = itsv->mCiPosLow;
        tmpai[1] = itsv->mCiPosHigh;
        bcf_update_info_int32(hdr, rec, "CIPOS", tmpai, 2);
        tmpai[0] = itsv->mCiEndLow;
        tmpai[1] = itsv->mCiEndHigh;
        bcf_update_info_int32(hdr, rec, "CIEND", tmpai, 2);
        // Precise SV specific tags
        if(itsv->mPrecise){
            bcf_update_info_int32(hdr, rec, "SR", &itsv->mSRSupport, 1);
            tmpi = itsv->mSRMapQuality;
            bcf_update_info_int32(hdr, rec, "SRMAPQ", &tmpi, 1);
            bcf_update_info_float(hdr, rec, "SRALNQ", &itsv->mSRAlignQuality, 1);
            bcf_update_info_int32(hdr, rec, "INSLEN", &itsv->mAlnInsLen, 1);
            bcf_update_info_int32(hdr, rec, "HOMLEN", &itsv->mHomLen, 1);
            bcf_update_info_string(hdr, rec, "CONSENSUS", itsv->mConsensus.c_str());
        }
        // Compute GLs
        if(itsv->mPrecise){
            if(mJctCnts[itsv->mID].mRefQualBeg.size() > mJctCnts[itsv->mID].mRefQualEnd.size()){
                svutil::computeGL(mJctCnts[itsv->mID].mRefQualBeg, mJctCnts[itsv->mID].mAltQual, gls, gts, gqval);
            }else{
                svutil::computeGL(mJctCnts[itsv->mID].mRefQualEnd, mJctCnts[itsv->mID].mAltQual, gls, gts, gqval);
            }
        }else{
            if(mSpnCnts[itsv->mID].mRefQualBeg.size() > mSpnCnts[itsv->mID].mRefQualEnd.size()){
                svutil::computeGL(mSpnCnts[itsv->mID].mRefQualBeg, mSpnCnts[itsv->mID].mAltQual, gls, gts, gqval);
            }else{
                svutil::computeGL(mSpnCnts[itsv->mID].mRefQualEnd, mSpnCnts[itsv->mID].mAltQual, gls, gts, gqval);
            }
        }
        bcf_update_genotypes(hdr, rec, gts, bcf_hdr_nsamples(hdr) * 2);
        bcf_update_format_float(hdr, rec, "GL", gls, bcf_hdr_nsamples(hdr) * 3);
        bcf_update_format_int32(hdr, rec, "GQ", gqval, bcf_hdr_nsamples(hdr));
        // Genotype filter
        if(gqval[0] < mOpt->passOpt->mMinGenoQual) ftarr = "LowQual";
        else ftarr = "PASS";
        bcf_update_format_string(hdr, rec, "FT", &ftarr, bcf_hdr_nsamples(hdr));
        // Add read/pair counts
        drcount[0] = mSpnCnts[itsv->mID].getRefDep();
        dvcount[0] = mSpnCnts[itsv->mID].getAltDep();
        bcf_update_format_int32(hdr, rec, "DR", drcount, bcf_hdr_nsamples(hdr));
        bcf_update_format_int32(hdr, rec, "DV", dvcount, bcf_hdr_nsamples(hdr));
        if(mOpt->libInfo->mIsHaploTagged){
            hp1drcount[0] = mSpnCnts[itsv->mID].mRefh1;
            hp2drcount[0] = mSpnCnts[itsv->mID].mRefh2;
            hp1dvcount[0] = mSpnCnts[itsv->mID].mAlth1;
            hp2dvcount[0] = mSpnCnts[itsv->mID].mAlth2;
            bcf_update_format_int32(hdr, rec, "HP1DR", hp1drcount, bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "HP2DR", hp2drcount, bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "HP1DV", hp1dvcount, bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "HP2DV", hp2dvcount, bcf_hdr_nsamples(hdr));
        }
        rrcount[0] = mJctCnts[itsv->mID].getRefDep();
        rvcount[0] = mJctCnts[itsv->mID].getAltDep();
        bcf_update_format_int32(hdr, rec, "RR", rrcount, bcf_hdr_nsamples(hdr));
        bcf_update_format_int32(hdr, rec, "RV", rvcount, bcf_hdr_nsamples(hdr));
        if(mOpt->libInfo->mIsHaploTagged){
            hp1rrcount[0] = mJctCnts[itsv->mID].mRefh1;
            hp2rrcount[0] = mJctCnts[itsv->mID].mRefh2;
            hp1rvcount[0] = mJctCnts[itsv->mID].mAlth1;
            hp2rvcount[0] = mJctCnts[itsv->mID].mAlth2;
            bcf_update_format_int32(hdr, rec, "HP1RR", hp1rrcount, bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "HP2RR", hp2rrcount, bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "HP1RV", hp1rvcount, bcf_hdr_nsamples(hdr));
            bcf_update_format_int32(hdr, rec, "HP2RV", hp2rvcount, bcf_hdr_nsamples(hdr));
        }
        bcf_write(fp, hdr, rec);
        bcf_clear(rec);
    }
    // Close BCF file
    bcf_close(fp);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    // Build BCF file index
    bcf_index_build(mOpt->bcfOut.c_str(), 14);
    // Close BAM file
    sam_close(samfp);
    bam_hdr_destroy(bamhdr);
    // Cleanup allocated memory
    free(gts);
    free(gls);
    free(drcount);
    free(dvcount);
    free(hp1drcount);
    free(hp2drcount);
    free(hp1dvcount);
    free(hp2dvcount);
    free(rrcount);
    free(rvcount);
    free(hp1rrcount);
    free(hp2rrcount);
    free(hp1rvcount);
    free(hp2rvcount);
    free(gqval);
}
