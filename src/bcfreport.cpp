#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "stats.h"
#include <util.h>

void Stats::reportSVBCF(const SVSet* svs){
    // Open file handler
    samFile* samfp = sam_open(mOpt->bamfile.c_str(), "r");
    hts_set_fai_filename(samfp, mOpt->alnref.c_str());
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
    bcf_hdr_append(hdr, "##INFO=<ID=CHREND,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SVEND,Number=1,Type=Integer,Description=\"End position of the structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PECNT,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">");
    bcf_hdr_append(hdr, "##INFO=<ID=PEMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRCNT,Number=1,Type=Integer,Description=\"Split-read support\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">");
    bcf_hdr_append(hdr, "##INFO=<ID=SRALNQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CONSENSUSSEQ,Number=1,Type=String,Description=\"Split-read consensus sequence\">");
    bcf_hdr_append(hdr, "##INFO=<ID=CATT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">");
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
    bcf_hdr_append(hdr, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">");
    // Add reference
    std::string refloc = "##reference=" + mOpt->alnref;
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
    std::string sinf = "##softwoare=svscan v" + mOpt->softEnv->ver;
    std::string scmd = "##command=" + mOpt->softEnv->cmd;
    std::string scwd = "##cwd=" + mOpt->softEnv->cwd;
    bcf_hdr_append(hdr, sinf.c_str());
    bcf_hdr_append(hdr, scmd.c_str());
    bcf_hdr_append(hdr, scwd.c_str());
    // Add Samples
    bcf_hdr_add_sample(hdr, "SAMPLE");
    assert(bcf_hdr_write(fp, hdr) >= 0);
    // Add Records
    if(svs->empty()){// if empty set, just return
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
    for(uint32_t i = 0; i < svs->size(); ++i){
        // Prepare Filter field
        int filter = bcf_hdr_id2int(hdr, BCF_DT_ID, "PASS");
        if(svs->at(i)->mChr1 == svs->at(i)->mChr2){// Intra-chromosomal sv
            if((svs->at(i)->mPESupport < mOpt->passOpt->mIntraChrSVMinDPCnt || svs->at(i)->mPEMapQuality < mOpt->passOpt->mIntraChrSVMinDPQual) ||
               (svs->at(i)->mSRSupport < mOpt->passOpt->mIntraChrSVMinSRCnt || svs->at(i)->mSRMapQuality < mOpt->passOpt->mIntraChrSVMinSRQual)){
                filter = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
            }
        }else{// inter-chromosomal sv
            if((svs->at(i)->mPESupport < mOpt->passOpt->mInterChrSVMinDPCnt || svs->at(i)->mPEMapQuality < mOpt->passOpt->mInterChrSVMinDPQual) ||
               (svs->at(i)->mSRSupport < mOpt->passOpt->mInterChrSVMinSRCnt || svs->at(i)->mSRMapQuality < mOpt->passOpt->mInterChrSVMinSRQual)){
                filter = bcf_hdr_id2int(hdr, BCF_DT_ID, "LowQual");
            }
        }
        rec->rid = bcf_hdr_name2id(hdr, bamhdr->target_name[svs->at(i)->mChr1]); // CHROM
        int32_t svStartPos = std::max(1, svs->at(i)->mSVStart - 1);
        int32_t svEndPos = std::max(1, svs->at(i)-> mSVEnd);
        if(svEndPos >= (int32_t) bamhdr->target_len[svs->at(i)->mChr2]) svEndPos = bamhdr->target_len[svs->at(i)->mChr2] - 1;
        rec->pos = svStartPos; // POS
        std::string id = svutil::addID(svs->at(i)->mSVT) + std::to_string(svs->at(i)->mID);
        bcf_update_id(hdr, rec, id.c_str()); // ID
        bcf_update_alleles_str(hdr, rec, svs->at(i)->mAlleles.c_str()); // REF and ALT
        bcf_update_filter(hdr, rec, &filter, 1); // FILTER
        // Add INFO fields
        if(svs->at(i)->mPrecise) bcf_update_info_flag(hdr, rec, "PRECISE", NULL, 1);
        else bcf_update_info_flag(hdr, rec, "IMPRECISE", NULL, 1);
        bcf_update_info_string(hdr, rec, "SVTYPE", svutil::addID(svs->at(i)->mSVT).c_str());
        bcf_update_info_string(hdr, rec, "CHREND", bamhdr->target_name[svs->at(i)->mChr2]);
        bcf_update_info_int32(hdr, rec, "SVEND", &svs->at(i)->mSVEnd, 1);
        bcf_update_info_int32(hdr, rec, "PECNT", &svs->at(i)->mPESupport, 1);
        tmpi = svs->at(i)->mPEMapQuality;
        bcf_update_info_int32(hdr, rec, "PEMAPQ", &tmpi, 1);
        bcf_update_info_string(hdr, rec, "CATT", svutil::addOrientation(svs->at(i)->mSVT).c_str());
        int32_t tmpai[2];
        tmpai[0] = svs->at(i)->mCiPosLow;
        tmpai[1] = svs->at(i)->mCiPosHigh;
        bcf_update_info_int32(hdr, rec, "CIPOS", tmpai, 2);
        tmpai[0] = svs->at(i)->mCiEndLow;
        tmpai[1] = svs->at(i)->mCiEndHigh;
        bcf_update_info_int32(hdr, rec, "CIEND", tmpai, 2);
        // Precise SV specific tags
        if(svs->at(i)->mPrecise){
            bcf_update_info_int32(hdr, rec, "SRCNT", &svs->at(i)->mSRSupport, 1);
            tmpi = svs->at(i)->mSRMapQuality;
            bcf_update_info_int32(hdr, rec, "SRMAPQ", &tmpi, 1);
            bcf_update_info_float(hdr, rec, "SRALNQ", &svs->at(i)->mSRAlignQuality, 1);
            bcf_update_info_int32(hdr, rec, "INSLEN", &svs->at(i)->mAlnInsLen, 1);
            bcf_update_info_int32(hdr, rec, "HOMLEN", &svs->at(i)->mHomLen, 1);
            bcf_update_info_string(hdr, rec, "CONSENSUSSEQ", svs->at(i)->mConsensus.c_str());
        }
        // Compute GLs
        if(svs->at(i)->mPrecise){
            if(mJctCnts[svs->at(i)->mID].mRefQualBeg.size() > mJctCnts[svs->at(i)->mID].mRefQualEnd.size()){
                svutil::computeGL(mJctCnts[svs->at(i)->mID].mRefQualBeg, mJctCnts[svs->at(i)->mID].mAltQual, gls, gts, gqval);
            }else{
                svutil::computeGL(mJctCnts[svs->at(i)->mID].mRefQualEnd, mJctCnts[svs->at(i)->mID].mAltQual, gls, gts, gqval);
            }
        }else{
            if(mSpnCnts[svs->at(i)->mID].mRefQualBeg.size() > mSpnCnts[svs->at(i)->mID].mRefQualEnd.size()){
                svutil::computeGL(mSpnCnts[svs->at(i)->mID].mRefQualBeg, mSpnCnts[svs->at(i)->mID].mAltQual, gls, gts, gqval);
            }else{
                svutil::computeGL(mSpnCnts[svs->at(i)->mID].mRefQualEnd, mSpnCnts[svs->at(i)->mID].mAltQual, gls, gts, gqval);
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
        drcount[0] = mSpnCnts[svs->at(i)->mID].getRefDep();
        dvcount[0] = mSpnCnts[svs->at(i)->mID].getAltDep();
        bcf_update_format_int32(hdr, rec, "DR", drcount, bcf_hdr_nsamples(hdr));
        bcf_update_format_int32(hdr, rec, "DV", dvcount, bcf_hdr_nsamples(hdr));
        rrcount[0] = mJctCnts[svs->at(i)->mID].getRefDep();
        rvcount[0] = mJctCnts[svs->at(i)->mID].getAltDep();
        bcf_update_format_int32(hdr, rec, "RR", rrcount, bcf_hdr_nsamples(hdr));
        bcf_update_format_int32(hdr, rec, "RV", rvcount, bcf_hdr_nsamples(hdr));
        assert(bcf_write(fp, hdr, rec) >= 0);
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
    free(rrcount);
    free(rvcount);
    free(hp1rrcount);
    free(hp2rrcount);
    free(hp1rvcount);
    free(hp2rvcount);
    free(gqval);
}
