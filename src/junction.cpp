#include "junction.h"

bool JunctionMap::insertJunction(const bam1_t* b){
    bool inserted = false;
    bool fw = !(b->core.flag & BAM_FREVERSE);
    size_t seed = svutil::hashString(bam_get_qname(b));
    int32_t refpos = b->core.pos;
    int32_t seqpos = 0, readpos = 0;
    int32_t readStart = refpos;
    if(b->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) readStart = -1;
    uint32_t* cigar = bam_get_cigar(b);
    int32_t seqlen = bam_cigar2qlen(b->core.n_cigar, cigar);
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        int opint = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        if(opint == BAM_CMATCH || opint == BAM_CEQUAL || opint == BAM_CDIFF){
            refpos += oplen;
            seqpos += oplen;
        }else if(opint == BAM_CDEL){
            readpos = (seqpos <= seqlen && !fw) ? seqlen - seqpos : seqpos;
            if(oplen > mOpt->filterOpt->mMinRefSep){
                mJunctionReads[seed].push_back(Junction(fw, false, b->core.tid, readStart, refpos, readpos));
                refpos += oplen;
                mJunctionReads[seed].push_back(Junction(fw, true, b->core.tid, readStart, refpos, readpos));
                inserted = true;
            }else refpos += oplen;
        }else if(opint == BAM_CINS){
            readpos = (seqpos <= seqlen && !fw) ? seqlen - seqpos : seqpos;
            if(oplen > mOpt->filterOpt->mMinRefSep){
                mJunctionReads[seed].push_back(Junction(fw, false, b->core.tid, readStart, refpos, readpos));
                mJunctionReads[seed].push_back(Junction(fw, true, b->core.tid, readStart, refpos, readpos + oplen));
                inserted = true;
            }
            seqpos += oplen;
        }else if(opint == BAM_CSOFT_CLIP || opint == BAM_CHARD_CLIP){
            int32_t lastSeqPos = seqpos;
            bool scleft = false;
            if(seqpos == 0){
                lastSeqPos += oplen;
                scleft = true;
            }
            seqpos += oplen;
            readpos = (lastSeqPos <= seqlen && !fw) ? seqlen - lastSeqPos : lastSeqPos;
            if(oplen > mOpt->filterOpt->minClipLen){
                mJunctionReads[seed].push_back(Junction(fw, scleft, b->core.tid, readStart, refpos, readpos));
                inserted = true;
            }
        }else if(opint == BAM_CREF_SKIP) refpos += oplen;
    }
    return inserted;
}

void JunctionMap::sortJunctions(){
    for(auto iter = mJunctionReads.begin(); iter != mJunctionReads.end(); ++iter){
        std::sort(iter->second.begin(), iter->second.end());
    }
    mSorted = true;
}
