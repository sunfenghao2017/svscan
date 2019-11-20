#include "junction.h"
#include "bamutil.h"

bool JunctionMap::insertJunction(const bam1_t* b, bam_hdr_t* h){
    // parse non-supplementary alignment record first
    bool inserted = false;
    bool fw = !(b->core.flag & BAM_FREVERSE);
    size_t seed = svutil::hashString(bam_get_qname(b));
    int32_t refpos = b->core.pos;
    int32_t seqpos = 0, readpos = 0, seqmatch = 0;
    int32_t readStart = refpos;
    uint32_t* cigar = bam_get_cigar(b);
    int32_t seqlen = bamutil::getSeqLen(b);
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        int opint = bam_cigar_op(cigar[i]);
        int oplen = bam_cigar_oplen(cigar[i]);
        if(opint == BAM_CMATCH || opint == BAM_CEQUAL || opint == BAM_CDIFF){
            refpos += oplen;
            seqpos += oplen;
        }else if(opint == BAM_CDEL){
            readpos = (seqpos <= seqlen && !fw) ? seqlen - seqpos : seqpos;
            if(oplen > mOpt->filterOpt->mMinRefSep){
                mJunctionReads[seed].push_back(Junction(fw, false, 0, b->core.tid, readStart, refpos, readpos, 0));
                refpos += oplen;
                mJunctionReads[seed].push_back(Junction(fw, true, 0, b->core.tid, readStart, refpos, readpos, 0));
                inserted = true;
            }else refpos += oplen;
        }else if(opint == BAM_CINS){
            readpos = (seqpos <= seqlen && !fw) ? seqlen - seqpos : seqpos;
            if(oplen > mOpt->filterOpt->mMinRefSep){
                mJunctionReads[seed].push_back(Junction(fw, false, 0, b->core.tid, readStart, refpos, readpos, 0));
                mJunctionReads[seed].push_back(Junction(fw, true, 0, b->core.tid, readStart, refpos, readpos + oplen, 0));
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
            seqmatch = seqlen - oplen;
            if(oplen > mOpt->filterOpt->minClipLen){
                mJunctionReads[seed].push_back(Junction(fw, scleft, oplen, b->core.tid, readStart, refpos, readpos, seqmatch));
                inserted = true;
            }
        }else if(opint == BAM_CREF_SKIP) refpos += oplen;
    }
    // parse supplenmentary alignment record then
    uint8_t* sa = bam_aux_get(b, "SA");
    if(sa){
        std::string sastr = bam_aux2Z(sa);
        if(sastr.find_first_of(";") != sastr.find_last_of(";")) return inserted;
        std::vector<std::string> vstr;
        util::split(sastr, vstr, ",");
        int32_t tid = bam_name2id(h, vstr[0].c_str());
        refpos = std::atoi(vstr[1].c_str());
        fw = (vstr[2][0] == '+');
        readStart = -1;
        seqpos = 0, readpos = 0, seqmatch = 0;
        std::vector<std::pair<int32_t, char>> pcigar;
        bamutil::parseCigar(vstr[3], pcigar);
        for(auto& e: pcigar){
            char opchr = e.second;
            int oplen = e.first;
            if(opchr == 'M' || opchr == '=' || opchr == 'X'){
                refpos += oplen;
                seqpos += oplen;
            }else if(opchr == 'D'){
                readpos = (seqpos <= seqlen && !fw) ? seqlen - seqpos : seqpos;
                if(oplen > mOpt->filterOpt->mMinRefSep){
                    mJunctionReads[seed].push_back(Junction(fw, false, 0, tid, readStart, refpos, readpos, 0));
                    refpos += oplen;
                    mJunctionReads[seed].push_back(Junction(fw, true, 0, tid, readStart, refpos, readpos, 0));
                    inserted = true;
                }else refpos += oplen;
            }else if(opchr == 'I'){
                readpos = (seqpos <= seqlen && !fw) ? seqlen - seqpos : seqpos;
                if(oplen > mOpt->filterOpt->mMinRefSep){
                    mJunctionReads[seed].push_back(Junction(fw, false, 0, tid, readStart, refpos, readpos, 0));
                    mJunctionReads[seed].push_back(Junction(fw, true, 0, tid, readStart, refpos, readpos + oplen, 0));
                    inserted = true;
                }
                seqpos += oplen;
            }else if(opchr == 'S' || opchr == 'H'){
                int32_t lastSeqPos = seqpos;
                bool scleft = false;
                if(seqpos == 0){
                    lastSeqPos += oplen;
                    scleft = true;
                }
                seqpos += oplen;
                readpos = (lastSeqPos <= seqlen && !fw) ? seqlen - lastSeqPos : lastSeqPos;
                seqmatch = seqlen - oplen;
                if(oplen > mOpt->filterOpt->minClipLen){
                    mJunctionReads[seed].push_back(Junction(fw, scleft, oplen, tid, readStart, refpos, readpos, seqmatch));
                    inserted = true;
                }
            }else if(opchr == BAM_CREF_SKIP) refpos += oplen;
        }
    }
    return inserted;
}

void JunctionMap::sortJunctions(){
    for(auto iter = mJunctionReads.begin(); iter != mJunctionReads.end(); ++iter){
        std::sort(iter->second.begin(), iter->second.end());
    }
    mSorted = true;
}
