#ifndef BAMUTIL_H
#define BAMUTIL_H

#include <vector>
#include <memory>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include <iostream>
#include <algorithm>
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"

/** some usefule functions to operate bam file */
namespace bamutil{
    /** get read name of an alignment record
     * @param b pointer to bam1_t struct
     * @return read name of the alignment
     */ 
    inline std::string getQName(const bam1_t* b){
        return bam_get_qname(b);
    }
    
    /** get barcode sequence of an alignment record\n
     * @param b pointer to bam1_t struct
     * @return OX tag of b if parsed properly else empty string
     */
    inline std::string getBarcode(const bam1_t* b){
        const char tagBC[2] = {'O', 'X'};
        uint8_t* dataBC = bam_aux_get(b, tagBC);
        if(!dataBC){
            return "";
        }
        return bam_aux2Z(dataBC);
    }

    /** get read sequence of an alignment record
     * @param b pointer to bam1_t struct
     * @return read seqence of the alignment
     */
    inline std::string getSeq(const bam1_t* b){
        uint8_t* data = bam_get_seq(b);
        std::string seq(b->core.l_qseq, '\0');
        for(int32_t i = 0; i < b->core.l_qseq; ++i){
            seq[i] = seq_nt16_str[bam_seqi(data, i)];
        }
        return seq;
    }
    
    /** get quality string of an alignment record
     * @param b pointer to bam1_t struct
     * @return quality string of the alignment(phred33 based)
     */
    inline std::string getQual(const bam1_t* b){
        uint8_t* data = bam_get_qual(b);
        std::string qseq(b->core.l_qseq, '\0');
        for(int32_t i = 0; i < b->core.l_qseq; ++i){
            qseq[i] = data[i] + 33;
        }
        return qseq;
    }
    
    /** get cigar string of an alignment record
     * @param b pointer to bam1_t struct
     * @return cigar string
     */
    inline std::string getCigar(const bam1_t* b){
        uint32_t* data = bam_get_cigar(b);
        std::stringstream cigarSeq;
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            cigarSeq << bam_cigar_oplen(data[i]);
            cigarSeq << bam_cigar_opchr(data[i]);
        }
        return cigarSeq.str();
    }
    
    /** output an alignment record to std::cerr */
    inline void dump(const bam1_t* b){
        std::cerr << "R:     " << b->core.tid << ":" << b->core.pos << "\n";
        std::cerr << "M:     " << b->core.mtid << ":" <<  b->core.mpos << "\n";
        std::cerr << "TLEN:  " << b->core.isize << "\n";
        std::cerr << "QName: " << getQName(b) << "\n";
        std::cerr << "Cigar: " << getCigar(b) << "\n";
        std::cerr << "Seq:   " << getSeq(b) << "\n";
        std::cerr << "Qual:  " << getQual(b) << std::endl;
    }

    /** test wheather an alignment record is part of another alignment record
     * @param part pointer to bam1_t struct
     * @param whole pointer to bam1_t struct
     * @param isLeft compare from left to right if true,\n
     */
    inline bool isPartOf(const bam1_t* part, const bam1_t* whole, bool isLeft){
        if(!part || !whole || part->core.tid != whole->core.tid){
            return false;
        }
        if(whole->core.n_cigar < part->core.n_cigar){
            return false;
        }
        uint32_t* cigarDataPart = bam_get_cigar(part);
        uint32_t* cigarDataWhole = bam_get_cigar(whole);
        for(uint32_t i = 0; i < part->core.n_cigar; ++i){
            uint32_t valPart = isLeft ? cigarDataPart[i] : cigarDataPart[part->core.n_cigar - 1 - i];
            char opchrPart = bam_cigar_opchr(valPart);
            uint32_t oplenPart = bam_cigar_oplen(valPart);
            uint32_t valWhole = isLeft ? cigarDataWhole[i] : cigarDataWhole[whole->core.n_cigar - 1 - i];
            char opchrWhole = bam_cigar_opchr(valWhole);
            uint32_t oplenWhole = bam_cigar_oplen(valWhole);
            if(opchrPart != opchrWhole || oplenPart > oplenWhole){
                return false;
            }
            if(oplenPart < oplenWhole){
                if(i == part->core.n_cigar - 1){
                    return true;
                }else if(i == part->core.n_cigar - 2){
                    int opPartNext = isLeft ? bam_cigar_op(cigarDataPart[i + 1]) : bam_cigar_op(cigarDataPart[part->core.n_cigar - 2 - i]);
                    if(opPartNext != BAM_CHARD_CLIP){
                        return false;
                    }
                }else{
                    return false;
                }
            }
        }
        return true;
    }
    
    /** return length of reference consumed in range [b->core.pos to b->core.pos + offset)
     * @param b pointer to bam1_t struct
     * @param offset offset from read mapping starting position of b
     * @return length of reference consumed in range [b->core.pos to b->core.pos + offset)\n
     * or -1 if offset is invalid or cigar is invalid(unmapped)
     */
    inline int getRefOffset(const bam1_t* b, int offset){
        uint32_t* data = bam_get_cigar(b);
        int r = 0; // reference consumed length
        int q = 0; // query consumed length
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            uint32_t oplen = bam_cigar_oplen(data[i]);
            int opmask = bam_cigar_op(data[i]);
            switch(bam_cigar_type(opmask)){
                case 1:
                    q += oplen;
                    break;
                case 2:
                    r += oplen;
                    break;
                case 3:
                    q += oplen;
                    r += oplen;
                    break;
                default:
                    break;
            }
            if(q >= offset){
                if(opmask == BAM_CINS || opmask == BAM_CSOFT_CLIP){
                    return -1;
                }else{
                    return r - (q - offset);
                }
            }
        }
        return -1;
    }

    /** return length of reference consumed in alignment record
     * @param b pointer to bam1_t struct
     * @return length of reference consumed in alignmen record
     */
    inline int getRefLen(const bam1_t* b){
        uint32_t* data = bam_get_cigar(b);
        int r = 0; // reference consumed length
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            uint32_t oplen = bam_cigar_oplen(data[i]);
            int opmask = bam_cigar_op(data[i]);
            switch(bam_cigar_type(opmask)){
                case 2:
                    r += oplen;
                    break;
                case 3:
                    r += oplen;
                    break;
                default:
                    break;
            }
        }
        return r;
    }

    /** return total length of query seq(include hard clips)
     * @param b pointer to bam1_t struct
     * @return length of query seq
     */
    inline int getSeqLen(const bam1_t* b){
        uint32_t* data = bam_get_cigar(b);
        int q = 0; // query consumed length
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            int opint = bam_cigar_op(data[i]);
            if(opint == BAM_CINS || opint == BAM_CDIFF || opint == BAM_CEQUAL || opint == BAM_CMATCH || 
               opint == BAM_CHARD_CLIP || opint == BAM_CSOFT_CLIP) q += bam_cigar_oplen(data[i]);
        }
        return q;
    }

    /** set read name of one alignment record 
     * @param b pointer to bam1_t struct
     * @param n new name to be used in this read
     */
    inline void setQName(bam1_t* b, const std::string qname){
        int nonQnameLen = b->l_data - b->core.l_qname;
        uint8_t* nonQnameData = (uint8_t*)calloc(nonQnameLen, sizeof(uint8_t));
        std::memcpy(nonQnameData, b->data + b->core.l_qname, nonQnameLen);
        std::free(b->data);
        b->data = (uint8_t*)calloc(nonQnameLen + qname.length() + 1, sizeof(uint8_t));
        std::memcpy(b->data, (uint8_t*)qname.c_str(), qname.length());
        b->data[qname.length()] = '\0';
        b->l_data = b->l_data - b->core.l_qname + qname.length() + 1;
        b->core.l_qname = qname.length() + 1;
        std::memcpy(b->data + b->core.l_qname, nonQnameData, nonQnameLen);
        std::free(nonQnameData);
        b->m_data = b->l_data;
        b->core.l_extranul = 0;
    }

    /** set read name of one alignment record to another
     * @param from pointer to bam1_t struct
     * @param to pointer to bam1_t struct
     */
    inline void setQName(const bam1_t* from, bam1_t* to){
        std::string qname = bam_get_qname(from);
        setQName(to, qname);
    }

    /** modify the ith base of one alignment record
     * @param b pointer to bam1_t struct
     * @param i position at which to change nucleotide(0 based)
     * @param ch nucleotide to change to
     */
    inline void setBase(bam1_t* b, int i, char ch){
        if(i < 0 || i > b->core.l_qseq - 1){
            std::cerr << "Invalid position to set base!" << std::endl;
            std::exit(1);
        }
        uint8_t* seq = bam_get_seq(b);
        if(i % 2 == 1){
            seq[i/2] = (seq[i/2] & 0xF0) | seq_nt16_table[ch - 0];
        }else{
            seq[i/2] = (seq[i/2] & 0x0F) | (seq_nt16_table[ch - 0] << 4);
        }
    }
    
    /** get the first matched status
     * @param b pointer to bam1_t struct
     * @param cigarMOffset length of sequence consumed in read befor the first cigar M
     * @param cigarMOpLen the first cigar M consumed length
     */
    inline void getFirstCigarM(const bam1_t* b, int& cigarMOffset, int& cigarMOpLen){
        uint32_t* data = bam_get_cigar(b);
        int len = 0;
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            if(bam_cigar_op(data[i]) == BAM_CMATCH){
                cigarMOffset = len;
                cigarMOpLen = bam_cigar_oplen(data[i]);
                return;
            }
            switch(bam_cigar_type(bam_cigar_op(data[i]))){
                case 1: case 3:
                    len += bam_cigar_oplen(data[i]);
                    break;
                default:
                    break;
            }
        }
        cigarMOffset = -1;
        cigarMOpLen = 0;
    }

    /** get the first softclip status
     * @param b pointer to bam1_t struct
     * @param cigarSOffset length of sequence consumed in read before the first cigar S
     * @param cigarSOpLen the first cigar S consumed length
     */
    inline void getFirstCigarS(const bam1_t* b, int& cigarSOffset, int& cigarSOpLen){
        uint32_t* data = bam_get_cigar(b);
        int len = 0;
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            if(bam_cigar_op(data[i]) == BAM_CSOFT_CLIP){
                cigarSOffset = len;
                cigarSOpLen = bam_cigar_oplen(data[i]);
                return;
            }
            switch(bam_cigar_type(bam_cigar_op(data[i]))){
                case 1: case 3:
                    len += bam_cigar_oplen(data[i]);
                    break;
                default:
                    break;
            }
        }
        cigarSOffset = -1;
        cigarSOpLen = 0;
    }

    /** get edit distance of this alignment record
     * @param b pointer to bam1_t struct
     * @return edit distance of this alignment record if NM tag parsed properly, else 0
     */
    inline int getEditDistance(const bam1_t* b){
        uint8_t* data = bam_aux_get(b, "NM");
        if(!data){
            return 0;
        }
        return bam_aux2i(data);
    }

    /** get quality sum of a read
     * @param b pointer to bam1_t struct
     * @return quality sum of read in b
     */
    inline int32_t getQualSum(const bam1_t* b){
       int32_t qsum = 0;
       uint8_t* qseq = bam_get_qual(b);
       for(int32_t i = 0; i < b->core.l_qseq; ++i){
           qsum += qseq[i];
       }
       return qsum;
    }

    /** get effective bases
     * @param b pointer to bam1_t struct
     * @return effective bases number
     */
    inline int32_t getEffecBases(const bam1_t* b){
        int32_t effBase = 0;
        uint32_t* data = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            switch(bam_cigar_op(data[i])){
                case BAM_CMATCH: case BAM_CINS: case BAM_CEQUAL: case BAM_CDIFF:
                    effBase += bam_cigar_oplen(data[i]);
                    break;
                default:
                    break;
            }
        }
        return effBase;
    }
    
    /** get mismatches position on reference(refered bam_md.c in samtools)
     * @param h pointer to bam_hdr_t
     * @param b pointer to bam1_t struct
     * @param f reference file
     * @return mismatches positions on reference(0 based)
     */
    inline std::vector<int32_t> getMismatchPos(const bam_hdr_t* h, const bam1_t* b, const faidx_t* f){
        int32_t len = 0;
        char* ref = faidx_fetch_seq(f, h->target_name[b->core.tid], b->core.pos, bam_endpos(b) - 1, &len);
        std::vector<int32_t> ret;
        int rl = 0; // reference consumed length
        int ql = 0; // query consumed length
        uint32_t* cigar = bam_get_cigar(b);
        uint8_t* seq = bam_get_seq(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            int oplen = bam_cigar_oplen(cigar[i]);
            int opmask = bam_cigar_op(cigar[i]);
            if(opmask == BAM_CMATCH || opmask == BAM_CEQUAL || opmask == BAM_CDIFF){
                for(int j = 0; j < oplen; ++j){
                    int32_t qb = bam_seqi(seq, ql + j);
                    int32_t rb = seq_nt16_table[ref[rl + j] - 0];
                    if((qb == rb && qb != 15 && rb != 15) || qb == 0){
                        continue;
                    }else{
                        ret.push_back(b->core.pos + rl + j);
                    }
                }
                rl += oplen;
                ql += oplen;
            }else if(opmask == BAM_CDEL){
                for(int j = 0; j < oplen; ++j){
                    ret.push_back(b->core.pos + rl + j);
                }
                rl += oplen;
            }else if(opmask == BAM_CINS || opmask == BAM_CSOFT_CLIP){
                ql += oplen;
            }else if(opmask == BAM_CREF_SKIP){
                rl += oplen;
            }
        }
        free(ref);
        return ret;
    }

    /** get indel position on read of an alignment record
     * @param p pointer to bam1_t struct
     * @param indelPosVec list of indel positions on read, 0 based
     */
    inline void getIndelPos(const bam1_t* b, std::vector<int>& indelPosVec){
        indelPosVec.clear();
        uint32_t offset = 0;
        uint32_t* data = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            char op = bam_cigar_opchr(data[i]);
            uint32_t len = bam_cigar_oplen(data[i]);
            // skip cigar which does not consume query except BAM_CDEL
            if(op == BAM_CPAD || op == BAM_CHARD_CLIP || op == BAM_CREF_SKIP){
                continue;
            }
            // increase offset only if cigar consume query except BAM_CINS
            if(op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF || op == BAM_CSOFT_CLIP){
                offset += len;
                continue;
            }
            // record BAM_CINS
            if(op == BAM_CINS){
                indelPosVec.push_back(offset);
                offset += len;
                indelPosVec.push_back(offset);
                continue;
            }
            // record BAM_CDEL
            if(op == BAM_CDEL){
                indelPosVec.push_back(offset);
                continue;
            }
        }
    }

    /** check if indels are around the current position
     * @param pos position on read, 0-based
     * @param range range to be checked around pos
     * @param ivec indel position on read of a alignment record
     * @return true if a indel occurs around pos±range inclusive[]
     */
    inline bool posIsNearIndel(int pos, int range, const std::vector<int>& ivec){
        for(uint32_t i = 0; i < ivec.size(); ++i){
            if(std::abs(pos - ivec[i]) <= range){
                return true;
            }
        }
        return false;
    }

    /** get compacted mismatches of a read against the reference in an alignment record
     * @param b pointer to bam1_t struct
     * @return compacted mismatches of b(each indel is counted as one mismatch, ambiguous N mismatch discounted)
     */
    inline int getCompactMismatch(const bam1_t* b){
        int inDelLen = 0;
        int inDelNum = 0;
        int nBaseNum = 0;
        int mismatchNum = 0;
        std::string seq = getSeq(b);
        // ambiguous N base would not be counted as mismatch
        for(uint32_t i = 0; i < seq.length(); ++i){
            if(seq[i] == 'N'){
                ++nBaseNum;
                --mismatchNum;
            }
        }
        uint32_t* cigar = bam_get_cigar(b);
        bool cigarMGot = false;
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            int32_t op = bam_cigar_op(cigar[i]);
            int32_t ol = bam_cigar_oplen(cigar[i]);
            switch(op){
                case BAM_CMATCH:
                    cigarMGot = true;
                    break;
                case BAM_CDIFF:
                    ++mismatchNum;
                    break;
                case BAM_CINS: case BAM_CDEL:
                    ++mismatchNum; // each indel counted as one mismatch
                    inDelLen += ol;
                    ++inDelNum;
                    break;
                default:
                    break;
            }
        }
        // adjust cigar M
        if(cigarMGot){
            uint8_t* nmTag = (uint8_t*)bam_aux_get(b, "NM");
            if(nmTag){
                mismatchNum = bam_aux2i(nmTag) - inDelLen + inDelNum - nBaseNum;
            }
        }
        return mismatchNum < 0 ? 0 : mismatchNum;
    }

    /** get max read length in a bam by iterate at most the first 10000 record
     * @param bamFile bam filename
     * @return max read length scanned
     */
    inline int getBamReadLength(const std::string& bamFile){
        samFile* fp = sam_open(bamFile.c_str(), "r");
        bam_hdr_t* h = sam_hdr_read(fp);
        bam1_t* b = bam_init1();
        int readLen = 0;
        int count = 0;
        while(count < 10000 && sam_read1(fp, h, b) >= 0){
            if(b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL)){
                continue;
            }
            if(b->core.n_cigar == 1){
                readLen = std::max(b->core.l_qseq, readLen);
                ++count;
            }
        }
        bam_destroy1(b);
        bam_hdr_destroy(h);
        sam_close(fp);
        return readLen;
    }

    /** check if a alignment record in a bam_pileup1_t pileup near some range of read 3' end
     * @param p pointer to bam_pileup1_t struct
     * @param range range to be checked near 3'end of read
     * @return true if pipeup position is in read 3'end-range, inclusive[]
     */
    inline bool pileupIsInReadRear(const bam_pileup1_t* p, int range){
        if(bam_is_rev(p->b)){
            if(p->qpos <= range){
                return true;
            }
        }else{
            if(p->qpos >= p->b->core.l_qseq - range){
                return true;
            }
        }
        return false;
    }

    /** get average base quality around predefined range of of pileup posision
     * @pram p pointer to bam_pileup1_t 
     * @param range max number of bases to counted around p->qpos
     * @return average base quality around qpos±range
     */
    inline float getAverageBaseQualityAroundPileupPos(const bam_pileup1_t* p, int range){
        uint8_t* qualStr = bam_get_qual(p->b);
        float totalQual = 0.0;
        int i = 0;
        for(i = std::max(0, p->qpos - range); i <= std::min(p->b->core.l_qseq, p->qpos + range); ++i){
            totalQual += qualStr[i];
        }
        return totalQual/i;
    }

    /** get at most insertion sequence at pileup position of a bam record
     * @param p pointer to bam_pileup1_t 
     * @param maxLen max insertion sequence to extract
     * @return insertion sequence extracted
     */
    inline std::string getInsertionSeq(const bam_pileup1_t* p, int maxLen){
        int insLenToGet = std::min(p->indel, maxLen);
        uint8_t* seq = bam_get_seq(p->b);
        std::string ret(insLenToGet, '\0');
        for(int i = 0; i < insLenToGet; ++i){
            ret[i] = seq_nt16_str[bam_seqi(seq, p->qpos + i)];
        }
        return ret;
    }

    /** get softclip length of left/right softclip
     * @param b pointer to bam1_t
     * @return <leftClipLength, rightClipLength>
     */
    inline std::pair<int, int> getSoftClipLength(const bam1_t* b){
        std::pair<int, int> ret = {0, 0};
        uint32_t* cigarStr = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            if(bam_cigar_op(cigarStr[i]) == BAM_CSOFT_CLIP){
                if(i == 0){
                    ret.first = bam_cigar_oplen(cigarStr[i]);
                }else{
                    ret.second = bam_cigar_oplen(cigarStr[i]);
                }
            }
        }
        return ret;
    }

    /** get minimal distance from pileup position to the end of read
     * @param p pointer to bam_pileup1_t
     * @return minimal distance from p->qpos to ends of read
     */
    inline int getMinDistToReadEnd(const bam_pileup1_t* p){
        std::pair<int, int> clipLen = getSoftClipLength(p->b);
        return std::min(p->qpos + 1 - clipLen.first, p->b->core.l_qseq - clipLen.second - p->qpos - 1);
    }

    /** check a bam is alignment result of pairend reads or not
     * @param bamFile bam file name
     * @return true if bam is alignment result of pairend reads
     */
    inline bool bamIsPE(const char* bamFile){
        samFile* fp = sam_open(bamFile, "r");
        bam_hdr_t* h = sam_hdr_read(fp);
        bam1_t* b = bam_init1();
        bool isPE = false;
        assert(sam_read1(fp, h, b) >= 0);
        if(b->core.flag & BAM_FPAIRED){
            isPE = true;
        }
        bam_hdr_destroy(h);
        bam_destroy1(b);
        sam_close(fp);
        return isPE;
    }

    /** check a bam alignment have soft clipped seq
     * @param b pointer to bam1_t
     * @return true if BAM_CSOFT_CLIP is in cigar
     */
    inline bool haveSoftClip(const bam1_t* b){
        uint32_t* cigar = bam_get_cigar(b);
        for(uint32_t i = 0; i < b->core.n_cigar; ++i){
            if(bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP){
                return true;
            }
        }
        return false;
    }

    /** parse a cigar string to opchr, oplen pair list
     * @param cigar cigar string
     * @param ret result vector
     */
    inline void parseCigar(const std::string& cigar, std::vector<std::pair<int32_t, char>>& ret){
        std::string::size_type cpos = 0, lpos = 0;
        while((cpos = cigar.find_first_of(BAM_CIGAR_STR, lpos)) != std::string::npos){
            ret.push_back({std::atoi(cigar.substr(lpos, cpos - lpos).c_str()), cigar[cpos]});
            lpos = cpos + 1;
        }
    }
}

#endif
