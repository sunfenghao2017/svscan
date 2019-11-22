#ifndef ONLINEBWA_H
#define ONLINEBWA_H

#define _set_pac(pac, l, c) ((pac)[(l)>>2] |= (c)<<((~(l)&3)<<1))
#define _get_pac(pac, l) ((pac)[(l)>>2]>>((~(l)&3)<<1)&3)

#include <string>
#include <memory>
#include <cstring>
#include <sstream>
#include <cassert>
#include <zlib.h>
#include <util.h>
#include <bwa/bwa.h>
#include <bwa/bwt.h>
#include <bwa/kseq.h>
#include <bwa/utils.h>
#include <bwa/bwamem.h>
#include <bwa/bntseq.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

extern "C"{
    int is_bwt(ubyte_t *T, int n);
}

/** Structure to hold unaligned sequence (name and bases) */
struct UnalignedSeq{
    std::string mName;    ///< name of the contig
    std::string mComment; ///< comment of the contig
    std::string mSeq;     ///< sequence of the contig (upper-case ACTGN)
    std::string mQual;    ///< quality scores
    char mStrand;         ///< strand of the sequence. Default is '*'

    /** Construct an empty sequence */
    UnalignedSeq(){}

    /** Construct an unaligned sequence with name and sequence
     * @param n name of the sequence
     * @param s sequence, stored as ACTG or N characters
     */
    UnalignedSeq(const std::string& n, const std::string& s){
        mName = n;
        mComment = "";
        mSeq = s;
        mQual = "";
        mStrand = '*';
    }

    /** Output an unaligned sequence to ostream
    * @param os ostream
    * @param us UnalignedSeq
    */
    inline friend std::ostream& operator<<(std::ostream& os, const UnalignedSeq& us){
        os << us.mName;
        if(!us.mComment.empty()) os << " " << us.mComment;
        os << "\n";
        os << us.mSeq << "\n";
        os << us.mStrand << "\n";;
        os << us.mQual << "\n";
        return os;
    }
};


class OnlineBWA{
    public:
        mem_opt_t* mMemOpt; ///< pointer to memory options for bwa alignment
        bool mCopyComment;  ///< copy fasta/q comment to SAM output if true
        bwaidx_t* mIndex;   ///< pointer to bwa index structure
        
        /** OnlineBWA constructor */ 
        OnlineBWA(){
            mIndex = 0;
            mCopyComment = false;
            mMemOpt = mem_opt_init();
            mMemOpt->flag |= MEM_F_SOFTCLIP;
        }
        
        /** OnlineBWA destructor */
        ~OnlineBWA(){
            if(mIndex){
                bwa_idx_destroy(mIndex);
            }
            if(mMemOpt){
                free(mMemOpt);
            }
        }

    public:
        /** convert a bns to string
         * @return pointer to bam_hdr_t
         */
        bam_hdr_t* getBamHeader(){
            if(!mIndex){
                return NULL;
            }
            bntseq_t* bns = mIndex->bns;
            std::stringstream ss;
            for(int i = 0; i < bns->n_seqs; ++i){
                ss << "@SQ\tSN:" << bns->anns[i].name << "\tLN:" << bns->anns[i].len << "\n";
            }
            kstring_t str;
            bam_hdr_t *hhh;
            str.l = str.m = 0; str.s = 0;
            std::istringstream iss(ss.str());
            std::string line;
            while(std::getline(iss, line, '\n')){
                if(line.length() == 0 || line.at(0) != '@'){
                    break;
                }
                kputsn(line.c_str(), line.length(), &str);
                kputc('\n', &str);
            }
            if(str.l == 0){
                kputsn("", 0, &str);
            }
            hhh = sam_hdr_parse(str.l, str.s);
            hhh->l_text = str.l; hhh->text = str.s; // hhh->text needs to be freed
            return hhh;
        }
        
        /** convert pac tos bwt_t
         * @param pac uint8_t array
         * @param bwt_seq_lenr bwt sequence length
         * @return pointer to bwt_t
         */
        bwt_t* pac2bwt(const uint8_t* pac, int bwt_seq_lenr){
            bwt_t *bwt;
            ubyte_t *buf;
            int i;
            // initialization
            bwt = (bwt_t*)calloc(1, sizeof(bwt_t));
            bwt->seq_len = bwt_seq_lenr; //bwa_seq_len(fn_pac); //dummy
            bwt->bwt_size = (bwt->seq_len + 15) >> 4;
            // prepare sequence
            memset(bwt->L2, 0, 5 * 4);
            buf = (ubyte_t*)calloc(bwt->seq_len + 1, 1);
            for(i = 0; i < (int)bwt->seq_len; ++i){
                buf[i] = pac[i>>2] >> ((3 - (i&3)) << 1) & 3;
                ++bwt->L2[1+buf[i]];
            }
            for(i = 2; i <= 4; ++i) bwt->L2[i] += bwt->L2[i-1];
            // Burrows-Wheeler Transform
            bwt->primary = is_bwt(buf, bwt->seq_len);
            bwt->bwt = (u_int32_t*)calloc(bwt->bwt_size, 4);
            for(i = 0; i < (int)bwt->seq_len; ++i) bwt->bwt[i>>4] |= buf[i] << ((15 - (i&15)) << 1);
            free(buf);
            return bwt;
        }
        
        /** add an annotation sequence
         * @param name reference of annotation name
         * @param seq reference of annotation sequence
         * @param ann pointer to bntann1_t struct
         * @param offset offset of bntann1_t to add annotation
         */
        bntann1_t* addAnns(const std::string& name, const std::string& seq, bntann1_t* ann, size_t offset){
            ann->offset = offset;
            ann->name = (char*)malloc(name.length()+1); // +1 for \0
            strncpy(ann->name, name.c_str(), name.length()+1);
            ann->anno = (char*)malloc(7);
            strcpy(ann->anno, "(null)\0");
            ann->len = seq.length();
            ann->n_ambs = 0; // number of "holes"
            ann->gi = 0; // gi?
            ann->is_alt = 0;
            return ann;
        }

        /** bwa mem add1 function */
        uint8_t* addSeq2pac(const kseq_t* seq, bntseq_t* bns, uint8_t* pac, int64_t *m_pac, int *m_seqs, int *m_holes, bntamb1_t **q){
            bntann1_t *p;
            int lasts;
            if(bns->n_seqs == *m_seqs){
                *m_seqs <<= 1;
                bns->anns = (bntann1_t*)realloc(bns->anns, *m_seqs * sizeof(bntann1_t));
            }
            p = bns->anns + bns->n_seqs;
            p->name = strdup((char*)seq->name.s);
            p->anno = seq->comment.l > 0? strdup((char*)seq->comment.s) : strdup("(null)");
            p->gi = 0; p->len = seq->seq.l;
            p->offset = (bns->n_seqs == 0)? 0 : (p-1)->offset + (p-1)->len;
            p->n_ambs = 0;
            for(size_t i = lasts = 0; i < seq->seq.l; ++i){
                int c = nst_nt4_table[(int)seq->seq.s[i]];
                if(c >= 4){ // N
                    if(lasts == seq->seq.s[i]){ // contiguous N
                        ++(*q)->len;
                    }else{
                        if(bns->n_holes == *m_holes){
                            (*m_holes) <<= 1;
                            bns->ambs = (bntamb1_t*)realloc(bns->ambs, (*m_holes) * sizeof(bntamb1_t));
                        }
                        *q = bns->ambs + bns->n_holes;
                        (*q)->len = 1;
                        (*q)->offset = p->offset + i;
                        (*q)->amb = seq->seq.s[i];
                        ++p->n_ambs;
                        ++bns->n_holes;
                    }
                }
                lasts = seq->seq.s[i];
                { // fill buffer
                    if(c >= 4) c = lrand48()&3;
                    if(bns->l_pac == *m_pac){ // double the pac size
                        *m_pac <<= 1;
                        pac = (uint8_t*)realloc(pac, *m_pac/4);
                        memset(pac + bns->l_pac/4, 0, (*m_pac - bns->l_pac)/4);
                    }
                    _set_pac(pac, bns->l_pac, c);
                    ++bns->l_pac;
                }
            }
            ++bns->n_seqs;
            return pac;
        }

        /** make the pac structure for a bunch of reads
         * @param v vector of UnalignedSeq 
         * @param addReverse add reverse sequence to pac if true
         */
        uint8_t* makePac(const std::vector<UnalignedSeq>& v, bool addReverse){
            bntseq_t * bns = (bntseq_t*)calloc(1, sizeof(bntseq_t));
            uint8_t *pac = 0;
            int32_t m_seqs, m_holes;
            int64_t m_pac, l;
            bntamb1_t *q;

            bns->seed = 11; // fixed seed for random generator
            m_seqs = m_holes = 8; m_pac = 0x10000;
            bns->anns = (bntann1_t*)calloc(m_seqs, sizeof(bntann1_t));
            bns->ambs = (bntamb1_t*)calloc(m_holes, sizeof(bntamb1_t));
            pac = (uint8_t*) calloc(m_pac/4, 1);
            q = bns->ambs;
            // move through the unaligned sequences
            for(size_t k = 0; k < v.size(); ++k) {
                // make the ref name kstring
                kstring_t * name = (kstring_t*)malloc(1 * sizeof(kstring_t));
                name->l = v[k].mName.length() + 1;
                name->m = v[k].mName.length() + 3;
                name->s = (char*)calloc(name->m, sizeof(char));
                memcpy(name->s, v[k].mName.c_str(), v[k].mName.length()+1);

                // make the sequence kstring
                kstring_t * t = (kstring_t*)malloc(sizeof(kstring_t));
                t->l = v[k].mSeq.length();
                t->m = v[k].mSeq.length() + 2;
                //t->s = (char*)calloc(v[k].mSeq.length(), sizeof(char));
                t->s = (char*)malloc(t->m);
                memcpy(t->s, v[k].mSeq.c_str(), v[k].mSeq.length());

                // put into a kstring
                kseq_t *ks = (kseq_t*)calloc(1, sizeof(kseq_t));
                ks->seq = *t;
                ks->name = *name;
                // make the forward only pac
                pac = addSeq2pac(ks, bns, pac, &m_pac, &m_seqs, &m_holes, &q);
                // clear it out
                free(name->s);
                free(name);
                free(t->s);
                free(t);
                free(ks);
            }
            if(!addReverse){
                // add the reverse complemented sequence
                m_pac = (bns->l_pac * 2 + 3) / 4 * 4;
                pac = (uint8_t*)realloc(pac, m_pac/4);
                memset(pac + (bns->l_pac+3)/4, 0, (m_pac - (bns->l_pac+3)/4*4) / 4);
                for(l = bns->l_pac - 1; l >= 0; --l, ++bns->l_pac) _set_pac(pac, bns->l_pac, 3-_get_pac(pac, l));
            }
            bns_destroy(bns);
            return pac;
        }

        /** write pac part of the index to file
         * @param file filename to write pac to
         */
        void writePacToFile(const std::string& file) const {
            FILE* fp;
            std::string fname = file + ".pac";
            fp = xopen(fname.c_str(), "wb");
            ubyte_t ct;
            err_fwrite(mIndex->pac, 1, (mIndex->bns->l_pac>>2) + ((mIndex->bns->l_pac&3) == 0? 0 : 1), fp);
            if(mIndex->bns->l_pac % 4 == 0){
                ct = 0;
                err_fwrite(&ct, 1, 1, fp);
            }
            ct = mIndex->bns->l_pac % 4;
            err_fwrite(&ct, 1, 1, fp);
            err_fflush(fp);
            err_fclose(fp);
        }

        /** align a sequence to reference
         * @name sequence name
         * @param seq nucleotide sequence to do align, upper/lower case both accepted
         * @result alignment result vector
         * @hardclip if true, output bamrecord should be hardcliped
         */
        void alignSeq(const std::string& name, const std::string& seq, std::vector<bam1_t*>& result, bool hardclip = false){
            if(!mIndex){
                return;
            }
            result.clear();
            mem_alnreg_v ar = mem_align1(mMemOpt, mIndex->bwt, mIndex->bns, mIndex->pac, seq.length(), seq.data());
            std::vector<mem_aln_t> a;
            for(size_t i = 0; i < ar.n; ++i){
                mem_aln_t a_aln = mem_reg2aln(mMemOpt, mIndex->bns, mIndex->pac, seq.length(), seq.data(), &ar.a[i]);
                a.push_back(a_aln);
            }
            for(size_t i = 0; i < a.size(); ++i){
                bam1_t* b = bam_init1();
                b->core.tid = a[i].rid;
                b->core.pos = a[i].pos;
                b->core.qual = a[i].mapq;
                b->core.flag = a[i].flag;
                b->core.n_cigar = a[i].n_cigar;
                // set dummy mate
                b->core.mtid = -1;
                b->core.mpos = -1;
                b->core.isize = 0;
                if(a[i].is_rev){
                    b->core.flag |= BAM_FREVERSE;
                }
                std::string new_seq = seq;
                if(hardclip){
                    int32_t tstart = 0;
                    int32_t len = 0;
                    for(int c = 0; c < a[i].n_cigar; ++c){
                        if((c == 0 && bam_cigar_op(a[i].cigar[c]) == BAM_CREF_SKIP)){
                            tstart = bam_cigar_oplen(a[i].cigar[c]);
                        }else if(bam_cigar_type(bam_cigar_op(a[i].cigar[c]))&1){
                            len += bam_cigar_oplen(a[i].cigar[c]);
                        }
                    }
                    assert(len > 0);
                    assert(tstart + len <= (int32_t)seq.length());
                    new_seq = seq.substr(tstart, len);
                }
                // allocate all the data
                b->core.l_qname = name.length() + 1;
                b->core.l_qseq = new_seq.length();
                b->l_data = b->core.l_qname + (a[i].n_cigar << 2) + ((b->core.l_qseq + 1) >> 1) + (b->core.l_qseq);
                b->data = (uint8_t*)std::malloc(b->l_data);
                // allocate the qname
                std::memcpy(b->data, name.c_str(), name.length() + 1);
                // allocate the cigar. Reverse if aligned to neg strand, since mem_aln_t stores
                // cigars relative to referemce string oreiatnion, not forward alignment
                std::memcpy(b->data + b->core.l_qname, (uint8_t*)a[i].cigar, a[i].n_cigar << 2);
                // convert N to S or H
                int new_val = hardclip ? BAM_CHARD_CLIP : BAM_CSOFT_CLIP;
                uint32_t* cigar = bam_get_cigar(b);
                for(uint32_t k = 0; k < b->core.n_cigar; ++k){
                    if((cigar[k] & BAM_CIGAR_MASK) == BAM_CREF_SKIP){
                        cigar[k] &= ~BAM_CIGAR_MASK;
                        cigar[k] |= new_val;
                    }
                }
                // allocate the sequence
                uint8_t* mbases = b->data + b->core.l_qname + (b->core.n_cigar << 2);
                std::string fseq = new_seq;
                if(a[i].is_rev) util::reverseComplement(fseq);
                for(uint32_t j = 0; j < fseq.length(); ++j){
                    uint8_t base = 15;
                    switch(fseq[j]){
                        case 'A':
                            base = 1;
                            break;
                        case 'C':
                            base = 2;
                            break;
                        case 'G':
                            base = 4;
                            break;
                        case 'T':
                            base = 8;
                            break;
                        default:
                            base = 15;
                    }
                    mbases[j >> 1] &= ~(0xF << ((~j & 1) << 2));
                    mbases[j >> 1] |= base << ((~j & 1) << 2);
                }
                uint8_t* quals = bam_get_qual(b);
                quals[0] = 0xff;
                size_t arn = ar.n;
                bam_aux_append(b, "NA", 'i', 4, (uint8_t*)(&arn));
                uint32_t irec;
                irec = a[i].NM;
                bam_aux_append(b, "NM", 'i', 4, (uint8_t*)(&irec));
                if(a[i].XA){
                    bam_aux_append(b, "XA", 'Z', std::strlen(a[i].XA), (uint8_t*)a[i].XA);
                }
                irec = a[i].score;
                bam_aux_append(b, "AS", 'i', 4, (uint8_t*)(&irec));
                result.push_back(b);
                free(a[i].cigar);
            }
            free(ar.a);
        }

        /** align a UnalignedSeq to reference
         * @param us reference of UnalignedSeq object
         * @param result alignment result vector
         * @hardclip if true, output bamrecord should be hardcliped
         */
        void alignSeq(const UnalignedSeq& us, std::vector<bam1_t*>& result, bool hardclip = false){
            alignSeq(us.mName, us.mSeq, result, hardclip);
            if(mCopyComment){
                for(uint32_t i = 0; i < result.size(); ++i){
                    bam_aux_append(result[i], "BC", 'Z', us.mComment.length(), (uint8_t*)us.mComment.c_str());
                }
            }
        }

        /** align a kseq_t read to reference
         * @param seq pointer to kseq_t struct
         * @param result alignment result vector
         * @hardclip if true, output bamrecord should be hardcliped
         */
        void alignSeq(const kseq_t* seq, std::vector<bam1_t*>& result, bool hardclip = false){
            alignSeq(seq->name.s, seq->seq.s, result, hardclip);
            if(mCopyComment){
                for(uint32_t i = 0; i < result.size(); ++i){
                    bam_aux_append(result[i], "BC", 'Z', seq->comment.l, (uint8_t*)seq->comment.s);
                }
            }
        }

        /** construct index from an UnalignedSeq
         * @param us reference of UnalignedSeq object
         */
        void constructIndex(const UnalignedSeq& us){
            std::vector<UnalignedSeq> uv = {us};
            constructIndex(uv);
        }

        /** construct index from an squence
         * @param n name of sequence
         * @param s sequence, upper/lower case both accepted
         */
        void constructIndex(const std::string& n, const std::string& s){
            std::vector<UnalignedSeq> uv = {{n, s}};
            constructIndex(uv);
        }

        /** construct index from a list of UnalignedSeq
         * @param uv vector of UnalignedSeq
         */
        void constructIndex(const std::vector<UnalignedSeq>& v){
            if(!v.size()) return;
            // check the integrity of the input data
            for(auto i = v.begin(); i != v.end(); ++i){
                if(i->mName.empty() || i->mSeq.empty()){
                    throw std::invalid_argument("OnlineBWA::constructIndex - Reference sequences must have non-empty name and seq");
                }
            }
            if(mIndex){
                std::cerr << "...clearing old index" << std::endl;
                bwa_idx_destroy(mIndex);
                mIndex = 0;
            }
            // allocate memory for idx
            mIndex = (bwaidx_t*)calloc(1, sizeof(bwaidx_t));;
            // construct the forward-only pac
            uint8_t* fwd_pac = makePac(v, true); //true->addReverse
            // construct the forward-reverse pac ("packed" 2 bit sequence)
            uint8_t* pac = makePac(v, false); // don't write, becasue only used to make BWT
            size_t tlen = 0;
            for(auto i = v.begin(); i != v.end(); ++i) tlen += i->mSeq.length();
            // make the bwt
            bwt_t *bwt;
            bwt = pac2bwt(pac, tlen*2); // *2 for fwd and rev
            bwt_bwtupdate_core(bwt);
            free(pac); // done with fwd-rev pac
            // construct sa from bwt and occ. adds it to bwt struct
            bwt_cal_sa(bwt, 32);
            bwt_gen_cnt_table(bwt);
            // make the bns
            bntseq_t * bns = (bntseq_t*) calloc(1, sizeof(bntseq_t));
            bns->l_pac = tlen;
            bns->n_seqs = v.size();
            bns->seed = 11;
            bns->n_holes = 0;
            // make the anns
            bns->anns = (bntann1_t*)calloc(v.size(), sizeof(bntann1_t));
            size_t offset = 0;
            for(size_t k = 0; k < v.size(); ++k){
                addAnns(v[k].mName, v[k].mSeq, &bns->anns[k], offset);
                offset += v[k].mSeq.length();
            }
            //ambs is "holes", like N bases
            bns->ambs = 0; //(bntamb1_t*)calloc(1, sizeof(bntamb1_t));
            // make the in-memory idx struct
            mIndex->bwt = bwt;
            mIndex->bns = bns;
            mIndex->pac = fwd_pac;
        }

        /** load external bwt index
         * @param file index file
         */
        void loadIndex(const std::string& file){
            bwaidx_t* newIndex = bwa_idx_load(file.c_str(), BWA_IDX_ALL);
            if(!newIndex){
                util::errorExit("error loading index");
            }
            if(newIndex){
                bwa_idx_destroy(mIndex);
            }
            mIndex = newIndex;
        }
        
        /** write index to file
         * @param file output file of index
         */
        void writeIndex(const std::string& file){
            if(!mIndex) return;
            std::string bwt_file = file + ".bwt";
            std::string sa_file = file + ".sa";
            bwt_dump_bwt(bwt_file.c_str(), mIndex->bwt);
            bwt_dump_sa(sa_file.c_str(), mIndex->bwt);
            bns_dump(mIndex->bns, file.c_str());
            writePacToFile(file);
        }

        /** set the gap open penalty
         * @param gapOpenPenalty gap open penalty, default 6
         */
        inline void setGapOpenPenalty(int gapOpenPenalty){
            assert(gapOpenPenalty > 0);
            mMemOpt->o_del = mMemOpt->o_ins = gapOpenPenalty;
        }

        /** set the gap extension penalty
         * @param gapExtPenalty gap extension penalty, default 1
         */
        inline void setGapExtendPenalty(int gapExtPenalty){
            assert(gapExtPenalty > 0);
            mMemOpt->e_del = mMemOpt->e_ins = gapExtPenalty;
        }

        /** set mismatch penalty
         * @param mismatchPenaly mismatch penalty, default 4
         */
        inline void setMismatchPenalty(int mismatchPenaly){
            assert(mismatchPenaly > 0);
            mMemOpt->b = mismatchPenaly;
            bwa_fill_scmat(mMemOpt->a, mMemOpt->b, mMemOpt->mat);
        }

        /** set the reseed trigger
         * @param reseed look for internal seeds inside a seed longer than seedlength * reseed, default 1.5
         */
        inline void setReseedTriger(float reseed){
            assert(reseed > 0);
            mMemOpt->split_factor = reseed;
        }

        /** set SW alignment bandwidth
         * @param width SW alignment bandwidth, default 100
         */
        inline void setBandWidth(int width){
            assert(width > 0);
            mMemOpt->w = width;
        }

        /** set the SW alignment off-diagonal X-dropoff
         * @param dropOff off-diagonal X-dropoff, default 100
         */
        inline void setXDropoff(int dropOff){
            assert(dropOff > 0);
            mMemOpt->zdrop = dropOff;
        }

        /** set the 3' clipping penalty
         * @param p3ClipPenalty penalty for 3'-end clipping, default 5
         */
        inline void set3PrimeClipPenalty(int p3ClipPenalty){
            assert(p3ClipPenalty > 0);
            mMemOpt->pen_clip3 = p3ClipPenalty;
        }

        /** set the 5' clipping penalty
         * @param p5ClipPenalty penalty for 5'-end clipping, default 5
         */
        inline void set5PrimeClipPenalty(int p5ClipPenalty){
            assert(p5ClipPenalty > 0);
            mMemOpt->pen_clip5 = p5ClipPenalty;
        }

       /** set the match score, this should be set first as it will scale penalty options 
        * @param matchScore score for a sequence match, which scales options -TdBOELU unless overridden
        */
        inline void setMatchScore(int matchScore){
            assert(matchScore > 0);
            mMemOpt->b *= matchScore;
            mMemOpt->T *= matchScore;
            mMemOpt->o_del *= matchScore;
            mMemOpt->o_ins *= matchScore;
            mMemOpt->e_del *= matchScore;
            mMemOpt->e_ins *= matchScore;
            mMemOpt->zdrop *= matchScore;
            mMemOpt->pen_clip3 *= matchScore;
            mMemOpt->pen_clip5 *= matchScore;
            mMemOpt->pen_unpaired *= matchScore;
            mMemOpt->a = matchScore;
        }
};

#endif
