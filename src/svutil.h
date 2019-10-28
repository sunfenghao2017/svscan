#ifndef SVUTIL_H
#define SVUTIL_H

#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cstdint>
#include "util.h"
#include <htslib/vcf.h>
#include <htslib/sam.h>

/* inversion
 * reference genome, two strands
 * 5'--------------aC--------------------Ag---------3' ref 5'->3'
 * 3'--------------tG--------------------Tc---------5' ref 3'->5'
 * sample genome, two strands
 * 5'--------------aT--------------------Gg---------3' smp 5'->3'
 * 3'--------------tA--------------------Cc---------5' smp 3'->5'
 * 
 * PE reads with orientation FR to smp 5'->3' strand 
 *     ----> ReadF
 * 5'--------------aT--------------------Gg---------3' smp 5'->3'
 * 3'--------------tA--------------------Cc---------5' smp 3'->5'
 *                       <---- ReadR
 * Mapping status on reference genome(3'->3' connection type, 3to3)
 *     ----> ReadF       ----> ReadR
 * 5'--------------aC--------------------Ag---------3' ref 5'->3'
 *
 * PE reads with orientation to smp 5'->3' strand
 *                       ----> ReadF
 * 5'--------------aT--------------------Gg---------3' smp 5'->3'
 * 3'--------------tA--------------------Cc---------5' smp 3'->5'
 *                                           <---- ReadR
 * Mapping status on reference genome(5'->5' connection type, 5to5)
 * 5'--------------aC--------------------Ag---------3' ref 5'->3'
 *                       <---- ReadF'        ----> ReadR'
 *
 * SR read
 *              |------>SRa          |--------->SRb
 * 5'--------------aT--------------------Gg---------3' smp 5'->3'
 *
 * Mapping status
 *              |----~~~SRa          ~~~~~----->SRb
 * 5'--------------aC--------------------Ag---------3' ref 5'->3'
 *             ~~~~~----|SRb'         <---~~~~SRa'
 */

/* deletion
 * reference genome, two strands
 * 5'--------------aC------------Ag-----------------3' ref 5'->3'
 * 3'--------------tG------------Tc-----------------5' ref 3'->5'
 * sample genome, two strands
 * 5'--------------ag---------------3' smp 5'->3'
 * 3'--------------tc---------------5' smp 3'->5'
 * 
 * PE reads with orientation to smp 5'->3' strand
 *       ----->ReadF
 * 5'--------------ag---------------3' smp 5'->3'
 * 3'--------------tc---------------5' smp 3'->5'
 *                     <------ReadR
 * 
 * Mapping status on reference
 *      ------>ReadF 
 * 5'--------------aC------------Ag-----------------3' ref 5'->3'
 *                                  ------->ReadR'
 * SR read with orientation to smp 5'->3' strand
 *             ---------->ReadF
 * 5'--------------ag---------------3' smp 5'->3'
 * 3'--------------tc---------------5' smp 3'->5'
 *             <----------ReadR
 * Mapping status on reference
 *        ReadF-----~~~~~~> ~~~~~~----->ReadF  
 * 5'--------------aC------------Ag-----------------3' ref 5'->3'
 *             ---~~~~~~~>ReadR'
 *                         ~~~~~~-------->ReadR'
 */

/* duplication
 * reference genome, two strands
 * 5'--------------aC------------Agg-----------------3' ref 5'->3'
 * 3'--------------tG------------Tcc-----------------5' ref 3'->5'
 * sample genome, two strands
 * 5'--------------aC------------AgaC------------Agg----------------3' smp 5'->3'
 * 3'--------------tG------------TctG------------Tcc----------------5' smp 3'->5'
 * 
 * PE reads with orientation to smp 5'->3' strand
 *                      ----->ReadF
 * 5'--------------aC------------AgaC------------Agg----------------3' smp 5'->3'
 * 3'--------------tG------------TctG------------Tcc----------------5' smp 3'->5'
 *                                    <------ReadR
 * Mapping status on reference
 *                      ----->ReadF
 * 5'--------------aC------------Agg-----------------3' ref 5'->3'
 *                    ----->ReadR'
 * 
 * SR read with orientation to smp 5'->3' strand
 *                           --------->ReadF
 * 5'--------------aC------------AgaC------------Agg----------------3' smp 5'->3'
 * 3'--------------tG------------TctG------------Tcc----------------5' smp 3'->5'
 *                            <-----------ReadR
 * 
 * Mapping status on reference
 *                           -----~~~~>ReadF
 * 5'--------------aC------------Agg-----------------3' ref 5'->3'
 *             ~~~~~----->ReadF
 *                             ---~~~~~~~~>ReadR'
 *             ~~~~---------->ReadR'
 */

/* insertion
 * reference genome, two strands
 * 5'-----------------aC-------------------3' ref 5'->3'
 * 3'-----------------tG-------------------5' ref 3'->5'
 * sample genome, two strands
 * 5'-----------------a--------C-----------------3' smp 5'->3'
 * 3'-----------------t--------G-----------------5' smp 3'->5'
 *
 * PE Reads with orientation to smp 5'->3' strand
 *            ------>ReadF
 * 5'-----------------a--------C-----------------3' smp 5'->3'
 * 3'-----------------t--------G-----------------5' smp 3'->5'
 *                               <------ReadR
 * Mapping status on reference
 *            ------>ReadF
 * 5'-----------------aC-----------------3' smp 5'->3'
 * 3'-----------------tG-----------------5' smp 3'->5'
 *                       <------ReadR
 *
 * SR Read with orientation to smp 5'->3' strand
 *                 --------------->Read
 * 5'-----------------a--------C-----------------3' smp 5'->3'
 * 3'-----------------t--------G-----------------5' smp 3'->5'
 * Mapping status on reference(hard clip)
 *                 ----~~~~~~~~>ReadF'  
 * 5'-----------------aC-------------------3' ref 5'->3'
 * 3'-----------------tG-------------------5' ref 3'->5'
 *              ~~~~~~~--->ReadR'
 */

/* translocation
 * chrBig
 * 5'-------------------------AC----------------------------3'
 * 3'-------------------------TG----------------------------5'
 * chrSmall
 * 5'-------------------------ac----------------------------3'
 * 3'-------------------------tg----------------------------5'
 * 
 * 5to5 catenation
 *          chrBig 5' part          chrSmall 5' part
 * 5'-------------------------At----------------------------3'
 * 3'-------------------------Ta----------------------------5'
 * 
 * 3to3 catenation
 *          chrBig 3' part          chrSmall 3' part
 * 5'-------------------------Gc----------------------------3'
 * 3'-------------------------Cg----------------------------5'
 *
 * 5to3 catenation
 *         chrBig 5' part           chrSmall 3' part
 * 5'-------------------------Ac----------------------------3'
 * 3'-------------------------Tg----------------------------5'
 *
 * 3to5 catenation
 *        chrSmall 5' part          chrBig 3' part
 * 5'-------------------------aC----------------------------3'
 * 3'-------------------------tG----------------------------5'
 */

#define FUSION_FALLGENE         0x1     ///< All partners are genes
#define FUSION_FNORMALCATDIRECT 0x2     ///< hgene 5' -> tgene 3' positive strand catenation
#define FUSION_FHOTGENE         0x4     ///< one of the partner is in hot fusion gene list
#define FUSION_FCOMMONHOTDIRECT 0x8     ///< hot gene partner is in the predefined direction in fusion gene list
#define FUSION_FINDB            0x10    ///< fusion gene is in database
#define FUSION_FMIRROR          0x20    ///< fusion gene has mirror fusion
#define FUSION_FBLACKPAIR       0x40    ///< fusion gene is in black list
#define FUSION_FFBG             0x80    ///< fusion gene is in background samples
#define FUSION_FBLACKGENE       0x100   ///< one partner is in black gene list
#define FUSION_FLOWCOMPLEX      0x200   ///< one partner has low complex concensus sequence
#define FUSION_FPRIMARY         0x400   ///< this fusion event should be reported primarily
#define FUSION_FSUPPLEMENTARY   0x800   ///< this fusion event should be reported as supplementary
#define FUSION_FTOOSMALLSIZE    0x1000  ///< this fusion event has too small size
#define FUSION_FINSAMEGENE      0x2000  ///< this fusion event occurs in the same gene
#define FUSION_FLOWAF           0x4000  ///< this fusion event has too low AF
#define FUSION_FLOWSUPPORT      0x8000  ///< this fusion event has too low support
#define FUSION_FLOWDEPTH        0x10000 ///< this fusion event has too low depth around breakpoint
#define FUSION_FPRECISE         0x20000 ///< this fusion event has precise breakpoint and concensus sequence
#define FUSION_FCALLFROMRNASEQ  0x40000 ///< this fusion event is called from rna seq
#define FUSION_FMIRRORINDB      0x80000 ///< this fusion event's mirror event is in public db

/** fusion flag type */
typedef uint32_t TFUSION_FLAG;

/** fuse gene struct */
struct FuseGene{
    std::string hgene;   ///< gene name of 5' part of fusion gene
    std::string hend;    ///< hgene 3' or 5' fused
    std::string hstrand; ///< hgene + or - strand fused
    std::string tgene;   ///< gene name of 3' part of fusion gene
    std::string tend;    ///< tgene 3' or 5' fused
    std::string tstrand; ///< hgene + or - strand fused
    int32_t hidx;        ///< hgene index
    bool hfrom1;         ///< hgene is from breakpoint1
    bool tfrom1;         ///< tgene is from breakpoint1
    int32_t tidx;        ///< tgene index
    TFUSION_FLAG status; ///< mask to show fusion status 1:gene,2:normal,4:hot,8:common,16:indb,32:mirror

    /** FuseGene constructor */
    FuseGene(){
        hgene = "-";
        hend = ".";
        hstrand = "-";
        tgene = "-";
        tend = "-";
        tstrand = ".";
        status = 0;
        hidx = -1;
        tidx = -1;
        hfrom1 = false;
        tfrom1 = false;
    }

    /** FuseGene destructor */
    ~FuseGene(){};

    /** get string representation of FuseGene
     * @return string representation of FuseGene
     */
    inline std::string toStr(){
        std::stringstream ss;
        ss << hgene << "->" << tgene << "(";
        ss << hgene << "," << hend << "," << hstrand << ";";
        ss << tgene << "," << tend << "," << tstrand << ")";
        return ss.str();
    }

    /** operator to output FuseGene to ostream
     * @param os reference of ostream
     * @param fg reference of FuseGene
     * @return ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const FuseGene& fg){
        os << fg.hgene << "->" << fg.tgene << "\t";
        os << fg.hgene << "\t" << fg.hend << "\t" << fg.hstrand << "\t";
        os << fg.tgene << "\t" << fg.tend << "\t" << fg.tstrand << "\t";
        return os;
    }
};

/** useful functions used in sv calling */
namespace svutil{
    /** convert an C string into its unsigned int hash value
     * @param s C string
     * @return unsigned hash value of s
     */
    inline size_t hashString(const char *s){
        size_t h = 37;
        while(*s){
            h = (h * 54059) ^ (s[0] * 76963);
            ++s;
        }
        return h;
    }

    /** combine an hash value with another number to generate a new hash value
     * @param seed input hash value
     * @param v number to be combined into seed
     */
    template<typename T>
    inline void hashCombine(size_t& seed, const T& v){
        seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }

    /** get hash value of current bam record
     * @param b pointer to bam1_t
     * @return hash value representation of b
     */
    inline size_t hashPairCurr(bam1_t* b){
        size_t seed = hashString(bam_get_qname(b));
        hashCombine(seed, b->core.tid);
        hashCombine(seed, b->core.pos);
        hashCombine(seed, b->core.mtid);
        hashCombine(seed, b->core.mpos);
        return seed;
    }

    /** get hash value of mate bam record
     * @param b pointer to bam1_t
     * @return hash value of mate of b in a pair
     */
    inline size_t hashPairMate(bam1_t* b){
        size_t seed = hashString(bam_get_qname(b));
        hashCombine(seed, b->core.mtid);
        hashCombine(seed, b->core.mpos);
        hashCombine(seed, b->core.tid);
        hashCombine(seed, b->core.pos);
        return seed;
    }

    /** transform an allele string into IUPAC format
     * @param alleles allele string
     * @return IUPAC format allele string
     */
    inline std::string replaceIUPAC(std::string const& alleles){
        std::vector<char> out(alleles.size());
        int32_t inTag = 0;
        const std::string validStr = "ACGTNacgtn<>[],";
        for(uint32_t i = 0; i<alleles.size(); ++i){
            // Output valid character of mutation and sv alt info
            if((validStr.find_first_of(alleles[i]) != std::string::npos) || (inTag)){
                out[i] = alleles[i];
                if(inTag == 1 && alleles[i] == '>') inTag = 0;
                else if(inTag == 2 && alleles[i] == ']') inTag = 0;
                else if(inTag == 3 && alleles[i] == '[') inTag = 0;
                else if(alleles[i] == '<') inTag = 1;
                else if(alleles[i] == ']') inTag = 2;
                else if(alleles[i] == '[') inTag = 3;
            }else{
                // Replace IUPAC
                if((alleles[i] == 'U') || (alleles[i] == 'u')) out[i] = 'T';
                else if((alleles[i] == 'R') || (alleles[i] == 'r')) out[i] = 'A';
                else if((alleles[i] == 'Y') || (alleles[i] == 'y')) out[i] = 'C';
                else if((alleles[i] == 'S') || (alleles[i] == 's')) out[i] = 'C';
                else if((alleles[i] == 'W') || (alleles[i] == 'w')) out[i] = 'A';
                else if((alleles[i] == 'K') || (alleles[i] == 'k')) out[i] = 'G';
                else if((alleles[i] == 'M') || (alleles[i] == 'm')) out[i] = 'A';
                else if((alleles[i] == 'B') || (alleles[i] == 'b')) out[i] = 'C';
                else if((alleles[i] == 'D') || (alleles[i] == 'd')) out[i] = 'A';
                else if((alleles[i] == 'H') || (alleles[i] == 'h')) out[i] = 'A';
                else if((alleles[i] == 'V') || (alleles[i] == 'v')) out[i] = 'A';
                else out[i] = 'N';
            }
        }
        return std::string(out.begin(), out.end());
    }

    /** get abbreviation of sv type
     * @param svt sv type
     * @return str rep of svt
     */ 
    inline std::string addID(int svt){
        switch(svt){
            case 0: case 1:
                return "INV";
                break;
            case 2:
                return "DEL";
                break;
            case 3:
                return "DUP";
                break;
            case 4:
                return "INS";
                break;
            default:
                return "BND";
        }
    }

    /** get bp description str of SV type
     * @param svt SV type
     * @return string description of bp
     */
    inline std::string getBpMark(int svt){
        switch(svt){
            case 0:
                return "L";
                break;
            case 1:
                return "R";
                break;
            case 2:
                return "-";
                break;
            case 3:
                return "-";
                break;
            case 4:
                return "-";
                break;
            case 5:
                return "5'->5'";
                break;
            case 6:
                return "3'->3'";
                break;
            case 7:
                return "5'->3'";
                break;
            case 8:
                return "3'->5'";
                break;
            default:
                return "-";
                break;
        }
    }

    /** get string representation of SV type
     * @param svt SV type
     * @return string representation of svt
     */
    inline std::string addOrientation(int svt){
        if(svt >= 5) svt -= 5;
        switch(svt){
            case 0:
                return "5to5";
                break;
            case 1:
                return "3to3";
                break;
            case 2:
                return "5to3";
                break;
            case 3:
                return "3to5";
                break;
            default:
                return "NtoN";
                break;
        }
    }

    /** get entropy of an string 
     * @param st
     * @return entropy of st
     */
    inline double entropy(const std::string& st){
        std::map<char, int32_t> countMap;
        for(uint32_t i = 0; i < st.size(); ++i){
            ++countMap[st[i]];
        }
        double ent = 0, freq = 0, len = st.size();
        for(auto& e: countMap){
            freq = e.second/len;
            ent += freq * std::log2(freq);
        }
        return -ent;
    }

    /** compute genotype likelihood based on mapping quality of REF and ALT at one site
     * @param mapqRef mapping quality of reads supporting REF at one site
     * @param mapqAlt mapping quality of reads supporting ALT at the site
     * @param gls phred-scaled genotype likelihoods of alleles at the site
     * @param gts genotypes of each alleles at the site
     * @param gqval genotype quality of each alleles at the site
     */
    inline void computeGL(const std::vector<uint8_t>& mapqRef, const std::vector<uint8_t>& mapqAlt, float* gls, int32_t* gts, int32_t* gqval){
        if(mapqRef.empty() || mapqAlt.empty()){
            gqval[0] = 0;
            gts[0] = bcf_gt_missing;
            gts[1] = bcf_gt_missing;
            gls[0] = 0;
            gls[1] = 0;
            gls[2] = 0;
            return;
        }
        const double minGL = -1000;
        double gl[3] = {0}; // gl[0] = log10(p(ALT|ReadsObserved)), gl[2] = log10(p(REF|ReadsObserved)), gl[1] = log10(p(RandomALT/REF| ReadsObserved))
        // Compute genotype likelihoods
        int32_t depth = mapqRef.size() + mapqAlt.size();
        for(uint32_t i = 0; i < mapqRef.size(); ++i){
            gl[0] += (double)mapqRef[i]/(-10); // log10(p(ALT|RefRead))
            gl[2] += std::log10(1 - std::pow(10, (double)mapqRef[i]/(-10))); // log10(p(REF|RefRead))
        }
        for(uint32_t i = 0; i < mapqAlt.size(); ++i){
            gl[0] += std::log10(1 - std::pow(10, (double)mapqAlt[i]/(-10))); // log10(p(ALT|AltRead))
            gl[2] += (double)mapqAlt[i]/(-10);// log10(p(REF|AltRead))
        }
        gl[1] += -(double)(depth) * std::log10(2.0);// log10(p(Random ALT/REF of each read));
        // Get largest genotype likelihood
        uint32_t glbesti = 0;
        double glbestv = gl[glbesti];
        for(int geno = 1; geno < 3; ++geno){
            if(gl[geno] > glbestv){
                glbestv = gl[geno];
                glbesti = geno;
            }
        }
        // Rescale genotype likelihoods by minus largest likelihood
        for(int geno = 0; geno < 3; ++geno){
            gl[geno] -= glbestv;
            gl[geno] = (gl[geno] > minGL ) ? gl[geno] : minGL;// cap at minGL
        }
        // Phred-scaled genotype liklihoods
        uint32_t pl[3] = {0};// pl[0] = phredQ of ALT, pl[2] = phredQ of REF, pl[1] = phredQ of REF/ALT
        pl[0] = (uint32_t)std::round(-10 * gl[0]);
        pl[1] = (uint32_t)std::round(-10 * gl[1]);
        pl[2] = (uint32_t)std::round(-10 * gl[2]);
        if(depth && (pl[0] + pl[1] + pl[2]) > 0){
            double likelihood = std::log10(1 - 1 / (std::pow(10, gl[0]) + std::pow(10, gl[1]) + std::pow(10, gl[2])));
            likelihood = likelihood > minGL ? likelihood : minGL;
            gqval[0] = std::round(-10 * likelihood);
            if(glbesti == 0){
                gts[0] = bcf_gt_unphased(1); ///< 110
                gts[1] = bcf_gt_unphased(1); ///< 110
            }else if(glbesti == 1){
                gts[0] = bcf_gt_unphased(0); ///< 010
                gts[1] = bcf_gt_unphased(1); ///< 110
            }else{
                gts[0] = bcf_gt_unphased(0); ///< 010
                gts[1] = bcf_gt_unphased(0); ///< 010
            }
        }else{
            gts[0] = bcf_gt_missing;
            gts[1] = bcf_gt_missing;
            gqval[0] = 0;
        }
        gls[0] = gl[2];
        gls[1] = gl[1];
        gls[2] = gl[0];
    }

    /** get fusion gene information at a breakpoint
     * @param gene1 gene name of first breakpoint(lower coordinate breakpoint if sv on same chr, else breakpoint on bigger chr)
     * @param gene2 gene name of second breakpoint(higher coordinate breakpoint if sv on same chr, else breakpoint on little chr)
     * @param strand1 strand of gene1
     * @param strand2 strand of gene2
     * @param svt SV type
     * @return FuseGene information
     */
    inline FuseGene getFusionGene(std::string gene1, std::string gene2, char strand1, char strand2, int32_t svt){
        FuseGene ret;
        if(gene1 != "-" && gene2 != "-") ret.status |= FUSION_FALLGENE;
        if(gene1 == "-" || gene2 == "-" || svt == 4) return ret;
        if(svt >= 5) svt -= 5;
        if(svt == 3){// convert 3to5 to 5to3
            std::swap(gene1, gene2);
            std::swap(strand1, strand2);
        }
        if(svt == 0){// 5to5
            if(strand1 == '+' && strand2 == '+'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "5";
                ret.hstrand = "+";
                ret.tstrand = "-";
            }else if(strand1 == '+' && strand2 == '-'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
            }else if(strand1 == '-' && strand2 == '+'){
                ret.hgene = gene2;
                ret.tgene = gene1;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
            }else if(strand1 == '-' && strand2 == '-'){
                ret.hgene = gene2;
                ret.tgene = gene1;
                ret.hend = "3";
                ret.tend = "3";
                ret.hstrand = "-";
                ret.tstrand = "+";
            }
        }else if(svt == 1){// 3to3
            if(strand1 == '+' && strand2 == '+'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "3";
                ret.tend = "3";
                ret.hstrand = "-";
                ret.tstrand = "+";
            }else if(strand1 == '+' && strand2 == '-'){
                ret.hgene = gene2;
                ret.tgene = gene1;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
            }else if(strand1 == '-' && strand2 == '+'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
            }else if(strand1 == '-' && strand2 == '-'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "5";
                ret.hstrand = "+";
                ret.tstrand = "-";
            }
        }else if(svt == 2 || svt == 3){// 5to3 and 3to5
            if(strand1 == '+' && strand2 == '+'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
            }else if(strand1 == '+' && strand2 == '-'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "5";
                ret.hstrand = "+";
                ret.tstrand = "-";
            }else if(strand1 == '-' && strand2 == '+'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "3";
                ret.tend = "3";
                ret.hstrand = "-";
                ret.tstrand = "+";
            }else if(strand1 == '-' && strand2 == '-'){
                ret.hgene = gene2;
                ret.tgene = gene1;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
            }
        }
        return ret;
    }

    /** convert string representation of SV event back to integer representation
     * @param ct catenation type of SV
     * @param svt string representation of SV
     * @return integer representation of SV
     */
    inline int32_t str2svt(const std::string& ct, const std::string& svt){
        if(svt == "INV"){
            if(ct == "5to5") return 0;
            if(ct == "3to3") return 1;
        }else if(svt == "DEL"){
            if(ct == "5to3") return 2;
        }else if(svt == "DUP"){
            if(ct == "3to5") return 3;
        }else if(svt == "INS"){
            if(ct == "NtoN") return 4;
        }else if(svt == "BND"){
            if(ct == "5to5") return 5;
            if(ct == "3to3") return 6;
            if(ct == "5to3") return 7;
            if(ct == "3to5") return 8;
        }
        return -1;
    }

    /** test a sequence is too simple
     * @param seq sequence to be tested
     * @return true if a sequence is too simple
     */
    inline bool simpleSeq(const std::string& seq){
        int32_t cnt[4] = {0, 0, 0, 0};
        for(auto& e: seq){
            switch(e){
                case 'A':
                    cnt[0]++;
                    break;
                case 'T':
                    cnt[1]++;
                    break;
                 case 'C':
                    cnt[2]++;
                    break;
                 case 'G':
                    cnt[3]++;
                 default:
                    break;
            }
        }
        int32_t maxAllowed = 0.8 * seq.length();
        for(int i = 0; i < 4; ++i){
            for(int j = i + 1; j < 4; ++j){
                if((cnt[i] + cnt[j]) > maxAllowed) return true;
            }
        }
        return false;
    }

    /** test whether two transcript units are near each other
     * @param trs1 transcript unit representation str vector, str format: name,unit,strand,no
     * @param trs2 transcript uint representation str vector, str format: name,unit,strand,no
     * @param maxoffset max intron(exon) offset to be near
     * @return true if trs1 and trs2 is near
     */
    inline bool trsUnitIsNear(const std::vector<std::string>& trs1, const std::vector<std::string>& trs2, int32_t maxoffset = 1){
        std::map<std::string, std::map<std::string, int32_t>> m1;
        std::map<std::string, std::map<std::string, int32_t>> m2;
        // parse trs1
        std::vector<std::string> vstr;
        for(auto& e: trs1){
            util::split(e, vstr, ",");
            if(util::startsWith(vstr[1], "utr")){
                m1[vstr[0]]["utr"] = 0;
            }else{
                m1[vstr[0]]["ie"] = std::atoi(vstr[3].c_str());
            }
        }
        // parse trs2
        for(auto& e: trs2){
            util::split(e, vstr, ",");
            if(util::startsWith(vstr[1], "utr")){
                m2[vstr[0]]["utr"] = 0;
            }else{
                m2[vstr[0]]["ie"] = std::atoi(vstr[3].c_str());
            }
        }
        // calculate distance
        for(auto& e: m1){
            auto iter1 = m2.find(e.first);
            if(iter1 != m2.end()){
                for(auto& f: e.second){
                    auto iter2 = iter1->second.find(f.first);
                    if(iter2 != iter1->second.end()){
                        return std::abs(iter2->second - f.second) <= maxoffset;
                    }
                }
            }
        }
        return false;
    }

    /** test whether two transcript units are near each other
     * @param trs1 transcript unit representation str, format: name,unit,strand,no[;name,unit,strand,no]..
     * @param trs2 transcript uint representation str, format: name,unit,strand,no[;name,unit,strand,no]..
     * @param maxoffset max intron(exon) offset to be near
     * @return true if trs1 and trs2 is near
     */
    inline bool trsUnitIsNear(const std::string& trs1, const std::string& trs2, int32_t maxoffset = 1){
        std::vector<std::string> vstr1;
        std::vector<std::string> vstr2;
        // parse trs1 and trs2
        util::split(trs1, vstr1, ";");
        util::split(trs2, vstr2, ";");
        return trsUnitIsNear(vstr1, vstr2, maxoffset);
    }

    /** convert transcript coordinate into genome coordinate
     * @param tpos transcript position
     * @param tipos transcriptome starting position containing this pos
     * @param tepos transcriptome ending position containing this pos
     * @param gipos genome starting posision containing this pos
     * @param strand genome strand this transcript comes from
     */
    inline int32_t trpos2gnpos(int32_t tpos, int32_t tipos, int32_t tepos, int32_t gipos, char strand){
        if(strand == '+'){
            return gipos + (tpos - tipos);
        }else{
            return gipos + (tepos - tpos);
        }
    }
}

#endif
