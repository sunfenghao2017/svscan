#ifndef SVUTIL_H
#define SVUTIL_H

#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <cstdint>
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
 * 5to5
 * 5'----------chrBigF-----------3'||5'--------chrSmallR------------3'
 * 3'----------chrBigR-----------5'||3'--------chrSmallF------------5'
 * 3to3
 * 5'----------chrBigR-----------3'||5'--------chrSmallR------------3'
 * 3'----------chrBigF-----------5'||3'--------chrSmallF------------5'
 * 5to3
 * 5'----------chrBigF-----------3'||5'--------chrSmallF------------3'
 * 3'----------chrBigR-----------3'||3'--------chrSmallR------------3'
 * 3to5
 * 5'----------chrSmallF---------3'||5'--------chrBigF--------------3'
 * 3'----------chrSmallR---------3'||3'--------chrBigR--------------3'
 */

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
                return "Inversion Left BP";
                break;
            case 1:
                return "Inversion Right BP";
                break;
            case 2:
                return "Deletion BP";
                break;
            case 3:
                return "Duplication BP";
                break;
            case 4:
                return "Insertion BP";
                break;
            case 5:
                return "5'->5' Across Chr BP";
                break;
            case 6:
                return "3'->3' Across Chr BP";
                break;
            case 7:
                return "5'->3' Across Chr BP";
                break;
            case 8:
                return "3'->5' Across Chr BP";
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
}

#endif
