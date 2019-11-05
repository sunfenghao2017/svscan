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
#include "trsrec.h"
#include "fusegene.h"
#include <htslib/vcf.h>
#include <htslib/sam.h>

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
                ret.status |= FUSION_FHTFLSWAPPED;
            }else if(strand1 == '-' && strand2 == '-'){
                ret.hgene = gene2;
                ret.tgene = gene1;
                ret.hend = "3";
                ret.tend = "3";
                ret.hstrand = "-";
                ret.tstrand = "+";
                ret.status |= FUSION_FHTFLSWAPPED;
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
                ret.status |= FUSION_FHTFLSWAPPED;
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
                if(svt == 3) ret.status |= FUSION_FHTFLSWAPPED;
            }else if(strand1 == '+' && strand2 == '-'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "5";
                ret.tend = "5";
                ret.hstrand = "+";
                ret.tstrand = "-";
                if(svt == 3) ret.status |= FUSION_FHTFLSWAPPED;
            }else if(strand1 == '-' && strand2 == '+'){
                ret.hgene = gene1;
                ret.tgene = gene2;
                ret.hend = "3";
                ret.tend = "3";
                ret.hstrand = "-";
                ret.tstrand = "+";
                if(svt == 3) ret.status |= FUSION_FHTFLSWAPPED;
            }else if(strand1 == '-' && strand2 == '-'){
                ret.hgene = gene2;
                ret.tgene = gene1;
                ret.hend = "5";
                ret.tend = "3";
                ret.hstrand = "+";
                ret.tstrand = "+";
                ret.status |= FUSION_FNORMALCATDIRECT;
                if(svt == 2) ret.status |= FUSION_FHTFLSWAPPED;
            }
        }else if(svt == 4){
            ret.hgene = gene1;
            ret.tgene = gene2;
            ret.hend = "5";
            ret.tend = "3";
            ret.hstrand = "+";
            ret.tstrand = "+";
            ret.status |= FUSION_FNORMALCATDIRECT;
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

    /** get cigar string representation of breakpoint
     * @param htrs 5' part transcript 
     * @param ttrs 3' part transcript
     * @param svt SV type
     * @return cigar string representation of breakpoint
     */
    inline std::string bp2cigar(const TrsRec& htrs, const TrsRec& ttrs, int32_t svt){
        std::stringstream css;
        switch(svt){
            case 5: // 5->5 catenation
                if(htrs.ioffset < 10){
                    css << htrs.unit << std::atoi(htrs.number.c_str()) - 1 << ",tI" << htrs.ioffset;
                }else if(htrs.eoffset < 10){
                    css << htrs.unit << htrs.number << ",tD" << htrs.eoffset;
                }else{
                    css << htrs.unit << std::atoi(htrs.number.c_str()) << ",hM" << htrs.ioffset;
                }
                if(ttrs.ioffset < 10){
                    css << "tI" << ttrs.ioffset << "," << ttrs.unit << std::atoi(ttrs.number.c_str()) - 1;
                }else if(ttrs.eoffset < 10){
                    css << "tD" << ttrs.eoffset << "," << ttrs.unit << ttrs.number;
                }else{
                    css << "hM" << ttrs.ioffset << "," << ttrs.unit << ttrs.number;
                }
                break;
            case 6: // 3->3 catenation
                if(htrs.ioffset < 10){
                    css << htrs.unit << htrs.number << ",hD" << htrs.ioffset;
                }else if(htrs.eoffset < 10){
                    css << htrs.unit << std::atoi(htrs.number.c_str()) + 1 << ",hI" << htrs.eoffset;
                }else{
                    css << htrs.unit << htrs.number << ",tM" << htrs.eoffset;
                }
                if(ttrs.ioffset < 10){
                    css << "hD" << ttrs.ioffset << "," << ttrs.unit << ttrs.number;
                }else if(ttrs.eoffset < 10){
                    css << "hI" << ttrs.eoffset << "," << ttrs.unit << std::atoi(ttrs.number.c_str()) + 1;
                }else{
                    css << "tM" << ttrs.eoffset << "," << ttrs.unit << ttrs.number;
                }
                break;
            case 7: case 8: // 5->3 catenation
                if(htrs.ioffset < 10){
                    css << htrs.unit << std::atoi(htrs.number.c_str()) - 1 << ",tI" << htrs.ioffset;
                }else if(htrs.eoffset < 10){
                    css << htrs.unit << htrs.number << ",tD" << htrs.eoffset;
                }else{
                    css << htrs.unit << htrs.number << ",hM" << htrs.ioffset;
                }
                if(ttrs.ioffset < 10){
                    css << "hD" << ttrs.ioffset << "," << ttrs.unit << ttrs.number;
                }else if(ttrs.eoffset < 10){
                    css << "hI" << ttrs.eoffset << "," << ttrs.unit << std::atoi(htrs.number.c_str()) + 1;
                }else{
                    css << "tM" << ttrs.ioffset << "," << ttrs.unit << ttrs.number;
                }
                break;
             default:
                break;
        }
        return css.str();
    }
}

#endif
