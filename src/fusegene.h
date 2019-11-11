#ifndef FUSE_GENE_H
#define FUSE_GENE_H

#include <string>
#include <sstream>
#include <iostream>
#include <cstdint>
#include "util.h"

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

#define FUSION_FALLGENE                 0x1                  ///< All partners are genes
#define FUSION_FNORMALCATDIRECT         0x2                  ///< hgene 5' -> tgene 3' positive strand catenation
#define FUSION_FHOTGENE                 0x4                  ///< one of the partner is in hot fusion gene list
#define FUSION_FCOMMONHOTDIRECT         0x8                  ///< hot gene partner is in the predefined direction in fusion gene list
#define FUSION_FINDB                    0x10                 ///< fusion gene is in database
#define FUSION_FMIRROR                  0x20                 ///< fusion gene has mirror fusion
#define FUSION_FBLACKPAIR               0x40                 ///< fusion gene is in black list
#define FUSION_FFBG                     0x80                 ///< fusion gene is in background samples
#define FUSION_FBLACKGENE               0x100                ///< one partner is in black gene list
#define FUSION_FLOWCOMPLEX              0x200                ///< one partner has low complex concensus sequence
#define FUSION_FPRIMARY                 0x400                ///< this fusion event should be reported primarily
#define FUSION_FSUPPLEMENTARY           0x800                ///< this fusion event should be reported as supplementary
#define FUSION_FTOOSMALLSIZE            0x1000               ///< this fusion event has too small size
#define FUSION_FINSAMEGENE              0x2000               ///< this fusion event occurs in the same gene
#define FUSION_FLOWAF                   0x4000               ///< this fusion event has too low AF
#define FUSION_FLOWSUPPORT              0x8000               ///< this fusion event has too low support
#define FUSION_FLOWDEPTH                0x10000              ///< this fusion event has too low depth around breakpoint
#define FUSION_FPRECISE                 0x20000              ///< this fusion event has precise breakpoint and concensus sequence
#define FUSION_FCALLFROMRNASEQ          0x40000              ///< this fusion event is called from rna seq
#define FUSION_FMIRRORINDB              0x80000              ///< this fusion event's mirror event is in public db
#define FUSION_FHTFLSWAPPED             0x100000             ///< this fusion event's h/t gene swapped against bp1/2 gene

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
    std::string cigar;   ///< cigar string to describe breakpoint(rna only)
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
        ss << hgene << "->" << tgene << "|";
        ss << hgene << "," << hend << "," << hstrand << "|";
        ss << tgene << "," << tend << "," << tstrand;
        if(!cigar.empty()) ss << "|" << cigar;
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

/** class to store a list of FuseGene */
typedef std::vector<FuseGene> FuseGeneList;

#endif