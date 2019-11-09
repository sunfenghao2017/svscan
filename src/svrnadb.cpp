#include "svrnadb.h"

void SVRNADBOpt::prepDB(){
    // parse main transcripts
    std::set<std::string> mainTrs;
    std::ifstream fr(cncTrsList);
    std::string tmpStr;
    while(std::getline(fr, tmpStr)) mainTrs.insert(tmpStr);
    // get gene trsncript map set
    std::map<std::string, std::set<std::string>> gene2TrsSet;
    BGZF* ifp = bgzf_open(refSeqDB.c_str(), "rb");
    kstring_t str = {0, 0, 0};
    bgzf_getline(ifp, '\n', &str);
    std::vector<std::string> vstr;
    std::map<std::string, std::string> outGeneCncMap;
    while(bgzf_getline(ifp, '\n', &str) > 0){
        util::split(str.s, vstr, "\t");
        std::string trs = vstr[1];
        std::string gene = vstr[12];
        gene2TrsSet[gene].insert(trs);
        if(mainTrs.find(trs) != mainTrs.end()){
            outGeneCncMap[gene] = trs;
        }
    }
    for(auto& e: gene2TrsSet){
        if(outGeneCncMap.find(e.first) == outGeneCncMap.end()){
            for(auto& f: e.second){
                outGeneCncMap[e.first] = f;
                break;
            }
        }
    }
    // write transcript used by gene to file
    std::ofstream oft(gene2trs);
    for(auto& e: outGeneCncMap){
        if(!keepNCRna && util::startsWith(e.second, "NR")) continue;
        oft << e.first << "\t" << e.second << "\n";
    }
    oft.close();
    // output selected gene with its main transcript database
    bgzf_close(ifp);
    ifp = bgzf_open(refSeqDB.c_str(), "rb");
    BGZF* ofa = bgzf_open(svAnnoDB.c_str(), "wb"); // annotation db
    BGZF* off = bgzf_open(refMrna.c_str(), "wb");  // refmrna sequence
    std::ofstream ofc(bedCDS); // cds bed region
    std::ofstream ofu(bedUnits); // unit bed region
    std::ofstream ofl(utr3len); // 3'utr length
    faidx_t* fai = fai_load(genome.c_str());
    bgzf_getline(ifp, '\n', &str);
    std::vector<std::string> istr, estr;
    std::vector<int32_t> iint, eint;
    std::stringstream aoss; // annotation db record buffer
    std::stringstream coss; // cds bed region record buffer
    std::stringstream uoss; // unit bed region buffer
    std::stringstream loss; // 3'utr length buffer
    std::set<std::string> processedTrs;
    // chr start end strand feature count trsname tgenename trsversion
    while(bgzf_getline(ifp, '\n', &str) > 0){
        aoss.str(""); coss.str(""); uoss.str(""); loss.str("");  
        aoss.clear(); coss.clear(); uoss.clear(); loss.clear();
        util::split(str.s, vstr, "\t");
        std::string trs = vstr[1];
        if(!keepNCRna && util::startsWith(trs, "NR")) continue; // drop ncrna if needed
        std::string chr = vstr[2];
        std::string gene = vstr[12];
        std::string version = vstr[16];
        std::string strand = vstr[3];
        int32_t trsStart = std::atoi(vstr[4].c_str());
        int32_t trsEnd = std::atoi(vstr[5].c_str()) - 1;
        int32_t cdsStart = std::atoi(vstr[6].c_str());
        int32_t cdsEnd = std::atoi(vstr[7].c_str()) - 1;
        if(!faidx_has_seq(fai, chr.c_str())){
            util::loginfo(chr + " does not exist in your reference !!!");
            continue; // only fetch chr in genome
        }
        if(outGeneCncMap[gene] != trs) continue; // only process computed canonical transcripts
        if(processedTrs.find(trs) != processedTrs.end()) continue;
        else processedTrs.insert(trs); // only process each trascript once
        // exon range (dna)
        util::split(util::rstrip(vstr[9], ","), istr, ",");
        util::split(util::rstrip(vstr[10], ","), estr, ",");
        util::strvec2intvec(istr, iint);
        util::strvec2intvec(estr, eint);
        std::vector<int32_t> exonlens(iint.size(), 0);
        // calculate exon length
        int32_t exonttl = 0;
        for(uint32_t eidx = 0; eidx < iint.size(); ++eidx){
            exonlens[eidx] = eint[eidx] - iint[eidx];
            exonttl += exonlens[eidx];

        }
        // calculate utr length
        int32_t lutrlen = 0; // left utr
        int32_t rutrlen = 0; // right utr
        for(uint32_t i = 0; i < iint.size(); ++i){
            if(cdsStart > iint[i]){
                if(cdsStart < eint[i]) lutrlen += cdsStart - iint[i];
                else lutrlen += exonlens[i];
            }
            if(cdsEnd < eint[i]){
                if(cdsEnd < iint[i]) rutrlen += exonlens[i];
                else rutrlen += eint[i] - cdsEnd - 1;
            }
        }
        // utr range (dna)
        int32_t utr5start = trsStart;
        int32_t utr5end = cdsStart - 1;
        int32_t utr3start = cdsEnd + 1;
        int32_t utr3end = trsEnd;
        // utr range (rna) strand +
        int32_t rutr5start = 0;
        int32_t rutr5end = lutrlen - 1;
        int32_t rutr3start = exonttl - rutrlen;
        int32_t rutr3end = exonttl - 1;
        if(strand[0] == '-'){ // strand -
            // dna
            utr5start = cdsEnd + 1;
            utr5end = trsEnd;
            utr3start = trsStart;
            utr3end = cdsStart - 1;
            // rna
            rutr5start = 0;
            rutr5end = rutrlen - 1;
            rutr3start = exonttl - lutrlen;
            rutr3end = exonttl - 1;
        }
        // generate annotation record for each unit
        if(utr5end > utr5start){ // 5'utr
            // transcriptname rna-start rna-end uint count genename 
            aoss << trs << "\t" << rutr5start << "\t" << rutr5end << "\tutr5\t0\t" << gene << "\t";
            // chr dna-start dna-end dna-strand trsversion
            aoss << chr << "\t" << utr5start << "\t" << utr5end << "\t" << strand << "\t" << version << "\n";
            // unit bed region
            uoss << trs << "\t" << rutr5start << "\t" << rutr5end + 1 << "\t" << trs << "\tutr5\n";
        }
        if(utr3start < utr3end){ // 3'utr
            // transcriptname rna-start rna-end uint count genename 
            aoss << trs << "\t" << rutr3start << "\t" << rutr3end << "\tutr3\t0\t" << gene << "\t";
            // chr dna-start dna-end dna-strand trsversion
            aoss << chr << "\t" << utr3start << "\t" << utr3end << "\t" << strand << "\t" << version << "\n";
            // unit bed region
            uoss << trs << "\t" << rutr3start << "\t" << rutr3end + 1 << "\t" << trs << "\tutr3\n";
            // 3'utr length
            loss << trs << "\t" << (rutr3end - rutr3start) << "\n"; 
        }
        if(strand == "+"){ // forward gene
            int32_t acculen = 0;
            for(uint32_t i = 0; i < exonlens.size(); ++i){
                aoss << trs << "\t" << acculen << "\t" << (acculen + exonlens[i] - 1) << "\texon\t" << (i + 1) << "\t" << gene << "\t";
                aoss << chr << "\t" << iint[i] << "\t" << (eint[i] - 1) << "\t" << strand << "\t" << version << "\n";
                // unit bed region
                uoss << trs << "\t" << acculen << "\t" << acculen + exonlens[i] << "\t" << trs << "\texon" << (i + 1) << "\n";
                acculen += exonlens[i];
            }
        }else{// reverse gene
            int32_t acculen = 0;
            for(int32_t i = exonlens.size() - 1; i >= 0; --i){
                aoss << trs << "\t" << acculen << "\t" << (acculen + exonlens[i] - 1) << "\texon\t" << (exonlens.size() - i) << "\t" << gene << "\t";
                aoss << chr << "\t" << iint[i] << "\t" << (eint[i] - 1) << "\t" << strand << "\t" << version << "\n";
                // unit bed region
                uoss << trs << "\t" << acculen << "\t" << acculen + exonlens[i] << "\t" << trs << "\texon" << (exonlens.size() - i) << "\n";
                acculen += exonlens[i];
            }
        }
        // cds bed region
        int32_t rcdsBeg = std::max(rutr5end, 0);
        int32_t rcdsEnd = rutr3start;
        if(rcdsBeg < rcdsEnd) coss << trs << "\t" << rcdsBeg << "\t" << rcdsEnd << "\t" << trs << "\n";
        // output annotation record
        std::string arec = aoss.str();
        assert(bgzf_write(ofa, arec.c_str(), arec.size()) >= 0);
        // output cds region bed
        std::string crec = coss.str();
        ofc << crec;
        // output uint region bed
        std::string urec = uoss.str();
        ofu << urec;
        // output 3utr length
        std::string lrec = loss.str();
        ofl << lrec;
        // fetch trascript sequence
        std::string seqName = ">" + trs + "\n";
        std::string rnaSeq = "";
        int32_t len = -1;
        char* refseq = NULL;
        for(uint32_t i = 0; i < iint.size(); ++i){
            refseq = faidx_fetch_seq(fai, chr.c_str(), iint[i], eint[i] - 1, &len);
            rnaSeq.append(refseq);
            free(refseq);
        }
        util::str2upper(rnaSeq);
        if(strand == "-") util::reverseComplement(rnaSeq);
        rnaSeq.append("\n");
        // output seq
        assert(bgzf_write(off, seqName.c_str(), seqName.length()) >= 0);
        assert(bgzf_write(off, rnaSeq.c_str(), rnaSeq.length()) >= 0);
    }
    bgzf_close(ofa);
    bgzf_close(off);
    ofc.close();
    ofu.close();
    ofl.close();
    bgzf_close(ifp);
    fai_destroy(fai);
}
