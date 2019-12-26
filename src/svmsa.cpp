#include "msa.h"
#include "util.h"
#include "aligner.h"

/** class to store MSA argument */
struct MSAOpt{
    int32_t mMinCovForCS = 3;          ///< minimum coverage needed for a position in msa result to be included in consensus sequence
    float mMinBaseRateForCS = 0.5;     ///< minimum base ratio needed for a position in msa result to be included in consensus sequence
    bool mAlignHorzEndGapFree = false; ///< use horizontal end gap penalty free strategy when get consensus sequence of SR
    bool mALignVertEndGapFree = false; ///< use vertical end gap penalty free strategy when get consensus sequence of SR

    /** MSAOpt consturcor */
    MSAOpt(){}

    /** MSAOpt destructor */
    ~MSAOpt(){}
};

int main(int argc, char** argv){
    if(argc < 2){
        std::cout << argv[0] << " <seqf> " << std::endl;
        return 0;
    }
    MSAOpt* msaOpt = new MSAOpt();
    AlignConfig alnCfg(5, -4, -10, -1, true, true);
    std::multiset<std::string> seqs;
    std::cout << "get input seqs" << std::endl;
    std::ifstream fr(argv[1]);
    std::string line;
    while(std::getline(fr, line)){
        seqs.insert(line);
    }
    fr.close();
    std::cout << "beg msa" << std::endl;
    MSA* msa = new MSA(&seqs, msaOpt->mMinCovForCS, msaOpt->mMinBaseRateForCS, &alnCfg);
    std::cout << "try 100" << std::endl;
    std::string ccseq;
    for(int i = 0; i < 100; ++i){
        ccseq = "";
        msa->msa(ccseq);
        std::cout << ">try" << i << "\n";
        std::cout << ccseq << "\n";
    }
}
