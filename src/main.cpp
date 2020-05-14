#include <CLI.hpp>
#include "options.h"
#include "svscanner.h"
#include <iostream>

int main(int argc, char** argv){
    if(argc == 1){
        std::string helpCMD = std::string(argv[0]) + " -h";
        std::system(helpCMD.c_str());
        return 0;
    }
    // parse commandline arguments
    Options* opt = new Options();
    CLI::App app("program: " + std::string(argv[0]) + "\n" + opt->softEnv->getSoftInfo());
    app.get_formatter()->column_width(50);
    // General Options options
    app.add_option("-b,--bam", opt->bamfile, "indexed input bam file sorted by coordinates")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-g,--alnref", opt->alnref, "indexed reference fasta file used in input bam")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-a,--anno", opt->annodb, "indexed annotation database file")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-r,--sreg", opt->reg, "file of bed regions to scan bam")->required(false)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-c,--creg", opt->creg, "file of bed regions sv must overlap")->required(false)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-p,--preg", opt->preg, "file of bed regions of hot fusion pairs")->required(false)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-o,--bcfout", opt->bcfOut, "output sv bcf file", true)->required(false)->group("General Options");
    app.add_option("-t,--tsvout", opt->tsvOut, "output sv tsv file", true)->required(false)->group("General Options");
    app.add_option("-v,--bamout", opt->bamout, "output sv bam file, non-providing will disable it", true)->group("General Options");
    app.add_option("-e,--b2excel", opt->bam2tb, "excel file to output fusion event supporting bam records", true)->required(false)->group("General Options");
    app.add_option("-f,--b2tsv", opt->bam2tt, "tsv file to output fusion event supporting bam records", true)->required(false)->group("General Options");
    app.add_option("-s,--svtype", opt->svtypes, "SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND")->check(CLI::Range(0, 4))->group("General Options");
    app.add_option("-n,--nthread", opt->nthread, "number of threads used", true)->check(CLI::Range(1, 100))->group("General Options");
    app.add_option("-d,--debug", opt->debug, "debug mask,1:call,2:cann,4:gann,8:out,16:realn,32:finsv:64:readrescue", true)->group("General Options");
    app.add_option("-q,--qname", opt->qndbg, "qname of read to debug")->group("General Options");
    CLI::Option* prna = app.add_flag("--rna", opt->rnamode, "discovery structural variants from rna data if set")->group("General Options");
    app.add_option("--genome", opt->genome, "indexed genome fasta corresponding to the rna transcriptome")->needs(prna)->group("General Options");
    app.add_option("--ganno", opt->gannodb, "indexed genome fasta annotation database")->needs(prna)->group("General Options");
    // Control options
    app.add_option("--min_ref_sep", opt->filterOpt->mMinRefSep, "min sv length to compute", true)->group("Threshold Options");
    app.add_option("--max_read_sep", opt->filterOpt->mMaxReadSep, "max read split mapping pos allowed to compute sv", true)->group("Threshold Options");
    app.add_option("--min_flk_len", opt->filterOpt->mMinFlankSize, "min flank length needed on each side of breakpoint", true)->group("Threshold Options");
    app.add_option("--min_map_qual", opt->filterOpt->minMapQual, "min mapping quality of reads used to compute sv", true)->group("Threshold Options");
    app.add_option("--min_clip_len", opt->filterOpt->minClipLen, "min clip length needed to compute sv", true)->group("Threshold Options");
    app.add_option("--min_tra_qual", opt->filterOpt->mMinTraQual, "min mapping quality of discordant pair used to compute sv", true)->group("Threshold Options");
    app.add_option("--min_dup_dp", opt->filterOpt->mMinDupSize, "min duplication size needed for a DP to compute sv", true)->group("Threshold Options");
    app.add_option("--min_inv_rpt", opt->filterOpt->mMinInversionRpt, "min inversion size to report", true)->group("Threshold Options");
    app.add_option("--min_del_rpt", opt->filterOpt->mMinDeletionRpt, "min deletion size to report", true)->group("Threshold Options");
    app.add_option("--min_dup_rpt", opt->filterOpt->mMinDupRpt, "min dup size to report", true)->group("Threshold Options");
    app.add_option("--min_seed_sr", opt->filterOpt->mMinSeedSR, "min seed split reads needed to compute SV", true)->check(CLI::Range(1, 1000))->group("Threshold Options");
    app.add_option("--min_seed_dp", opt->filterOpt->mMinSeedDP, "min seed discordant pair needed to computr SV", true)->check(CLI::Range(1, 1000))->group("Threshold Options");
    // Fusion report options
    app.add_option("--whitemindep", opt->fuseOpt->mWhiteFilter.mMinDepth, "min depth for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmindep", opt->fuseOpt->mUsualFilter.mMinDepth, "min depth for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminsrs", opt->fuseOpt->mWhiteFilter.mMinSRSeed, "min sr seeds for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminsrs", opt->fuseOpt->mUsualFilter.mMinSRSeed, "min sr seeds for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whitemindps", opt->fuseOpt->mWhiteFilter.mMinDPSeed, "min dp seeds for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmindps", opt->fuseOpt->mUsualFilter.mMinDPSeed, "min dp seeds for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminsrr", opt->fuseOpt->mWhiteFilter.mMinSRSupport, "min sr reads for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminsrr", opt->fuseOpt->mUsualFilter.mMinSRSupport, "min sr reads for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whitemindpr", opt->fuseOpt->mWhiteFilter.mMinDPSupport, "min dp reads for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmindpr", opt->fuseOpt->mUsualFilter.mMinDPSupport, "min dp reads for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminttr", opt->fuseOpt->mWhiteFilter.mMinSupport, "min molecules for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminttr", opt->fuseOpt->mUsualFilter.mMinSupport, "min molecules for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminaf", opt->fuseOpt->mWhiteFilter.mMinVAF, "min VAF for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminaf", opt->fuseOpt->mUsualFilter.mMinVAF, "min VAF for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminigs", opt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize, "min intra-gene(gene-nonegene) sv size for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminigs", opt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize, "min intra-gene(gene-nonegene) sv size for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whitemaxrph", opt->fuseOpt->mWhiteFilter.mMaxRepHit, "max hits on reference allowed for partial seqs of fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmaxrph", opt->fuseOpt->mUsualFilter.mMaxRepHit, "max hits on reference allowed for partial seqs of fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--maxbpoffset", opt->fuseOpt->mMaxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true)->group("Fusion Options");
    app.add_option("--bgbcf", opt->fuseOpt->mBgBCF, "indexed background sv vents BCF file")->group("Fusion Options");
    app.add_option("--whitelist", opt->fuseOpt->mWhiteList, "white list of fusion events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--blacklist", opt->fuseOpt->mBlackList, "black list of fusion events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--samegenel", opt->fuseOpt->mSameGeneSVList, "white list of gene with inner sv events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--extraanno", opt->fuseOpt->mExtraAnnoList, "extra annoation gene list")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--fsrptlist", opt->fuseOpt->mFsRptList, "report range list")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--genecrdlist", opt->fuseOpt->mGeneCrdList, "gene coordinate list")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--fusionrpt", opt->fuseOpt->mOutFile, "fusion report file path", true)->group("Fusion Options");
    app.add_option("--idpdropmask", opt->fuseOpt->mIDBDropMask, "fusion in db drop bit mask", true)->group("Fusion Options");
    app.add_option("--ndbdropmask", opt->fuseOpt->mNDBDropMask, "fusion not in db drop bit mask", true)->group("Fusion Options");
    app.add_option("--keepmasks", opt->fuseOpt->mKeepMasks, "fusion keep masks", true)->group("Fusion Options");
    app.add_option("--primarymask", opt->fuseOpt->mPrimaryMask, "primary fusion mask", true)->group("Fusion Options");
    app.add_option("--hotfspartner", opt->fuseOpt->mHotPartnerList, "hot fusion partner exon list")->check(CLI::ExistingFile)->group("Fusion Options");
    // parse arguments
    CLI11_PARSE(app, argc, argv);
    // validate arguments
    util::loginfo("Command line arguments parsed");
    opt->validate();
    util::loginfo("Command line arguments validated");
    // update some arguments
    opt->update(argc, argv);
    util::loginfo("Options updated");
    // print command line arguments
    util::loginfo("CMD: " + opt->softEnv->cmd);
    util::loginfo("CWD: " + opt->softEnv->cwd);
    util::loginfo("VER: " + opt->softEnv->ver);
    // output library basic information
    util::loginfo("Library basic information:\n" + opt->libInfo->toStr());
    // calling SV
    SVScanner* svScanner = new SVScanner(opt);
    util::loginfo("Beg scanning SV");
    svScanner->scanDPandSR();
    util::loginfo("End scanning SV");
    // write execution time
    util::loginfo("Time consumed: " + opt->softEnv->getExecutionTime());
}
