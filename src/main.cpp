#include "CLI.hpp"
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
    app.add_option("-b,--bam", opt->bamfile, "bam file")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-g,--alnref", opt->alnref, "reference file used in input bam")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-a,--anno", opt->annodb, "annotation database file")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-r,--sreg", opt->reg, "file of regions to scan bam")->required(false)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-c,--creg", opt->creg, "file of region sv must capture")->required(false)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-o,--bcfout", opt->bcfOut, "output sv bcf file", true)->required(false)->group("General Options");
    app.add_option("-t,--tsvout", opt->tsvOut, "output sv tsv file", true)->required(false)->group("General Options");
    app.add_option("-v,--bamout", opt->bamout, "output sv bam file, non-providing will disable it", true)->required(true)->group("General Options");
    app.add_option("-e,--b2excel", opt->bam2tb, "output ss/fs event supportint bam records excel", true)->required(false)->group("General Options");
    app.add_option("-s,--svtype", opt->svtypes, "SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND")->check(CLI::Range(0, 4))->group("General Options");
    app.add_option("-n,--nthread", opt->nthread, "number of threads used to process bam", true)->check(CLI::Range(1, 100))->group("General Options");
    app.add_option("--batchsvn", opt->batchsvn, "number of svs to annotate coverage in one thread", true)->check(CLI::Range(1000, 1000000))->group("General Options");
    CLI::Option* prna = app.add_flag("--rna", opt->rnamode, "discovery structural variants from rna data")->group("General Options");
    app.add_option("--genome", opt->genome, "genome file corresponding to the rna transcriptome")->needs(prna)->group("General Options");
    app.add_option("--ganno", opt->gannodb, "genome file annotation database")->needs(prna)->group("General Options");
    app.add_option("-d,--debug", opt->debug, "debug mask,1:call,2:cann,4:gann,8:out,16:realn,32:finsv", true)->group("General Options");
    // Control options
    app.add_option("--min_ref_sep", opt->filterOpt->mMinRefSep, "min sv length to compute", true)->group("Threshold Options");
    app.add_option("--max_read_sep", opt->filterOpt->mMaxReadSep, "max read split mapping pos allowed to compute sv", true)->group("Threshold Options");
    app.add_option("--min_flk_len", opt->filterOpt->mMinFlankSize, "min flank length needed on each side of breakpoint", true)->group("Threshold Options");
    app.add_option("--min_map_qual", opt->filterOpt->minMapQual, "min mapping quality of read used to compute sv", true)->group("Threshold Options");
    app.add_option("--min_clip_len", opt->filterOpt->minClipLen, "min clip length needed to compute sv", true)->group("Threshold Options");
    app.add_option("--min_tra_qual", opt->filterOpt->mMinTraQual, "min mapping quality of read pair used to compute sv", true)->group("Threshold Options");
    app.add_option("--min_dup_dp", opt->filterOpt->mMinDupSize, "min duplication size needed for an DP to compute sv", true)->group("Threshold Options");
    app.add_option("--min_inv_rpt", opt->filterOpt->mMinInversionRpt, "min inversion size to report", true)->group("Threshold Options");
    app.add_option("--min_del_rpt", opt->filterOpt->mMinDeletionRpt, "min deletion size to report", true)->group("Threshold Options");
    app.add_option("--min_dup_rpt", opt->filterOpt->mMinDupRpt, "min dup size to report", true)->group("Threshold Options");
    app.add_option("--min_seed_sr", opt->filterOpt->mMinSeedSR, "min seed split reads needed to compute SV", true)->check(CLI::Range(1, 100000))->group("Threshold Options");
    app.add_option("--min_seed_dp", opt->filterOpt->mMinSeedDP, "min seed discordant reads needed to computr SV", true)->check(CLI::Range(1, 100000))->group("Threshold Options");
    // Fusion report options
    app.add_option("--whitemindep", opt->fuseOpt->mWhiteFilter.mMinDepth, "min depth for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmindep", opt->fuseOpt->mUsualFilter.mMinDepth, "min depth for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminsrs", opt->fuseOpt->mWhiteFilter.mMinSRSeed, "min sr seeds for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminsrs", opt->fuseOpt->mUsualFilter.mMinSRSeed, "min sr seeds for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whitemindps", opt->fuseOpt->mWhiteFilter.mMinDPSeed, "min dp seeds for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmindps", opt->fuseOpt->mUsualFilter.mMinDPSeed, "min dp seeds for an valid fusion break point not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminttr", opt->fuseOpt->mWhiteFilter.mMinSupport, "min total reads for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminttr", opt->fuseOpt->mUsualFilter.mMinSupport, "min total reads for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminaf", opt->fuseOpt->mWhiteFilter.mMinVAF, "min VAF for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminaf", opt->fuseOpt->mUsualFilter.mMinVAF, "min VAF for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminigs", opt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminigs", opt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whitemaxrph", opt->fuseOpt->mWhiteFilter.mMaxRepHit, "max hits on reference allowed for partial seqs of fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmaxrph", opt->fuseOpt->mUsualFilter.mMaxRepHit, "max hits on reference allowed for partial seqs of fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--maxbpoffset", opt->fuseOpt->mMaxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true)->group("Fusion Options");
    app.add_option("--bgbcf", opt->fuseOpt->mBgBCF, "background events BCF file")->group("Fusion Options");
    app.add_option("--whitelist", opt->fuseOpt->mWhiteList, "white list of fusion events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--blacklist", opt->fuseOpt->mBlackList, "black list of fusion events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--samegenel", opt->fuseOpt->mSameGeneSVList, "white list of gene with inner sv events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--extraanno", opt->fuseOpt->mExtraAnnoList, "extra annoation gene list")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--fsrptlist", opt->fuseOpt->mFsRptList, "report range list")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--genecrdlist", opt->fuseOpt->mGeneCrdList, "gene coord list")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--fusionrpt", opt->fuseOpt->mOutFile, "primary fusion report file path", true)->group("Fusion Options");
    app.add_option("--supplerpt", opt->fuseOpt->mSupFile, "supplementary fusion report file path", true)->group("Fusion Options");
    // parse arguments
    CLI_PARSE(app, argc, argv);
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
