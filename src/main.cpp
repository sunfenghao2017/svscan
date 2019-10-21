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
    CLI::App app("program: " + std::string(argv[0]) + "\n" + opt->softEnv->cmp);
    app.get_formatter()->column_width(36);
    // General Options options
    app.add_option("-b,--bam", opt->bamfile, "bam file")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-g,--genome", opt->genome, "reference genome/transcriptome")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-a,--anno", opt->annodb, "annotation database file")->required(true)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-r,--reg", opt->reg, "valid region to disvover SV")->required(false)->check(CLI::ExistingFile)->group("General Options");
    app.add_option("-o,--bcfout", opt->bcfOut, "output sv bcf file", true)->required(false)->group("General Options");
    app.add_option("-t,--tsvout", opt->tsvOut, "output sv tsv file", true)->required(false)->group("General Options");
    app.add_option("-v,--bamout", opt->bamout, "output sv bam file, non-providing will disable it", true)->required(false)->group("General Options");
    app.add_option("-s,--svtype", opt->svtypes, "SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND")->check(CLI::Range(0, 4))->group("General Options");
    app.add_option("-n,--nthread", opt->nthread, "number of threads used to process bam", true)->check(CLI::Range(1, 100))->group("General Options");
    app.add_flag("--rna", opt->rnamode, "discovery structural variants from rna data")->group("General Options");
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
    // Fusion report options
    app.add_option("--whitemindep", opt->fuseOpt->mWhiteFilter.mMinDepth, "min depth for an valid fusion break point in whitelist", true)->group("Fusion Options");
    app.add_option("--usualmindep", opt->fuseOpt->mUsualFilter.mMinDepth, "min depth for an valid fusion break point ont in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminr", opt->fuseOpt->mWhiteFilter.mMinSupport, "min reads support for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminr", opt->fuseOpt->mUsualFilter.mMinSupport, "min reads support for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminaf", opt->fuseOpt->mWhiteFilter.mMinVAF, "min VAF for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminaf", opt->fuseOpt->mUsualFilter.mMinVAF, "min VAF for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--whiteminigs", opt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion in whitelist", true)->group("Fusion Options");
    app.add_option("--usualminigs", opt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion not in whitelist", true)->group("Fusion Options");
    app.add_option("--maxbpoffset", opt->fuseOpt->mMaxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true)->group("Fusion Options");
    app.add_option("--bgbcf", opt->fuseOpt->mBgBCF, "background events BCF file")->group("Fusion Options");
    app.add_option("--whitelist", opt->fuseOpt->mWhiteList, "white list of fusion events")->check(CLI::ExistingFile)->group("Fusion Options");
    app.add_option("--blacklist", opt->fuseOpt->mBlackList, "black list of fusion events")->check(CLI::ExistingFile)->group("Fusion Options");
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
    // output library basic information
    util::loginfo("Library basic information:\n" + opt->libInfo->toStr());
    // get Valid regions
    opt->getValidRegion();
    util::loginfo("Valid region parsed");
    // calling SV
    SVScanner* svScanner = new SVScanner(opt);
    util::loginfo("Beg scanning SV");
    svScanner->scanDPandSR();
    util::loginfo("End scanning SV");
    // write execution time
    util::loginfo("Time consumed: " + opt->softEnv->getExecutionTime());
}
