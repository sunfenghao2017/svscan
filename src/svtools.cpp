#include "CLI.hpp"
#include "svtools.h"

int main(int argc, char** argv){
    if(argc == 1){
        std::string helpCMD = std::string(argv[0]) + " -h";
        std::system(helpCMD.c_str());
        return 0;
    }
    SVToolOpts* svtopt = new SVToolOpts();
    svtopt->update(argc, argv);
    CLI::App app("program: " + std::string(argv[0]) + "\n" + svtopt->softEnv->getSoftInfo());
    app.get_formatter()->column_width(36);
    // merge
    VCFMerger* mergeOpt = new VCFMerger();
    CLI::App* pmerge = app.add_subcommand("merge", "merge sv bcfs generated by svscan");
    pmerge->add_option("-i,--in", mergeOpt->infilelist, "input bcf file list to be merged")->required(true)->check(CLI::ExistingFile);
    pmerge->add_option("-o,--out", mergeOpt->outfile, "merged bcf file output path", true);
    pmerge->add_option("-t,--tmp", mergeOpt->tmpdir, "temp directory to store intermediate files", true);
    pmerge->add_option("-v,--vaf", mergeOpt->vaf, "min VAF for SV to be merged", true);
    pmerge->add_option("-c,--cov", mergeOpt->coverage, "min depth for SV to be merged", true);
    pmerge->add_option("--covratio", mergeOpt->recoverlap, "min overlap ratio needed for two dup SVs", true);
    pmerge->add_option("--bpoffset", mergeOpt->bpoffset, "max position offset allowed for two dup SV", true);
    pmerge->add_option("--minsize", mergeOpt->minsize, "min size of SV to be merged", true);
    pmerge->add_option("--maxsize", mergeOpt->maxsize, "max size of SV to be merged", true);
    pmerge->add_option("--chunksize", mergeOpt->chunksize, "max input files to be merged at one batch", true);
    pmerge->add_flag("--passfilter", mergeOpt->filterForPass, "keep PASS records only if set");
    pmerge->add_flag("--preciseFilter", mergeOpt->filterForPrecise, "keep PRECISE records only if set");
    // filter
    SVFilter* filterOpt = new SVFilter();
    CLI::App* pfilter = app.add_subcommand("filter", "filter sv tsvs generated by svscan");
    pfilter->add_option("-i,--in", filterOpt->infile, "input tsv sv result of sver")->required(true)->check(CLI::ExistingFile);
    pfilter->add_option("-o,--out", filterOpt->outfile, "output tsv fusion result", true);
    pfilter->add_option("--minsr", filterOpt->minsr, "min SR support for an valid fusion", true);
    pfilter->add_option("--mindp", filterOpt->mindp, "min DP support for an valid fusion", true);
    pfilter->add_option("--minsu", filterOpt->minsu, "min SU support for an valid fusion", true);
    pfilter->add_option("--minaf", filterOpt->minaf, "min VAF for an valid fusion", true);
    pfilter->add_option("--maxbpoffset", filterOpt->maxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true);
    pfilter->add_option("--bgbcf", filterOpt->mBgBCF, "background events BCF file");
    // fuserpt
    FusionReporter* fuserptOpt = new FusionReporter();
    CLI::App* pfsrpt = app.add_subcommand("fuserpt", "report fusion event from sv tsv generate by svscan");
    pfsrpt->add_option("-i,--in", fuserptOpt->fuseOpt->mInfile, "input tsv sv result of sver")->required(true)->check(CLI::ExistingFile);
    pfsrpt->add_option("-o,--out", fuserptOpt->fuseOpt->mOutFile, "primary fusion result", true);
    pfsrpt->add_option("-s,--sup", fuserptOpt->fuseOpt->mSupFile, "supplementary fusion result", true);
    pfsrpt->add_option("-v,--svf", fuserptOpt->fuseOpt->mSVModFile, "updated sv tsv file, if empty, no update will done", true);
    pfsrpt->add_option("-f,--iflag", fuserptOpt->fuseOpt->mFsMaskInclude, "fusion flags which must be set", true);
    pfsrpt->add_option("-F,--eflag", fuserptOpt->fuseOpt->mFsMaskExclude, "fusion flags which must be unset", true);
    pfsrpt->add_option("--whitemindep", fuserptOpt->fuseOpt->mWhiteFilter.mMinDepth, "min depth for an valid fusion break point in whitelist", true);
    pfsrpt->add_option("--usualmindep", fuserptOpt->fuseOpt->mUsualFilter.mMinDepth, "min depth for an valid fusion break point ont in whitelist", true);
    pfsrpt->add_option("--whiteminr", fuserptOpt->fuseOpt->mWhiteFilter.mMinSupport, "min reads support for an valid fusion in whitelist", true);
    pfsrpt->add_option("--usualminr", fuserptOpt->fuseOpt->mUsualFilter.mMinSupport, "min reads support for an valid fusion not in whitelist", true);
    pfsrpt->add_option("--whiteminaf", fuserptOpt->fuseOpt->mWhiteFilter.mMinVAF, "min VAF for an valid fusion in whitelist", true);
    pfsrpt->add_option("--usualminaf", fuserptOpt->fuseOpt->mUsualFilter.mMinVAF, "min VAF for an valid fusion not in whitelist", true);
    pfsrpt->add_option("--whiteminigs", fuserptOpt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion in whitelist", true);
    pfsrpt->add_option("--usualminigs", fuserptOpt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion not in whitelist", true);
    pfsrpt->add_option("--maxbpoffset", fuserptOpt->fuseOpt->mMaxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true);
    pfsrpt->add_option("--bgbcf", fuserptOpt->fuseOpt->mBgBCF, "background events BCF file");
    pfsrpt->add_option("--whitelist", fuserptOpt->fuseOpt->mWhiteList, "white list of fusion events")->check(CLI::ExistingFile);
    pfsrpt->add_option("--blacklist", fuserptOpt->fuseOpt->mBlackList, "black list of fusion events")->check(CLI::ExistingFile);
    // dna annodb
    SVDNADBOpt* annDBOpt = new SVDNADBOpt();
    CLI::App* panndb = app.add_subcommand("dnadb", "prepare DNA sv annotation database for svscan");
    panndb->add_option("-i,--in", annDBOpt->refSeqDB, "refseq database with transcript version added")->required(true)->check(CLI::ExistingFile);
    panndb->add_option("-m,--mtrs", annDBOpt->gene2cnc, "gene to canonical transcript list")->required(true)->check(CLI::ExistingFile);
    panndb->add_option("-a,--anno", annDBOpt->svAnnoDB, "output file path of sv annotation db", true);
    panndb->add_option("-r,--r2g", annDBOpt->ref2gene, "output file path of refseq2gene tsv", true);
    panndb->add_option("-c,--cdsbed", annDBOpt->bedCDS, "output file path of cds region", true);
    panndb->add_option("-u,--unitbed", annDBOpt->bedUnits, "output file path of exon/utr units", true);
    // rna annodb
    SVRNADBOpt* rnaDBOpt = new SVRNADBOpt();
    CLI::App* prnadb = app.add_subcommand("rnadb", "prepare RNA sv annotation database for svscan");
    prnadb->add_option("-i,--in", rnaDBOpt->refSeqDB, "refseq database with transcript version added")->required(true)->check(CLI::ExistingFile);
    prnadb->add_option("-m,--mtrs", rnaDBOpt->cncTrsList, "canonical transcript name list")->required(true)->check(CLI::ExistingFile);
    prnadb->add_option("-g,--genome", rnaDBOpt->genome, "genome reference fasta file path")->required(true)->check(CLI::ExistingFile);
    prnadb->add_option("-a,--anno", rnaDBOpt->svAnnoDB, "output file path of sv annotation db", true);
    prnadb->add_option("-t,--gene2trs", rnaDBOpt->gene2trs, "output file path of transcript used by each gene", true);
    prnadb->add_option("-r,--refmrna", rnaDBOpt->refMrna, "output file path of refmrna fasta", true);
    prnadb->add_option("-c,--cdsbed", rnaDBOpt->bedCDS, "output file path of cds bed region", true);
    prnadb->add_option("-u,--unitbed", rnaDBOpt->bedUnits, "output file path of exon/utr units", true);
    prnadb->add_option("-l,--utr3len", rnaDBOpt->utr3len, "outut file path of 3'utr length", true);
    prnadb->add_flag("-n,--keepncrna", rnaDBOpt->keepNCRna, "keep ncrna gene in final results");
    //wlist
    FuseWOpt* fwOpt = new FuseWOpt();
    CLI::App* pfsw = app.add_subcommand("wlist", "prepare whitelist for svscan");
    pfsw->add_option("-g,--glist", fwOpt->genelist, "gene and common fuse direction list")->required(true)->check(CLI::ExistingFile);
    pfsw->add_option("-f,--fsdb", fwOpt->fusedb, "fusion pair database from public database")->required(true)->check(CLI::ExistingFile);
    pfsw->add_option("-o,--wlist", fwOpt->whitelist, "fusion white list usef for svscan", true);
    // flags
    FuseFlags* fuseflagOpt = new FuseFlags();
    CLI::App* pflags = app.add_subcommand("flags", "show fusion flags defined by svscan");
    // svbam
    SVBAMOpt* svbamOpt = new SVBAMOpt();
    CLI::App* psvbam = app.add_subcommand("svbam", "get bam records supporting an sv");
    psvbam->add_option("-i,--ibam", svbamOpt->ibam, "bam generated by svscan")->required(true)->check(CLI::ExistingFile);
    psvbam->add_option("-s,--svid", svbamOpt->svid, "SV ID")->required(true);
    psvbam->add_option("-o,--obam", svbamOpt->obam, "output bam path", true);
    // debug
    SVDebug* debugOpt = new SVDebug();
    CLI::App* pdebug = app.add_subcommand("debug", "debug fusion candidate");
    pdebug->add_option("-a,--anndb", debugOpt->annodb, "annodb generated by svtools rnadb")->check(CLI::ExistingFile);
    pdebug->add_option("-g,--g2tf", debugOpt->gene2trans, "gene to transcript table generated by svtools rnadb")->required(true)->check(CLI::ExistingFile);
    pdebug->add_option("-i,--ibam", debugOpt->inbam, "input bam to debug")->required(true)->check(CLI::ExistingFile);
    pdebug->add_option("-f,--hgene", debugOpt->hgene, "hgene in fusion")->required(true);
    pdebug->add_option("-l,--tgene", debugOpt->tgene, "tgene in fusion")->required(true);
    pdebug->add_option("-s,--srbam", debugOpt->srbam, "split read bam", true);
    pdebug->add_option("-o,--osbam", debugOpt->osbam, "other split read bam", true);
    pdebug->add_option("-d,--dpbam", debugOpt->dpbam, "discordant read pair bam", true);
    pdebug->add_option("-v,--tsv", debugOpt->table, "stat table", true);
    pdebug->add_flag("-r,--rna", debugOpt->rnamode, "sv event is called from rnaseq bam");
    // parse arguments
    CLI_PARSE(app, argc, argv);
    // merge
    if(pmerge->parsed()){
        mergeOpt->update(argc, argv);
        mergeOpt->mergeAllType();
        delete mergeOpt;
    }
    // filter
    if(pfilter->parsed()){
        filterOpt->update(argc, argv);
        filterOpt->filter();
        delete filterOpt;
    }
    // flags
    if(pflags->parsed()){
        fuseflagOpt->showFuseFlags();
        delete fuseflagOpt;
    }
    // fuserpt
    if(pfsrpt->parsed()){
        fuserptOpt->update(argc, argv);
        fuserptOpt->report();
        delete fuserptOpt;
    }
    // dna annodb
    if(panndb->parsed()){
        annDBOpt->update(argc, argv);
        annDBOpt->prepDB();
        delete annDBOpt;
    }
    // rna annodb
    if(prnadb->parsed()){
        rnaDBOpt->update(argc, argv);
        rnaDBOpt->prepDB();
        delete rnaDBOpt;
    }
    // wlist
    if(pfsw->parsed()){
        fwOpt->update(argc, argv);
        fwOpt->prepWlist();
        delete fwOpt;
    }
    // svbam
    if(psvbam->parsed()){
        svbamOpt->getBam();
        delete svbamOpt;
    }
    // debug
    if(pdebug->parsed()){
        if(debugOpt->rnamode){
            debugOpt->debugRNA();
        }else{
            debugOpt->debugDNA();
        }
        delete debugOpt;
    }
}
