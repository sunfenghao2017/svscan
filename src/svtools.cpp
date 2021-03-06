#include <CLI.hpp>
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
    app.get_formatter()->column_width(50);
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
    pfilter->add_option("-i,--in", filterOpt->infile, "input tsv sv result of svscan")->required(true)->check(CLI::ExistingFile);
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
    pfsrpt->add_option("-i,--in", fuserptOpt->fuseOpt->mInfile, "input tsv sv result of svscan")->required(true)->check(CLI::ExistingFile);
    pfsrpt->add_option("-o,--out", fuserptOpt->fuseOpt->mOutFile, "fusion result tsv file", true);
    pfsrpt->add_option("-v,--svf", fuserptOpt->fuseOpt->mSVModFile, "updated sv tsv file, if empty, no update will done", true);
    pfsrpt->add_option("-f,--iflag", fuserptOpt->fuseOpt->mFsMaskInclude, "fusion flags which must be set", true);
    pfsrpt->add_option("-F,--eflag", fuserptOpt->fuseOpt->mFsMaskExclude, "fusion flags which must be unset", true);
    pfsrpt->add_option("-g,--genome", fuserptOpt->fuseOpt->mRef, "reference used in alignment of BAM", true);
    pfsrpt->add_flag("-d,--debug", fuserptOpt->debug, "turn on debug mode");
    pfsrpt->add_option("--whitemindep", fuserptOpt->fuseOpt->mWhiteFilter.mMinDepth, "min depth for an valid fusion break point in whitelist", true);
    pfsrpt->add_option("--usualmindep", fuserptOpt->fuseOpt->mUsualFilter.mMinDepth, "min depth for an valid fusion break point ont in whitelist", true);
    pfsrpt->add_option("--whiteminsrs", fuserptOpt->fuseOpt->mWhiteFilter.mMinSRSeed, "min sr seeds for an valid fusion break point in whitelist", true);
    pfsrpt->add_option("--usualminsrs", fuserptOpt->fuseOpt->mUsualFilter.mMinSRSeed, "min sr seeds for an valid fusion break point not in whitelist", true);
    pfsrpt->add_option("--whitemindps", fuserptOpt->fuseOpt->mWhiteFilter.mMinDPSeed, "min dp seeds for an valid fusion break point in whitelist", true);
    pfsrpt->add_option("--usualmindps", fuserptOpt->fuseOpt->mUsualFilter.mMinDPSeed, "min dp seeds for an valid fusion break point not in whitelist", true);
    pfsrpt->add_option("--whiteminsrr", fuserptOpt->fuseOpt->mWhiteFilter.mMinSRSupport, "min sr reads for an valid fusion break point in whitelist", true);
    pfsrpt->add_option("--usualminsrr", fuserptOpt->fuseOpt->mUsualFilter.mMinSRSupport, "min sr reads for an valid fusion break point not in whitelist", true);
    pfsrpt->add_option("--whitemindpr", fuserptOpt->fuseOpt->mWhiteFilter.mMinDPSupport, "min dp reads for an valid fusion break point in whitelist", true);
    pfsrpt->add_option("--usualmindpr", fuserptOpt->fuseOpt->mUsualFilter.mMinDPSupport, "min dp reads for an valid fusion break point not in whitelist", true);
    pfsrpt->add_option("--whiteminttr", fuserptOpt->fuseOpt->mWhiteFilter.mMinSupport, "min molecule for an valid fusion in whitelist", true);
    pfsrpt->add_option("--usualminttr", fuserptOpt->fuseOpt->mUsualFilter.mMinSupport, "min molecule for an valid fusion not in whitelist", true);
    pfsrpt->add_option("--whiteminaf", fuserptOpt->fuseOpt->mWhiteFilter.mMinVAF, "min VAF for an valid fusion in whitelist", true);
    pfsrpt->add_option("--usualminaf", fuserptOpt->fuseOpt->mUsualFilter.mMinVAF, "min VAF for an valid fusion not in whitelist", true);
    pfsrpt->add_option("--whiteminigs", fuserptOpt->fuseOpt->mWhiteFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion in whitelist", true);
    pfsrpt->add_option("--usualminigs", fuserptOpt->fuseOpt->mUsualFilter.mMinIntraGeneSVSize, "min intra-gene sv size for an valid fusion not in whitelist", true);
    pfsrpt->add_option("--whitemaxrph", fuserptOpt->fuseOpt->mWhiteFilter.mMaxRepHit, "max hits on reference allowed for partial seqs of fusion in whitelist", true);
    pfsrpt->add_option("--usualmaxrph", fuserptOpt->fuseOpt->mUsualFilter.mMaxRepHit, "max hits on reference allowed for partial seqs of fusion not in whitelist", true);
    pfsrpt->add_option("--maxbpoffset", fuserptOpt->fuseOpt->mMaxBpOffset, "max breakpoint offset allowed for an SV excluded from background SVs", true);
    pfsrpt->add_option("--bgbcf", fuserptOpt->fuseOpt->mBgBCF, "background events BCF file");
    pfsrpt->add_option("--whitelist", fuserptOpt->fuseOpt->mWhiteList, "white list of fusion events")->check(CLI::ExistingFile);
    pfsrpt->add_option("--blacklist", fuserptOpt->fuseOpt->mBlackList, "black list of fusion events")->check(CLI::ExistingFile);
    pfsrpt->add_option("--samegenel", fuserptOpt->fuseOpt->mSameGeneSVList, "white list of gene with inner sv events")->check(CLI::ExistingFile);
    pfsrpt->add_option("--fsrptlist", fuserptOpt->fuseOpt->mFsRptList, "report range list")->check(CLI::ExistingFile);
    pfsrpt->add_option("--genecrdlist", fuserptOpt->fuseOpt->mGeneCrdList, "gene coord list")->check(CLI::ExistingFile);
    pfsrpt->add_option("--idpdropmask", fuserptOpt->fuseOpt->mIDBDropMask, "fusion in db drop bit mask", true);
    pfsrpt->add_option("--ndbdropmask", fuserptOpt->fuseOpt->mNDBDropMask, "fusion not in db drop bit mask", true);
    pfsrpt->add_option("--keepmasks", fuserptOpt->fuseOpt->mKeepMasks, "fusion keep masks", true);
    pfsrpt->add_option("--primarymask", fuserptOpt->fuseOpt->mPrimaryMask, "primary fusion mask", true);
    pfsrpt->add_option("--hotfspartner", fuserptOpt->fuseOpt->mHotPartnerList, "hot fusion partner exon list")->check(CLI::ExistingFile);
    // dna annodb
    SVDNADBOpt* annDBOpt = new SVDNADBOpt();
    CLI::App* panndb = app.add_subcommand("dnadb", "prepare DNA sv annotation database for svscan");
    panndb->add_option("-i,--in", annDBOpt->refSeqDB, "refseq database with transcript version added")->required(true)->check(CLI::ExistingFile);
    panndb->add_option("-m,--mtrs", annDBOpt->gene2cnc, "gene to canonical transcript list")->required(true)->check(CLI::ExistingFile);
    panndb->add_option("-a,--anno", annDBOpt->svAnnoDB, "output file path of sv annotation db", true);
    panndb->add_option("-r,--r2g", annDBOpt->ref2gene, "output file path of refseq2gene tsv", true);
    panndb->add_option("-c,--cdsbed", annDBOpt->bedCDS, "output file path of cds region", true);
    panndb->add_option("-u,--unitbed", annDBOpt->bedUnits, "output file path of exon/utr units", true);
    panndb->add_option("-g,--generng", annDBOpt->geneRng, "output file path of gene cnctrs coordinates", true);
    // rna annodb
    SVRNADBOpt* rnaDBOpt = new SVRNADBOpt();
    CLI::App* prnadb = app.add_subcommand("rnadb", "prepare RNA sv annotation database for svscan");
    prnadb->add_option("-i,--in", rnaDBOpt->refSeqDB, "refseq database with transcript version added")->required(true)->check(CLI::ExistingFile);
    prnadb->add_option("-m,--mtrs", rnaDBOpt->cncTrsList, "canonical transcript name list")->required(true)->check(CLI::ExistingFile);
    prnadb->add_option("-g,--genome", rnaDBOpt->alnref, "genome reference fasta file path")->required(true)->check(CLI::ExistingFile);
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
    pfsw->add_option("-O,--allgs", fwOpt->allglist, "all gene in fusion white list", true);
    // flags
    FuseFlags* fuseflagOpt = new FuseFlags();
    CLI::App* pflags = app.add_subcommand("flags", "show fusion flags defined by svscan");
    // svbam
    SVBAMOpt* svbamOpt = new SVBAMOpt();
    CLI::App* psvbam = app.add_subcommand("svbam", "get bam records supporting an sv");
    psvbam->add_option("-i,--ibam", svbamOpt->ibam, "bam generated by svscan")->required(true)->check(CLI::ExistingFile);
    psvbam->add_option("-s,--svid", svbamOpt->svid, "SV ID")->required(true);
    psvbam->add_option("-o,--obam", svbamOpt->obam, "output bam path", true);
    psvbam->add_option("-O,--obam2", svbamOpt->obam2, "output realign bam path");
    psvbam->add_option("-r,--ref", svbamOpt->ref, "reference used to realign");
    // bam2tb
    BamToTable* bam2tb = new BamToTable();
    CLI::App* pb2t = app.add_subcommand("bam2tb", "get bam table supportint fusion results");
    pb2t->add_option("-i,--svbam", bam2tb->svbam, "bam generated by svscan")->required(true)->check(CLI::ExistingFile);
    pb2t->add_option("-b,--obam", bam2tb->newbam, "fusion support bam", true)->required(false);
    pb2t->add_option("-f,--fstsv", bam2tb->fstsv, "fusion result tsv")->required(false)->check(CLI::ExistingFile);
    pb2t->add_option("-u,--usrid", bam2tb->usrid, "sv ids user intereseted")->required(false);
    pb2t->add_option("-g,--fsgene", bam2tb->fsgene, "fusion gene of each sv id")->required(false);
    pb2t->add_option("-c,--column", bam2tb->svidf, "svid column index in tsv[0-based]", true)->required(false);
    pb2t->add_option("-e,--outexcel", bam2tb->bamtb, "output excel path(ignore output if not provided)", true)->required(false);
    pb2t->add_flag("-r,--refinedp", bam2tb->refinedp, "refine dp if true");
    pb2t->add_option("-v,--outtsv", bam2tb->bamtt, "output tsv path(ignore output if not provided)", true)->required(false);
    // svcreg
    SVCreg* crgOpt = new SVCreg();
    CLI::App* psvcrg = app.add_subcommand("svcreg", "get creg for svscan");
    psvcrg->add_option("-i,--ibed", crgOpt->ibed, "input bed file")->required(true)->check(CLI::ExistingFile);
    psvcrg->add_option("-g,--glist", crgOpt->glist, "gene list file")->required(true)->check(CLI::ExistingFile);
    psvcrg->add_option("-o,--obed", crgOpt->obed, "output bed file", true);
    // debug
    SVDebug* debugOpt = new SVDebug();
    CLI::App* pdebug = app.add_subcommand("debug", "debug fusion candidate");
    pdebug->add_option("-a,--anndb", debugOpt->annodb, "annodb generated by svtools rnadb")->check(CLI::ExistingFile);
    pdebug->add_option("-g,--g2tf", debugOpt->gene2trans, "gene to transcript table generated by svtools rnadb")->required(true)->check(CLI::ExistingFile);
    pdebug->add_option("-i,--ibam", debugOpt->inbam, "input bam to debug")->required(true)->check(CLI::ExistingFile);
    CLI::Option* pdgl = pdebug->add_option("-s,--fglist", debugOpt->fglist, "fusion pair list")->required(false)->check(CLI::ExistingFile);
    pdebug->add_option("-f,--hgs", debugOpt->hgl, "hgenes in fusion")->required(false)->excludes(pdgl);
    pdebug->add_option("-l,--tgs", debugOpt->tgl, "tgenes in fusion")->required(false)->excludes(pdgl);
    pdebug->add_option("-b,--svbam", debugOpt->svbam, "sv support bam", true);
    pdebug->add_option("-v,--tsv", debugOpt->table, "stat table", true);
    pdebug->add_option("-t,--thread", debugOpt->nthread, "thread number", true);
    pdebug->add_flag("-r,--rna", debugOpt->rnamode, "sv event is called from rnaseq bam");
    // rsvg
    RSVGen* rsvgOpt = new RSVGen();
    CLI::App* prsvg = app.add_subcommand("rsvg", "generate rna fuse events");
    prsvg->add_option("-a,--anndb", rsvgOpt->annrna, "annodb generated by svtools rnadb")->required(true)->check(CLI::ExistingFile);
    prsvg->add_option("-g,--g2tf", rsvgOpt->g2tf, "gene to transcript table generated by svtools rnadb")->required(true)->check(CLI::ExistingFile);
    prsvg->add_option("-r,--refrna", rsvgOpt->refrna, "reference seq fa of rna corresponding to anndb")->required(true)->check(CLI::ExistingFile);
    prsvg->add_option("-i,--incfg", rsvgOpt->incfg, "input fusion gene info")->required(true)->check(CLI::ExistingFile);
    prsvg->add_option("-o,--outcfg", rsvgOpt->outcfg, "output fusion construct info", true);
    prsvg->add_option("-f,--outfa", rsvgOpt->outfa, "output fasta file of fusion seqs", true);
    // scope rpt
    ScopeRptOpt* scRptOpt = new ScopeRptOpt();
    CLI::App* pscrpt = app.add_subcommand("scrpt", "generate scope report");
    pscrpt->add_option("-i,--infs", scRptOpt->infs, "input fusion result")->required(true)->check(CLI::ExistingFile);
    pscrpt->add_option("-o,--outfs", scRptOpt->outfs, "output fusion result")->required(true);
    pscrpt->add_flag("-r,--rnamode", scRptOpt->fromrna, "from rna seq if set");
    // prod rpt
    ProdRptOpt* prRptOpt = new ProdRptOpt();
    CLI::App* pprrpt = app.add_subcommand("prrpt", "generate prod report");
    pprrpt->add_option("-i,--infs", prRptOpt->infs, "input fusion result")->required(true)->check(CLI::ExistingFile);
    pprrpt->add_option("-p,--outfp", prRptOpt->outfp, "output fusion result for product department")->required(true);
    pprrpt->add_option("-m,--outfm", prRptOpt->outfm, "output fusion result for medcine department")->required(true);
    pprrpt->add_option("-g,--hgl", prRptOpt->hlist, "hot gene list")->required(true);
    pprrpt->add_flag("-r,--rnamode", prRptOpt->fromrna, "from rna seq if set");
    // parse arguments
    CLI11_PARSE(app, argc, argv);
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
        svbamOpt->update(argc, argv);
        svbamOpt->getBam();
        delete svbamOpt;
    }
    // bam2tb
    if(pb2t->parsed()){
        bam2tb->b2t();
        delete bam2tb;
    }
    // svcreg
    if(psvcrg->parsed()){
        crgOpt->getCreg();
        delete crgOpt;
    }
    // debug
    if(pdebug->parsed()){
        if(debugOpt->hgl.size() != debugOpt->tgl.size()){
            util::errorExit("hgene list and tgene list lengths differ");
        }
        if(pdgl->count() == 0 && debugOpt->hgl.size() == 0){
            util::errorExit("either fuselist or h/t gene pairs should be provided");
        }
        debugOpt->init();
        if(debugOpt->rnamode){
            debugOpt->debugRNA();
        }else{
            if(debugOpt->annodb.empty()){
                util::errorExit("annodb must be provided in DNA SV debugging");
            }
            debugOpt->debugDNA();
        }
        delete debugOpt;
    }
    // rsvg
    if(prsvg->parsed()){
        rsvgOpt->init();
        rsvgOpt->gensv();
    }
    // scope rpt
    if(pscrpt->parsed()){
        scRptOpt->out();
    }
    // prod rpt
    if(pprrpt->parsed()){
        prRptOpt->get_hot();
        prRptOpt->out();
    }
}
