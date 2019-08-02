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
    app.add_option("-b,--bam", opt->bamfile, "bam file")->required(true)->check(CLI::ExistingFile)->group("General");
    app.add_option("-g,--genome", opt->genome, "reference genome")->required(true)->check(CLI::ExistingFile)->group("General");
    app.add_option("-o,--outfile", opt->bcfOut, "output bcf file", true)->required(false)->group("General");
    app.add_option("-c,--madcutoff", opt->madCutoff, "isize mad cutoff for Del deletion", true)->group("General");
    app.add_option("-e,--excreg", opt->excludeReg, "invalid region bed file", true)->check(CLI::ExistingFile)->group("General");
    app.add_option("-s,--svtype", opt->svtypes, "SV types to discover,0:INV,1:DEL,2:DUP,3:INS,4:BND")->check(CLI::Range(0, 4))->group("General");
    app.add_option("-n,--nthread", opt->nthread, "number of threads used to process bam", true)->check(CLI::Range(1, 20))->group("General");
    CLI_PARSE(app, argc, argv);
    // validate arguments
    util::loginfo("Command line arguments parsed");
    opt->validate();
    util::loginfo("Command line arguments validated");
    // update some arguments
    opt->update(argc, argv);
    util::loginfo("Options updated");
    // output library basic information
    util::loginfo("Library basic information:\n" + opt->libInfo->toStr());
    // get Valid regions
    opt->getValidRegion();
    util::loginfo("Valid region parsed");
    // calling SV
    SVScanner* svScanner = new SVScanner(opt);
    util::loginfo("Start scanning SV");
    svScanner->scanDPandSR();
    util::loginfo("Finish scanning SV");
}
