#ifndef REGION_H
#define REGION_H

#include <htslib/sam.h>
#include <string>
#include <vector>
#include <stack>
#include <set>
#include "util.h"

/** class used to store regions of genome */
class Region{
    public:
        std::vector<std::set<std::pair<int32_t, int32_t>>> mRegs; ///< regions in each contig
    public:
        /** Region constructor */
        Region(){}

        /** Region constructor
         * @param regfile region bed file
         * @param bamfile bam file covering this region
         */
        Region(const std::string& regfile, const std::string& bamfile);

        /** merge and sort regions of one contig
         * @param regs regions of one contig
         * @param merged and sorted regions
         */
        static std::vector<std::pair<int32_t, int32_t>> mergeAndSortRegions(std::vector<std::pair<int32_t, int32_t>>& regs);

        /** get first overlaped region of provided region
         * @param tid region contig id
         * @param p region to be found
         * @return iter to the first overlapped region or end
         */
        std::set<std::pair<int32_t, int32_t>>::iterator getFirstOverlap(int32_t tid, const std::pair<int32_t, int32_t>& p);
};

#endif
