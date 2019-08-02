#ifndef ANNOTATOR_H
#define ANNOTATOR_H

#include "stats.h"

/** SV breakpoint coverage annotator */
class Annotator{
    public:
        Options* mOpt; ///< pointer to Options

    public:
        /** constructor of Annotator
         * @param opt pointer to Options object
         */
        Annotator(Options* opt) : mOpt(opt) {}
        
        /** destructor of Annotator */
        ~Annotator(){}

        /** annotate SV coverage
         * @param svs reference of SVRecords
         */
        Stats* covAnnotate(SVSet& svs);

        /** get first overlap of an region against an region set
         * @param s region sets
         * @param p region to query
         * @return iterator to an region overlap with p or end
         */
        static std::set<std::pair<int32_t, int32_t>>::iterator getFirstOverlap(const std::set<std::pair<int32_t, int32_t>>& s, const std::pair<int32_t, int32_t>& p);
};

#endif
