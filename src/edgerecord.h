#ifndef EDGERECORD_H
#define EDGERECORD_H

#include <map>
#include <cstdint>
#include <iostream>

/** class to store cluster relationship bewteen nodes\n
 * a node can be an SRBamrecord or an DPBamrecord
 */
struct EdgeRecord{
    int32_t mSource; ///< index in node vector
    int32_t mTarget; ///< index in node vector
    int32_t mWeight; ///< difference of source and target node

    /** EdgeRecord constructor
     * @param source index of source node in vector
     * @param target index of target node in vector
     * @param weight weight of target and source node
     */
    EdgeRecord(int32_t source, int32_t target, int32_t weight){
        mSource = source;
        mTarget = target;
        mWeight = weight;
    }

    /** EdgeRecord destructor */
    ~EdgeRecord(){}

    /** operator to output an EdgeRecord to ostream
     * @param os reference of ostream
     * @param e reference of EdgeRecord
     * @return reference of ostream
     */
    inline friend std::ostream& operator<<(std::ostream& os, const EdgeRecord& e){
        os << "Source: " << e.mSource << "\n";
        os << "Target: " << e.mTarget << "\n";
        os << "Weight: " << e.mWeight << "\n";
        return os;
    }

    /** operator to compare two EdgeRecord
     * @param other reference of another EdgeRecord
     * @return true if this < that
     */
    inline bool operator<(const EdgeRecord& other) const {
        return mWeight < other.mWeight ||
               (mWeight == other.mWeight && mSource < other.mSource) ||
               (mWeight == other.mWeight && mSource == other.mSource && mTarget < other.mTarget);
    }
};

/** type to store an raw cluster of SR/DP */
typedef std::map<int32_t, std::vector<EdgeRecord>> Cluster;

#endif
