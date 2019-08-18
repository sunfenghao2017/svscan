#ifndef MSA_H
#define MSA_H

#include <set>
#include "aligner.h"
#include "matrix2d.h"
#include "aligncfg.h"

/** class to do multiple sequence alignment */
class MSA{
    public:
        std::multiset<std::string>* mSeqs = NULL; ///< sequences to do msa
        AlignConfig* mAlignConfig = NULL;         ///< alignment strategy used
        bool mDefaultConfigCreated = false;       ///< default mAlignConfig constructed if true
        int32_t mMinCovForCS = 3;                 ///< minimum coverage needed for a position in msa result to be included in consensus sequence
        float mMinBaseRatioForCS = 0.5;           ///< minimum base ratio needed for a position in msa result to be included in consensus sequence

    public:
        /** MSA constructor
         * @param seqs pointer to a set of sequences to do msa
         * @param minCovForCS minimum coverage needed for a position in msa result to be included in consensus sequence
         * @param minBaseRatioForCS minimum base ratio needed for a position in msa result to be included in consensus sequence
         * @param alignCfg alignment strategy used
         */
        MSA(std::multiset<std::string>* seqs, int32_t minCovForCS = 3, float minBaseRatioForCS = 0.5, AlignConfig* alignCfg = NULL){
            mSeqs = seqs;
            if(alignCfg) mAlignConfig = alignCfg;
            else{
                mAlignConfig = new AlignConfig();
                mDefaultConfigCreated = true;
            }
            mMinCovForCS = minCovForCS;
            mMinBaseRatioForCS = minBaseRatioForCS;
        }

        /** MSA destructor */
        ~MSA(){
            if(mDefaultConfigCreated){
                delete mAlignConfig;
            }
        }

    public:
        /** get longest common sequence of two sequence s1 and s2(gaps allowed)
         * @param s1 sequence
         * @param s2 sequence
         * @return longest common sequence length of s1 and s2
         */
        static int lcs(const std::string& s1, const std::string& s2);
        
        /** construct distance matrix of mSeqs, d is square matrix with (2 * mSeqs.size() + 1) rows\n
         * initially d[i, j] = 0 for i >= j, d[i, j] = -1 otherwise\n
         * in distanceMatrix function d[i, j] = 100 * lcs(seqi, seqj)/ min(seqi.size(), seqj.size()) for i < j( i < mSeqs.size(), j < mSeqs.size())\n
         * @param d matrix to store pair-wise distance of mSeqs
         */
        void distanceMatrix(Matrix2D<int>* d);

        /** get a pair of (clustered) sequences who has the minimal distance in d[0..num-1, 0...num-1]\n
         * minimal distance is the max d[i, j] with i < num and j < num\n
         * @param d distance matrix of mSeqs and all clustered pairs
         * @param num max columns/rows below which to get the closest pair
         * @param di store i of closest pair(i, j) in d
         * @param dj store j of closest pair(i, j) in d
         * @return max value of d[i, j] with i < num and j < num
         */
        int closestPair(Matrix2D<int>* d, int num, int& di, int& dj);

        /** update distance matrix with a closest pair (di, dj) found in d[0...num-1, 0...num-1]\n
         * set any seq k that haven't been clustered in mSeqs distance to seqnum equals 0.5 * (distance(di, k) + distance(dj, k))\n
         * that is d[k, num] = 0.5 * (distance(di, k) + distance(dj, k)) if p[k, 0] = -1\n
         * set column and row in d contains of di and dj to -1\n
         * @param d distance matrix whose closest pair in d[0...num-1, 0...num-1] is (di, dj)
         * @param p track matrix with dimension (2 * mSeqs.size() + 1, 3), initialized to -1 in all fields\n
         * p[i, 0] and p[j, 0] is set to num if closest pair (di, dj) found in d[0...num-1, 0...num-1]\n
         * p[num, 1] = i and p[num, 2] = j if closest pair (di, dj) found in d[0...num-1, 0...num-1]\n
         * @param num a number wich closest pair (di, dj) found in d[0...num-1, 0...num-1]
         * @param di i of closest pair(i, j) in d[0...num-1, 0...num-1]
         * @param dj j of closest pair(i, j) in d[0...num-1, 0...num-1]
         */
        void updateDistanceMatrix(Matrix2D<int>* d, const Matrix2D<int>* p, int num, int di, int dj);

        /** UPGMA (unweighted pair group method with arithmetic mean) \n
         * is a simple agglomerative (bottom-up) hierarchical clustering method.\n
         * The method is generally attributed to Sokal and Michener.\n
         * n is set to mSeqs.size() firstly and increased to 2 * num\n
         * for each n, find closest pair (di, dj) in d[0...n-1, 0...n-1]\n
         * after found, update distance matrix d use updateDistanceMatrix\n
         * which is mainly clear all distance with di and dj\n
         * and set new distance of all seqs(expcept di, dj) to clustered seq and store in column n in d\n
         * p[di, 0] = n, p[dj, 0] = n, to store clustered seqs and their clustered seq number marker\n
         * p[n, 1] = di, p[n, 2] = dj, to store clustered seq n comes from seq/cluster di and dj\n
         * @param d distance matrix
         * @param p phylogenetic matrix used in UPGMA
         * @param num find closest pair (di, dj) found in d[0...num-1, 0...num-1]
         * @return root of phylogenetic matrix p, just the last column number in p to store the final clustered seq distances
         */
        int upgma(Matrix2D<int>* d, Matrix2D<int>* p, int num);

        /** prograssive alignment of sequences in mSeqs base on UPGMA matrix and root(last clustered seq num in d)
         * from root of phylogenetic matrix p and track back to leaf of p[root, 1] and p[root, 2]  recursively\n
         * use gotoh global pairwise alignment to alignment p[root, 1] and p[root, 2] from bottom up to get final msa result\n
         * @param p phylogenetic matrix used in UPGMA
         * @param root root of phylogenetic matrix p, just the last column number in p to store the final clustered seq distances
         * @param aln matrix to store msa result
         */
        void palign(const Matrix2D<int>* p, int root, Matrix2D<char>* aln);

        /** get consensus string from msa result that satisfy coverage and majority base rate limits\n
         * @param aln matrix to store msa alignment result
         * @param cs consensus string got from msa alignment result
         */
        void consensus(Matrix2D<char>* aln, std::string& cs);
        
        /** do multiple seqeuence alignment of mSeqs\n
         * @param cs consensus sequence got from msa
         * @return consensus sequence supporting seq number
         */
        int msa(std::string& cs);

};

#endif
