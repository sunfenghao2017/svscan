#ifndef ALIGNER_H
#define ALIGNER_H

#include "util.h"
#include "matrix2d.h"
#include "aligncfg.h"

/** class to do global alignment of two sequnces */
class Aligner{
    public:
        Matrix2D<char>* mSeqVertical = NULL;   ///< vertical sequence matrix2d representation
        Matrix2D<char>* mSeqHorizontal = NULL; ///< horizoncal sequence matrix2d representation
        AlignConfig* mAlignConfig = NULL;      ///< alignment configuration
        bool mSeqMatrixCreated = false;        ///< mSeqVertical and mSeqHorizontal are created automatically
        bool mDefaultConfigCreated = false;    ///< default AlignConfig object constructed if true
    
    public:
        /** default Aligner constructor */
        Aligner(){}

        /** Aligner constructor
         * @param seqVertical vertical sequence matrix2d representation
         * @param seqHorizontal horizoncal sequence matrix2d representation
         * @param alignConfig alignment configuration object
         */
        Aligner(Matrix2D<char>* seqVertical, Matrix2D<char>* seqHorizontal, AlignConfig* alignConfig = NULL){
            mSeqVertical = seqVertical;
            mSeqHorizontal = seqHorizontal;
            if(alignConfig) mAlignConfig = alignConfig;
            else{
                mAlignConfig = new AlignConfig();
                mDefaultConfigCreated = true;
            }
        }
        
        /** Aligner constructor
         * @param seqVertical vertical sequence
         * @param seqHorizontal horizoncal sequence
         * @param alignConfig alignment configuration object
         */
        Aligner(const std::string& seqVertical, const std::string& seqHorizontal, AlignConfig* alignConfig = NULL){
            mSeqVertical = new Matrix2D<char>();
            mSeqHorizontal = new Matrix2D<char>();
            mSeqVertical->createProfile(seqVertical);
            mSeqHorizontal->createProfile(seqHorizontal);
            mSeqMatrixCreated = true;
            if(alignConfig) mAlignConfig = alignConfig;
            else{
                mAlignConfig = new AlignConfig();
                mDefaultConfigCreated = true;
            }
        }

        /** Aligner destructor */
        ~Aligner(){
            if(mAlignConfig && mDefaultConfigCreated){
                delete mAlignConfig;
                mAlignConfig = NULL;
            }
            if(mSeqHorizontal && mSeqMatrixCreated){
                delete mSeqHorizontal;
                mSeqHorizontal = NULL;
            }
            if(mSeqVertical && mSeqMatrixCreated){
                delete mSeqVertical;
                mSeqVertical = NULL;
            }
        }

    public:
        /** get match score of two (set of) sequence at some position
         * @param p1 profile of mSeqVertical
         * @param p2 profile of mSeqHorizontal
         * @param row index of mSeqVertical at which to match
         * @param col index of mSeqHorizontal at which to match
         * @return match score of a1[row] and a2[col]
         */
        template<typename T>
        int score(const Matrix2D<T>* p1, const Matrix2D<T>* p2, int row, int col);
        /******************************************intro*******************************************************
         * The algorithm by Osamu Gotoh (1982) computes the optimal global alignment of two sequences when    *
         * using an affine gap scoring. Here, the scoring of a long consecutive gap (insertion/deletion) is   *
         * favored over a collection of small gaps with the same combined length. This incorporates the       *
         * assumption that a single large insertion/deletion event is biologically more likely to happen      *
         * compared to many small insertions/deletions. While sophisticated gap scoring models can be applied *
         * in the generic algorithm by Waterman-Smith-Beyer (1976), affine gap scoring used in Gotoh's        *
         * algorithm enables a reasonable gap model with reduced runtime.                                     *
         ******************************************DP matrix***************************************************
         * gotoh algorithm DP variable explanations to align string a and b                                   *
         * A[row,col] := similarity of best alignment of a[1..row], b[1..col] with gap aligned against a[row] *
         * B[row,col] := similarity of best alignment of a[1..row], b[1..col] with gap aligned against b[col] *
         * D[row,col] := similarity of best alignment of a[1..row], b[1..col]                                 *
         * A[row,col] = max(A[row-1,col] + ge), D[row-1,col] + go)                                            *
         * B[row,col] = max(B[row,col-1] + ge), D[row,col-1] + go)                                            *
         * D[row,col] = max(D[row-1,col-1] + w(a[row], b[col]), A[row,col], B[row,col])                       *
         ******************************************************************************************************
         * DP matrix is 2D matrix of (a.size() + 1) * (b.size() + 1)                                          *
         *          b[0] b[1]...b[n]                                                                          *
         *       0                                                                                            *
         * a[0]                                                                                               *
         * a[1]                                                                                               *
         * ...                                                                                                *
         * a[m]                                                                                               *
         * ****************************************************************************************************
         * @param alnResult result of alignment
         * @return alignment score
         */
        int gotoh(Matrix2D<char>* alnResult);

        /******************************************needle******************************************************
         * Saul B. Needleman and Christian D. Wunsch introduced 1970 an approach to compute the optimal global*
         * alignment of two sequences. A minimizing variant was introduced 1974 by Peter H. Sellers.          *
         ******************************************DP matrix***************************************************
         * Needleman algorithm DP variable explanations to align string a and b                               *
         * D[row,col] := similarity of best alignment of a[1..row], b[1..col]                                 *
         * D[row,col] = max(D[row-1,col-1] + w(a[row], b[col]), D[row-1,col] + go, D[row,col-1] + go)         *
         ******************************************************************************************************
         * DP matrix is 2D matrix of (a.size() + 1) * (b.size() + 1)                                          *
         *          b[0] b[1]...b[n]                                                                          *
         *       0                                                                                            *
         * a[0]                                                                                               *
         * a[1]                                                                                               *
         * ...                                                                                                *
         * a[m]                                                                                               *
         * ****************************************************************************************************
         * @param alnResult result of alignment
         * @return alignment score
         */
        int needle(Matrix2D<char>* alnResult);

        /** do needle alignment of two sequences and their reverse complements to get a better split alignment result
         * @param s1 vertical sequence in DP, which can be manually gapped in middle to get a better alignment
         * @param s2 horizontal sequence in DP, which can not be manually gapped in middle to get a better alignment
         * @param alnResult Matrix2D to store better alignment result
         * @return true if a better split alignment found
         */
        bool splitAligner(const std::string& s1, const std::string& s2, Matrix2D<char>* alnResult);

        /** get homology sequence length of two sequences with max gaps limited
         * @param s1 sequence vertical sequence in DP, aligned to at most s1.length()
         * @param s2 sequence horizontal sequence in DP
         * @param maxGap max accumulated gaps allowed in horizontal/vertical direction
         * @return length of homolog stared from s1[0]
         */
        static int longestHomology(const std::string& s1, const std::string& s2, int maxGap);

        /** create alignment matrix of two set of sequences based on alignment trace
         * @param trace alignment trace of mSeqHorizontal and mSeqVertical
         * @param alnResult result of alignment
         */
        void createAlignment(const std::vector<char>& trace, Matrix2D<char>* alnResult);

        /** create alignment matrix of two sequences based on alignment trace
         * @param trace alignment trace of verSeq and hozSeq
         * @param verSeq vertical sequence in DP matrix
         * @param hozSeq horizontal sequence in DP matrix
         * @param alnResult matrix to store alignment result
         */
        void createAlignment(const std::vector<char>& trace, const std::string& verSeq, const std::string& hozSeq, Matrix2D<char>* alnResult);

        /** get horizontal gap penalty accumulated
         * @param pos position at which to open/extend gap
         * @param end position at which the sequence ends
         * @param len gap length until pos
         */
        inline int horizontalGapSum(int pos, int end, int len){
            return mAlignConfig->horizontalGapSum(pos, end, len);
        }
        
        /** get vertical gap penalty accumulated
         * @param pos position at which to open/extend gap
         * @param end position at which the sequence ends
         * @param len gap length until pos
         */
        inline int verticalGapSum(int pos, int end, int len){
            return mAlignConfig->verticalGapSum(pos, end, len);
        }
        
        /** get horizontal gap penalty extend cost
         * @param pos position at which to extend gap
         * @param end position at which the sequence ends
         */
        inline int horizontalGapExtend(int pos, int end){
            return mAlignConfig->horizontalGapExtend(pos, end);
        }

        /** get vertical gap penalty extend cost
         * @param pos position at which to extend gap
         * @param end position at which the sequence ends
         */
        inline int verticalGapExtend(int pos, int end){
            return mAlignConfig->verticalGapExtend(pos, end);
        }
};

#endif
