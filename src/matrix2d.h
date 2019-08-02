#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <memory>
#include <iostream>

/** a simple 2d matrix */
template<typename T>
class Matrix2D{
    T** mMatrix; ///< C 2d array of T
    int mRow;    ///< row number
    int mCol;    ///< column number
    
    public:
    /** construct an empty 2D matrix */
    Matrix2D(){
        mMatrix = NULL;
        mRow = 0;
        mCol = 0;
    }
    
    /** construct  an 2D matrix
     * @param row number of rows
     * @param col number of columns
     */
    Matrix2D(int row, int col){
        mRow = row;
        mCol = col;
        initMatrix(row, col);
    }

    /** destroy an 2D matrix */
    ~Matrix2D(){
        if(mMatrix){
            for(int i = 0; i < mRow; ++i){
                free(mMatrix[i]);
                mMatrix[i] = NULL;
            }
            free(mMatrix);
            mMatrix = NULL;
        }
    }
        
    public:
    /** operator to get value
     * @param r row index 
     * @param c column index
     * @return mMatrix[r][c]
     */
    inline T operator()(int r, int c) const {
        return mMatrix[r][c];
    }

    /** operator to set value
     * @param r row index
     * @param c column index
     * @return reference of mMatrix[r][c]
     */
    inline T& operator()(int r, int c){
        return mMatrix[r][c];
    }

    /** get value
     * @param r row index 
     * @param c column index
     * @return mMatrix[r][c]
     */
    inline T get(int r, int c) const {
        return mMatrix[r][c];
    }

    /** set value
     * @param r row index
     * @param c column index
     * @return reference of mMatrix[r][c]
     */
    inline T& set(int r, int c){
        return mMatrix[r][c];
    } 

    /** resize array to specified dimension
     * @param r row number
     * @param c column number
     */
    inline void resize(int r, int c){
        freeMatrix();
        mRow = r;
        mCol = c;
        initMatrix(r, c);
    }
    /** get row number
     * @return row number
     */
    inline int nrow() const {
        return mRow;
    }

    /** get column number
     * @return column number
     */
    inline int ncol() const {
        return mCol;
    }

    /** initialize matrix to predefiend dimension
     * @param r row dimension
     * @param c column dimension
     */
    inline void initMatrix(int r, int c){
        mMatrix = (T**)std::calloc(sizeof(T*), r);
        for(int i = 0; i < r; ++i){
            mMatrix[i] = (T*)std::calloc(sizeof(T), c);
        }
    }

    /** free memory used by array */
    inline void freeMatrix(){
        if(mMatrix){
            for(int i = 0; i < mRow; ++i){
                free(mMatrix[i]);
                mMatrix[i] = NULL;
            }
            free(mMatrix);
            mMatrix = NULL;
        }
    }

    /** operator to output matrix */
    friend std::ostream& operator<<(std::ostream& os, Matrix2D* m){
        os << ' ' << "\t";
        for(int i = 0; i < m->ncol(); ++i){
            os << '[' << i << ']' << "\t";
        }
        os << std::endl;
        for(int i = 0; i < m->nrow(); ++i){
            os << '[' << i << ']' << "\t";
            for(int j = 0; j < m->ncol(); ++j){
                os << m->get(i, j) << "\t";
            }
            os << std::endl;
        }
        return os;
    }

    /** operator to output matrix */
    friend std::ostream& operator<<(std::ostream& os, const Matrix2D& m){
        os << "\n";
        for(int i = 0; i < m.nrow(); ++i){
            for(int j = 0; j < m.ncol(); ++j){
                os << m.get(i, j);
            }
            os << std::endl;
        }
        return os;
    }

    /** create profile matrix of sequence
     * @param seq sequence
     */
    inline void createProfile(const std::string& seq){
        resize(1, seq.size());
        for(size_t j = 0; j < seq.size(); ++j){
            set(0, j) = seq[j];
        }
    }
    
    /** create profile of multiple sequences
     * @param am multiple sequcne matrix
     */
    inline void createProfile(Matrix2D<char>* am){
        resize(6, am->ncol());
        for(int j = 0; j < am->ncol(); ++j){
            int sum = 0;
            for(int i = 0; i < am->nrow(); ++i){
                ++sum;
                switch(std::toupper(am->get(i, j))){
                    case 'A':
                        set(0, j) += 1;
                        break;
                    case 'C':
                        set(1, j) += 1;
                        break;
                    case 'G':
                        set(2, j) += 1;
                        break;
                    case 'T':
                        set(3, j) += 1;
                        break;
                    case 'N':
                        set(4, j) += 1;
                        break;
                    case '-':
                        set(5, j) += 1;
                        break;
                    default:
                        --sum;
                        break;
                }
            }
            for(int k = 0; k < 6; ++k){
                set(k, j) /= sum;
            }
        }
    }

    /** get identity percent of a gapped alignment result\n
     * skip leading and taling continuous gaps as well as predefined inner gap\n
     * @param gs skipped inner gap starting index
     * @param ge skipped inner gap ending index
     * @return ma / (ma + mm), mm = (internal gapped mismatch + non gapped mismatch), ma = bases matched
     */
    inline float identityPercent(int32_t gs, int32_t ge){
        bool s1NonGapMet = false;
        bool s2NonGapMet = false;
        int32_t gapMM = 0; // internal gap mismatch
        int32_t mm = 0; // total mismatch
        int32_t ma = 0; // base match
        bool inGap = false;
        for(int32_t j = 0; j < ncol(); ++j){
            if(j < gs || j > ge){
                if(get(0, j) != '-') s1NonGapMet = true;
                if(get(1, j) != '-') s2NonGapMet = true;
                if(get(0, j) == '-' || get(1, j) == '-'){
                    if(s1NonGapMet && s2NonGapMet){
                        if(!inGap){
                            inGap = true;
                            gapMM = 0;
                        }
                        gapMM += 1;
                    }
                }else{
                    if(inGap){
                        mm += gapMM;
                        inGap = false;
                    }
                    if(get(0, j) == get(1, j)) ma += 1;
                    else mm += 1;
                }
            }
        }
        return (float)ma/(float)(ma + mm);
    }
};

#endif
