#include "aligner.h"

template<typename T>
int Aligner::score(const Matrix2D<T>* p1, const Matrix2D<T>* p2, int row, int col){
    if(mSeqVertical->nrow() == 1 && mSeqHorizontal->nrow() == 1){
        if(mSeqVertical->get(0, row) == mSeqHorizontal->get(0, col)){
            return mAlignConfig->mMatch;
        }else{
            return mAlignConfig->mMisMatch;
        }
    }else{
        float score = 0;
        for(int k1 = 0; k1 < 5; ++k1){
            for(int k2 = 0; k2 < 5; ++k2){
                score += p1->get(k1, row) * p2->get(k2, col) * ( k1 == k2 ? mAlignConfig->mMatch : mAlignConfig->mMisMatch);
            }
        }
        return score;
    }
}

int Aligner::gotoh(Matrix2D<char>* alnResult){
    // DP variables
    int m = mSeqVertical->ncol();
    int n = mSeqHorizontal->ncol();
    std::vector<int> s(n + 1, 0); // D[row, col]
    std::vector<int> v(n + 1, 0); // A[row, col]
    int newHoz = 0; // B[row, col]
    int preSub = 0; // D[row - 1, col - 1]
    // Trace Matrix
    int mf = n + 1;
    int totalBits = (m + 1) * (n + 1);
    bool* bitHor = (bool*)std::calloc(totalBits, sizeof(bool)); // bitHor[row*mf+col] = true if DP[row,col] on horizontal gap route
    bool* bitVer = (bool*)std::calloc(totalBits, sizeof(bool)); // bitVer[row*mf+col] = true if DP[row,col] on vertical gap route
    // Create profile
    Matrix2D<float>* p1 = new Matrix2D<float>();
    Matrix2D<float>* p2 = new Matrix2D<float>();
    if(mSeqVertical->nrow() != 1 || mSeqHorizontal->nrow() != 1){
        p1->createProfile(mSeqVertical);
        p2->createProfile(mSeqHorizontal);
    }
    // DP
    for(int row = 0; row <= m; ++row){
        for(int col = 0; col <= n; ++col){
            // Initialization
            if(row == 0 && col == 0){
                s[0] = 0; // D[0,0]
                v[0] = -mAlignConfig->mInf; // A[0, 0]
                newHoz = -mAlignConfig->mInf; // B[0, 0]
            }else if(row == 0){
                s[col] = mAlignConfig->horizontalGapSum(0, m, col); // D[0, col]
                v[col] = -mAlignConfig->mInf; // A[0, col]
                newHoz = s[col]; // B[0, col]
                bitHor[col] = true;
            }else if(col == 0){
                s[0] = mAlignConfig->verticalGapSum(0, n, row); // D[row, 0]
                v[0] = s[0]; // A[row, 0] 
                newHoz = -mAlignConfig->mInf; // B[row, 0]
                if(row == 1){
                    preSub = 0; // D[0, -1]
                }else{
                    preSub = mAlignConfig->verticalGapSum(0, n, row - 1);// D[row -1, col - 1]
                }
                bitVer[row * mf] = true;
            }else{
                // Recursion
                int preHoz = newHoz; // B[row, col - 1]
                int preVer = v[col]; // A[row - 1, col]
                int prePreSub = preSub; // D[row - 1,col - 1]
                preSub = s[col]; // D[row - 1, col]
                newHoz = std::max(s[col - 1] + horizontalGapSum(row, m, 1),  preHoz + horizontalGapExtend(row, m)); // B[row, col]
                v[col] = std::max(preSub + mAlignConfig->verticalGapSum(col, n, 1), preVer + mAlignConfig->verticalGapExtend(col, n)); // A[row, col]
                s[col] = std::max(std::max(prePreSub + score(p1, p2, row - 1, col - 1), newHoz), v[col]); // D[row, col]
                // Trace
                if(s[col] == newHoz){
                    bitHor[row * mf + col] = true; // D[row, col] == B[row, col]
                }else if(s[col] == v[col]){
                    bitVer[row * mf + col] = true; // D[row, col] == A[row, col]
                }
            }
        }
    }
    delete p1;
    p1 = NULL;
    delete p2;
    p2 = NULL;
    // Trace-back using pointers
    // 's' : mSeqVertical, mSeqHorizontal both consumed
    // 'v' : vertical gap, mSeqVertical consumed
    // 'h' : horizontal gap, mSeqHorizontal consumed
    int32_t row = m;
    int32_t col = n;
    std::vector<char> btr;
    while((row > 0) || (col > 0)){
        if(col > 0 && bitHor[row * mf + col]){
            --col;
            btr.push_back('h');
        }else if(row > 0 && bitVer[row * mf + col]){
            --row;
            btr.push_back('v');
        }else{
            --row;
            --col;
            btr.push_back('s');
        }
    }
    free(bitHor);
    bitHor = NULL;
    free(bitVer);
    bitVer = NULL;
    // Create alignment
    createAlignment(btr, alnResult);
    // Return sm
    return s[n];
}

int Aligner::needle(Matrix2D<char>* alnResult){
    // DP Matrix
    int m = mSeqVertical->ncol();
    int n = mSeqHorizontal->ncol();
    std::vector<int> s(n + 1, 0); // D[row, col]
    int preSub = 0; // D[row - 1, col - 1]
    // Trace Matrix
    int mf = n + 1;
    int totalBits = (m + 1) * (n + 1);
    bool* bitHor = (bool*)std::calloc(totalBits, sizeof(bool));
    bool* bitVer = (bool*)std::calloc(totalBits, sizeof(bool));
    // DP
    for(int row = 0; row <= m; ++row){
        for(int col = 0; col <= n; ++col){
            // Initialization
            if(row == 0 && col == 0){
                s[0] = 0;
                preSub = 0;
            }else if(row == 0){
                s[col] = horizontalGapExtend(0, m) * col;
                bitHor[col] = true;
            }else if(col == 0){
                s[0] = verticalGapExtend(0, n) * row;
                if(row == 1){
                    preSub = 0;
                }else{
                    preSub = verticalGapExtend(0, n) * (row - 1);
                }
                bitVer[row * mf] = 0;
            }else{
                // Recursion
                int prePreSub = preSub;
                preSub = s[col];
                s[col] = std::max(std::max(prePreSub + score(mSeqVertical, mSeqHorizontal, row  - 1, col -1), s[col - 1] + horizontalGapExtend(row, m)),
                                  preSub + verticalGapExtend(col, n));
                // Trace
                if(s[col] == s[col - 1] + horizontalGapExtend(row, m)) bitHor[row * mf + col] = true;
                else if(s[col] == preSub + verticalGapExtend(col, n)) bitVer[row * mf + col] = true;
            }
        }
    }
    // Trace-back using pointers
    // 's' : mSeqVertical, mSeqHorizontal both consumed
    // 'v' : vertical gap, mSeqVertical consumed
    // 'h' : horizontal gap, mSeqHorizontal consumed
    int32_t row = m;
    int32_t col = n;
    std::vector<char> btr;
    while((row > 0) || (col > 0)){
        if(col > 0 && bitHor[row * mf + col]){
            --col;
            btr.push_back('h');
        }else if(row > 0 && bitVer[row * mf + col]){
            --row;
            btr.push_back('v');
        }else{
            --row;
            --col;
            btr.push_back('s');
        }
    }
    free(bitHor);
    bitHor = NULL;
    free(bitVer);
    bitVer = NULL;
    // Create alignment
    createAlignment(btr, alnResult);
    // Return sm
    return s[n];
}

void Aligner::createAlignment(const std::vector<char>& trace, Matrix2D<char>* alnResult){
    int r1 = mSeqVertical->nrow();
    int r2 = mSeqHorizontal->nrow();
    alnResult->resize(r1 + r2, trace.size());
    int row = 0;
    int col = 0;
    int chn = 0; // index of aligned result trace string
    for(int c = trace.size() - 1; c >= 0; --c, ++chn){
        if(trace[c] == 's'){
            for(int i = 0; i < r1; ++i){
                alnResult->set(i, chn) = mSeqVertical->get(i, row);
            }
            for(int i = 0; i < r2; ++i){
                alnResult->set(r1 + i, chn) = mSeqHorizontal->get(i, col);
            }
            ++row;
            ++col;
        }else if(trace[c] == 'h'){
            for(int i = 0; i < r1; ++i){
                alnResult->set(i, chn) = '-';
            }
            for(int i = 0; i < r2; ++i){
                alnResult->set(r1 + i, chn) = mSeqHorizontal->get(i, col);
            }
            ++col;
        }else{
            for(int i = 0; i < r1; ++i){
                alnResult->set(i, chn) = mSeqVertical->get(i, row);
            }
            for(int i = 0; i < r2; ++i){
                alnResult->set(r1 + i, chn) = '-';
            }
            ++row;
        }
    }
}

void Aligner::createAlignment(const std::vector<char>& trace, const std::string& verSeq, const std::string& hozSeq, Matrix2D<char>* alnResult){
    alnResult->resize(2, trace.size());
    int row = 0;
    int col = 0;
    int chn = 0; // index of aligned result trace string
    for(int c = trace.size() - 1; c >= 0; --c, ++chn){
        if(trace[c] == 's'){
            alnResult->set(0, chn) = verSeq[row++];
            alnResult->set(1, chn) = hozSeq[col++];
        }else if(trace[c] == 'h'){
            alnResult->set(0, chn) = '-';
            alnResult->set(1, chn) = hozSeq[col++];
        }else{
            alnResult->set(0, chn) = verSeq[row++];
            alnResult->set(1, chn) = '-';
        }
    }
}

bool Aligner::splitAligner(const std::string& s1, const std::string& s2, Matrix2D<char>* alnResult){
    // DP Matrix
    int m = s1.size();
    int n = s2.size();
    Matrix2D<int>* mat = new Matrix2D<int>(m + 1, n + 1);
    // Initialization
    for(int col = 1; col <= n; ++col) mat->set(0, col) = mat->get(0, col - 1) + horizontalGapExtend(0, m);
    for(int row = 1; row <= m; ++row) mat->set(row, 0) = mat->get(row - 1, 0) + verticalGapExtend(0, n);
    // Forward alignment
    for(int row = 1; row <= m; ++row){
        for(int col = 1; col <= n; ++col){
            mat->set(row, col) = std::max(std::max(mat->get(row, col - 1) + horizontalGapExtend(row, m), mat->get(row - 1, col) + verticalGapExtend(col, n)),
                                          mat->get(row - 1, col - 1) + (s1[row - 1] == s2[col - 1] ? mAlignConfig->mMatch : mAlignConfig->mMisMatch));
        }
    }
    // Reverse alignment
    std::string sRev1 = util::reverseComplement(s1);
    std::string sRev2 = util::reverseComplement(s2);
    Matrix2D<int>* rev = new Matrix2D<int>(m + 1, n + 1);
    for(int col = 1; col <= n; ++col) rev->set(0, col) = rev->get(0, col - 1) + horizontalGapExtend(0, m);
    for(int row = 1; row <= m; ++row) rev->set(row, 0) = rev->get(row - 1, 0) + verticalGapExtend(0, n);
    for(int row = 1; row <= m; ++row){
        for(int col = 1; col <= n; ++col){
            rev->set(row, col) = std::max(std::max(rev->get(row, col - 1) + horizontalGapExtend(row, m), rev->get(row - 1, col) + verticalGapExtend(col, n)),
                                          rev->get(row - 1, col - 1) + (sRev1[row - 1] == sRev2[col - 1] ? mAlignConfig->mMatch : mAlignConfig->mMisMatch));
        }
    }
    if(mat->get(m, n) != rev->get(m, n)) return false;
    // Find best join
    Matrix2D<int>* bestMat = new Matrix2D<int>(m + 1, n + 1);
    for(int row = 0; row <= m; ++row){
        bestMat->set(row, 0) = mat->get(row, 0);
        for(int col = 1; col <= n; ++col){
            if(mat->get(row, col) > bestMat->get(row, col - 1)){
                bestMat->set(row, col) = mat->get(row, col);
            }else{
                bestMat->set(row, col) = bestMat->get(row, col -1);
            }
        }
    }
    Matrix2D<int>* bestRev = new Matrix2D<int>(m + 1, n + 1);
    for(int row = 0; row <= m; ++row){
        bestRev->set(row, 0) = rev->get(row, 0);
        for(int col = 1; col <= n; ++col){
            if(rev->get(row, col) > bestRev->get(row, col - 1)){
               bestRev->set(row, col) = rev->get(row, col);
            }else{
               bestRev->set(row, col) = bestRev->get(row, col - 1);
            }
        }
    }
    int bestScore = mat->get(m, n);
    int s1Left = 0;
    int s2Left = 0;
    // Find s1 right bound
    for(int row = 0; row <= m; ++row){
        for(int col = 0; col <= n; ++col){
            if(bestMat->get(row, col) + bestRev->get(m - row, n - col) > bestScore){
                bestScore = bestMat->get(row, col) + bestRev->get(m - row, n - col);
                s1Left = row;
                s2Left = col;
            }
        }
    }
    delete bestMat;
    bestMat = NULL;
    delete bestRev;
    bestRev = NULL;
    // Check if better split found
    if(bestScore < mat->get(m, n)){
        delete mat;
        mat = NULL;
        delete rev;
        rev = NULL;
        return false;
    }
    int s1Right = m - s1Left;
    int s2Right = 0;
    // Find s2 right bound
    for(int right = 0; right <= n - s2Left; ++right){
        if(mat->get(s1Left, s2Left) + rev->get(s1Right, right) == bestScore){
            s2Right = right;
        }
    }
    // Trace-back foward
    int rr = s1Left;
    int cc = s2Left;
    std::vector<char> trace;
    while(rr > 0 || cc > 0){
        if(rr > 0 && mat->get(rr, cc) == mat->get(rr - 1, cc) + verticalGapExtend(cc, n)){
            --rr;
            trace.push_back('v');
        }else if(cc > 0 && mat->get(rr, cc) == mat->get(rr, cc - 1) + horizontalGapExtend(rr, m)){
            --cc;
            trace.push_back('h');
        }else{
            --rr;
            --cc;
            trace.push_back('s');
        }
    }
    delete mat;
    mat = NULL;
    Matrix2D<char>* fwd = new Matrix2D<char>();
    createAlignment(trace, s1.substr(0, s1Left), s2.substr(0, s2Left), fwd);
    // Trace-back rev
    rr = s1Right;
    cc = s2Right;
    std::vector<char> rtrace;
    while(rr > 0 || cc > 0){
        if(rr > 0 && rev->get(rr, cc) == rev->get(rr - 1, cc) + verticalGapExtend(cc, n)){
            --rr;
            rtrace.push_back('v');
        }else if(cc > 0 && rev->get(rr, cc) == rev->get(rr, cc - 1) + horizontalGapExtend(rr, m)){
            --cc;
            rtrace.push_back('h');
        }else{
            --rr;
            --cc;
            rtrace.push_back('s');
        }
    }
    delete rev;
    rev = NULL;
    Matrix2D<char>* rvs = new Matrix2D<char>();
    createAlignment(rtrace, sRev1.substr(0, s1Right), sRev2.substr(0, s2Right), rvs);
    // Concat alignments
    int gaps2 = n - s2Right - s2Left;
    int alnlen = fwd->ncol() + rvs->ncol() + gaps2;
    alnResult->resize(2, alnlen);
    for(int i = 0; i < 2; ++i){
        int alncol = 0;
        for(;alncol < fwd->ncol(); ++alncol){
            alnResult->set(i, alncol) = fwd->get(i, alncol);
        }
        for(int j = s2Left; j < n - s2Right; ++j, ++alncol){
            if(i == 0){
                alnResult->set(i, alncol) = '-';
            }else{
                alnResult->set(i, alncol) = s2[j];
            }
        }
        for(int j = rvs->ncol() - 1; j >= 0; --j, ++alncol){
            switch(rvs->get(i, j)){
                case 'A':
                    alnResult->set(i, alncol) = 'T';
                    break;
                case 'C':
                    alnResult->set(i, alncol) = 'G';
                    break;
                case 'G':
                    alnResult->set(i, alncol) = 'C';
                    break;
                case 'T':
                    alnResult->set(i, alncol) = 'A';
                    break;
                case 'N':
                    alnResult->set(i, alncol) = 'N';
                    break;
                case '-':
                    alnResult->set(i, alncol) = '-';
                    break;
                default:
                    break;
            }
        }
    }
    delete fwd;
    fwd = NULL;
    delete rvs;
    rvs = NULL;
    return true;
}

int Aligner::longestHomology(const std::string& s1, const std::string& s2, int maxGap){
    // DP Matrix
    int m = s1.size();
    int n = s2.size();
    Matrix2D<int>* mat = new Matrix2D<int>(m + 1, n + 1);
    // Initialization
    int k = std::abs(maxGap);
    int c = std::min(k, std::min(m, n));
    for(int i = 1; i <= c; ++i){
        mat->set(0, i) = mat->get(0, i - 1) - 1;
        mat->set(i, 0) = mat->get(i - 1, 0) - 1;
    }
    // Edit distance
    for(int row = 1; row <= m; ++row){
        int bestCol = -maxGap - 1;
        for(int h = -k; h <= k; ++h){
            int32_t col = row + h;
            if((col >= 1) && (col <= n)){
                mat->set(row, col) = mat->get(row - 1, col - 1) + (s1[row - 1] == s2[col - 1] ? 0 : -1);
                if((row - 1 - col >= -k) && (row - 1 - col <= k)){// vertical gap
                    mat->set(row, col) = std::max(mat->get(row, col), mat->get(row - 1, col) - 1);
                }
                if((row - col + 1 >= -k) && (row - col + 1 <= k)){// horizontal gap
                    mat->set(row, col) = std::max(mat->get(row, col), mat->get(row, col - 1) - 1);
                }
                if(mat->get(row, col) > bestCol){
                    bestCol = mat->get(row, col);
                }
            }
        }
        if(bestCol < -maxGap){// longest homolog found
            delete mat;
            mat = NULL;
            return row - 1;
        }
    }
    delete mat;
    mat = NULL;
    return s1.length();
}
