#include "msa.h"

int MSA::lcs(const std::string& s1, const std::string& s2){
    int m = s1.size();
    int n = s2.size();
    int preDiag = 0;
    int prePreDiag = 0;
    std::vector<int> s(n + 1, 0);
    for(int i = 0; i <= m; ++i){
        for(int j = 0; j <= n; ++j){
            if(i == 0 || j == 0){
                s[j] = 0;
                preDiag = 0;
                prePreDiag = 0;
            }else{
                prePreDiag = preDiag;
                preDiag = s[j];
                if(s1[i - 1] == s2[i - 1]){
                    s[j] = prePreDiag + 1;
                }else{
                    s[j] = (s[j] > s[j-1] ? s[j] : s[j-1]);
                }
            }
        }
    }
    return s[n];
}

void MSA::distanceMatrix(Matrix2D<int>* d){
    auto iter1 = mSeqs->cbegin();
    auto iter2 = iter1;
    for(int i = 0; iter1 != mSeqs->cend(); ++iter1, ++i){
        iter2 = iter1;
        ++iter2;
        for(int j = i + 1; iter2 != mSeqs->cend(); ++iter2, ++j){
            d->set(i, j) = lcs(*iter1, *iter2) * 100 / (std::min(iter1->length(), iter2->length()));
        }
    }
}

int MSA::closestPair(Matrix2D<int>* d, int num, int& di, int& dj){
    int dMax = -1;
    for(int i = 0; i < num; ++i){
        for(int j = i + 1; j < num; ++j){
            if(d->get(i, j) > dMax){
                di = i;
                dj = j;
                dMax = d->get(i, j);
            }
        }
    }
    return dMax;
}

void MSA::updateDistanceMatrix(Matrix2D<int>* d, const Matrix2D<int>* p, int num, int di, int dj){
    // set d[i, num] = (distance(di, i) + distance(di, j))/2
    for(int i = 0; i < num; ++i){
        if(p->get(i, 0) == -1){
            d->set(i, num) = ((di < i ? d->get(di, i) : d->get(i, di)) + (dj < i ? d->get(dj, i) : d->get(i, dj))) / 2;
        }
    }
    // set all distance to di to -1(column)
    for(int i = 0; i < di; ++i){
        d->set(i, di) = -1;
    }
    // set all distance to di to -1(row)
    for(int i = di + 1; i < num + 1; ++i){
        d->set(di, i) = -1;
    }
    // set all distance to dj to -1(column)
    for(int j = 0; j < dj; ++j){
        d->set(j, dj) = -1;
    }
    // set all distance to dj to -1(row)
    for(int j = dj + 1; j < num + 1; ++j){
        d->set(dj, j) = -1;
    }
}

int MSA::upgma(Matrix2D<int>* d, Matrix2D<int>* p, int num){
    int nn = num;
    for(; nn < 2 * num + 1; ++nn){
        int di = 0;
        int dj = 0;
        if(closestPair(d, nn, di, dj) == -1){
            break;
        }
        p->set(di, 0) = nn;
        p->set(dj, 0) = nn;
        p->set(nn, 1) = di;
        p->set(nn, 2) = dj;
        updateDistanceMatrix(d, p, nn, di, dj);
    }
    return nn > 0 ? (nn - 1) : 0;
}

void MSA::palign(const Matrix2D<int>* p, int root, Matrix2D<char>* aln){
    if(p->get(root, 1) == -1 && p->get(root, 2) == -1){
        auto iter = mSeqs->begin();
        if(root){
            std::advance(iter, root);
        }
        aln->resize(1, iter->size());
        int ind = 0;
        for(auto siter = iter->begin(); siter != iter->end(); ++siter){
            aln->set(0, ind++) = *siter;
        }
    }else{
        Matrix2D<char>* aln1 = new Matrix2D<char>();
        palign(p, p->get(root, 1), aln1);
        Matrix2D<char>* aln2 = new Matrix2D<char>();
        palign(p, p->get(root, 2), aln2);
        Aligner* aligner = new Aligner(aln1, aln2, mAlignConfig);
        aligner->gotoh(aln);
        delete aln1;
        delete aln2;
        delete aligner;
    }
}

void MSA::consensus(Matrix2D<char>* aln, std::string& cs){
    // Calculate coverage of non-gaps
    Matrix2D<bool>* fl = new Matrix2D<bool>(aln->nrow(), aln->ncol());
    std::vector<int> cov(aln->ncol(), 0);
    for(int i = 0; i < aln->nrow(); ++i){
        int beg = 0;
        int end = -1;
        // Skipping leading/ending continuous gaps
        for(int j = 0; j < aln->ncol(); ++j){
            fl->set(i, j) = false;
            if(aln->get(i, j) != '-'){
                end = j;
            }else if(end == -1){
                beg = j + 1;
            }
        }
        for(int j = beg; j <= end; ++j){
            ++cov[j];
            fl->set(i, j) = true;
        }
    }
    // Get consensus sequence
    int j = 0;
    for(auto itCov = cov.begin(); itCov != cov.end(); ++itCov, ++j){
        if(*itCov >= mMinCovForCS){
            // Get consensus letter
            const std::string base = "ACGT";
            std::vector<int> countBase(4, 0);
            for(int i = 0; i < aln->nrow(); ++i){
                if(fl->get(i, j)){
                    switch(aln->get(i, j)){
                        case 'A': case 'a':
                            ++countBase[0];
                            break;
                        case 'C': case 'c':
                            ++countBase[1];
                            break;
                        case 'G': case 'g':
                            ++countBase[2];
                            break;
                        case 'T': case 't':
                            ++countBase[3];
                            break;
                        default:
                            break;
                    }
                }
            }
            int countAligned = std::accumulate(countBase.begin(), countBase.end(), 0);
            if(countAligned > mMinBaseRatioForCS * (*itCov)){
                int maxInd = 0;
                int maxCnt = countBase[0];
                for(int ind = 1; ind < 4; ++ind){
                    if(countBase[ind] > maxCnt){
                        maxCnt = countBase[ind];
                        maxInd = ind;
                    }
                }
                cs.push_back(base[maxInd]);
            }
        }
    }
    delete fl;
}



int MSA::msa(std::string& cs){
    // Compute distance matrix
    int num = mSeqs->size();
    Matrix2D<int>* d = new Matrix2D<int>(2 * num + 1, 2 * num + 1);
    for(int i = 0; i < 2 * num + 1; ++i){
        for(int j = i + 1; j < 2 * num + 1; ++j){
            d->set(i, j) = -1;
        }
    }
    distanceMatrix(d);
    // UPGMA
    Matrix2D<int>* p = new Matrix2D<int>(2 * num + 1, 3);
    for(int i = 0; i < 2 * num + 1; ++i){
        for(int j = 0; j < 3; ++j){
            p->set(i, j) = -1;
        }
    }
    int root = upgma(d, p, num);
    // Progressive Alignment
    Matrix2D<char>* aln = new Matrix2D<char>();
    palign(p, root, aln);
    // Consensus calling
    consensus(aln, cs);
    // Cleanup Resources
    delete d;
    delete p;
    int support = aln->nrow();
    delete aln;
    // Return split-read support;
    return support;
}
