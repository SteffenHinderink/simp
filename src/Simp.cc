#include <iostream>
#include <iomanip>
#include <vector>

struct SimplicialComplex {
    int n;
    std::vector<std::vector<int>> facets;
};

SimplicialComplex exampleSC() {
    SimplicialComplex sc;
    sc.n = 7;
    std::vector<std::vector<int>> facets(7);
    facets[0] = {2, 3, 7};
    facets[1] = {4, 5, 7};
    facets[2] = {5, 6, 7};
    facets[3] = {1, 2};
    facets[4] = {3, 4};
    facets[5] = {6, 1};
    facets[6] = {1, 7};
    sc.facets = facets;
    return sc;
}

std::vector<int> betti(SimplicialComplex sc) {

    auto equal = [](std::vector<int> a, std::vector<int> b) {
        int n = a.size();
        for (int i = 0; i < n; i++) {
            bool contained = false;
            for (int j = 0; j < n && !contained; j++) {
                if (a[i] == b[j]) {
                    contained = true;
                }
            }
            if (!contained) {
                return false;
            }
        }
        return true;
    };

    auto signum = [](std::vector<int> a, std::vector<int> b) {
        int n = a.size();
        std::vector<int> pi(n);
        for (int i = 0; i < n; i++) {
            bool found = false;
            for (int j = 0; j < n && !found; j++) {
                if (b[j] == a[i]) {
                    pi[j] = i;
                    found = true;
                }
            }
        }
        float signum = 1;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                signum *= (float) (pi[j] - pi[i]) / (j - i);
            }
        }
        return (int) signum;
    };

    auto rank = [](std::vector<std::vector<int>> m) {
        std::vector<std::vector<float>> a(m.size());
        for (int i = 0; i < a.size(); i++) {
            a[i] = std::vector<float>(m[i].size());
            for (int j = 0; j < a[i].size(); j++) {
                a[i][j] = m[i][j];
            }
        }

        int cnt = 0;
        for (int col = 0; col < a[0].size(); col++) {
            int k = -1;
            for (int row = 0; row < a.size(); row++) {
                if (k < 0) {
                    if (a[row][col] != 0) {
                        k = row;
                        float value = a[k][col];
                        for (int i = col; i < a[k].size(); i++) {
                            a[k][i] /= value;
                        }
                        cnt++;
                    }
                } else {
                    float value = a[row][col];
                    for (int i = col; i < a[row].size(); i++) {
                        a[row][i] -= value * a[k][i];
                    }
                }
            }
            if (k >= 0) {
                a.erase(a.begin() + k);
            }
        }
        return cnt;
    };

    int maxFacetDim = 0;
    for (int i = 0; i < sc.facets.size(); i++) {
        if (sc.facets[i].size() > maxFacetDim) {
            maxFacetDim = sc.facets[i].size();
        }
    }
    std::vector<int> b(maxFacetDim);

    std::vector<std::vector<int>> cCurr;
    for (int i = 0; i < sc.facets.size(); i++) {
        if (sc.facets[i].size() == maxFacetDim) {
            cCurr.push_back(sc.facets[i]);
        }
    }

    int preImgDim = 0;
    for (int idx = maxFacetDim - 1; idx > 0; idx--) {
        std::vector<std::vector<int>> cNext;
        for (int i = 0; i < sc.facets.size(); i++) {
            if (sc.facets[i].size() == idx) {
                cNext.push_back(sc.facets[i]);
            }
        }
        for (int i = 0; i < cCurr.size(); i++) {
            for (int j = 0; j < cCurr[i].size(); j++) {
                std::vector<int> facet;
                for (int k = 0; k < cCurr[i].size(); k++) {
                    if (k != j) {
                        facet.push_back(cCurr[i][k]);
                    }
                }
                bool newFacet = true;
                for (int k = 0; k < cNext.size() && newFacet; k++) {
                    if (equal(facet, cNext[k])) {
                        newFacet = false;
                    }
                }
                if (newFacet) {
                    cNext.push_back(facet);
                }
            }
        }

        std::vector<std::vector<int>> delta(cNext.size());
        for (int i = 0; i < cNext.size(); i++) {
            delta[i] = std::vector<int>(cCurr.size());
            for (int j = 0; j < cCurr.size(); j++) {
                delta[i][j] = 0;
            }
        }
        for (int i = 0; i < cCurr.size(); i++) {
            int sign = 1;
            for (int j = 0; j < cCurr[i].size(); j++) {
                std::vector<int> facet;
                for (int k = 0; k < cCurr[i].size(); k++) {
                    if (k != j) {
                        facet.push_back(cCurr[i][k]);
                    }
                }
                bool inserted = false;
                for (int k = 0; k < cNext.size() && !inserted; k++) {
                    if (equal(facet, cNext[k])) {
                        delta[k][i] = sign * signum(facet, cNext[k]);
                        inserted = true;
                    }
                }
                sign *= -1;
            }
        }

        int r = rank(delta);
        b[idx] = cCurr.size() - r - preImgDim;
        preImgDim = r;
        cCurr = cNext;
    }
    b[0] = sc.n - 1 - preImgDim;

    return b;
}

int main() {
    SimplicialComplex sc = exampleSC();
    std::vector<int> b = betti(sc);
    for (int i = 0; i < b.size(); i++) {
        std::cout << "b_" << i << ": " << b[i] << std::endl;
    }
}
