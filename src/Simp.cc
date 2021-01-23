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

void printMatrix(std::vector<std::vector<int>> m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m[i].size(); j++) {
            std::cout << std::setw(3) << m[i][j];
        }
        std::cout << std::endl;
    }
}

std::vector<int> b(SimplicialComplex sc) {

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

    bool firstImg = true;
    auto imgDim = [firstImg](std::vector<std::vector<int>> m) mutable {

        // TODO: Bilddimension von Tarek

        if (firstImg) {
            firstImg = false;
            return 3;
        }
        return 6;
    };

    bool firstKer = true;
    auto kerDim = [firstKer](std::vector<std::vector<int>> m) mutable {

        // TODO: Kerndimension von Tarek

        if (firstKer) {
            firstKer = false;
            return 0;
        }
        return 6;
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
        std::cout << "--------" << std::endl << "delta_" << idx << ": C_" << idx << " -> C_" << idx - 1 << std::endl;

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

        std::cout << "C_" << idx << ":" << std::endl;
        printMatrix(cCurr);
        std::cout << "C_" << idx - 1 << ":" << std::endl;
        printMatrix(cNext);

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

        std::cout << "delta_" << idx << ":" << std::endl;
        printMatrix(delta);

        b[idx] = kerDim(delta) - preImgDim;
        preImgDim = imgDim(delta);
        cCurr = cNext;
    }
    b[0] = sc.n - 1 - preImgDim;

    return b;
}

int main() {
    SimplicialComplex sc = exampleSC();
    printMatrix(sc.facets);
    std::cout << "----------------" << std::endl;
    std::vector<int> betas = b(sc);
    std::cout << "----------------" << std::endl;
    for (int i = 0; i < betas.size(); i++) {
        std::cout << "b_" << i << ": " << betas[i] << std::endl;
    }
}
