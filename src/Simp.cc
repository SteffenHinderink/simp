#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>

struct SimplicialComplex {
    int n;
    std::vector<std::vector<int>> facets;
};

bool permuted(std::vector<int> a, std::vector<int> b) {
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

std::vector<int> betti(SimplicialComplex sc) {

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
    if (maxFacetDim == 0) {
        return b;
    }

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
                    if (permuted(facet, cNext[k])) {
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
                    if (permuted(facet, cNext[k])) {
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

std::vector<float> sigma(SimplicialComplex sc) {

    auto intersect = [](SimplicialComplex sc, int mask) {
        SimplicialComplex intersection;
        intersection.n = 0;
        for (int i = 0; i < sc.n; i++) {
            if ((mask >> i) & 1) {
                intersection.n++;
            }
        }
        for (int i = 0; i < sc.facets.size(); i++) {
            std::vector<int> facet;
            for (int j = 0; j < sc.facets[i].size(); j++) {
                if ((mask >> (sc.facets[i][j] - 1)) & 1) { // Requires the vertices to be 1 to n
                    facet.push_back(sc.facets[i][j]);
                }
            }
            if (facet.size() > 0) {
                bool newFacet = true;
                for (int j = 0; j < intersection.facets.size() && newFacet; j++) {
                    if (permuted(facet, intersection.facets[j])) {
                        newFacet = false;
                    }
                }
                if (newFacet) {
                    intersection.facets.push_back(facet);
                }
            }
        }
        return intersection;
    };

    auto binomialCoefficient = [](int n, int k) {
        float bc = 1;
        for (int j = 1; j <= k; j++) {
            bc *= (float) (n + 1 - j) / j;
        }
        return (int) bc;
    };

    std::vector<int> b = betti(sc);
    std::vector<float> sigma(b.size());
    for (int i = 0; i < sigma.size(); i++) {
        sigma[i] = (float) b[i];
    }
    for (int i = 0; i < (1 << sc.n) - 1; i++) {
        SimplicialComplex intersection = intersect(sc, i);
        std::vector<int> b = betti(intersection);
        int bc = binomialCoefficient(sc.n, intersection.n);
        for (int j = 0; j < sigma.size(); j++) {
            sigma[j] += (float) (j < b.size() ? b[j] : 0) / bc;
        }
    }
    return sigma;
}

std::vector<float> mu(SimplicialComplex sc) {

    auto link = [](SimplicialComplex sc, int x) {
        int idx = 0;
        std::vector<int> map(sc.n);
        for (int j = 0; j < map.size(); j++) {
            map[j] = 0;
        }
        SimplicialComplex link;
        for (int i = 0; i < sc.facets.size(); i++) {
            bool xContained = false;
            for (int j = 0; j < sc.facets[i].size() && !xContained; j++) {
                if (sc.facets[i][j] == x) {
                    xContained = true;
                }
            }
            if (xContained) {
                std::vector<int> facet;
                for (int j = 0; j < sc.facets[i].size(); j++) {
                    if (sc.facets[i][j] != x) {
                        if (map[sc.facets[i][j] - 1] == 0) {
                            idx++;
                            map[sc.facets[i][j] - 1] = idx;
                        }
                        facet.push_back(map[sc.facets[i][j] - 1]);
                    }
                }
                if (facet.size() > 0) {
                    link.facets.push_back(facet);
                }
            }
        }
        link.n = idx;
        return link;
    };

    std::vector<int> b = betti(sc);
    std::vector<float> mu(b.size());
    for (int i = 0; i < mu.size(); i++) {
        mu[i] = 0;
    }
    for (int i = 1; i <= sc.n; i++) {
        SimplicialComplex l = link(sc, i);
        std::vector<float> s = sigma(l);
        for (int j = 0; j < mu.size(); j++) {
            mu[j] += (j < s.size() ? s[j] : 0) / (1 + l.n);
        }
    }
    return mu;
}

int main(int argc, char** argv) {
    SimplicialComplex sc;
    std::istream is(std::cin.rdbuf());
    std::ifstream file;
    if (argc == 2) {
        file = std::ifstream(argv[1]);
        is.rdbuf(file.rdbuf());
    } else {
        std::cout << "Please enter the vertex indices of the facets of a simplicial complex" << std::endl;
        std::cout << "separated by spaces, stop with empty line:" << std::endl;
    }
    sc.n = 0;
    std::string line;
    while (std::getline(is, line)) {
        if (line.empty()) {
            break;
        }
        std::istringstream iss(line);
        std::vector<int> facet;
        int x;
        while (iss >> x) {
            if (x > sc.n) {
                sc.n = x;
            }
            facet.push_back(x);
        }
        sc.facets.push_back(facet);
    }

    std::vector<int> b = betti(sc);
    std::vector<float> s = sigma(sc);
    std::vector<float> u = mu(sc);
    
    std::cout << "   |    Betti |    sigma |       mu" << std::endl << "---+----------+----------+----------" << std::endl;
    std::cout << " 0 | " << std::setw(8) << b[0] << " | " << std::setw(8) << s[0] << " |        -" << std::endl;
    for (int i = 1; i < b.size(); i++) {
        std::cout << std::setw(2) << i << " | " << std::setw(8) << b[i] << " | " << std::setw(8) << s[i] << " | " << std::setw(8) << u[i - 1] << std::endl;
    }
    std::cout << std::setw(2) << b.size() << " |        0 |        0 | " << std::setw(8) << u[b.size() - 1] << std::endl;
    return 0;
}
