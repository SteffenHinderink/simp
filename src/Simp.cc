#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

/// A struct for a simplicial complex.
struct SimplicialComplex {
    /// Number of vertices.
    int n;
    /// Coloring of the vertices.
    std::vector<int> coloring;
    /// Facets of the simplicial complex given as vertex indices.
    std::vector<std::vector<int>> facets;
};

/**
 * Checks if the values in two vectors are the same, they may be permuted.
 * For facets that means that they are the same facet.
 * @param a First vector that is compared.
 * @param b Second vector that is compared.
 * @return True if the values are the same, false otherwise.
 */
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

/**
 * Calculates the signum of a permutation between two vectors.
 * The vectors need to have the same values and do not need to be sorted.
 * @param a First vector.
 * @param b Second vector with the same possibly permutated values as the first vector.
 * @return Signum of the permutation.
 */
int signum(std::vector<int> a, std::vector<int> b) {
    int n = a.size();
    // Get permutation in relation to ascending sequence
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
    // sgn(pi) = prod_1<=i<j<=n (pi(j) - pi(i)) / (j - i)
    float signum = 1;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            signum *= (float) (pi[j] - pi[i]) / (j - i);
        }
    }
    return (int) signum;
};

/**
 * Calculates the rank of a matrix.
 * @param m Matrix of which the rank is calculated.
 * @return Rank of the matrix.
 */
int rank(std::vector<std::vector<int>> m) {
    // Change type to float
    std::vector<std::vector<float>> a(m.size());
    for (int i = 0; i < a.size(); i++) {
        a[i] = std::vector<float>(m[i].size());
        for (int j = 0; j < a[i].size(); j++) {
            a[i][j] = m[i][j];
        }
    }

    int rank = 0;
    for (int col = 0; col < a[0].size(); col++) {
        int k = -1;
        for (int row = 0; row < a.size(); row++) {
            if (a[row][col] != 0) {
                if (k < 0) {
                    // Scale first non zero row (row k) in column to have 1 in that column
                    k = row;
                    float value = a[k][col];
                    for (int i = col; i < a[k].size(); i++) {
                        a[k][i] /= value;
                    }
                    // A non zero row increases the rank
                    rank++;
                } else {
                    // Subtract multiples of row k from any further non zero rows to have 0 in that column
                    float value = a[row][col];
                    for (int i = col; i < a[row].size(); i++) {
                        a[row][i] -= value * a[k][i];
                    }
                }
            }
        }
        // Delete row k as that index cannot increase the rank anymore
        if (k >= 0) {
            a.erase(a.begin() + k);
        }
    }
    return rank;
};

/**
 * Intersects a simplicial complex with a subset of its vertices.
 * The vertices of the intersection are not guaranteed to be 1 to n.
 * The coloring of the vertices of the intersection is defective.
 * @param sc Simplicial complex that is intersected.
 * @param mask Vertices with which the simplicial complex is intersected.
 *             The vertices are given as a number.
 *             The lowest bits of its binary representation stand for the vertices.
 *             1 means that the subset contains the vertex and 0 not.
 * @return Intersection of the simplicial complex with the subset of vertices.
 */
SimplicialComplex intersect(SimplicialComplex sc, int mask) {
    SimplicialComplex intersection;
    // Get number of vertices of the intersection
    intersection.n = 0;
    for (int i = 0; i < sc.n; i++) {
        if ((mask >> i) & 1) {
            intersection.n++;
        }
    }
    // Get facets of the intersection from the facets of the input
    for (int i = 0; i < sc.facets.size(); i++) {
        std::vector<int> facet;
        for (int j = 0; j < sc.facets[i].size(); j++) {
            if ((mask >> (sc.facets[i][j] - 1)) & 1) { // Requires the vertices of the input to be 1 to n
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

/**
 * Calculates the binomial coefficient (n choose k).
 * It depicts the number of subsets with size k of a set with size n.
 * @param n Number of elements in the set.
 * @param k Number of elements in the subsets.
 * @return Binomial coefficient.
 */
int binomialCoefficient(int n, int k) {
    // Multiplicative formula: n choose k = prod_i=1^k (n + 1 - i) / i
    float bc = 1;
    for (int i = 1; i <= k; i++) {
        bc *= (float) (n + 1 - i) / i;
    }
    return (int) bc;
};

/**
 * Returns the link of a vertex in a colored simplicial complex.
 * @param sc Colored simplicial complex in which the link exists.
 * @param x Vertex of which the link is returned.
 * @return Link of the vertex.
 */
SimplicialComplex link(SimplicialComplex sc, int x) {
    // Map from vertex indices in the simplicial complex -1 to vertex indices in the link
    int idx = 0;
    std::vector<int> map(sc.n);
    for (int i = 0; i < map.size(); i++) {
        map[i] = 0;
    }
    SimplicialComplex link;
    // Get facets of the link from the facets of the input
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
                        // Create new map entry
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

    // Get the coloring of the link respecting the mapped vertices
    int colorIdx = 0;
    std::vector<int> colorMap(sc.n);
    for (int i = 0; i < colorMap.size(); i++) {
        colorMap[i] = 0;
    }
    for (int i = 0; i < link.n; i++) {
        int color;
        for (int j = 0; j < sc.n; j++) {
            if (map[j] == i + 1) {
                if (colorMap[sc.coloring[j] - 1] == 0) {
                    colorIdx++;
                    colorMap[sc.coloring[j] - 1] = colorIdx;
                }
                link.coloring.push_back(colorMap[sc.coloring[j] - 1]);
                break;
            }
        }
    }
    
    return link;
};

/**
 * Returns the Betti numbers of a simplicial complex.
 * @param sc Simplicial complex of which the Betti numbers are calculated.
 * @return Betti numbers.
 */
std::vector<int> betti(SimplicialComplex sc) {
    // Get maximum facet dimension
    int maxFacetDim = 0;
    for (int i = 0; i < sc.facets.size(); i++) {
        if (sc.facets[i].size() > maxFacetDim) {
            maxFacetDim = sc.facets[i].size();
        }
    }
    std::vector<int> b(maxFacetDim);
    // Empty simplicial complex
    if (maxFacetDim == 0) {
        return b;
    }

    // Fill first C with facets of maximum dimension
    std::vector<std::vector<int>> cCurr;
    for (int i = 0; i < sc.facets.size(); i++) {
        if (sc.facets[i].size() == maxFacetDim) {
            cCurr.push_back(sc.facets[i]);
        }
    }

    /* 
     *           C_m          C_m-1         C_m-2                        C_1           C_0
     *            0    ----> Z^f_m-1 ----> Z^f_m-2 ---->   ...   ---->  Z^f_1  ---->   Z^n   ---->    Z
     *                delta_m      delta_m-1     delta_m-2      delta_2       delta_1       delta_0
     * 
     * dim(img)          0 ----\       x ----\       x --  ...     x ----\       x ----\
     *                          \             \                           \             \
     * dim(ker)                 |\---- x      |\---- x     ...  -- x      |\---- x      |\--- n-1
     *                          |             |                           |             |
     *                          V             V                           V             V
     * Betti number           b_m-1         b_m-2          ...           b_1           b_0
     *                       \                                                     /
     *                        \------------------------\ /------------------------/
     *                                                  V
     *                            Calculate image and kernel dimensions in loop
     */
    int preImgDim = 0;
    for (int idx = maxFacetDim - 1; idx > 0; idx--) {
        // C_idx-1
        std::vector<std::vector<int>> cNext;
        // New facets from the simplicial complex
        for (int i = 0; i < sc.facets.size(); i++) {
            if (sc.facets[i].size() == idx) {
                cNext.push_back(sc.facets[i]);
            }
        }
        // Facets that are subfacets of facets in C_idx
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

        // Calculate matrix delta_idx : C_idx -> C_idx-1
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
                        /*
                         * C_idx = (+)_i f_i Z
                         * f_i |-> f_i\(f_i)_1 - f_i\(f_i)_2 + f_i\(f_i)_2 - ... +- f_i\(f_i)_idx+1
                         */
                        delta[k][i] = sign * signum(facet, cNext[k]);
                        inserted = true;
                    }
                }
                sign *= -1;
            }
        }

        int r = rank(delta);
        // dim(ker(delta_idx)) = dim(C_idx) - rank(delta_idx)
        b[idx] = cCurr.size() - r - preImgDim;
        // dim(img(delta_idx)) = rank(delta_idx)
        preImgDim = r;
        cCurr = cNext;
    }
    b[0] = sc.n - 1 - preImgDim;

    return b;
}

/**
 * Returns the sigma numbers of a simplicial complex.
 * @param sc Simplicial complex of which the sigma numbers are calculated.
 * @return Sigma numbers.
 */
std::vector<float> sigma(SimplicialComplex sc) {
    // sigma_i(S) = sum_(W subset of V(S)) b_i(intersection of S with W) / (#V(S) choose #W)
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
            sigma[j] += j < b.size() ? (float) b[j] / bc : 0;
        }
    }
    return sigma;
}

/**
 * Returns the mu numbers of a simplicial complex.
 * @param sc Simplicial complex of which the mu numbers are calculated.
 * @return Mu numbers.
 */
std::vector<float> mu(SimplicialComplex sc) {
    // mu_i(S) = sum_(x in V(S)) sigma_i-1(link(x)) / (1 + #V(link(x)))
    std::vector<int> b = betti(sc);
    std::vector<float> mu(b.size() + 1);
    for (int i = 0; i < mu.size(); i++) {
        mu[i] = 0;
    }
    for (int i = 1; i <= sc.n; i++) {
        SimplicialComplex l = link(sc, i);
        std::vector<float> s = sigma(l);
        mu[0] += 1.0 / (1 + l.n);
        for (int j = 1; j < mu.size(); j++) {
            mu[j] += j - 1 < s.size() ? s[j - 1] / (1 + l.n) : 0;
        }
    }
    return mu;
}

/**
 * Returns the sigma' numbers of a colored simplicial complex.
 * @param sc Colored simplicial complex of which the sigma' numbers are calculated.
 * @return Sigma' numbers.
 */
std::vector<float> sigmaColored(SimplicialComplex sc) {
    // sigma'_i(S) = sum_(s subset of C(V(S))) b_i(intersection of S with C^-1(s)) / (#C(V(S)) choose #s)
    std::vector<int> b = betti(sc);
    std::vector<float> sigma(b.size());
    for (int i = 0; i < sigma.size(); i++) {
        sigma[i] = (float) b[i];
    }
    int numC = 0;
    for (int i = 0; i < sc.coloring.size(); i++) {
        if (sc.coloring[i] > numC) {
            numC = sc.coloring[i];
        }
    }
    for (int s = 0; s < (1 << numC) - 1; s++) {
        int mask = 0;
        int numS = 0;
        for (int i = 0; i < numC; i++) {
            if ((s >> i) & 1) {
                // Add vertices of this color to mask
                for (int j = 0; j < sc.n; j++) {
                    if (sc.coloring[j] == i + 1) {  // Requires the vertices of the input to be 1 to n
                        mask |= (1 << j);
                    }
                }
                numS++;
            }
        }
        SimplicialComplex intersection = intersect(sc, mask);
        std::vector<int> b = betti(intersection);
        int bc = binomialCoefficient(numC, numS);
        for (int j = 0; j < sigma.size(); j++) {
            sigma[j] += j < b.size() ? (float) b[j] / bc : 0;
        }
    }
    return sigma;
}

/**
 * Returns the tau numbers of a colored simplicial complex.
 * @param sc Colored simplicial complex of which the tau numbers are calculated.
 * @return Tau numbers.
 */
std::vector<float> tau(SimplicialComplex sc) {
    // tau_i(S) = sum_(x in V(S)) sigma'_i-1(link(x)) / (1 + #C(V(link(x))))
    std::vector<int> b = betti(sc);
    std::vector<float> tau(b.size() + 1);
    for (int i = 0; i < tau.size(); i++) {
        tau[i] = 0;
    }
    for (int i = 1; i <= sc.n; i++) {
        SimplicialComplex l = link(sc, i);
        std::vector<float> s = sigmaColored(l);
        int numC = 0;
        for (int j = 0; j < l.coloring.size(); j++) {
            if (l.coloring[j] > numC) {
                numC = l.coloring[j];
            }
        }
        tau[0] += 1.0 / (1 + numC);
        for (int j = 1; j < tau.size(); j++) {
            tau[j] += j - 1 < s.size() ? s[j - 1] / (1 + numC) : 0;
        }
    }
    return tau;
}

/**
 * Reads a colored simplicial complex from a file or the command line,
 * calculates the Betti, sigma, mu, sigma' and tau numbers of that simplicial complex and outputs them.
 * @param argc Number of command line arguments.
 * @param argv Command line arguments.
 *             If a second argument is given, it is the file containing a simplicial complex.
 * @return 0 if the program ran successfully.
 */
int main(int argc, char** argv) {
    // Read simplicial complex
    SimplicialComplex sc;
    std::istream is(std::cin.rdbuf());
    std::ifstream file;
    if (argc == 2) {
        file = std::ifstream(argv[1]);
        is.rdbuf(file.rdbuf());
    } else {
        std::cout << "Please enter the colors of the vertices of a simplicial complex" << std::endl;
        std::cout << "separated by spaces:" << std::endl;
    }
    sc.n = 0;
    std::string line;
    std::getline(is, line);
    std::istringstream iss(line);
    int x;
    while (iss >> x) {
        sc.coloring.push_back(x);
        sc.n++;
    }
    if (argc != 2) {
        std::cout << "Please enter the vertex indices of the facets of the simplicial complex" << std::endl;
        std::cout << "separated by spaces, stop with empty line:" << std::endl;
    }
    while (std::getline(is, line)) {
        if (line.empty()) {
            break;
        }
        iss = std::istringstream(line);
        std::vector<int> facet;
        while (iss >> x) {
            facet.push_back(x);
        }
        sc.facets.push_back(facet);
    }

    // Calculate Betti, sigma, mu, sigma' and tau numbers
    std::vector<int> b = betti(sc);
    std::vector<float> s = sigma(sc);
    std::vector<float> u = mu(sc);
    std::vector<float> swift = sigmaColored(sc);
    std::vector<float> taylor = tau(sc);
    
    // Output results
    std::cout << "   |    Betti |    sigma |       mu |   sigma' |      tau " << std::endl
              << "---+----------+----------+----------+----------+----------" << std::endl;
    for (int i = 0; i < b.size(); i++) {
        std::cout << std::setw(2) << i << " | " << std::setw(8) << b[i] << " | " << std::setw(8) << s[i] << " | " << std::setw(8) << u[i]
                  << " | " << std::setw(8) << swift[i] << " | " << std::setw(8) << taylor[i] << std::endl;
    }
    std::cout << std::setw(2) << b.size() << " |        0 |        0 | " << std::setw(8) << u[b.size()]
              << " |        0 | " << std::setw(8) << taylor[b.size()] << std::endl;
    return 0;
}
