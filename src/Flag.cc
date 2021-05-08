#include <iostream>
#include <vector>

int binomial_coefficient(int n, int k) {
    // Multiplicative formula: n choose k = prod_i=1^k (n + 1 - i) / i
    float bc = 1;
    for (int i = 1; i <= k; i++) {
        bc *= (float) (n + 1 - i) / i;
    }
    return (int) bc;
};

int main() {
    int n;
    std::cout << "Dimension?" << std::endl;
    std::cin >> n;
    std::vector<int> f(n + 1);
    for (int i = 0; i < n + 1; i++) {
        std::cout << i << "? ";
        std::cin >> f[i];
    }

    std::vector<std::vector<int>> s(1 << n + 1);

    std::vector<int> fs(1 << n + 1);
    for (int i = 0; i < 1 << n + 1; i++) {
        for (int j = 0; j < n + 1; j++) {
            if (i >> j & 1) {
                s[i].push_back(j + 1);
            }
        }

        fs[i] = 1;
        if (s[i].size() > 0) {
            fs[i] = f[s[i][s[i].size() - 1] - 1];
            for (int j = 0; j < s[i].size() - 1; j++) {
                fs[i] *= binomial_coefficient(s[i][j + 1], s[i][j]);
            }
        }
    }

    std::vector<int> h(1 << n + 1);
    for (int i = 0; i < 1 << n + 1; i++) {
        h[i] = 0;
        for (int j = 0; j < 1 << s[i].size(); j++) {
            int index = 0;
            int length = 0;
            for (int k = 0; k < s[i].size(); k++) {
                if (j >> k & 1) {
                    index += 1 << (s[i][k] - 1);
                    length++;
                }
            }
            h[i] += ((s[i].size() - length) % 2 == 0 ? 1 : -1) * fs[index];
        }
    }

    for (int i = 0; i < 1 << n + 1; i++) {
        std::cout << std::endl
                  << "S = {";
        for (int j = 0; j < s[i].size(); j++) {
            std::cout << s[i][j];
            if (j < s[i].size() - 1) {
                std::cout << " ";
            }
        }
        std::cout << "}" << std::endl;
        std::cout << "f_S = " << fs[i] << std::endl
                  << "h^S = " << h[i] << std::endl;
    }
}
