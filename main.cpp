#include <bits/stdc++.h>
#include "problem.hpp"

using namespace std;

using ld = long double;

struct LogShiftProblem : Problem {
    vector<ld> shifts;
    ld coef1, coef2;

    LogShiftProblem(
        size_t n_, const vector<ld>& shifts_, const vector<ld>& l_, const vector<ld>& r_, ld a_
    ) : Problem(n_, l_, r_, a_), shifts(shifts_), coef1(0), coef2(0)  {}

    ~LogShiftProblem() {}

    ld f(size_t k, ld x) const {
        return log(x + shifts[k]);
    };

    ld g(size_t k, ld x) const {
        return 1. / (x + shifts[k]);
    };

    ld q(size_t k, ld x) const {
        return 1. / x - shifts[k];
    };

    void prepare_h0() {
        coef1 = coef2 = 0;
        for (auto k : h0_ids) {
            ++coef1;
            coef2 -= shifts[k];
        }
    }

    ld h0(ld x) const {
        return coef1 / x + coef2;
    }
};

mt19937 rng(random_device{}());

LogShiftProblem gen(size_t n, ld shift_l, ld shift_r, ld ineq_l, ld ineq_r) {
    uniform_real_distribution<ld> shift_rd(shift_l, shift_r);
    vector<ld> shifts;
    for (size_t k = 0; k < n; ++k) {
        shifts.push_back(shift_rd(rng));
    }
    vector<ld> l, r;
    uniform_real_distribution<ld> ineq_rd(ineq_l, ineq_r);
    for (size_t k = 0; k < n; ++k) {
        l.push_back(ineq_rd(rng));
        r.push_back(ineq_rd(rng));
        if (l[k] > r[k]) {
            swap(l[k], r[k]);
        }
    }
    auto l_sum = accumulate(l.begin(), l.end(), .0L);
    auto r_sum = accumulate(r.begin(), r.end(), .0L);
    uniform_real_distribution<ld> a_rd(l_sum, r_sum);
    ld a = a_rd(rng);
    return LogShiftProblem(n, shifts, l, r, a);
}

int main() {
    size_t iterations = 20;
    size_t n;
    cin >> n;
    ld shift_l = 1, shift_r = 2;
    ld ineq_l = 0, ineq_r = 3;

    ld time_total = 0;
    for (int i = 0; i < iterations; ++i) {
        auto P = gen(n, shift_l, shift_r, ineq_l, ineq_r);
        auto start = clock();
        P.solve();
        auto finish = clock();
        time_total += ld(finish - start) / CLOCKS_PER_SEC;
        cout << (i + 1) << endl;
    }

    cout << setprecision(10) << fixed;
    cout << "Average Time: " << (time_total / iterations) << endl;
}
