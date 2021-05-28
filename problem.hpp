#include <bits/stdc++.h>

using namespace std;

using ld = long double;

struct Problem {
    static constexpr ld eps = 1e-12;  // for comparison
    static constexpr ld tol = 1e-12;  // biscetion tolerance

    size_t n;
    vector<ld> l, r;
    ld a;

    vector<size_t> h0_ids;  // Iₒᵥₑᵣ
    ld b;  // H₀(λ) = b

    Problem(size_t n_, const vector<ld>& l_, const vector<ld>& r_, ld a_) :
        n(n_), l(l_), r(r_), a(a_), b(0) {}

    virtual ~Problem() {};

    virtual ld f(size_t k, ld x) const = 0;  // fₖ(x)
    virtual ld g(size_t k, ld x) const = 0;  // gₖ(x) = f'ₖ(x)
    virtual ld q(size_t k, ld x) const = 0;  // g⁻¹ₖ(x)

    ld h(size_t k, ld x) const {
        if (x <= g(k, r[k])) return r[k];
        if (x >= g(k, l[k])) return l[k];
        return q(k, x);
    }

    ld h(ld x) const {
        ld res = 0;
        for (size_t k = 0; k < n; ++k) {
            res += h(k, x);
        }
        return res;
    }

    void prepare_h0() {}

    ld h0(ld x) const {
        ld res = 0;
        for (auto k : h0_ids) {
            res += h(k, x);
        }
        return res;
    }

    bool nonempty() const {
        auto l_sum = accumulate(l.begin(), l.end(), .0L);
        auto r_sum = accumulate(r.begin(), r.end(), .0L);
        return l_sum - eps <= a && a <= r_sum + eps;
    }

    vector<ld> breakpoints() const {
        vector<ld> res;
        for (size_t k = 0; k < n; ++k) {
            res.push_back(g(k, l[k]));
            res.push_back(g(k, r[k]));
        }
        return res;
    }

    pair<ld, ld> segment(vector<ld> bps) const {
        sort(bps.begin(), bps.end());
        int lo = -1, hi = 2 * n - 1;
        while (hi - lo > 1) {
            int mid = (lo + hi) / 2;
            if (h(bps[mid]) <= a) {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        if (hi == 0) {
            return {bps[0], bps[0]};
        }
        return {bps[hi - 1], bps[hi]};
    }

    void prepare_eq(ld seg_l, ld seg_r) {
        b = a;
        for (size_t k = 0; k < n; ++k) {
            if (g(k, r[k]) >= seg_r) {
                b -= r[k];
                continue;
            }
            if (g(k, l[k]) <= seg_l) {
                b -= l[k];
                continue;
            }
            h0_ids.push_back(k);
        }
    }

    ld solve_eq(ld seg_l, ld seg_r) const {
        ld lo = seg_l, hi = seg_r;
        while (hi - lo > tol) {
            ld mid = (lo + hi) / 2;
            if (h0(mid) > b) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        return (lo + hi) / 2;
    }

    vector<ld> solution(ld lambda) const {
        vector<ld> res(n);
        for (size_t k = 0; k < n; ++k) {
            res[k] = h(k, lambda);
        }
        return res;
    }

    vector<ld> solve() {
        assert(nonempty());
        auto bps = breakpoints();
        auto [seg_l, seg_r] = segment(bps);
        prepare_eq(seg_l, seg_r);
        prepare_h0();
        auto sol = solve_eq(seg_l, seg_r);
        return solution(sol);
    }
};
