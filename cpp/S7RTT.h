// ==============================================================================
// File Name:    S7RTT.h
// Author:       feecat
// Version:      V1.8.1
// Description:  Simple 7seg Real-Time Trajectory Generator
// Website:      https://github.com/feecat/S7RTT
// License:      Apache License Version 2.0
// ==============================================================================
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ==============================================================================

#ifndef S7RTT_H
#define S7RTT_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>
#include <array>

namespace S7RTT_Lib {

// ==============================================================================
// MotionState
// ==============================================================================
struct MotionState {
    double dt; // Duration of this segment
    double p;  // Position
    double v;  // Velocity
    double a;  // Acceleration
    double j;  // Jerk

    MotionState(double _dt = 0.0, double _p = 0.0, double _v = 0.0, double _a = 0.0, double _j = 0.0)
        : dt(_dt), p(_p), v(_v), a(_a), j(_j) {}
};

// ==============================================================================
// Internal Helpers
// ==============================================================================
struct Segment {
    double dt;
    double j;
};

// Stack-allocated container to avoid heap allocation during profile construction
template<int N=3>
struct TinyProfile {
    std::array<Segment, N> segs;
    int count = 0;

    inline void push(double dt, double j) {
        if (count < N) {
            segs[count].dt = dt;
            segs[count].j = j;
            count++;
        }
    }

    inline void clear() { count = 0; }
    const Segment* begin() const { return segs.data(); }
    const Segment* end() const { return segs.data() + count; }
    inline Segment& operator[](int i) { return segs[i]; }
    inline const Segment& operator[](int i) const { return segs[i]; }
};

// ==============================================================================
// S7RTT Class
// ==============================================================================
class S7RTT {
private:
    // --- Constants ---
    static constexpr double EPS_TIME   = 1e-10;
    static constexpr double EPS_DIST   = 1e-10;
    static constexpr double EPS_VEL    = 1e-10;
    static constexpr double EPS_ACC    = 1e-10;
    static constexpr double MATH_EPS   = 1e-10;

    static constexpr double EPS_SOLVER = 1e-4;
    static constexpr int    SOLVER_ITER = 50;
    static constexpr double SOLVER_TOL  = 1e-10;

    static constexpr double ONE_SIXTH  = 1.0 / 6.0;
    static constexpr double ONE_HALF   = 0.5;

public:
    S7RTT() = default;

    // ==========================================================================
    // 1. Core Integrator (Inlined & Optimized)
    // ==========================================================================

    // Optimized inplace integrator
    static inline void _integrate_state_inplace(MotionState& s, double dt, double j) {
        double dt2 = dt * dt;
        double dt3 = dt2 * dt;

        s.p += s.v * dt + s.a * dt2 * ONE_HALF + j * dt3 * ONE_SIXTH;
        s.v += s.a * dt + j * dt2 * ONE_HALF;
        s.a += j * dt;
    }

    // Helper to return a new state
    static inline MotionState _integrate_step(const MotionState& s, double dt, double j) {
        if (dt <= EPS_TIME) return s;
        MotionState next = s;
        _integrate_state_inplace(next, dt, j);
        next.dt = 0.0;
        next.j = j;
        return next;
    }

private:
    // ==========================================================================
    // 2. Simulation & Append Helpers (Memory Optimized)
    // ==========================================================================

    // Optimization: Append directly to output vector
    template <int N>
    MotionState _append_from_profile(std::vector<MotionState>& nodes, const MotionState& start_s, const TinyProfile<N>& shapes) {
        MotionState curr = start_s;
        for (int i = 0; i < shapes.count; ++i) {
            double dt = shapes[i].dt;
            double j = shapes[i].j;

            if (dt < EPS_TIME) continue;

            MotionState seg_start = curr;
            seg_start.dt = dt;
            seg_start.j = j;
            nodes.push_back(seg_start);

            _integrate_state_inplace(curr, dt, j);
        }
        return curr;
    }

    // Solver helper: Only simulates the endpoint (Fast path, no vector ops)
    template <int N>
    static inline void _simulate_endpoint_inplace(MotionState& curr, const TinyProfile<N>& shapes) {
        for (int i = 0; i < shapes.count; ++i) {
            if (shapes[i].dt >= EPS_TIME) {
                _integrate_state_inplace(curr, shapes[i].dt, shapes[i].j);
            }
        }
    }

    // --- Optimized Saturation Logic ---

    // Calculates the saturated state for solver use (Math only)
    static inline void _integrate_saturated_state_only(MotionState& curr, double t, double j_apply, double a_max) {
        double limit_a = (j_apply > 0) ? a_max : -a_max;
        double dist_to_lim = limit_a - curr.a;

        if (std::abs(j_apply) < MATH_EPS) {
            // If Jerk is 0, we are just moving with constant acceleration (which might be limit_a)
            if (t > EPS_TIME) _integrate_state_inplace(curr, t, j_apply);
            return;
        }

        // Check if we are moving towards the limit
        bool is_same_dir = (j_apply > 0 && dist_to_lim > -MATH_EPS) ||
                           (j_apply < 0 && dist_to_lim < MATH_EPS);

        double t_ramp = is_same_dir ? (dist_to_lim / j_apply) : 0.0;

        if (t <= t_ramp) {
            if (t > EPS_TIME) _integrate_state_inplace(curr, t, j_apply);
        } else {
            if (t_ramp > EPS_TIME) _integrate_state_inplace(curr, t_ramp, j_apply);

            curr.a = limit_a; // Clamp a to limit to avoid numerical drift

            double t_hold = t - t_ramp;
            if (t_hold > EPS_TIME) _integrate_state_inplace(curr, t_hold, 0.0);
        }
    }

    // Generates and appends saturated profile nodes (Used for final generation)
    MotionState _append_saturated_profile(std::vector<MotionState>& nodes, const MotionState& s, double t, double j_apply, double a_max) {
        MotionState curr = s;
        if (t <= EPS_TIME) return curr;

        double limit_a = (j_apply > 0) ? a_max : -a_max;
        double dist_to_lim = limit_a - s.a;

        double t_ramp = 0.0;
        if (std::abs(j_apply) < MATH_EPS) {
            t_ramp = std::numeric_limits<double>::infinity();
        } else {
            bool is_same_dir = (j_apply > 0 && dist_to_lim > -MATH_EPS) ||
                               (j_apply < 0 && dist_to_lim < MATH_EPS);
            t_ramp = is_same_dir ? (dist_to_lim / j_apply) : 0.0;
        }

        if (t <= t_ramp) {
            if (t > EPS_TIME) {
                MotionState node = curr;
                node.dt = t; node.j = j_apply;
                nodes.push_back(node);
                _integrate_state_inplace(curr, t, j_apply);
            }
        } else {
            if (t_ramp > EPS_TIME) {
                MotionState node = curr;
                node.dt = t_ramp; node.j = j_apply;
                nodes.push_back(node);
                _integrate_state_inplace(curr, t_ramp, j_apply);
            }

            double t_hold = t - t_ramp;
            if (t_hold > EPS_TIME) {
                curr.a = limit_a; // Precise clamp
                MotionState node = curr;
                node.dt = t_hold; node.j = 0.0;
                nodes.push_back(node);
                _integrate_state_inplace(curr, t_hold, 0.0);
            }
        }
        return curr;
    }

    // ==========================================================================
    // 3. Math & Solvers
    // ==========================================================================

    struct VelChangeTimes {
        double t1, t2, t3;
        double dir;
    };

    static inline VelChangeTimes _calc_vel_change_times(double v0, double a0, double v1, double a_max, double j_max) {
        double _a0 = std::max(-a_max, std::min(a_max, a0));

        double t_to_zero = std::abs(_a0) / j_max;
        double j_restore = (std::abs(_a0) > MATH_EPS) ? -std::copysign(j_max, _a0) : 0.0;
        double v_min_feasible = v0 + _a0 * t_to_zero + 0.5 * j_restore * t_to_zero * t_to_zero;

        double direction = 1.0;
        if (v1 < v_min_feasible - MATH_EPS) {
            direction = -1.0;
        }

        double _v0 = v0 * direction;
        double _a0_scaled = _a0 * direction;
        double _v1 = v1 * direction;

        double t1 = 0.0, t2 = 0.0, t3 = 0.0;

        double t1_max = (a_max - _a0_scaled) / j_max;
        if (t1_max < 0) t1_max = 0.0;

        double t3_max = a_max / j_max;

        double dv_inflection = (_a0_scaled * t1_max + 0.5 * j_max * t1_max * t1_max) +
                               (a_max * t3_max - 0.5 * j_max * t3_max * t3_max);

        double dv_req = _v1 - _v0;

        if (dv_req > dv_inflection) {
            double dv_missing = dv_req - dv_inflection;
            t2 = dv_missing / a_max;
            t1 = t1_max;
            t3 = t3_max;
        } else {
            double term = j_max * dv_req + 0.5 * _a0_scaled * _a0_scaled;
            if (term < 0) term = 0.0;
            double a_peak = std::sqrt(term);

            t1 = (a_peak - _a0_scaled) / j_max;
            t3 = a_peak / j_max;
            if (t1 < 0) t1 = 0.0;
            if (t3 < 0) t3 = 0.0;
        }

        return {t1, t2, t3, direction};
    }

    static inline void _build_vel_profile(TinyProfile<3>& nodes, const MotionState& curr, double v_target, double a_max, double j_max) {
        nodes.clear();
        VelChangeTimes res = _calc_vel_change_times(curr.v, curr.a, v_target, a_max, j_max);

        if (res.t1 > EPS_TIME) nodes.push(res.t1, res.dir * j_max);
        if (res.t2 > EPS_TIME) nodes.push(res.t2, 0.0);
        if (res.t3 > EPS_TIME) nodes.push(res.t3, -res.dir * j_max);
    }

    template <typename Func>
    double _solve_brent(Func&& func, double lower, double upper) {
        double a = lower, b = upper;
        double fa = func(a), fb = func(b);

        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b); std::swap(fa, fb);
        }

        double c = a, fc = fa;
        double d = b - a, e = b - a;

        for (int i = 0; i < SOLVER_ITER; ++i) {
            if (std::abs(fb) < SOLVER_TOL) return b;

            if (std::abs(fc) < std::abs(fb)) {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }

            double xm = 0.5 * (c - b);
            if (std::abs(xm) < SOLVER_TOL) return b;

            if (std::abs(e) >= SOLVER_TOL && std::abs(fa) > std::abs(fb)) {
                double s = fb / fa;
                double p, q;
                if (a == c) {
                    p = 2.0 * xm * s; q = 1.0 - s;
                } else {
                    q = fa / fc; double r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }
                if (p > 0) q = -q;
                p = std::abs(p);
                double min_term = std::min(3.0 * xm * q - std::abs(SOLVER_TOL * q), std::abs(e * q));
                if (2.0 * p < min_term) { e = d; d = p / q; }
                else { d = xm; e = d; }
            } else {
                d = xm; e = d;
            }
            a = b; fa = fb;
            if (std::abs(d) > SOLVER_TOL) b += d;
            else b += (xm > 0 ? SOLVER_TOL : -SOLVER_TOL);
            fb = func(b);
            if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
                c = a; fc = fa; d = e = b - a;
            }
        }
        return b;
    }

    double _solve_via_bisection(const MotionState& curr, double target_p, double target_v, double v_max, double a_max, double j_max) {
        TinyProfile<3> shapes_1, shapes_2;

        auto get_error = [&](double v_mid) -> double {
            MotionState s_sim = curr;
            _build_vel_profile(shapes_1, s_sim, v_mid, a_max, j_max);
            _simulate_endpoint_inplace(s_sim, shapes_1);
            _build_vel_profile(shapes_2, s_sim, target_v, a_max, j_max);
            _simulate_endpoint_inplace(s_sim, shapes_2);
            return s_sim.p - target_p;
        };

        double low = -v_max;
        double high = v_max;

        if (get_error(low) > 0) return -v_max;
        if (get_error(high) < 0) return v_max;

        return _solve_brent(get_error, low, high);
    }

    double _calc_max_reach(const MotionState& curr, double v_limit, double target_v, double a_max, double j_max) {
        TinyProfile<3> shapes;
        MotionState s_sim = curr;
        _build_vel_profile(shapes, s_sim, v_limit, a_max, j_max);
        _simulate_endpoint_inplace(s_sim, shapes);
        _build_vel_profile(shapes, s_sim, target_v, a_max, j_max);
        _simulate_endpoint_inplace(s_sim, shapes);
        return s_sim.p - curr.p;
    }

    // ==========================================================================
    // 4. Trajectory Planning Logic
    // ==========================================================================

    struct CandidateResult {
        bool valid = false;
        double total_duration = std::numeric_limits<double>::infinity();
        double switch_time = 0.0;
    };

    // Corresponds to Python's inner `solve_for_jerk`
    CandidateResult _solve_jerk_specific(const MotionState& curr, double target_p, double target_v,
                                         double a_max, double j_max, double v_max, double j_apply) {
        CandidateResult res;

        // 1. Estimate search horizon
        double t_est = (a_max > 0) ? (std::abs(curr.v) + v_max) / a_max : 1.0;
        double search_horizon = t_est * 2.0 + 5.0;

        // Shared calculation object to avoid reallocation
        TinyProfile<3> shapes_rem;

        // 2. Error function for Brent's method
        auto get_pos_error = [&](double t) -> double {
            if (t < 0) t = 0;
            MotionState s_sim = curr;
            _integrate_saturated_state_only(s_sim, t, j_apply, a_max);
            _build_vel_profile(shapes_rem, s_sim, target_v, a_max, j_max);
            _simulate_endpoint_inplace(s_sim, shapes_rem);
            return s_sim.p - target_p;
        };

        // --- Boundary & Feasibility Checks ---
        double err_0 = get_pos_error(0.0);

        // Check 1: Zero Drift
        if (std::abs(err_0) < EPS_SOLVER) {
            res.switch_time = 0.0;
            res.valid = true;
        }
        // Check 2: Unreachable (Fail Fast)
        else if (err_0 * get_pos_error(search_horizon) > 0) {
            return res; // Invalid
        }
        else {
            // Check 3: Solve Root
            res.switch_time = _solve_brent(get_pos_error, 0.0, search_horizon);

            // Verify precision
            if (std::abs(get_pos_error(res.switch_time)) > EPS_SOLVER) {
                return res; // Invalid
            }
            res.valid = true;
        }

        // --- Calculate Total Duration (for competition) ---
        // Need to reconstruct the profile steps mathematically to sum dt
        MotionState s_switch = curr;
        _integrate_saturated_state_only(s_switch, res.switch_time, j_apply, a_max);

        // Note: s_switch.dt is not used here, we use the solved time.
        // Get the velocity profile duration
        VelChangeTimes vct = _calc_vel_change_times(s_switch.v, s_switch.a, target_v, a_max, j_max);

        // Total duration = saturated_time + t1 + t2 + t3
        res.total_duration = res.switch_time + vct.t1 + vct.t2 + vct.t3;

        return res;
    }

    void _append_fallback_cruise(std::vector<MotionState>& nodes, MotionState curr, double target_p, double target_v, double v_max, double a_max, double j_max) {
        // 1. Find optimal cruise velocity
        double best_v = _solve_via_bisection(curr, target_p, target_v, v_max, a_max, j_max);

        // 2. Accel to best_v
        TinyProfile<3> shapes;
        _build_vel_profile(shapes, curr, best_v, a_max, j_max);
        curr = _append_from_profile(nodes, curr, shapes);

        // Numerical cleanup before cruise
        curr.a = 0.0;

        // 3. Calculate Cruise Duration
        _build_vel_profile(shapes, curr, target_v, a_max, j_max);
        MotionState s_dec_sim = curr;
        _simulate_endpoint_inplace(s_dec_sim, shapes);

        double dist_gap = target_p - s_dec_sim.p;
        double effective_v = (std::abs(curr.v) < MATH_EPS) ? std::copysign(MATH_EPS, dist_gap) : curr.v;

        if (std::abs(dist_gap) > EPS_DIST) {
            double cruise_time = dist_gap / effective_v;
            if (cruise_time > EPS_TIME) {
                MotionState n = curr;
                n.dt = cruise_time;
                n.j = 0.0;
                nodes.push_back(n);
                _integrate_state_inplace(curr, cruise_time, 0.0);
            }
        }

        // 4. Decel to target
        _build_vel_profile(shapes, curr, target_v, a_max, j_max);
        _append_from_profile(nodes, curr, shapes);
    }

    MotionState _append_safety_decel(std::vector<MotionState>& nodes, MotionState curr, double a_max, double j_max) {
        if (std::abs(curr.a) > a_max + EPS_ACC) {
            double j_rec = -std::copysign(j_max, curr.a);
            double tgt_a = std::copysign(a_max, curr.a);
            double t_rec = (curr.a - tgt_a) / (-j_rec);

            if (t_rec > EPS_TIME) {
                MotionState n = curr;
                n.dt = t_rec;
                n.j = j_rec;
                nodes.push_back(n);
                _integrate_state_inplace(curr, t_rec, j_rec);
                curr.a = tgt_a; // Hard clamp
            }
        }
        return curr;
    }

    void _refine_trajectory_precision(std::vector<MotionState>& nodes, const MotionState& start_state, double target_p) {
        if (nodes.empty()) return;

        MotionState sim_s = start_state;
        int correction_idx = -1;
        double max_cruise_dt = -1.0;

        // Phase 1: Simulation and Identification
        for (size_t i = 0; i < nodes.size(); ++i) {
            if (std::abs(nodes[i].j) < MATH_EPS && std::abs(nodes[i].a) < EPS_ACC) {
                if (nodes[i].dt > max_cruise_dt) {
                    max_cruise_dt = nodes[i].dt;
                    correction_idx = (int)i;
                }
            }
            _integrate_state_inplace(sim_s, nodes[i].dt, nodes[i].j);
        }

        double pos_error = target_p - sim_s.p;

        // Phase 2: Correction
        if (std::abs(pos_error) > EPS_DIST && correction_idx != -1) {
            double v_cruise = nodes[correction_idx].v;

            if (std::abs(v_cruise) > EPS_VEL) {
                double dt_fix = pos_error / v_cruise;
                double new_dt = nodes[correction_idx].dt + dt_fix;

                if (new_dt < EPS_TIME) new_dt = EPS_TIME;

                nodes[correction_idx].dt = new_dt;

                // Propagate changes by full re-integration
                MotionState curr = start_state;
                for (auto& node : nodes) {
                    node.p = curr.p;
                    node.v = curr.v;
                    node.a = curr.a;
                    _integrate_state_inplace(curr, node.dt, node.j);
                }
            }
        }
    }

public:
    std::vector<MotionState> plan(const MotionState& start_state, double target_p, double target_v, double v_max, double a_max, double j_max) {
        if (v_max <= 0 || a_max <= 0 || j_max <= 0) return {};

        std::vector<MotionState> final_nodes;
        final_nodes.reserve(32);

        MotionState curr = start_state;

        // 1. Safety Decel
        curr = _append_safety_decel(final_nodes, curr, a_max, j_max);

        // 2. Capacity Check
        double dist_req = target_p - curr.p;
        double d_pos_limit = _calc_max_reach(curr, v_max, target_v, a_max, j_max);
        double d_neg_limit = _calc_max_reach(curr, -v_max, target_v, a_max, j_max);

        bool use_optimal_solver = true;
        if (dist_req > d_pos_limit + EPS_DIST) use_optimal_solver = false;
        if (dist_req < d_neg_limit - EPS_DIST) use_optimal_solver = false;

        // 3. Execution Strategy (Bidirectional Competition)
        bool opt_found = false;

        if (use_optimal_solver) {
            CandidateResult best_res;
            double best_j = 0.0;

            // Competition: Try both +j_max and -j_max
            std::array<double, 2> candidates_j = { j_max, -j_max };

            for (double j_try : candidates_j) {
                CandidateResult res = _solve_jerk_specific(curr, target_p, target_v, a_max, j_max, v_max, j_try);

                if (res.valid) {
                    if (!opt_found || res.total_duration < best_res.total_duration) {
                        best_res = res;
                        best_j = j_try;
                        opt_found = true;
                    }
                }
            }

            // Reconstruction
            if (opt_found) {
                // Phase 1: Variable acceleration (Switching phase)
                curr = _append_saturated_profile(final_nodes, curr, best_res.switch_time, best_j, a_max);

                // Phase 2: Velocity profile to target (Remaining phase)
                TinyProfile<3> shapes_rem;
                _build_vel_profile(shapes_rem, curr, target_v, a_max, j_max);
                _append_from_profile(final_nodes, curr, shapes_rem);
            }
        }

        // Fallback
        if (!opt_found) {
            _append_fallback_cruise(final_nodes, curr, target_p, target_v, v_max, a_max, j_max);
        }

        // 4. Precision Refinement
        _refine_trajectory_precision(final_nodes, start_state, target_p);

        if (final_nodes.empty()) {
            curr.dt = 0.0;
            curr.j = 0.0;
            curr.a = 0.0;
            final_nodes.push_back(curr);
        }
		
        return final_nodes;
    }

    std::vector<MotionState> plan_velocity(const MotionState& start_state, double target_v, double v_max, double a_max, double j_max) {
        if (v_max <= 0 || a_max <= 0 || j_max <= 0) return {};
        std::vector<MotionState> final_nodes;
        final_nodes.reserve(32);
        MotionState curr = start_state;
        curr = _append_safety_decel(final_nodes, curr, a_max, j_max);
        double safe_target_v = std::max(-v_max, std::min(v_max, target_v));
        TinyProfile<3> shapes;
        _build_vel_profile(shapes, curr, safe_target_v, a_max, j_max);
        _append_from_profile(final_nodes, curr, shapes);
        if (final_nodes.empty()) {
            curr.dt = 0.0;
            curr.j = 0.0;
            curr.a = 0.0; 
            final_nodes.push_back(curr);
        }
        return final_nodes;
    }

    MotionState at_time(const std::vector<MotionState>& trajectory, double t) {
        if (trajectory.empty()) return MotionState();

        double elapsed = 0.0;
        for (const auto& node : trajectory) {
            if (t <= elapsed + node.dt) {
                return _integrate_step(node, t - elapsed, node.j);
            }
            elapsed += node.dt;
        }

        const MotionState& last = trajectory.back();
        MotionState final_s = _integrate_step(last, last.dt, last.j);
        return _integrate_step(final_s, t - elapsed, 0.0);
    }
};

} // namespace S7RTT_Lib

#endif // S7RTT_H
