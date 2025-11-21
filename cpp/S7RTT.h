// ==============================================================================
// File Name:    S7RTT.h
// Author:       feecat
// Version:      V1.3
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
#include <limits>
#include <algorithm>
#include <iostream>
#include <tuple>
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
        : dt((_dt > 0.0) ? _dt : 0.0), p(_p), v(_v), a(_a), j(_j) {}

    friend std::ostream& operator<<(std::ostream& os, const MotionState& s) {
        os << "State(dt=" << s.dt << ", P=" << s.p << ", V=" << s.v
           << ", A=" << s.a << ", J=" << s.j << ")";
        return os;
    }
};

// ==============================================================================
// Internal Lightweight Structures (Stack Allocated)
// ==============================================================================
// Represents a single simplified segment for calculation (time + jerk)
struct SegmentLite {
    double dt;
    double j;
};

// Represents an acceleration or deceleration phase (max 3 segments: ramp up, flat acc, ramp down)
struct ProfilePlan {
    int count = 0;
    std::array<SegmentLite, 3> segs; // Fixed size, no heap allocation

    void add(double t, double jerk) {
        if (t > 1e-9) {
            segs[count++] = {t, jerk};
        }
    }
};

// ==============================================================================
// Solver (Templated Brent's Method)
// ==============================================================================
class Solver {
public:
    static constexpr double TOL = 1e-12;
    static constexpr int MAX_ITER = 100;

    // Templated to allow inlining of the lambda/functor
    template <typename Func>
    static double solve_brent(Func&& func, double low, double high) {
        double a = low;
        double b = high;
        double fa = func(a);
        double fb = func(b);

        const double MACH_EPS = std::numeric_limits<double>::epsilon();

        if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) {
            return (std::abs(fa) < std::abs(fb)) ? a : b;
        }

        if (fa == 0) return a;
        if (fb == 0) return b;

        double c = a;
        double fc = fa;
        double d = b - a;
        double e = b - a;

        bool mflag = true;

        for (int i = 0; i < MAX_ITER; ++i) {
            if (std::abs(fc) < std::abs(fb)) {
                a = b; b = c; c = a;
                fa = fb; fb = fc; fc = fa;
            }

            double tol1 = 2.0 * MACH_EPS * std::abs(b) + 0.5 * TOL;
            double xm = 0.5 * (c - b);

            if (std::abs(xm) <= tol1 || fb == 0) {
                return b;
            }

            if (std::abs(e) >= tol1 && std::abs(fa) > std::abs(fb)) {
                double s = fb / fa;
                double p, q;

                if (a == c) {
                    p = 2.0 * xm * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    double r = fb / fc;
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0) q = -q;
                p = std::abs(p);

                double min1 = 3.0 * xm * q - std::abs(tol1 * q);
                double min2 = std::abs(e * q);

                if (2.0 * p < min1) {
                    if (mflag) {
                        if (2.0 * p < min2) { e = d; d = p / q; mflag = false; }
                        else { d = xm; e = d; mflag = true; }
                    } else {
                        if (2.0 * p < std::abs(d * q)) { e = d; d = p / q; mflag = false; }
                        else { d = xm; e = d; mflag = true; }
                    }
                } else { d = xm; e = d; mflag = true; }
            } else { d = xm; e = d; mflag = true; }

            a = b; fa = fb;

            if (std::abs(d) > tol1) b += d;
            else b += std::copysign(tol1, xm);

            fb = func(b);

            if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
                c = a; fc = fa;
                d = e = b - a;
                mflag = true;
            }
        }
        return b;
    }
};

// ==============================================================================
// S7RTT Class
// ==============================================================================
class S7RTT {
private:
    static constexpr double EPS_TIME = 1e-9;
    static constexpr double EPS_VAL  = 1e-6;
    static constexpr double EPS_DIST = 1e-5;
    static constexpr double ONE_SIXTH = 1.0 / 6.0;

    // Helper to integrate state forward without creating objects
    static void _integrate_inplace(double& p, double& v, double& a, double dt, double j) {
        double dt2 = dt * dt;
        double dt3 = dt2 * dt;
        p += v * dt + 0.5 * a * dt2 + ONE_SIXTH * j * dt3;
        v += a * dt + 0.5 * j * dt2;
        a += j * dt;
    }

    // Calculates the shape parameters (times and jerks) without allocating vectors.
    static ProfilePlan _compute_profile_params(double v_start, double a_start, double v_target, double a_max, double j_max) {
        ProfilePlan plan;

        // 1. Clamp start acceleration
        double acc_clamped = std::clamp(a_start, -a_max, a_max);

        // 2. Calc velocity reached if we immediately reduce accel to zero
        double t_to_zero = std::abs(acc_clamped) / j_max;
        double j_to_zero = (acc_clamped > 0) ? -j_max : j_max;
        double dv_base = acc_clamped * t_to_zero + 0.5 * j_to_zero * t_to_zero * t_to_zero;
        double v_base = v_start + dv_base;

        // 3. Determine direction
        double direction = (v_target < v_base) ? -1.0 : 1.0;

        double j_up = (direction > 0) ? j_max : -j_max;
        double j_down = (direction > 0) ? -j_max : j_max;

        // 4. Constraints
        double v_req_total = v_target - v_start;
        double a_limit = a_max * direction;

        double t1_max = (a_limit - acc_clamped) / j_up;
        double t3_max = std::abs(a_limit) / j_max;

        // Safety for precision issues
        if (t1_max < 0) t1_max = 0;

        double dv_trapezoid = (acc_clamped * t1_max + 0.5 * j_up * t1_max * t1_max) +
                              (a_limit * t3_max + 0.5 * j_down * t3_max * t3_max);

        bool needs_flat = (direction > 0) ? (v_req_total > dv_trapezoid) : (v_req_total < dv_trapezoid);

        if (needs_flat) {
            // Trapezoidal
            double v_missing = v_req_total - dv_trapezoid;
            double t_flat = v_missing / a_limit;
            if (t_flat < 0) t_flat = 0;

            plan.add(t1_max, j_up);
            plan.add(t_flat, 0.0);
            plan.add(t3_max, j_down);
        } else {
            // Triangular
            double term = v_req_total * j_up + 0.5 * acc_clamped * acc_clamped;
            if (term < 0) term = 0.0; // Numerical noise safety
            double a_peak_mag = std::sqrt(term);
            double a_peak = (direction > 0) ? a_peak_mag : -a_peak_mag;

            double t1 = (a_peak - acc_clamped) / j_up;
            double t3 = (0.0 - a_peak) / j_down;

            if (t1 < 0) t1 = 0; // Safety

            plan.add(t1, j_up);
            plan.add(t3, j_down);
        }
        return plan;
    }

    // Returns total distance covered by a plan phase, updates final v and a
    static double _integrate_plan_dist(double& v, double& a, const ProfilePlan& plan) {
        double dist = 0.0;
        // Local integrate to get distance
        for (int i = 0; i < plan.count; ++i) {
            double t = plan.segs[i].dt;
            double j = plan.segs[i].j;
            double dt2 = t * t;

            dist += v * t + 0.5 * a * dt2 + ONE_SIXTH * j * dt2 * t;
            v += a * t + 0.5 * j * dt2;
            a += j * t;
        }
        return dist;
    }

    // Calculates total distance for specific v_peak without ANY heap allocation
    static double _calc_distance_for_v_peak(const MotionState& current_state, double v_peak, double target_v, double a_max, double j_max) {
        double v_curr = current_state.v;
        double a_curr = current_state.a;

        // Phase 1: Start -> V_Peak
        ProfilePlan acc_plan = _compute_profile_params(v_curr, a_curr, v_peak, a_max, j_max);
        double d1 = _integrate_plan_dist(v_curr, a_curr, acc_plan);

        // Phase 2: V_Peak -> Target
        ProfilePlan dec_plan = _compute_profile_params(v_curr, a_curr, target_v, a_max, j_max);
        double d2 = _integrate_plan_dist(v_curr, a_curr, dec_plan);

        return d1 + d2;
    }

    // Helper to actually construct the vectors only when needed
    void _append_to_trajectory(std::vector<MotionState>& traj, MotionState& curr, const ProfilePlan& plan) {
        for (int i = 0; i < plan.count; ++i) {
            MotionState s = curr;
            s.dt = plan.segs[i].dt;
            s.j = plan.segs[i].j;
            traj.push_back(s);
            // Update current state for next node
            curr = _integrate_state(curr, s.dt, s.j);
        }
    }

    MotionState _integrate_state(const MotionState& state, double dt, double j) {
        double p = state.p;
        double v = state.v;
        double a = state.a;
        _integrate_inplace(p, v, a, dt, j);
        return MotionState(0.0, p, v, a, 0.0);
    }

public:
    S7RTT() {}

    std::vector<MotionState> plan(const MotionState& start_state, double target_p, double target_v, double v_max, double a_max, double j_max) {
        if (v_max <= 0 || a_max <= 0 || j_max <= 0) return {};

        // Pre-reserve to avoid re-allocations. 7 is typical max segments (recover + acc(3) + cruise + dec(3))
        std::vector<MotionState> final_trajectory;
        final_trajectory.reserve(16);

        MotionState current_state = start_state;

        // 2. Acceleration Recovery
        double t_recover = 0.0;
        double j_recover = 0.0;
        bool recovering = false;

        if (current_state.a > a_max + EPS_VAL) {
            t_recover = (current_state.a - a_max) / j_max;
            j_recover = -j_max;
            recovering = true;
        } else if (current_state.a < -a_max - EPS_VAL) {
            t_recover = (-a_max - current_state.a) / j_max;
            j_recover = j_max;
            recovering = true;
        }

        if (recovering) {
            MotionState rec_state = current_state;
            rec_state.dt = t_recover;
            rec_state.j = j_recover;
            final_trajectory.push_back(rec_state);
            current_state = _integrate_state(current_state, t_recover, j_recover);
            // Force snap to limit to remove float noise
            current_state.a = (j_recover < 0) ? a_max : -a_max;
        }

        double dist_req = target_p - current_state.p;

        // 3. Inertial Reference (Velocity reached if accel goes to 0)
        double t_to_zero = std::abs(current_state.a) / j_max;
        double j_to_zero = (current_state.a > 0) ? -j_max : j_max;
        double v_inertial = current_state.v + current_state.a * t_to_zero + 0.5 * j_to_zero * t_to_zero * t_to_zero;

        // 4. Direct Profile Check (Can we just go straight to target V?)
        double d_direct = _calc_distance_for_v_peak(current_state, current_state.v, target_v, a_max, j_max);

        // Variables for result construction
        std::vector<ProfilePlan> plans_to_stitch;
        double cruise_time = 0.0;

        if (std::abs(d_direct - dist_req) <= EPS_DIST) {
            // Case A: Go directly to target velocity
            ProfilePlan p = _compute_profile_params(current_state.v, current_state.a, target_v, a_max, j_max);
            plans_to_stitch.push_back(p);
        } else {
            // 5. Boundary Calc
            double dist_upper = _calc_distance_for_v_peak(current_state, v_max, target_v, a_max, j_max);
            double dist_lower = _calc_distance_for_v_peak(current_state, -v_max, target_v, a_max, j_max);

            if (dist_req > dist_upper + EPS_DIST) {
                // Case B: Cruise at +V_max
                double gap = dist_req - dist_upper;
                cruise_time = gap / v_max;
                plans_to_stitch.push_back(_compute_profile_params(current_state.v, current_state.a, v_max, a_max, j_max));
                plans_to_stitch.push_back(_compute_profile_params(v_max, 0.0, target_v, a_max, j_max));
            } else if (dist_req < dist_lower - EPS_DIST) {
                // Case C: Cruise at -V_max
                double gap = dist_req - dist_lower;
                cruise_time = gap / -v_max; // gap is negative, -v_max is positive divisor result
                plans_to_stitch.push_back(_compute_profile_params(current_state.v, current_state.a, -v_max, a_max, j_max));
                plans_to_stitch.push_back(_compute_profile_params(-v_max, 0.0, target_v, a_max, j_max));
            } else {
                // Case D: Within bounds, use Brent's Solver
                double dist_inertial = _calc_distance_for_v_peak(current_state, v_inertial, target_v, a_max, j_max);

                double s_low, s_high;
                double limit_safe = std::max({v_max, std::abs(current_state.v), std::abs(target_v)}) * 1.5 + 10.0;
                double overlap_v = EPS_VAL * 10.0;

                if (dist_req > dist_inertial) {
                    s_low = v_inertial - overlap_v;
                    s_high = limit_safe;
                } else {
                    s_high = v_inertial + overlap_v;
                    s_low = -limit_safe;
                }

                if (s_low >= s_high) s_low = s_high - EPS_VAL;

                // Capture params by copy to avoid reference issues, lightweight enough
                auto cost_func = [current_state, target_v, a_max, j_max, dist_req](double v_p) -> double {
                    return _calc_distance_for_v_peak(current_state, v_p, target_v, a_max, j_max) - dist_req;
                };

                double best_v = Solver::solve_brent(cost_func, s_low, s_high);

                // Fallback check (sometimes roots are outside predicted bracket if j makes weird shapes)
                bool hit_low = std::abs(best_v - s_low) < EPS_VAL;
                bool hit_high = std::abs(best_v - s_high) < EPS_VAL;

                if (hit_low || hit_high) {
                    best_v = Solver::solve_brent(cost_func, -limit_safe, limit_safe);
                }

                // Construct final shape based on best_v
                double v_mid_final = best_v;
                double a_mid_final = 0.0; // Approximation, need precise state

                // Recalculate parameters for the best_v
                ProfilePlan p1 = _compute_profile_params(current_state.v, current_state.a, v_mid_final, a_max, j_max);

                // We need the state at the peak to calculate the second phase precisely
                double v_temp = current_state.v;
                double a_temp = current_state.a;
                _integrate_plan_dist(v_temp, a_temp, p1); // updates v_temp, a_temp

                ProfilePlan p2 = _compute_profile_params(v_temp, a_temp, target_v, a_max, j_max);

                plans_to_stitch.push_back(p1);
                plans_to_stitch.push_back(p2);
            }
        }

        // 6. Stitching
        // First profile (Accel)
        if (!plans_to_stitch.empty()) {
            _append_to_trajectory(final_trajectory, current_state, plans_to_stitch[0]);
        }

        // Cruise
        if (cruise_time > EPS_TIME) {
            MotionState cruise = current_state;
            cruise.dt = cruise_time;
            cruise.j = 0.0;
            // Clean up small acceleration noise
            cruise.a = 0.0;
            final_trajectory.push_back(cruise);
            current_state = _integrate_state(current_state, cruise_time, 0.0);
        }

        // Second profile (Decel) - if exists
        if (plans_to_stitch.size() > 1) {
            _append_to_trajectory(final_trajectory, current_state, plans_to_stitch[1]);
        }

        return final_trajectory;
    }

    // Velocity planner (simpler, no solver needed)
    std::vector<MotionState> plan_velocity(const MotionState& start_state, double target_v, double v_max, double a_max, double j_max) {
        if (v_max <= 0 || a_max <= 0 || j_max <= 0) return {};

        std::vector<MotionState> final_trajectory;
        final_trajectory.reserve(8);

        double safe_target_v = std::clamp(target_v, -v_max, v_max);
        MotionState current_state = start_state;

        // 1. Recovery
        if (current_state.a > a_max + EPS_VAL) {
            double t_recover = (current_state.a - a_max) / j_max;
            MotionState rec_state = current_state;
            rec_state.dt = t_recover;
            rec_state.j = -j_max;
            final_trajectory.push_back(rec_state);
            current_state = _integrate_state(current_state, t_recover, -j_max);
            current_state.a = a_max;
        } else if (current_state.a < -a_max - EPS_VAL) {
            double t_recover = (-a_max - current_state.a) / j_max;
            MotionState rec_state = current_state;
            rec_state.dt = t_recover;
            rec_state.j = j_max;
            final_trajectory.push_back(rec_state);
            current_state = _integrate_state(current_state, t_recover, j_max);
            current_state.a = -a_max;
        }

        // 2. Build Profile
        ProfilePlan plan = _compute_profile_params(current_state.v, current_state.a, safe_target_v, a_max, j_max);
        _append_to_trajectory(final_trajectory, current_state, plan);

        return final_trajectory;
    }

    MotionState at_time(const std::vector<MotionState>& trajectory, double dt) {
        if (trajectory.empty()) return MotionState();

        double t_remaining = dt;

        // Optimized loop: avoid recalculating constants if possible, but logic is simple enough
        for (const auto& node : trajectory) {
            if (t_remaining <= node.dt + 1e-12) {
                // Ensure negative time is handled as 0
                double t = (t_remaining < 0) ? 0 : t_remaining;
                double j = node.j;
                double dt2 = t * t;

                double p = node.p + node.v * t + 0.5 * node.a * dt2 + ONE_SIXTH * j * dt2 * t;
                double v = node.v + node.a * t + 0.5 * j * dt2;
                double a = node.a + j * t;
                return MotionState(0.0, p, v, a, j);
            } else {
                t_remaining -= node.dt;
            }
        }

        // Extrapolate
        const MotionState& last_node = trajectory.back();
        double t = last_node.dt;
        double j = last_node.j;
        double dt2 = t * t;

        double p_end = last_node.p + last_node.v * t + 0.5 * last_node.a * dt2 + ONE_SIXTH * j * dt2 * t;
        double v_end = last_node.v + last_node.a * t + 0.5 * j * dt2;

        double p_final = p_end + v_end * t_remaining;
        return MotionState(0.0, p_final, v_end, 0.0, 0.0);
    }
};

} // namespace S7RTT_Lib

#endif // S7RTT_H
