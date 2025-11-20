#ifndef S7RTT_H
#define S7RTT_H

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

namespace S7RTT_Lib {

// ==========================================
// 1. Data Structures
// ==========================================

/**
 * @brief Represents the instantaneous kinematic state of the system.
 */
struct MotionState {
    double p; // Position
    double v; // Velocity
    double a; // Acceleration
    double j; // Jerk

    MotionState(double _p = 0.0, double _v = 0.0, double _a = 0.0, double _j = 0.0)
        : p(_p), v(_v), a(_a), j(_j) {}
};

/**
 * @brief Represents a constant jerk segment (a standard block of an S-Curve).
 */
struct Segment {
    double jerk;
    double duration;

    Segment(double _jerk, double _duration) : jerk(_jerk) {
        duration = (_duration > 0.0) ? _duration : 0.0;
    }
};

// ==========================================
// 2. Mathematical Solver (Templated for Performance)
// ==========================================

class Solver {
public:
    static constexpr double EPS_TOL = 1e-8;
    static constexpr double EPS_CONVERGE = 1e-12;
    static constexpr int MAX_ITER = 50;

    /**
     * @brief Finds the root of a monotonic function using a hybrid method (Brent's method variant).
     *
     * @tparam Func Template type for the callable (Lambda/Functor).
     *              Using templates avoids std::function overhead and allows inlining.
     * @param func The function to solve (returns error value).
     * @param low Lower bound of the search interval.
     * @param high Upper bound of the search interval.
     * @return The input value that results in a function output close to 0.
     */
    template <typename Func>
    static double solve_monotonic_brent(Func&& func, double low, double high) {
        double f_low = func(low);
        double f_high = func(high);

        // Basic bracket check
        if (f_low * f_high > 0) {
            // Root is not bracketed or multiple roots exist. Return closest edge.
            if (std::abs(f_low) < std::abs(f_high)) return low;
            else return high;
        }

        double a = low;
        double b = high;
        double fa = f_low;
        double fb = f_high;

        double root_est = b;

        for (int i = 0; i < MAX_ITER; ++i) {
            // Ensure 'b' is the better guess
            if (std::abs(fa) < std::abs(fb)) {
                std::swap(a, b);
                std::swap(fa, fb);
            }

            // Convergence check
            if (std::abs(fb - fa) < EPS_CONVERGE) {
                root_est = b;
                break;
            }

            // Interpolation (Secant method)
            if (std::abs(fa - fb) > 1e-14) {
                root_est = b - fb * (b - a) / (fb - fa);
            } else {
                root_est = b;
            }

            // Bisection safeguards
            double mid_val = 0.5 * (a + b);
            bool cond1 = (root_est < mid_val && b < mid_val);
            bool cond2 = (root_est > mid_val && b > mid_val);

            if (cond1 || cond2) {
                root_est = mid_val;
            }

            double min_val = (a < b) ? a : b;
            double max_val = (a > b) ? a : b;

            if (root_est < min_val || root_est > max_val) {
                root_est = mid_val;
            }

            double f_root = func(root_est);

            if (std::abs(f_root) < EPS_TOL || std::abs(b - a) < EPS_TOL) {
                return root_est;
            }

            // Update brackets
            if (fa * f_root < 0) {
                b = root_est;
                fb = f_root;
            } else {
                a = root_est;
                fa = f_root;
            }
        }
        return root_est;
    }
};

// ==========================================
// 3. Planner Class (S-Curve / 7-Segment)
// ==========================================

class S7RTT {
public:
    static constexpr double EPS_TIME = 1e-9;
    static constexpr double EPS_VAL = 1e-6;
    static constexpr double EPS_DIST = 1e-5;

    S7RTT() {}

    /**
     * @brief Inline helper to perform kinematic integration for a single step.
     *        Update p, v, a based on jerk j and time t.
     */
    inline static void step_physics(double j, double t, double& p, double& v, double& a) {
        double t2 = t * t;
        double t3 = t2 * t;
        p += v * t + 0.5 * a * t2 + (1.0 / 6.0) * j * t3;
        v += a * t + 0.5 * j * t2;
        a += j * t;
    }

    /**
     * @brief Forward integrates a list of segments from a start state to find the end state.
     */
    MotionState integrate_path(const MotionState& start_state, const std::vector<Segment>& segments) const {
        MotionState curr = start_state;
        for (const auto& seg : segments) {
            double t = seg.duration;
            double j = seg.jerk;
            double t2 = t * t;
            double t3 = t2 * t;

            curr.p += curr.v * t + 0.5 * curr.a * t2 + (1.0 / 6.0) * j * t3;
            curr.v += curr.a * t + 0.5 * j * t2;
            curr.a += j * t;
            curr.j = j;
        }
        return curr;
    }

    /**
     * @brief Main planning function.
     *
     * @param start_state Current P, V, A.
     * @param target_p Target Position.
     * @param target_v Target Velocity (usually 0).
     * @param v_max Maximum Velocity limit.
     * @param a_max Maximum Acceleration limit.
     * @param j_max Maximum Jerk limit.
     * @return std::vector<Segment> Sequence of jerk segments to execute.
     */
    std::vector<Segment> plan(const MotionState& start_state, double target_p, double target_v, double v_max, double a_max, double j_max) {
        if (v_max <= 0 || a_max <= 0 || j_max <= 0) return {};

        std::vector<Segment> pre_segments;
        pre_segments.reserve(2); // Reserve for potential safety deceleration segments

        MotionState current_state = start_state;

        // 1. Initial State Sanitization
        // If current acceleration exceeds limits, generate a segment to bring it back within bounds.
        if (current_state.a > a_max + EPS_VAL) {
            double t_recover = (current_state.a - a_max) / j_max;
            Segment seg(-j_max, t_recover);
            pre_segments.push_back(seg);
            step_physics(-j_max, t_recover, current_state.p, current_state.v, current_state.a);
            current_state.a = a_max; // Clamp strictly
        }
        else if (current_state.a < -a_max - EPS_VAL) {
            double t_recover = (-a_max - current_state.a) / j_max;
            Segment seg(j_max, t_recover);
            pre_segments.push_back(seg);
            step_physics(j_max, t_recover, current_state.p, current_state.v, current_state.a);
            current_state.a = -a_max; // Clamp strictly
        }

        double dist_req = target_p - current_state.p;

        // 2. Compute Physical Boundaries (Full Calculation)
        // Calculate the distance traveled if we accelerate to +V_max and -V_max immediately.
        double d_upper;
        std::vector<Segment> acc_up, dec_up;
        compute_trajectory_full(current_state, v_max, target_v, a_max, j_max, d_upper, acc_up, dec_up);

        double d_lower;
        std::vector<Segment> acc_lo, dec_lo;
        compute_trajectory_full(current_state, -v_max, target_v, a_max, j_max, d_lower, acc_lo, dec_lo);

        std::vector<Segment> plan_segments;
        plan_segments.reserve(7); // S-Curve typically has max 7 segments

        // 3. Three-Branch Decision Logic
        if (dist_req > d_upper + EPS_DIST) {
            // Case 1: Positive Saturation (Cruising at +V_max)
            double gap = dist_req - d_upper;
            double t_cruise = gap / v_max;

            plan_segments.insert(plan_segments.end(), acc_up.begin(), acc_up.end());
            if (t_cruise > EPS_TIME) plan_segments.emplace_back(0.0, t_cruise);
            plan_segments.insert(plan_segments.end(), dec_up.begin(), dec_up.end());
        }
        else if (dist_req < d_lower - EPS_DIST) {
            // Case 2: Negative Saturation (Cruising at -V_max)
            double gap = dist_req - d_lower;
            double t_cruise = gap / (-v_max);

            plan_segments.insert(plan_segments.end(), acc_lo.begin(), acc_lo.end());
            if (t_cruise > EPS_TIME) plan_segments.emplace_back(0.0, t_cruise);
            plan_segments.insert(plan_segments.end(), dec_lo.begin(), dec_lo.end());
        }
        else {
            // Case 3: Peak Search (Performance Bottleneck)
            // The target is reachable without hitting V_max cruise, but we need to find the specific Peak Velocity.

            // Lambda for the solver: Returns (Calculated_Distance - Required_Distance)
            auto get_error = [&](double v_p) -> double {
                // FAST PATH: Only calculate scalar distance, no vector allocation
                double d = calc_trajectory_dist_only(current_state, v_p, target_v, a_max, j_max);
                return d - dist_req;
            };

            double max_abs_v = std::max({ v_max, std::abs(current_state.v), std::abs(target_v) }) * 2.0;

            // Solve for the exact peak velocity
            double best_v = Solver::solve_monotonic_brent(get_error, -max_abs_v, max_abs_v);

            // Generate actual segments using the found velocity (SLOW PATH)
            double d_final;
            std::vector<Segment> acc_fin, dec_fin;
            compute_trajectory_full(current_state, best_v, target_v, a_max, j_max, d_final, acc_fin, dec_fin);

            plan_segments.insert(plan_segments.end(), acc_fin.begin(), acc_fin.end());
            plan_segments.insert(plan_segments.end(), dec_fin.begin(), dec_fin.end());
        }

        // Merge any safety recovery segments from the start
        pre_segments.insert(pre_segments.end(), plan_segments.begin(), plan_segments.end());
        return pre_segments;
    }

private:
    // ==========================================================
    // FAST PATH: Scalar calculations only, no heap allocation.
    // Used repeatedly inside the Solver loop.
    // ==========================================================

    /**
     * @brief Calculates the distance traveled for a single velocity change profile (e.g., Start -> Peak).
     *
     * @return Accumulated distance.
     */
    double calc_profile_dist_only(double v_start, double a_start, double v_target,
                                  double a_max, double j_max,
                                  double& v_end_out, double& a_end_out)
    {
        // 1. Clamp Initial Acceleration
        double acc_clamped = a_start;
        if (acc_clamped > a_max) acc_clamped = a_max;
        else if (acc_clamped < -a_max) acc_clamped = -a_max;

        // 2. Calculate Base Velocity (velocity reached if we reduce accel to 0 immediately)
        double t_to_zero = std::abs(acc_clamped) / j_max;
        double j_to_zero = (acc_clamped > 0) ? -j_max : j_max;
        double dv_base = acc_clamped * t_to_zero + 0.5 * j_to_zero * t_to_zero * t_to_zero;
        double v_base = v_start + dv_base;

        // Determine Direction
        double direction = (v_target < v_base) ? -1.0 : 1.0;
        double j_up = (direction > 0) ? j_max : -j_max;
        double j_down = (direction > 0) ? -j_max : j_max;

        // 3. Solve Physics Constraints
        double v_req_total = v_target - v_start;
        double a_limit = a_max * direction;

        // Time to ramp acceleration from current to limit, and limit to zero
        double t1_max = (a_limit - acc_clamped) / j_up;
        double t3_max = std::abs(a_limit) / j_max;

        // Velocity change provided by a full trapezoidal acceleration profile
        double dv_trapezoid = (acc_clamped * t1_max + 0.5 * j_up * t1_max * t1_max) +
                              (a_limit * t3_max + 0.5 * j_down * t3_max * t3_max);

        // Check if we hit the acceleration ceiling (Trapezoidal vs Triangular profile)
        bool needs_flat = (direction > 0) ? (v_req_total > dv_trapezoid) : (v_req_total < dv_trapezoid);

        // Accumulators
        double dist_accum = 0.0;
        double v_curr = v_start;
        double a_curr = acc_clamped;
        // Note: We ignore the tiny displacement of clamping a_start here, assuming it's handled externally
        // or negligible for the root finder's purpose.

        // Lambda to step scalar physics
        auto fast_step = [&](double j, double t) {
            if (t <= EPS_TIME) return;
            double t2 = t * t;
            dist_accum += v_curr * t + 0.5 * a_curr * t2 + (1.0 / 6.0) * j * t2 * t;
            v_curr += a_curr * t + 0.5 * j * t2;
            a_curr += j * t;
        };

        if (needs_flat) {
            // Case: Acceleration saturates (Trapezoidal Accel Profile)
            double v_missing = v_req_total - dv_trapezoid;
            double t_flat = v_missing / a_limit;

            fast_step(j_up, t1_max);    // Ramp up accel
            fast_step(0.0, t_flat);     // Constant accel
            fast_step(j_down, t3_max);  // Ramp down accel
        } else {
            // Case: Acceleration peaks before saturation (Triangular Accel Profile)
            // Solve: v_req = Integral of accel (area under triangle)
            double term = v_req_total * j_up + 0.5 * acc_clamped * acc_clamped;
            if (term < 0) term = 0.0;
            double a_peak_mag = std::sqrt(term);
            double a_peak = (direction > 0) ? a_peak_mag : -a_peak_mag;

            double t1 = (a_peak - acc_clamped) / j_up;
            double t3 = (0.0 - a_peak) / j_down;

            fast_step(j_up, t1);
            fast_step(j_down, t3);
        }

        v_end_out = v_curr;
        a_end_out = a_curr;
        return dist_accum;
    }

    /**
     * @brief Calculates total distance for a full trajectory (Start -> V_peak -> Target) without generating segments.
     */
    double calc_trajectory_dist_only(const MotionState& start_state, double v_peak, double target_v,
                                     double a_max, double j_max)
    {
        double v_mid, a_mid;
        // Phase 1: Accelerate/Decelerate to Peak Velocity
        double d1 = calc_profile_dist_only(start_state.v, start_state.a, v_peak, a_max, j_max, v_mid, a_mid);

        // Phase 2: Accelerate/Decelerate from Peak Velocity to Target Velocity
        // Note: Use the end state of Phase 1 as start of Phase 2
        double v_end, a_end;
        double d2 = calc_profile_dist_only(v_mid, a_mid, target_v, a_max, j_max, v_end, a_end);

        return d1 + d2;
    }


    // ==========================================================
    // SLOW PATH: Generates actual Segment objects.
    // Called only once per plan when the optimal solution is found.
    // ==========================================================

    std::vector<Segment> _build_profile_segs(double v_start, double a_start, double v_target, double a_max, double j_max) {
        std::vector<Segment> segments;
        segments.reserve(3);

        double acc_clamped = a_start;
        if (acc_clamped > a_max) acc_clamped = a_max;
        else if (acc_clamped < -a_max) acc_clamped = -a_max;

        // Calculate theoretical velocity if we drop accel to zero now
        double t_to_zero = std::abs(acc_clamped) / j_max;
        double j_to_zero = (acc_clamped > 0) ? -j_max : j_max;
        double dv_base = acc_clamped * t_to_zero + 0.5 * j_to_zero * t_to_zero * t_to_zero;
        double v_base = v_start + dv_base;

        double direction = (v_target < v_base) ? -1.0 : 1.0;
        double j_up = (direction > 0) ? j_max : -j_max;
        double j_down = (direction > 0) ? -j_max : j_max;

        double v_req_total = v_target - v_start;
        double a_limit = a_max * direction;

        double t1_max = (a_limit - acc_clamped) / j_up;
        double t3_max = std::abs(a_limit) / j_max;

        double dv_trapezoid = (acc_clamped * t1_max + 0.5 * j_up * t1_max * t1_max) +
                              (a_limit * t3_max + 0.5 * j_down * t3_max * t3_max);

        bool needs_flat = (direction > 0) ? (v_req_total > dv_trapezoid) : (v_req_total < dv_trapezoid);

        if (needs_flat) {
            double v_missing = v_req_total - dv_trapezoid;
            double t_flat = v_missing / a_limit;
            if (t1_max > EPS_TIME) segments.emplace_back(j_up, t1_max);
            if (t_flat > EPS_TIME) segments.emplace_back(0.0, t_flat);
            if (t3_max > EPS_TIME) segments.emplace_back(j_down, t3_max);
        } else {
            double term = v_req_total * j_up + 0.5 * acc_clamped * acc_clamped;
            if (term < 0) term = 0.0;
            double a_peak_mag = std::sqrt(term);
            double a_peak = (direction > 0) ? a_peak_mag : -a_peak_mag;
            double t1 = (a_peak - acc_clamped) / j_up;
            double t3 = (0.0 - a_peak) / j_down;
            if (t1 > EPS_TIME) segments.emplace_back(j_up, t1);
            if (t3 > EPS_TIME) segments.emplace_back(j_down, t3);
        }
        return segments;
    }

    void compute_trajectory_full(const MotionState& start_state, double v_peak, double target_v, double a_max, double j_max,
                                 double& out_dist, std::vector<Segment>& out_acc_segs, std::vector<Segment>& out_dec_segs) {

        // Build acceleration phase
        out_acc_segs = _build_profile_segs(start_state.v, start_state.a, v_peak, a_max, j_max);

        // Integrate to find state at the peak
        MotionState state_mid = integrate_path(start_state, out_acc_segs);

        // Build deceleration phase
        out_dec_segs = _build_profile_segs(state_mid.v, state_mid.a, target_v, a_max, j_max);

        // Calculate total distance by summing segments (Double verification)
        double d1 = 0;
        double v = start_state.v, a = start_state.a;
        for(auto& s : out_acc_segs) {
            double t = s.duration; double j = s.jerk;
            d1 += v*t + 0.5*a*t*t + (1.0/6.0)*j*t*t*t;
            v += a*t + 0.5*j*t*t; a += j*t;
        }

        double d2 = 0;
        v = state_mid.v; a = state_mid.a;
        for(auto& s : out_dec_segs) {
            double t = s.duration; double j = s.jerk;
            d2 += v*t + 0.5*a*t*t + (1.0/6.0)*j*t*t*t;
            v += a*t + 0.5*j*t*t; a += j*t;
        }

        out_dist = d1 + d2;
    }
};

} // namespace S7RTT_Lib

#endif // S7RTT_H
