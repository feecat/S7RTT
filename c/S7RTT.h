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
// Usage:
// 1. In ONE .c file, define S7RTT_IMPLEMENTATION before including this header:
//    #define S7RTT_IMPLEMENTATION
//    #include "S7RTT.h"
//
// 2. In other files, just #include "S7RTT.h"
// ==============================================================================


#ifndef S7RTT_H
#define S7RTT_H

#ifdef __cplusplus
extern "C" {
#endif

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <float.h>
#include <stdbool.h>

// ==============================================================================
// Configuration & Typedefs
// ==============================================================================

// Define S7RTT_REAL to float if your MCU lacks double FPU
#ifndef S7RTT_REAL
#define S7RTT_REAL double
#endif

// Maximum segments in a trajectory (Recovery + Acc(3) + Cruise + Dec(3) = ~8)
// 16 is a safe upper bound.
#define S7RTT_MAX_SEGS 16

typedef struct {
    S7RTT_REAL dt; // Duration of this segment
    S7RTT_REAL p;  // Position
    S7RTT_REAL v;  // Velocity
    S7RTT_REAL a;  // Acceleration
    S7RTT_REAL j;  // Jerk
} S7RTT_State;

// ==============================================================================
// API Function Prototypes
// ==============================================================================

/**
 * @brief Plan a position trajectory (S-Curve)
 *
 * @param start     Current state (p, v, a)
 * @param target_p  Target position
 * @param target_v  Target final velocity
 * @param v_max     Maximum velocity limit (>0)
 * @param a_max     Maximum acceleration limit (>0)
 * @param j_max     Maximum jerk limit (>0)
 * @param out_buf   Pointer to an array of S7RTT_State to store segments
 * @param max_len   Size of the out_buf array
 * @return int      Number of segments generated, or -1 on error/buffer overflow
 */
int S7RTT_Plan(const S7RTT_State* start,
               S7RTT_REAL target_p, S7RTT_REAL target_v,
               S7RTT_REAL v_max, S7RTT_REAL a_max, S7RTT_REAL j_max,
               S7RTT_State* out_buf, int max_len);

/**
 * @brief Plan a velocity only trajectory
 */
int S7RTT_PlanVelocity(const S7RTT_State* start,
                       S7RTT_REAL target_v,
                       S7RTT_REAL v_max, S7RTT_REAL a_max, S7RTT_REAL j_max,
                       S7RTT_State* out_buf, int max_len);

/**
 * @brief Sample the trajectory at a specific time delta from start
 *
 * @param trajectory Array of segments returned by Plan
 * @param count      Number of segments
 * @param t          Time to sample
 * @param out_state  Resulting state
 */
void S7RTT_AtTime(const S7RTT_State* trajectory, int count, S7RTT_REAL t, S7RTT_State* out_state);

// ==============================================================================
// Implementation Section
// ==============================================================================
#ifdef S7RTT_IMPLEMENTATION

// Constants
#define S7_EPS_TIME 1e-9
#define S7_EPS_VAL  1e-6
#define S7_EPS_DIST 1e-5
#define S7_ONE_SIXTH (1.0 / 6.0)

// Utils
static inline S7RTT_REAL s7_fmax(S7RTT_REAL a, S7RTT_REAL b) { return (a > b) ? a : b; }
static inline S7RTT_REAL s7_fmin(S7RTT_REAL a, S7RTT_REAL b) { return (a < b) ? a : b; }
static inline S7RTT_REAL s7_clamp(S7RTT_REAL v, S7RTT_REAL min, S7RTT_REAL max) {
    if (v < min) return min;
    if (v > max) return max;
    return v;
}
static inline S7RTT_REAL s7_sign(S7RTT_REAL x) {
    return (x < 0.0) ? -1.0 : 1.0;
}

// Internal lightweight structures
typedef struct {
    S7RTT_REAL dt;
    S7RTT_REAL j;
} S7_SegLite;

typedef struct {
    int count;
    S7_SegLite segs[3];
} S7_ProfilePlan;

// ------------------------------------------------------------------------------
// Solver (Brent's Method) - Refactored for C
// ------------------------------------------------------------------------------
typedef S7RTT_REAL (*S7_CostFunc)(S7RTT_REAL x, void* user_data);

static S7RTT_REAL s7_solve_brent(S7_CostFunc func, S7RTT_REAL low, S7RTT_REAL high, void* user_data) {
    S7RTT_REAL a = low, b = high;
    S7RTT_REAL fa = func(a, user_data);
    S7RTT_REAL fb = func(b, user_data);

    const int MAX_ITER = 50; // Reduced for embedded safety
    const S7RTT_REAL TOL = 1e-9;
    const S7RTT_REAL MACH_EPS = DBL_EPSILON; // from float.h

    if ((fa > 0 && fb > 0) || (fa < 0 && fb < 0)) {
        return (fabs(fa) < fabs(fb)) ? a : b;
    }
    if (fa == 0) return a;
    if (fb == 0) return b;

    S7RTT_REAL c = a, fc = fa, d = b - a, e = b - a;
    bool mflag = true;

    for (int i = 0; i < MAX_ITER; ++i) {
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }

        S7RTT_REAL tol1 = 2.0 * MACH_EPS * fabs(b) + 0.5 * TOL;
        S7RTT_REAL xm = 0.5 * (c - b);

        if (fabs(xm) <= tol1 || fb == 0) return b;

        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            S7RTT_REAL s = fb / fa;
            S7RTT_REAL p, q;
            if (a == c) {
                p = 2.0 * xm * s; q = 1.0 - s;
            } else {
                q = fa / fc; S7RTT_REAL r = fb / fc;
                p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if (p > 0) q = -q;
            p = fabs(p);
            S7RTT_REAL min1 = 3.0 * xm * q - fabs(tol1 * q);
            S7RTT_REAL min2 = fabs(e * q);
            if (2.0 * p < min1) {
                if (2.0 * p < min2) { e = d; d = p / q; }
                else { d = xm; e = d; }
            } else { d = xm; e = d; }
        } else { d = xm; e = d; }

        a = b; fa = fb;
        if (fabs(d) > tol1) b += d;
        else b += copysign(tol1, xm);
        fb = func(b, user_data);

        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a; fc = fa; d = e = b - a; mflag = true;
        }
    }
    return b;
}

// ------------------------------------------------------------------------------
// Private Logic
// ------------------------------------------------------------------------------

static void _integrate_inplace(S7RTT_REAL* p, S7RTT_REAL* v, S7RTT_REAL* a, S7RTT_REAL dt, S7RTT_REAL j) {
    S7RTT_REAL dt2 = dt * dt;
    S7RTT_REAL dt3 = dt2 * dt;
    *p += (*v) * dt + 0.5 * (*a) * dt2 + S7_ONE_SIXTH * j * dt3;
    *v += (*a) * dt + 0.5 * j * dt2;
    *a += j * dt;
}

static S7RTT_State _integrate_state(S7RTT_State s, S7RTT_REAL dt, S7RTT_REAL j) {
    _integrate_inplace(&s.p, &s.v, &s.a, dt, j);
    s.dt = 0; s.j = 0; // Reset purely for return object cleanliness
    return s;
}

static S7_ProfilePlan _compute_profile_params(S7RTT_REAL v_start, S7RTT_REAL a_start, S7RTT_REAL v_target, S7RTT_REAL a_max, S7RTT_REAL j_max) {
    S7_ProfilePlan plan;
    plan.count = 0;

    S7RTT_REAL acc_clamped = s7_clamp(a_start, -a_max, a_max);

    // Velocity reached if we immediately reduce acc to zero
    S7RTT_REAL t_to_zero = fabs(acc_clamped) / j_max;
    S7RTT_REAL j_to_zero = (acc_clamped > 0) ? -j_max : j_max;
    S7RTT_REAL dv_base = acc_clamped * t_to_zero + 0.5 * j_to_zero * t_to_zero * t_to_zero;
    S7RTT_REAL v_base = v_start + dv_base;

    S7RTT_REAL direction = (v_target < v_base) ? -1.0 : 1.0;
    S7RTT_REAL j_up = (direction > 0) ? j_max : -j_max;
    S7RTT_REAL j_down = (direction > 0) ? -j_max : j_max;

    S7RTT_REAL v_req_total = v_target - v_start;
    S7RTT_REAL a_limit = a_max * direction;

    // Calculate times for Trapezoidal profile
    S7RTT_REAL t1_max = (a_limit - acc_clamped) / j_up;
    S7RTT_REAL t3_max = fabs(a_limit) / j_max;
    if (t1_max < 0) t1_max = 0;

    S7RTT_REAL dv_trapezoid = (acc_clamped * t1_max + 0.5 * j_up * t1_max * t1_max) +
                              (a_limit * t3_max + 0.5 * j_down * t3_max * t3_max);

    bool needs_flat = (direction > 0) ? (v_req_total > dv_trapezoid) : (v_req_total < dv_trapezoid);

    if (needs_flat) {
        S7RTT_REAL v_missing = v_req_total - dv_trapezoid;
        S7RTT_REAL t_flat = v_missing / a_limit;
        if (t_flat < 0) t_flat = 0;

        if (t1_max > S7_EPS_TIME) { plan.segs[plan.count++] = (S7_SegLite){t1_max, j_up}; }
        if (t_flat > S7_EPS_TIME) { plan.segs[plan.count++] = (S7_SegLite){t_flat, 0.0}; }
        if (t3_max > S7_EPS_TIME) { plan.segs[plan.count++] = (S7_SegLite){t3_max, j_down}; }
    } else {
        // Triangular
        S7RTT_REAL term = v_req_total * j_up + 0.5 * acc_clamped * acc_clamped;
        if (term < 0) term = 0.0;
        S7RTT_REAL a_peak_mag = sqrt(term);
        S7RTT_REAL a_peak = (direction > 0) ? a_peak_mag : -a_peak_mag;

        S7RTT_REAL t1 = (a_peak - acc_clamped) / j_up;
        S7RTT_REAL t3 = (0.0 - a_peak) / j_down;
        if (t1 < 0) t1 = 0;

        if (t1 > S7_EPS_TIME) { plan.segs[plan.count++] = (S7_SegLite){t1, j_up}; }
        if (t3 > S7_EPS_TIME) { plan.segs[plan.count++] = (S7_SegLite){t3, j_down}; }
    }
    return plan;
}

static S7RTT_REAL _integrate_plan_dist(S7RTT_REAL* v, S7RTT_REAL* a, const S7_ProfilePlan* plan) {
    S7RTT_REAL dist = 0.0;
    for (int i = 0; i < plan->count; ++i) {
        S7RTT_REAL t = plan->segs[i].dt;
        S7RTT_REAL j = plan->segs[i].j;
        S7RTT_REAL dt2 = t * t;
        dist += (*v) * t + 0.5 * (*a) * dt2 + S7_ONE_SIXTH * j * dt2 * t;
        *v += (*a) * t + 0.5 * j * dt2;
        *a += j * t;
    }
    return dist;
}

static S7RTT_REAL _calc_distance_for_v_peak(const S7RTT_State* root_state, S7RTT_REAL v_peak, S7RTT_REAL target_v, S7RTT_REAL a_max, S7RTT_REAL j_max) {
    S7RTT_REAL v_curr = root_state->v;
    S7RTT_REAL a_curr = root_state->a;

    S7_ProfilePlan plan1 = _compute_profile_params(v_curr, a_curr, v_peak, a_max, j_max);
    S7RTT_REAL d1 = _integrate_plan_dist(&v_curr, &a_curr, &plan1);

    S7_ProfilePlan plan2 = _compute_profile_params(v_curr, a_curr, target_v, a_max, j_max);
    S7RTT_REAL d2 = _integrate_plan_dist(&v_curr, &a_curr, &plan2);

    return d1 + d2;
}

static void _append_to_buf(S7RTT_State* buf, int* count, int max_len, S7RTT_State* curr, const S7_ProfilePlan* plan) {
    for (int i = 0; i < plan->count; ++i) {
        if (*count >= max_len) return;

        S7RTT_State s = *curr;
        s.dt = plan->segs[i].dt;
        s.j = plan->segs[i].j;

        buf[*count] = s;
        (*count)++;

        *curr = _integrate_state(*curr, s.dt, s.j);
    }
}

// Context for solver callback
typedef struct {
    S7RTT_State state;
    S7RTT_REAL target_v;
    S7RTT_REAL a_max;
    S7RTT_REAL j_max;
    S7RTT_REAL dist_req;
} _SolverCtx;

static S7RTT_REAL _cost_brent_cb(S7RTT_REAL v_p, void* data) {
    _SolverCtx* ctx = (_SolverCtx*)data;
    S7RTT_REAL dist = _calc_distance_for_v_peak(&ctx->state, v_p, ctx->target_v, ctx->a_max, ctx->j_max);
    return dist - ctx->dist_req;
}

// ------------------------------------------------------------------------------
// Public API Implementation
// ------------------------------------------------------------------------------

int S7RTT_Plan(const S7RTT_State* start,
               S7RTT_REAL target_p, S7RTT_REAL target_v,
               S7RTT_REAL v_max, S7RTT_REAL a_max, S7RTT_REAL j_max,
               S7RTT_State* out_buf, int max_len)
{
    if (v_max <= 0 || a_max <= 0 || j_max <= 0 || !out_buf || max_len < 2) return -1;

    int count = 0;
    S7RTT_State current_state = *start;

    // 1. Recovery
    if (current_state.a > a_max + S7_EPS_VAL) {
        S7RTT_REAL t_rec = (current_state.a - a_max) / j_max;
        S7RTT_State rec = current_state; rec.dt = t_rec; rec.j = -j_max;
        out_buf[count++] = rec;
        current_state = _integrate_state(current_state, t_rec, -j_max);
        current_state.a = a_max; // Snap
    } else if (current_state.a < -a_max - S7_EPS_VAL) {
        S7RTT_REAL t_rec = (-a_max - current_state.a) / j_max;
        S7RTT_State rec = current_state; rec.dt = t_rec; rec.j = j_max;
        out_buf[count++] = rec;
        current_state = _integrate_state(current_state, t_rec, j_max);
        current_state.a = -a_max; // Snap
    }

    S7RTT_REAL dist_req = target_p - current_state.p;

    // Inertial check
    S7RTT_REAL t_z = fabs(current_state.a) / j_max;
    S7RTT_REAL j_z = (current_state.a > 0) ? -j_max : j_max;
    S7RTT_REAL v_inertial = current_state.v + current_state.a * t_z + 0.5 * j_z * t_z * t_z;

    S7RTT_REAL d_direct = _calc_distance_for_v_peak(&current_state, current_state.v, target_v, a_max, j_max);

    S7_ProfilePlan plans[2];
    int plans_cnt = 0;
    S7RTT_REAL cruise_time = 0.0;

    if (fabs(d_direct - dist_req) <= S7_EPS_DIST) {
        // Case A: Direct
        plans[plans_cnt++] = _compute_profile_params(current_state.v, current_state.a, target_v, a_max, j_max);
    } else {
        S7RTT_REAL d_upper = _calc_distance_for_v_peak(&current_state, v_max, target_v, a_max, j_max);
        S7RTT_REAL d_lower = _calc_distance_for_v_peak(&current_state, -v_max, target_v, a_max, j_max);

        if (dist_req > d_upper + S7_EPS_DIST) {
            // Case B: Cruise +Vmax
            cruise_time = (dist_req - d_upper) / v_max;
            plans[plans_cnt++] = _compute_profile_params(current_state.v, current_state.a, v_max, a_max, j_max);
            plans[plans_cnt++] = _compute_profile_params(v_max, 0.0, target_v, a_max, j_max);
        } else if (dist_req < d_lower - S7_EPS_DIST) {
            // Case C: Cruise -Vmax
            cruise_time = (dist_req - d_lower) / -v_max;
            plans[plans_cnt++] = _compute_profile_params(current_state.v, current_state.a, -v_max, a_max, j_max);
            plans[plans_cnt++] = _compute_profile_params(-v_max, 0.0, target_v, a_max, j_max);
        } else {
            // Case D: Solve
            S7RTT_REAL s_low, s_high;
            S7RTT_REAL lim = s7_fmax(v_max, s7_fmax(fabs(current_state.v), fabs(target_v))) * 1.5 + 10.0;

            if (dist_req > _calc_distance_for_v_peak(&current_state, v_inertial, target_v, a_max, j_max)) {
                s_low = v_inertial; s_high = lim;
            } else {
                s_low = -lim; s_high = v_inertial;
            }

            _SolverCtx ctx = {current_state, target_v, a_max, j_max, dist_req};
            S7RTT_REAL best_v = s7_solve_brent(_cost_brent_cb, s_low, s_high, &ctx);

            // Construct
            S7_ProfilePlan p1 = _compute_profile_params(current_state.v, current_state.a, best_v, a_max, j_max);

            S7RTT_REAL v_temp = current_state.v;
            S7RTT_REAL a_temp = current_state.a;
            _integrate_plan_dist(&v_temp, &a_temp, &p1);

            S7_ProfilePlan p2 = _compute_profile_params(v_temp, a_temp, target_v, a_max, j_max);

            plans[plans_cnt++] = p1;
            plans[plans_cnt++] = p2;
        }
    }

    // Stitch
    if (plans_cnt > 0) {
        _append_to_buf(out_buf, &count, max_len, &current_state, &plans[0]);
    }

    if (cruise_time > S7_EPS_TIME && count < max_len) {
        S7RTT_State cruise = current_state;
        cruise.dt = cruise_time;
        cruise.j = 0.0; cruise.a = 0.0;
        out_buf[count++] = cruise;
        current_state = _integrate_state(current_state, cruise_time, 0.0);
    }

    if (plans_cnt > 1) {
        _append_to_buf(out_buf, &count, max_len, &current_state, &plans[1]);
    }

    return count;
}

int S7RTT_PlanVelocity(const S7RTT_State* start,
                       S7RTT_REAL target_v,
                       S7RTT_REAL v_max, S7RTT_REAL a_max, S7RTT_REAL j_max,
                       S7RTT_State* out_buf, int max_len)
{
    if (v_max <= 0 || a_max <= 0 || j_max <= 0 || !out_buf || max_len < 2) return -1;

    int count = 0;
    S7RTT_State current_state = *start;
    S7RTT_REAL safe_target = s7_clamp(target_v, -v_max, v_max);

    // Recovery logic same as position
    if (current_state.a > a_max + S7_EPS_VAL) {
        S7RTT_REAL t_rec = (current_state.a - a_max) / j_max;
        S7RTT_State rec = current_state; rec.dt = t_rec; rec.j = -j_max;
        out_buf[count++] = rec;
        current_state = _integrate_state(current_state, t_rec, -j_max);
        current_state.a = a_max;
    } else if (current_state.a < -a_max - S7_EPS_VAL) {
        S7RTT_REAL t_rec = (-a_max - current_state.a) / j_max;
        S7RTT_State rec = current_state; rec.dt = t_rec; rec.j = j_max;
        out_buf[count++] = rec;
        current_state = _integrate_state(current_state, t_rec, j_max);
        current_state.a = -a_max;
    }

    S7_ProfilePlan p = _compute_profile_params(current_state.v, current_state.a, safe_target, a_max, j_max);
    _append_to_buf(out_buf, &count, max_len, &current_state, &p);

    return count;
}

void S7RTT_AtTime(const S7RTT_State* trajectory, int count, S7RTT_REAL t_req, S7RTT_State* out_state) {
    if (count <= 0) return;

    S7RTT_REAL t_rem = t_req;

    for (int i = 0; i < count; ++i) {
        const S7RTT_State* node = &trajectory[i];

        if (t_rem <= node->dt + 1e-12) {
            S7RTT_REAL t = (t_rem < 0) ? 0 : t_rem;
            S7RTT_REAL dt2 = t * t;

            out_state->p = node->p + node->v * t + 0.5 * node->a * dt2 + S7_ONE_SIXTH * node->j * dt2 * t;
            out_state->v = node->v + node->a * t + 0.5 * node->j * dt2;
            out_state->a = node->a + node->j * t;
            out_state->j = node->j;
            out_state->dt = 0; // irrelevant in result
            return;
        }
        t_rem -= node->dt;
    }

    // Extrapolate
    const S7RTT_State* last = &trajectory[count-1];
    S7RTT_REAL t = last->dt;
    S7RTT_REAL dt2 = t * t;
    S7RTT_REAL p_end = last->p + last->v * t + 0.5 * last->a * dt2 + S7_ONE_SIXTH * last->j * dt2 * t;
    S7RTT_REAL v_end = last->v + last->a * t + 0.5 * last->j * dt2;

    out_state->p = p_end + v_end * t_rem; // Accessing v with 0 accel
    out_state->v = v_end;
    out_state->a = 0;
    out_state->j = 0;
}

#endif // S7RTT_IMPLEMENTATION

#ifdef __cplusplus
}
#endif

#endif // S7RTT_H
