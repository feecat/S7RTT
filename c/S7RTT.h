// ==============================================================================
// File Name:    S7RTT.h
// Author:       feecat
// Version:      V1.5.2
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
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdbool.h>

/* --- Data Structures --- */

typedef struct {
    double dt; /* Duration of this segment */
    double p;  /* Position */
    double v;  /* Velocity */
    double a;  /* Acceleration */
    double j;  /* Jerk */
} S7RTT_MotionState;

/* Dynamic array wrapper for the resulting trajectory */
typedef struct {
    S7RTT_MotionState* nodes;
    int count;
    int capacity;
} S7RTT_Path;

/* --- API Functions --- */

/* Main Planning Function
 * Returns a path object. User must call s7rtt_path_free() on the result.
 */
S7RTT_Path s7rtt_plan(S7RTT_MotionState start, double target_p, double target_v,
                      double v_max, double a_max, double j_max);

/* Helper to sample the trajectory at a specific time t */
S7RTT_MotionState s7rtt_at_time(const S7RTT_Path* path, double t);

/* Helper to free path memory */
void s7rtt_path_free(S7RTT_Path* path);

/* Helper: Initialize a state */
static inline S7RTT_MotionState s7rtt_state_init(double dt, double p, double v, double a, double j) {
    S7RTT_MotionState s;
    s.dt = dt; s.p = p; s.v = v; s.a = a; s.j = j;
    return s;
}

#ifdef __cplusplus
}
#endif

#endif /* S7RTT_H */

/* ========================================================================== */
/* IMPLEMENTATION                                                             */
/* ========================================================================== */

#ifdef S7RTT_IMPLEMENTATION

/* --- Constants --- */
#define S7_EPS_TIME   1e-9
#define S7_EPS_DIST   1e-8
#define S7_EPS_VEL    1e-7
#define S7_EPS_ACC    1e-6
#define S7_EPS_SOLVER 1e-3
#define S7_MATH_EPS   1e-9
#define S7_SOLVER_TOL 1e-8
#define S7_SOLVER_ITER 30

#define S7_ONE_SIXTH  (1.0 / 6.0)
#define S7_ONE_HALF   0.5

/* --- Internal Helpers & Math --- */

static inline double s7_fclamp(double x, double lower, double upper) {
    return fmax(lower, fmin(x, upper));
}

static inline double s7_fsign(double x) {
    if (x > 0.0) return 1.0;
    if (x < 0.0) return -1.0;
    return 0.0;
}

static inline double s7_copysign(double x, double y) {
    return copysign(x, y);
}

/* --- Path (Vector) Management --- */

static void s7_path_init(S7RTT_Path* path, int initial_cap) {
    path->count = 0;
    path->capacity = initial_cap;
    path->nodes = (S7RTT_MotionState*)malloc(sizeof(S7RTT_MotionState) * initial_cap);
    if (path->nodes == NULL) {
        path->capacity = 0;
    }
}

static void s7_path_push(S7RTT_Path* path, S7RTT_MotionState s) {
    if (path->nodes == NULL) return;
    if (path->count >= path->capacity) return;
    path->nodes[path->count++] = s;
}

void s7rtt_path_free(S7RTT_Path* path) {
    if (path->nodes) {
        free(path->nodes);
        path->nodes = NULL;
    }
    path->count = 0;
    path->capacity = 0;
}

/* --- Core Integrator --- */

static inline void s7_integrate_state_inplace(S7RTT_MotionState* s, double dt, double j) {
    double dt2 = dt * dt;
    double dt3 = dt2 * dt;

    s->p += s->v * dt + s->a * dt2 * S7_ONE_HALF + j * dt3 * S7_ONE_SIXTH;
    s->v += s->a * dt + j * dt2 * S7_ONE_HALF;
    s->a += j * dt;
}

static inline S7RTT_MotionState s7_integrate_step(S7RTT_MotionState s, double dt, double j) {
    if (dt <= S7_EPS_TIME) return s;
    S7RTT_MotionState next = s;
    s7_integrate_state_inplace(&next, dt, j);
    next.dt = 0.0;
    next.j = j;
    return next;
}

/* --- TinyProfile Replacement (Stack Array) --- */
#define S7_MAX_SEGS 8
typedef struct {
    struct { double dt; double j; } segs[S7_MAX_SEGS];
    int count;
} S7_Profile;

static inline void s7_profile_push(S7_Profile* p, double dt, double j) {
    if (p->count < S7_MAX_SEGS) {
        p->segs[p->count].dt = dt;
        p->segs[p->count].j = j;
        p->count++;
    }
}

/* --- Simulation Helpers --- */

/* Append profile to path and return final state */
static S7RTT_MotionState s7_append_from_profile(S7RTT_Path* nodes, S7RTT_MotionState start_s, const S7_Profile* shapes) {
    S7RTT_MotionState curr = start_s;
    for (int i = 0; i < shapes->count; ++i) {
        double dt = shapes->segs[i].dt;
        double j = shapes->segs[i].j;

        if (dt < S7_EPS_TIME) continue;

        S7RTT_MotionState seg_start = curr;
        seg_start.dt = dt;
        seg_start.j = j;
        s7_path_push(nodes, seg_start);

        s7_integrate_state_inplace(&curr, dt, j);
    }
    return curr;
}

/* Fast endpoint simulation */
static inline void s7_simulate_endpoint_inplace(S7RTT_MotionState* curr, const S7_Profile* shapes) {
    for (int i = 0; i < shapes->count; ++i) {
        if (shapes->segs[i].dt >= S7_EPS_TIME) {
            s7_integrate_state_inplace(curr, shapes->segs[i].dt, shapes->segs[i].j);
        }
    }
}

static inline void s7_integrate_saturated_state_only(S7RTT_MotionState* curr, double t, double j_apply, double a_max) {
    if (t <= S7_EPS_TIME) return;

    double limit_a = (j_apply > 0) ? a_max : -a_max;
    double dist_to_lim = limit_a - curr->a;

    if (fabs(j_apply) < S7_MATH_EPS) {
        s7_integrate_state_inplace(curr, t, j_apply);
        return;
    }

    bool same_dir = (j_apply > 0) ? (dist_to_lim > -S7_MATH_EPS) : (dist_to_lim < S7_MATH_EPS);
    double t_ramp = same_dir ? (dist_to_lim / j_apply) : 0.0;

    if (t <= t_ramp) {
        s7_integrate_state_inplace(curr, t, j_apply);
    } else {
        if (t_ramp > S7_EPS_TIME) {
            s7_integrate_state_inplace(curr, t_ramp, j_apply);
        }
        curr->a = limit_a;
        double t_hold = t - t_ramp;
        if (t_hold > S7_EPS_TIME) {
            s7_integrate_state_inplace(curr, t_hold, 0.0);
        }
    }
}

static S7RTT_MotionState s7_append_saturated_profile(S7RTT_Path* nodes, S7RTT_MotionState s, double t, double j_apply, double a_max) {
    S7RTT_MotionState curr = s;
    if (t <= S7_EPS_TIME) return curr;

    double limit_a = (j_apply > 0) ? a_max : -a_max;
    double dist_to_lim = limit_a - s.a;

    double t_ramp = 0.0;
    if (fabs(j_apply) < S7_MATH_EPS) {
        t_ramp = INFINITY;
    } else {
        bool same_dir = (j_apply > 0) ? (dist_to_lim > -S7_MATH_EPS) : (dist_to_lim < S7_MATH_EPS);
        t_ramp = same_dir ? (dist_to_lim / j_apply) : 0.0;
    }

    if (t <= t_ramp) {
        S7RTT_MotionState node = curr;
        node.dt = t; node.j = j_apply;
        s7_path_push(nodes, node);
        s7_integrate_state_inplace(&curr, t, j_apply);
    } else {
        if (t_ramp > S7_EPS_TIME) {
            S7RTT_MotionState node = curr;
            node.dt = t_ramp; node.j = j_apply;
            s7_path_push(nodes, node);
            s7_integrate_state_inplace(&curr, t_ramp, j_apply);
        }
        curr.a = limit_a;
        double t_hold = t - t_ramp;
        if (t_hold > S7_EPS_TIME) {
            S7RTT_MotionState node = curr;
            node.dt = t_hold; node.j = 0.0;
            s7_path_push(nodes, node);
            s7_integrate_state_inplace(&curr, t_hold, 0.0);
        }
    }
    return curr;
}

/* --- Solvers --- */

typedef struct {
    double t1, t2, t3, dir;
} S7_VelChange;

static inline S7_VelChange s7_calc_vel_change_times(double v0, double a0, double v1, double a_max, double j_max) {
    double _a0 = s7_fclamp(a0, -a_max, a_max);
    double abs_a0 = fabs(_a0);
    double t_to_zero = abs_a0 / j_max;
    double j_restore = (abs_a0 > S7_MATH_EPS) ? -s7_copysign(j_max, _a0) : 0.0;
    double v_min_feasible = v0 + _a0 * t_to_zero + 0.5 * j_restore * t_to_zero * t_to_zero;

    double direction = 1.0;
    if (v1 < v_min_feasible - S7_MATH_EPS) {
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
        double a_peak = (term > 0) ? sqrt(term) : 0.0;
        t1 = (a_peak - _a0_scaled) / j_max;
        t3 = a_peak / j_max;
        if (t1 < 0) t1 = 0.0;
        if (t3 < 0) t3 = 0.0;
    }
    S7_VelChange res = {t1, t2, t3, direction};
    return res;
}

static inline void s7_build_vel_profile(S7_Profile* nodes, S7RTT_MotionState curr, double v_target, double a_max, double j_max) {
    nodes->count = 0;
    S7_VelChange res = s7_calc_vel_change_times(curr.v, curr.a, v_target, a_max, j_max);
    if (res.t1 > S7_EPS_TIME) s7_profile_push(nodes, res.t1, res.dir * j_max);
    if (res.t2 > S7_EPS_TIME) s7_profile_push(nodes, res.t2, 0.0);
    if (res.t3 > S7_EPS_TIME) s7_profile_push(nodes, res.t3, -res.dir * j_max);
}

/* Callback definition for Brent's method */
typedef double (*S7_SolverFunc)(double t, void* ctx);

static double s7_solve_brent(S7_SolverFunc func, void* ctx, double lower, double upper) {
    double a = lower, b = upper;
    double fa = func(a, ctx), fb = func(b, ctx);

    if (fabs(fa) < fabs(fb)) {
        double tmp = a; a = b; b = tmp;
        tmp = fa; fa = fb; fb = tmp;
    }

    double c = a, fc = fa;
    double d = b - a, e = b - a;

    for (int i = 0; i < S7_SOLVER_ITER; ++i) {
        if (fabs(fb) < S7_SOLVER_TOL) return b;
        if (fabs(fc) < fabs(fb)) {
            a = b; b = c; c = a;
            fa = fb; fb = fc; fc = fa;
        }
        double xm = 0.5 * (c - b);
        if (fabs(xm) < S7_SOLVER_TOL) return b;

        if (fabs(e) >= S7_SOLVER_TOL && fabs(fa) > fabs(fb)) {
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
            p = fabs(p);
            double min_term = fmin(3.0 * xm * q - fabs(S7_SOLVER_TOL * q), fabs(e * q));
            if (2.0 * p < min_term) { e = d; d = p / q; }
            else { d = xm; e = d; }
        } else {
            d = xm; e = d;
        }
        a = b; fa = fb;
        if (fabs(d) > S7_SOLVER_TOL) b += d;
        else b += (xm > 0 ? S7_SOLVER_TOL : -S7_SOLVER_TOL);
        fb = func(b, ctx);
        if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
            c = a; fc = fa; d = e = b - a;
        }
    }
    return b;
}

/* Contexts for solvers */
typedef struct {
    S7RTT_MotionState curr;
    double target_p, target_v, a_max, j_max;
} S7_BisectCtx;

typedef struct {
    S7RTT_MotionState curr;
    double target_p, target_v, a_max, j_max, j_action;
} S7_TimeCtx;

static double s7_bisect_cost(double v_mid, void* data) {
    S7_BisectCtx* c = (S7_BisectCtx*)data;
    S7_Profile p1, p2;
    S7RTT_MotionState s = c->curr;

    s7_build_vel_profile(&p1, s, v_mid, c->a_max, c->j_max);
    s7_simulate_endpoint_inplace(&s, &p1);

    s7_build_vel_profile(&p2, s, c->target_v, c->a_max, c->j_max);
    s7_simulate_endpoint_inplace(&s, &p2);

    return s.p - c->target_p;
}

static double s7_solve_via_bisection(S7RTT_MotionState curr, double target_p, double target_v,
                                     double v_max, double a_max, double j_max) {
    S7_BisectCtx ctx = {curr, target_p, target_v, a_max, j_max};

    double low = -v_max;
    double high = v_max;

    if (s7_bisect_cost(low, &ctx) > 0) return -v_max;
    if (s7_bisect_cost(high, &ctx) < 0) return v_max;

    return s7_solve_brent(s7_bisect_cost, &ctx, low, high);
}

static double s7_calc_max_reach(S7RTT_MotionState curr, double v_limit, double target_v, double a_max, double j_max) {
    S7_Profile p;
    S7RTT_MotionState s = curr;
    s7_build_vel_profile(&p, s, v_limit, a_max, j_max);
    s7_simulate_endpoint_inplace(&s, &p);
    s7_build_vel_profile(&p, s, target_v, a_max, j_max);
    s7_simulate_endpoint_inplace(&s, &p);
    return s.p - curr.p;
}

static double s7_time_opt_cost(double t, void* data) {
    if (t < 0) t = 0;
    S7_TimeCtx* c = (S7_TimeCtx*)data;
    S7RTT_MotionState s = c->curr;
    S7_Profile rem;

    s7_integrate_saturated_state_only(&s, t, c->j_action, c->a_max);
    s7_build_vel_profile(&rem, s, c->target_v, c->a_max, c->j_max);
    s7_simulate_endpoint_inplace(&s, &rem);

    return s.p - c->target_p;
}

static double s7_solve_time_optimal(S7RTT_MotionState curr, double target_p, double target_v,
                                    double a_max, double j_max, double v_max, double j_action) {
    S7_TimeCtx ctx = {curr, target_p, target_v, a_max, j_max, j_action};

    double t_est = (a_max > 0) ? (fabs(curr.v) + v_max) / a_max : 1.0;
    double t_search_max = t_est * 3.0 + 5.0;

    double e0 = s7_time_opt_cost(0.0, &ctx);
    double e_max = s7_time_opt_cost(t_search_max, &ctx);

    if (e0 * e_max > 0) return -1.0;

    double best_t = s7_solve_brent(s7_time_opt_cost, &ctx, 0.0, t_search_max);

    if (fabs(s7_time_opt_cost(best_t, &ctx)) > S7_EPS_SOLVER) return -1.0;
    return best_t;
}

/* --- Planning Steps --- */

static S7RTT_MotionState s7_append_safety_decel(S7RTT_Path* nodes, S7RTT_MotionState curr, double a_max, double j_max) {
    if (fabs(curr.a) > a_max + S7_EPS_ACC) {
        double j_rec = -s7_copysign(j_max, curr.a);
        double tgt_a = s7_copysign(a_max, curr.a);
        double t_rec = (curr.a - tgt_a) / (-j_rec);

        if (t_rec > S7_EPS_TIME) {
            S7RTT_MotionState n = curr;
            n.dt = t_rec; n.j = j_rec;
            s7_path_push(nodes, n);
            s7_integrate_state_inplace(&curr, t_rec, j_rec);
            curr.a = tgt_a;
        }
    }
    return curr;
}

static void s7_append_fallback_cruise(S7RTT_Path* nodes, S7RTT_MotionState curr, double target_p, double target_v,
                                      double v_max, double a_max, double j_max) {
    double best_v = s7_solve_via_bisection(curr, target_p, target_v, v_max, a_max, j_max);

    S7_Profile shapes;
    s7_build_vel_profile(&shapes, curr, best_v, a_max, j_max);
    curr = s7_append_from_profile(nodes, curr, &shapes);

    curr.a = 0.0;

    s7_build_vel_profile(&shapes, curr, target_v, a_max, j_max);
    S7RTT_MotionState s_dec = curr;
    s7_simulate_endpoint_inplace(&s_dec, &shapes);

    double dist_gap = target_p - s_dec.p;
    double effective_v = (fabs(curr.v) < S7_MATH_EPS) ? s7_copysign(S7_MATH_EPS, dist_gap) : curr.v;

    if (fabs(dist_gap) > S7_EPS_DIST) {
        double cruise_time = dist_gap / effective_v;
        if (cruise_time > S7_EPS_TIME) {
            S7RTT_MotionState n = curr;
            n.dt = cruise_time; n.j = 0.0;
            s7_path_push(nodes, n);
            s7_integrate_state_inplace(&curr, cruise_time, 0.0);
        }
    }

    s7_build_vel_profile(&shapes, curr, target_v, a_max, j_max);
    s7_append_from_profile(nodes, curr, &shapes);
}

static void s7_refine_trajectory(S7RTT_Path* nodes, S7RTT_MotionState start, double target_p) {
    if (nodes->count == 0) return;

    S7RTT_MotionState sim_s = start;
    int correction_idx = -1;
    double max_cruise_dt = -1.0;

    for (int i = 0; i < nodes->count; ++i) {
        if (fabs(nodes->nodes[i].j) < S7_MATH_EPS && fabs(nodes->nodes[i].a) < S7_EPS_ACC) {
            if (nodes->nodes[i].dt > max_cruise_dt) {
                max_cruise_dt = nodes->nodes[i].dt;
                correction_idx = i;
            }
        }
        s7_integrate_state_inplace(&sim_s, nodes->nodes[i].dt, nodes->nodes[i].j);
    }

    double pos_error = target_p - sim_s.p;
    if (fabs(pos_error) <= S7_EPS_DIST || correction_idx == -1) return;

    double v_cruise = nodes->nodes[correction_idx].v;
    if (fabs(v_cruise) > S7_EPS_VEL) {
        double dt_fix = pos_error / v_cruise;
        double new_dt = nodes->nodes[correction_idx].dt + dt_fix;
        if (new_dt < S7_EPS_TIME) new_dt = S7_EPS_TIME;

        nodes->nodes[correction_idx].dt = new_dt;

        S7RTT_MotionState curr = start;
        for (int k = 0; k < nodes->count; ++k) {
            nodes->nodes[k].p = curr.p;
            nodes->nodes[k].v = curr.v;
            nodes->nodes[k].a = curr.a;
            s7_integrate_state_inplace(&curr, nodes->nodes[k].dt, nodes->nodes[k].j);
        }
    }
}

/* --- Main Plan Implementation --- */

S7RTT_Path s7rtt_plan(S7RTT_MotionState start_state, double target_p, double target_v,
                      double v_max, double a_max, double j_max) {
    S7RTT_Path path;
    s7_path_init(&path, 32);

    if (v_max <= 0 || a_max <= 0 || j_max <= 0) return path;

    S7RTT_MotionState curr = start_state;

    /* 1. Safety Decel */
    curr = s7_append_safety_decel(&path, curr, a_max, j_max);

    /* 2. Feasibility */
    double dist_req = target_p - curr.p;
    double d_pos_lim = s7_calc_max_reach(curr, v_max, target_v, a_max, j_max);
    double d_neg_lim = s7_calc_max_reach(curr, -v_max, target_v, a_max, j_max);

    bool use_optimal = true;
    if (dist_req > d_pos_lim + S7_EPS_DIST || dist_req < d_neg_lim - S7_EPS_DIST) {
        use_optimal = false;
    }

    bool opt_success = false;
    if (use_optimal) {
        S7_Profile stop_shapes;
        s7_build_vel_profile(&stop_shapes, curr, target_v, a_max, j_max);
        S7RTT_MotionState s_stop = curr;
        s7_simulate_endpoint_inplace(&s_stop, &stop_shapes);

        bool is_pos_move = (fabs(curr.v) > S7_EPS_VEL) ? (curr.v > 0) : (target_p > curr.p);
        double gap = target_p - s_stop.p;
        double j_action = (is_pos_move) ? ((gap < -S7_EPS_DIST) ? -j_max : j_max)
                                        : ((gap > S7_EPS_DIST) ? j_max : -j_max);

        double best_t = s7_solve_time_optimal(curr, target_p, target_v, a_max, j_max, v_max, j_action);

        if (best_t >= 0.0) {
            curr = s7_append_saturated_profile(&path, curr, best_t, j_action, a_max);
            S7_Profile rem;
            s7_build_vel_profile(&rem, curr, target_v, a_max, j_max);
            s7_append_from_profile(&path, curr, &rem);
            opt_success = true;
        }
    }

    if (!opt_success) {
        s7_append_fallback_cruise(&path, curr, target_p, target_v, v_max, a_max, j_max);
    }

    s7_refine_trajectory(&path, start_state, target_p);

    return path;
}

S7RTT_MotionState s7rtt_at_time(const S7RTT_Path* path, double t) {
    if (!path || path->count == 0) {
        S7RTT_MotionState s = {0};
        return s;
    }

    double elapsed = 0.0;
    for (int i = 0; i < path->count; ++i) {
        double dt = path->nodes[i].dt;
        if (t <= elapsed + dt) {
            return s7_integrate_step(path->nodes[i], t - elapsed, path->nodes[i].j);
        }
        elapsed += dt;
    }

    /* Extrapolate last */
    S7RTT_MotionState last = path->nodes[path->count - 1];
    S7RTT_MotionState final_s = s7_integrate_step(last, last.dt, last.j);
    return s7_integrate_step(final_s, t - elapsed, 0.0);
}

#endif /* S7RTT_IMPLEMENTATION */
