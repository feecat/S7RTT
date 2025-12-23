// ==============================================================================
// File Name:    S3RTT.h
// Author:       feecat
// Version:      V1.0.1
// Description:  Trapezoidal Velocity Profile Generator (C Implementation)
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
// Description:
//   This library generates motion trajectories using a Trapezoidal Velocity 
//   Profile (T-Curve). It is designed as a lightweight alternative to 
//   S-Curve (S7RTT) algorithms, prioritizing execution speed and low 
//   computational overhead for real-time embedded systems.
//
// Performance Benchmark (Tested on ESP32-S3):
//   Benchmarks reveal a massive reduction in calculation time compared to
//   S-Curve implementations:
//     - S7RTT (S-Curve, double): ~4000 us (Max)
//     - S7RTT (S-Curve, float):  ~800 us (Max)
//     - S3RTT (T-Curve):         ~30 us  (Max)
//
// Application Recommendation:
//   Due to its extremely low latency (~30us), S3RTT is highly recommended 
//   for Servo or FOC (Field Oriented Control) drivers that require high-frequency
//   PID adjustments. It is ideally suited for microcontrollers such as the 
//   ESP32, STM32, or similar platforms where CPU resources are critical.
// ==============================================================================
#ifndef S3RTT_H
#define S3RTT_H

#include <math.h>
#include <string.h> // for memset
#include <stdint.h>

// ==========================================
// Configuration
// ==========================================

#ifndef S3_MAX_NODES
#define S3_MAX_NODES 16
#endif

#ifndef S3_FLOAT_TYPE
#define S3_FLOAT_TYPE float
#endif

#define S3_EPS ((S3_FLOAT_TYPE)1e-5)

typedef S3_FLOAT_TYPE s3_float;

// ==========================================
// Data Structures
// ==========================================

typedef struct {
    s3_float dt;
    s3_float p;
    s3_float v;
    s3_float a;
} S3MotionState;

typedef struct {
    S3MotionState nodes[S3_MAX_NODES];
    int count;
} S3Path;

// ==========================================
// API Prototypes
// ==========================================

static inline void s3_path_init(S3Path* path);

static inline s3_float s3_path_total_time(const S3Path* path);

void s3_plan(S3MotionState start, s3_float target_p, s3_float target_v, 
             s3_float v_max, s3_float a_max, S3Path* path);

void s3_plan_velocity(S3MotionState start, s3_float target_v, 
                      s3_float v_max, s3_float a_max, S3Path* path);

S3MotionState s3_at_time(const S3Path* path, s3_float t);

// ==========================================
// Implementation
// ==========================================

static inline s3_float s3_sign(s3_float x) {
    return (x >= 0.0f) ? 1.0f : -1.0f;
}

static inline s3_float s3_sq(s3_float x) {
    return x * x;
}

static inline void s3_path_init(S3Path* path) {
    path->count = 0;
    // Optional: clear memory for debug safety
    // memset(path->nodes, 0, sizeof(path->nodes)); 
}

static inline void s3_path_push(S3Path* path, s3_float dt, s3_float p, s3_float v, s3_float a) {
    if (path->count < S3_MAX_NODES && dt > S3_EPS) {
        S3MotionState* node = &path->nodes[path->count++];
        node->dt = dt;
        node->p = p;
        node->v = v;
        node->a = a;
    }
}

static inline s3_float s3_path_total_time(const S3Path* path) {
    s3_float total = 0.0f;
    for (int i = 0; i < path->count; i++) {
        total += path->nodes[i].dt;
    }
    return total;
}

void s3_plan(S3MotionState start, s3_float target_p, s3_float target_v, 
             s3_float v_max, s3_float a_max, S3Path* path) 
{
    // Extract parameters
    s3_float vi = start.v;
    s3_float pi = start.p;
    s3_float pf = target_p;
    s3_float vf = target_v;
    s3_float acc = (s3_float)fabs(a_max);
    s3_float v_cap = (s3_float)fabs(v_max);

    // Clamp target velocity
    if (vf > v_cap) vf = v_cap;
    if (vf < -v_cap) vf = -v_cap;

    // Negligible acceleration check
    if (acc < S3_EPS) {
        s3_path_push(path, 0.0f, pi, vi, 0.0f);
        return;
    }

    s3_float dist = pf - pi;
    s3_float s = ((s3_float)fabs(dist) > S3_EPS) ? s3_sign(dist) : 1.0f;

    // ==========================================
    // Phase 1: Kinematic Violation Handling
    // ==========================================
    
    int is_wrong_dir = (vi * s < -S3_EPS);
    int is_overshoot = 0;

    if (!is_wrong_dir && (vi * s > S3_EPS)) {
        // Case A: Target velocity requires moving backwards
        if (vf * s < -S3_EPS) {
            is_overshoot = 1;
        }
        // Case B: Kinetic energy too high to stop
        else if (s3_sq(vi) > s3_sq(vf) + 2.0f * acc * (s3_float)fabs(dist) + S3_EPS) {
            is_overshoot = 1;
        }
    }

    if (is_wrong_dir || is_overshoot) {
        // Strategy: Full stop -> Recurse
        s3_float d_stop = s3_sq(vi) / (2.0f * acc) * s3_sign(vi);
        s3_float stop_pos = pi + d_stop;
        s3_float stop_dt = (s3_float)fabs(vi) / acc;
        
        // Add braking segment
        s3_path_push(path, stop_dt, pi, vi, -s3_sign(vi) * acc);
        
        // Recurse from stop state
        S3MotionState stopped_state = {0, stop_pos, 0.0f, 0.0f};
        s3_plan(stopped_state, pf, vf, v_max, a_max, path);
        return;
    }

    // Check 2: Insufficient Run-up Distance
    if ((s3_float)fabs(vf) > (s3_float)fabs(vi) && (vf * s > 0)) {
        s3_float max_reachable_sq = s3_sq(vi) + 2.0f * acc * (s3_float)fabs(dist);
        if (s3_sq(vf) > max_reachable_sq + S3_EPS) {
            // Strategy: Move backward (Run-up) -> Recurse
            s3_float p_turn = pf - s * (s3_sq(vf) / (2.0f * acc));
            
            // Plan to turnaround point
            s3_plan((S3MotionState){0, pi, vi, 0}, p_turn, 0.0f, v_max, a_max, path);
            
            // Plan from turnaround to target
            // Note: velocity at p_turn is 0
            s3_plan((S3MotionState){0, p_turn, 0, 0}, pf, vf, v_max, a_max, path);
            return;
        }
    }

    // ==========================================
    // Phase 2: Standard Profile Generation
    // ==========================================
    
    // Theoretical peak velocity
    s3_float vp_sq = (2.0f * acc * (s3_float)fabs(dist) + s3_sq(vi) + s3_sq(vf)) / 2.0f;
    s3_float vp = (s3_float)sqrtf((vp_sq > 0.0f) ? vp_sq : 0.0f);
    
    s3_float v_peak = ((vp < v_cap) ? vp : v_cap) * s;

    // 1. Acceleration Phase
    s3_float t1 = (s3_float)fabs(v_peak - vi) / acc;
    s3_float a1 = s3_sign(v_peak - vi) * acc;
    s3_path_push(path, t1, pi, vi, a1);

    s3_float p_curr = pi + (vi + v_peak) * t1 * 0.5f;

    // 2. Cruising Phase
    s3_float t3 = (s3_float)fabs(vf - v_peak) / acc;
    s3_float d3 = (v_peak + vf) * t3 * 0.5f;
    s3_float d_cruise = (pf - p_curr) - d3;

    if (d_cruise * s > S3_EPS) {
        s3_float t2 = (s3_float)fabs(d_cruise) / (s3_float)fabs(v_peak);
        s3_path_push(path, t2, p_curr, v_peak, 0.0f);
        p_curr += d_cruise;
    }

    // 3. Deceleration Phase
    s3_float a3 = s3_sign(vf - v_peak) * acc;
    s3_path_push(path, t3, p_curr, v_peak, a3);
}

void s3_plan_velocity(S3MotionState start, s3_float target_v, 
                      s3_float v_max, s3_float a_max, S3Path* path)
{
    s3_float vi = start.v;
    s3_float pi = start.p;
    s3_float acc = (s3_float)fabs(a_max);
    s3_float v_cap = (s3_float)fabs(v_max);

    s3_float vf = target_v;
    if (vf > v_cap) vf = v_cap;
    if (vf < -v_cap) vf = -v_cap;

    if (acc < S3_EPS || (s3_float)fabs(vf - vi) < S3_EPS) {
        s3_path_push(path, 0.0f, pi, vi, 0.0f);
        return;
    }

    s3_float duration = (s3_float)fabs(vf - vi) / acc;
    s3_float sign = (vf > vi) ? 1.0f : -1.0f;
    
    s3_path_push(path, duration, pi, vi, sign * acc);
}

S3MotionState s3_at_time(const S3Path* path, s3_float t) {
    S3MotionState res = {0};
    
    if (path->count == 0) return res;
    if (t < 0.0f) t = 0.0f;

    s3_float elapsed = 0.0f;
    
    for (int i = 0; i < path->count; i++) {
        const S3MotionState* node = &path->nodes[i];
        if (t < elapsed + node->dt) {
            s3_float dt = t - elapsed;
            res.a = node->a;
            res.v = node->v + node->a * dt;
            res.p = node->p + node->v * dt + 0.5f * node->a * dt * dt;
            return res;
        }
        elapsed += node->dt;
    }

    // Extrapolate using last state
    const S3MotionState* last = &path->nodes[path->count - 1];
    s3_float dt_seg = last->dt;
    s3_float v_end = last->v + last->a * dt_seg;
    s3_float p_end = last->p + last->v * dt_seg + 0.5f * last->a * dt_seg * dt_seg;
    
    s3_float dt_ex = t - elapsed;
    
    res.a = 0.0f;
    res.v = v_end;
    res.p = p_end + v_end * dt_ex;
    
    return res;
}

#endif // S3RTT_H