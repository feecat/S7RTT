# ==============================================================================
# File Name:    S7RTT.py
# Author:       feecat
# Version:      V1.7
# Description:  Simple 7seg Real-Time Trajectory Generator
# Website:      https://github.com/feecat/S7RTT
# License:      Apache License Version 2.0
# ==============================================================================
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

import math

class MotionState:
    """
    Represents a kinematic state at a specific point in time.
    Contains delta time, position, velocity, acceleration, and jerk.
    """
    __slots__ = ('dt', 'p', 'v', 'a', 'j')

    def __init__(self, dt=0.0, p=0.0, v=0.0, a=0.0, j=0.0):
        self.dt = float(dt)
        self.p = float(p)
        self.v = float(v)
        self.a = float(a)
        self.j = float(j)

    def copy(self):
        """Creates a deep copy of the current state."""
        return MotionState(self.dt, self.p, self.v, self.a, self.j)

    def __repr__(self):
        return (f"State(dt={self.dt:.9f}, p={self.p:.4f}, "
                f"v={self.v:.4f}, a={self.a:.4f}, j={self.j:.1f})")


class S7RTT:
    """
    S-Curve Trajectory Planner (Double S-Velocity Profile).
    Calculates time-optimal or continuous trajectories with jerk constraints.
    """
    
    # --- Numerical Tolerances ---
    EPS_TIME = 1e-10     # Minimum significant time duration
    EPS_DIST = 1e-10     # Position tolerance
    EPS_VEL  = 1e-10    # Velocity tolerance
    EPS_ACC  = 1e-10    # Acceleration tolerance
    MATH_EPS = 1e-10    # General epsilon
    EPS_SOLVER = 1e-4   # Solver tolerance
    SOLVER_ITER = 50    # Max iterations for Brent's method
    SOLVER_TOL  = 1e-10  # Solver precision

    def __init__(self):
        pass

    # ==========================================================================
    # 1. Core Discrete Integrator
    # ==========================================================================

    def _integrate_step(self, s, dt, j):
        """
        Integrates the kinematic state forward by 'dt' with constant jerk 'j'.
        Returns a new MotionState object representing the end state.
        """
        if dt <= 0.0: 
            return s.copy()
        
        dt2 = dt * dt
        dt3 = dt2 * dt
        
        p = s.p + s.v * dt + 0.5 * s.a * dt2 + (1.0/6.0) * j * dt3
        v = s.v + s.a * dt + 0.5 * j * dt2
        a = s.a + j * dt
        
        return MotionState(0.0, p, v, a, j)

    def _simulate_profile(self, start_s, shapes):
        """
        Simulates a sequence of (dt, jerk) segments starting from 'start_s'.
        Returns (final_state, list_of_node_states).
        """
        curr = start_s.copy()
        node_states = []
        
        for dt, j in shapes:
            if dt < self.EPS_TIME: 
                continue
            
            seg_start = curr.copy()
            seg_start.dt = dt
            seg_start.j = j
            node_states.append(seg_start)
            
            curr = self._integrate_step(curr, dt, j)
            
        return curr, node_states

    def _integrate_saturated(self, s, t, j_apply, a_max):
        """
        Integrates state 's' for time 't' using jerk 'j_apply', clamping accel at a_max.
        Used primarily for time-optimal calculations.
        """
        limit_a = a_max if j_apply > 0 else -a_max
        dist_to_lim = limit_a - s.a
        
        if abs(j_apply) < self.MATH_EPS:
            t_ramp = float('inf')
        else:
            # Check if we are moving towards the limit
            is_same_dir = (j_apply > 0 and dist_to_lim > -self.MATH_EPS) or \
                          (j_apply < 0 and dist_to_lim < self.MATH_EPS)
            t_ramp = dist_to_lim / j_apply if is_same_dir else 0.0
            
        segments = []
        curr = s.copy()
        
        if t <= t_ramp:
            if t > self.EPS_TIME:
                curr = self._integrate_step(curr, t, j_apply)
                segments.append([t, j_apply])
        else:
            if t_ramp > self.EPS_TIME:
                curr = self._integrate_step(curr, t_ramp, j_apply)
                segments.append([t_ramp, j_apply])
                
            t_hold = t - t_ramp
            if t_hold > self.EPS_TIME:
                curr.a = limit_a
                curr = self._integrate_step(curr, t_hold, 0.0)
                segments.append([t_hold, 0.0])
                
        return curr, segments

    # ==========================================================================
    # 2. Profile Math & Solvers
    # ==========================================================================

    def _calc_vel_change_times(self, v0, a0, v1, a_max, j_max):
        _a0 = max(-a_max, min(a_max, a0))
        
        t_to_zero = abs(_a0) / j_max
        j_restore = -math.copysign(j_max, _a0) if abs(_a0) > self.MATH_EPS else 0.0
        v_min_feasible = v0 + _a0 * t_to_zero + 0.5 * j_restore * t_to_zero**2
        
        direction = 1.0
        if v1 < v_min_feasible - self.MATH_EPS: 
            direction = -1.0
            
        _v0 = v0 * direction
        _a0 = _a0 * direction
        _v1 = v1 * direction
        
        t1, t2, t3 = 0.0, 0.0, 0.0
        
        t1_max = max(0.0, (a_max - _a0) / j_max)
        t3_max = a_max / j_max
        
        dv_inflection = (_a0 * t1_max + 0.5 * j_max * t1_max**2) + \
                        (a_max * t3_max - 0.5 * j_max * t3_max**2)
        dv_req = _v1 - _v0
        
        if dv_req > dv_inflection:
            dv_missing = dv_req - dv_inflection
            t2 = dv_missing / a_max
            t1 = t1_max
            t3 = t3_max
        else:
            term = j_max * dv_req + 0.5 * _a0**2
            if term < 0: term = 0.0
            a_peak = math.sqrt(term)
            t1 = max(0.0, (a_peak - _a0) / j_max)
            t3 = max(0.0, a_peak / j_max)
            
        return t1, t2, t3, direction

    def _build_vel_profile(self, curr, v_target, a_max, j_max):
        t1, t2, t3, direction = self._calc_vel_change_times(curr.v, curr.a, v_target, a_max, j_max)
        nodes = []
        if t1 > self.EPS_TIME: nodes.append([t1, direction * j_max])
        if t2 > self.EPS_TIME: nodes.append([t2, 0.0])
        if t3 > self.EPS_TIME: nodes.append([t3, -direction * j_max])
        return nodes

    def _solve_brent(self, func, lower, upper):
        a, b = lower, upper
        fa, fb = func(a), func(b)
        if abs(fa) < abs(fb): 
            a, b, fa, fb = b, a, fb, fa
        
        c, fc, d, e = a, fa, b - a, b - a
        
        for _ in range(self.SOLVER_ITER):
            if abs(fc) < abs(fb):
                a, b, c = b, c, b
                fa, fb, fc = fb, fc, fb
            
            xm = 0.5 * (c - b)
            
            if abs(xm) < self.SOLVER_TOL or fb == 0: 
                return b
            
            if abs(e) >= self.SOLVER_TOL and abs(fa) > abs(fb):
                s = fb / fa
                if a == c: 
                    p = 2.0 * xm * s; q = 1.0 - s
                else: 
                    q = fa / fc; r = fb / fc
                    p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0))
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0)
                
                if p > 0: q = -q
                p = abs(p)
                
                if 2.0 * p < min(3.0 * xm * q - abs(self.SOLVER_TOL * q), abs(e * q)): 
                    e = d; d = p / q
                else: 
                    d = xm; e = d
            else: 
                d = xm; e = d
            
            a = b; fa = fb
            if abs(d) > self.SOLVER_TOL: 
                b += d
            else: 
                b += math.copysign(self.SOLVER_TOL, xm)
            fb = func(b)
            
            if (fb > 0 and fc > 0) or (fb < 0 and fc < 0): 
                c = a; fc = fa; d = e = b - a
        return b

    def _solve_via_bisection(self, curr, target_p, target_v, v_max, a_max, j_max):
        """
        Finds the optimal intermediate velocity to reach target_p using bisection.
        """
        def get_error(v_mid):
            shapes_1 = self._build_vel_profile(curr, v_mid, a_max, j_max)
            s_mid, _ = self._simulate_profile(curr, shapes_1)
            shapes_2 = self._build_vel_profile(s_mid, target_v, a_max, j_max)
            s_end, _ = self._simulate_profile(s_mid, shapes_2)
            return s_end.p - target_p

        low, high = -v_max, v_max
        if get_error(low) > 0: return -v_max
        if get_error(high) < 0: return v_max
        
        return self._solve_brent(get_error, low, high)

    def _calc_max_reach(self, curr, v_limit, target_v, a_max, j_max):
        shapes_1 = self._build_vel_profile(curr, v_limit, a_max, j_max)
        s_peak, _ = self._simulate_profile(curr, shapes_1)
        shapes_2 = self._build_vel_profile(s_peak, target_v, a_max, j_max)
        s_end, _ = self._simulate_profile(s_peak, shapes_2)
        return s_end.p - curr.p

    # ==========================================================================
    # 3. Trajectory Planning Sub-routines
    # ==========================================================================

    def _handle_safety_decel(self, curr, a_max, j_max):
        """
        Generates a deceleration segment if the initial acceleration exceeds limits.
        Returns (list_of_nodes, updated_current_state).
        """
        nodes = []
        if abs(curr.a) > a_max + self.EPS_ACC:
            j_rec = -math.copysign(j_max, curr.a)
            tgt_a = math.copysign(a_max, curr.a)
            t_rec = (curr.a - tgt_a) / (-j_rec) if abs(j_rec) > 0 else 0
            
            if t_rec > self.EPS_TIME:
                n = curr.copy(); n.dt = t_rec; n.j = j_rec
                nodes.append(n)
                curr = self._integrate_step(curr, t_rec, j_rec)
                curr.a = tgt_a # Force precise value
        
        return nodes, curr

    def _try_time_optimal_plan(self, curr, target_p, target_v, a_max, j_max, v_max):
        """
        Attempts to find a time-optimal trajectory using a bidirectional competition strategy.
        Solves for the optimal switching time 't' using Brent's method.
        
        Returns:
            (nodes, final_state) if successful, else (None, None).
        """
        
        # 1. Estimate search horizon
        t_est = (abs(curr.v) + v_max) / a_max if a_max > 0 else 1.0
        search_horizon = t_est * 2.0 + 5.0
        
        # 2. Internal Solver Function
        def solve_for_jerk(j_apply):
            
            # Calculates position error after integrating for time 't' and planning the rest.
            def get_pos_error(t):
                t = max(0.0, t) # Clamp time
                s1, _ = self._integrate_saturated(curr, t, j_apply, a_max)
                shapes_rem = self._build_vel_profile(s1, target_v, a_max, j_max)
                s_final, _ = self._simulate_profile(s1, shapes_rem)
                return s_final.p - target_p

            # --- Boundary & Feasibility Checks ---
            err_0 = get_pos_error(0.0)
            
            # Check 1: Zero Drift Handling
            if abs(err_0) < self.EPS_SOLVER:
                best_t = 0.0
            # Check 2: Unreachable Target (Fail Fast)
            elif err_0 * get_pos_error(search_horizon) > 0:
                return None
            else:
                # Check 3: Solve Root (Brent's Method)
                best_t = self._solve_brent(get_pos_error, 0.0, search_horizon)
                
                # Verify precision strictness
                if abs(get_pos_error(best_t)) > self.EPS_SOLVER:
                    return None
            
            # --- Trajectory Reconstruction ---
            # Re-integrate to generate the actual node list.
            nodes = []
            _, switch_segments = self._integrate_saturated(curr, best_t, j_apply, a_max)
            
            # Use a temporary state variable to prevent side-effects on 'curr'
            running_state = curr.copy()
            
            # Phase 1: Variable acceleration (Switching phase)
            for dt, j in switch_segments:
                if dt < self.EPS_TIME: continue
                
                # Snapshot state at the start of the segment
                n = running_state.copy()
                n.dt, n.j = dt, j
                nodes.append(n)
                
                # Advance state
                running_state = self._integrate_step(running_state, dt, j)
            
            # Phase 2: Velocity profile to target (Remaining phase)
            shapes_rem = self._build_vel_profile(running_state, target_v, a_max, j_max)
            final_state, rem_nodes = self._simulate_profile(running_state, shapes_rem)
            nodes.extend(rem_nodes)
            
            # Calculate total duration for competition comparison
            total_duration = sum(n.dt for n in nodes)
            return (total_duration, nodes, final_state)

        # 3. Bidirectional Competition
        #    Try both Positive Jerk (+J) and Negative Jerk (-J).
        #    Collect valid results into a candidate list.
        candidates = []
        for j in [j_max, -j_max]:
            res = solve_for_jerk(j)
            if res: candidates.append(res)
        
        # 4. Selection
        #    If no valid trajectory found, return None (triggers fallback).
        if not candidates:
            return None, None
            
        #    Select the trajectory with the minimum total duration (Time Optimal).
        best_res = min(candidates, key=lambda x: x[0])
        
        return best_res[1], best_res[2]

    def _plan_fallback_cruise(self, curr, target_p, target_v, v_max, a_max, j_max):
        """
        Fallback strategy:
        1. Solve for ideal peak velocity using bisection.
        2. Cruise if there is a spatial gap.
        3. Decelerate to target.
        """
        nodes = []
        
        # 1. Reach optimal intermediate velocity
        best_v = self._solve_via_bisection(curr, target_p, target_v, v_max, a_max, j_max)
        shapes_a = self._build_vel_profile(curr, best_v, a_max, j_max)
        curr, nodes_a = self._simulate_profile(curr, shapes_a)
        nodes.extend(nodes_a)
        
        # Force acceleration to zero for pure cruise calculation
        curr.a = 0.0
        
        # 2. Calculate gap and add cruise segment if needed
        # Simulate deceleration to see where we land
        s_dec_sim, _ = self._simulate_profile(curr, self._build_vel_profile(curr, target_v, a_max, j_max))
        dist_gap = target_p - s_dec_sim.p
        
        effective_v = curr.v
        if abs(effective_v) < self.MATH_EPS:
            effective_v = self.MATH_EPS * math.copysign(1, dist_gap)
        
        if abs(dist_gap) > self.EPS_DIST: 
             cruise_time = dist_gap / effective_v
             if cruise_time > self.EPS_TIME:
                 n = curr.copy(); n.dt = cruise_time; n.j = 0.0
                 nodes.append(n)
                 curr = self._integrate_step(curr, cruise_time, 0.0)

        # 3. Decelerate to final velocity
        shapes_b = self._build_vel_profile(curr, target_v, a_max, j_max)
        curr, nodes_b = self._simulate_profile(curr, shapes_b)
        nodes.extend(nodes_b)
        
        return nodes

    def _refine_trajectory_precision(self, nodes, start_state, target_p):
        """
        Simulates the generated nodes to verify the final position.
        If a position error exists (> 1e-8) and a cruise segment is available,
        adjusts the cruise duration and propagates the correction to subsequent nodes.
        
        Optimized to avoid repetitive manual integration loops.
        """
        if not nodes:
            return

        # --- Phase 1: Simulation and Identification ---
        sim_s = start_state.copy()
        correction_idx = -1
        max_cruise_dt = -1.0
        
        # We simulate strictly using the integrator to capture accumulated floating point errors
        for i, node in enumerate(nodes):
            # Check if this node is a valid cruise candidate (a=0, j=0)
            # Note: node.a/j describes the state at the START of the segment
            if abs(node.j) < self.MATH_EPS and abs(node.a) < self.EPS_ACC:
                if node.dt > max_cruise_dt:
                    max_cruise_dt = node.dt
                    correction_idx = i
            
            # Advance simulation state
            sim_s = self._integrate_step(sim_s, node.dt, node.j)

        pos_error = target_p - sim_s.p

        # --- Phase 2: Correction and Propagation ---
        if abs(pos_error) > self.EPS_DIST and correction_idx != -1:
            v_cruise = nodes[correction_idx].v
            
            if abs(v_cruise) > self.EPS_VEL:
                dt_fix = pos_error / v_cruise
                new_dt = nodes[correction_idx].dt + dt_fix
                
                if new_dt < self.EPS_TIME:
                    new_dt = self.EPS_TIME
                
                # Apply time correction
                nodes[correction_idx].dt = new_dt
                
                # Propagate changes:
                # We must re-calculate the start state (p, v, a) for all nodes 
                # AFTER the modified one to ensure continuity.
                
                # 1. Get the state *after* the modified cruise segment
                #    We need the state at the start of the cruise segment to integrate forward
                curr_s = start_state.copy()
                for k in range(correction_idx):
                    curr_s = self._integrate_step(curr_s, nodes[k].dt, nodes[k].j)
                
                # Now integrate the modified cruise segment
                curr_s = self._integrate_step(curr_s, nodes[correction_idx].dt, nodes[correction_idx].j)
                
                # 2. Update subsequent nodes
                for k in range(correction_idx + 1, len(nodes)):
                    # Update the node's start state properties
                    nodes[k].p = curr_s.p
                    nodes[k].v = curr_s.v
                    nodes[k].a = curr_s.a
                    
                    # Integrate to find the start state for the next node
                    curr_s = self._integrate_step(curr_s, nodes[k].dt, nodes[k].j)

    # ==========================================================================
    # 4. Main Entry Point
    # ==========================================================================

    def plan(self, start_state, target_p, target_v, v_max, a_max, j_max):
        if v_max <= 0 or a_max <= 0 or j_max <= 0: return []
        
        final_nodes = []
        curr = start_state.copy()
        
        # 1. Safety Deceleration
        #    Handle cases where initial acceleration is already violating limits
        decel_nodes, curr = self._handle_safety_decel(curr, a_max, j_max)
        final_nodes.extend(decel_nodes)

        # 2. Capacity Check
        #    Determine if we can reach the target using max profiles or if we need optimal time solving
        dist_req = target_p - curr.p
        d_pos_limit = self._calc_max_reach(curr, v_max, target_v, a_max, j_max)
        d_neg_limit = self._calc_max_reach(curr, -v_max, target_v, a_max, j_max)
        
        use_optimal_solver = True
        if dist_req > d_pos_limit + self.EPS_DIST: use_optimal_solver = False
        if dist_req < d_neg_limit - self.EPS_DIST: use_optimal_solver = False
        
        # 3. Execution Strategy
        opt_nodes = None
        if use_optimal_solver:
            opt_nodes, _ = self._try_time_optimal_plan(curr, target_p, target_v, a_max, j_max, v_max)
        
        if opt_nodes is not None:
            final_nodes.extend(opt_nodes)
        else:
            # Fallback: Bisection solve for velocity + Cruise segment
            fallback_nodes = self._plan_fallback_cruise(curr, target_p, target_v, v_max, a_max, j_max)
            final_nodes.extend(fallback_nodes)

        # 4. Precision Refinement
        #    Fix small floating point errors by adjusting cruise segments if available
        self._refine_trajectory_precision(final_nodes, start_state, target_p)

        return final_nodes
    
    def plan_velocity(self, start_state, target_v, v_max, a_max, j_max):
        """
        Calculates a time-optimal velocity profile to reach target_v from start_state.
        Does not consider position constraints.
        """
        # 0. Basic parameter validation
        if v_max <= 0 or a_max <= 0 or j_max <= 0: 
            return []

        nodes = []
        curr = start_state.copy()

        # 1. Safety Deceleration
        # If initial acceleration is exceeding a_max, first ramp it down safely.
        safety_nodes, curr = self._handle_safety_decel(curr, a_max, j_max)
        nodes.extend(safety_nodes)

        # 2. Clamp Target Velocity
        # Ensure the target does not exceed the physical maximum velocity
        safe_target_v = max(-v_max, min(v_max, target_v))

        # 3. Build Velocity Profile
        # Uses internal math to find t1 (jerk), t2 (const accel), t3 (jerk)
        # to reach safe_target_v with final acceleration = 0.
        shapes = self._build_vel_profile(curr, safe_target_v, a_max, j_max)

        # 4. Generate Trajectory Nodes
        # Simulates the shapes to create the MotionState objects
        _, profile_nodes = self._simulate_profile(curr, shapes)
        nodes.extend(profile_nodes)

        return nodes

    def at_time(self, trajectory, t):
        if not trajectory: return MotionState()
        elapsed = 0.0
        for node in trajectory:
            if t <= elapsed + node.dt: 
                return self._integrate_step(node, t - elapsed, node.j)
            elapsed += node.dt
        last = trajectory[-1]
        final_state = self._integrate_step(last, last.dt, last.j)
        return self._integrate_step(final_state, t - elapsed, 0.0)
