# ==============================================================================
# File Name:    S7RTT.py
# Author:       feecat
# Version:      V1.1
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

class Solver:
    """
    A numerical solver utility using Brent's Method.
    Primarily used to find the exact peak velocity required to align the
    trajectory's total distance with the target position.
    """
    
    EPS_TOL = 1e-8       # Tolerance for function value convergence
    EPS_CONVERGE = 1e-12 # Tolerance for interval convergence
    MAX_ITER = 50        # Maximum iterations to prevent infinite loops

    @staticmethod
    def solve_monotonic_brent(func, low, high):
        """
        Finds the root of a monotonic function 'func' within the interval [low, high].
        """
        f_low = func(low)
        f_high = func(high)
        
        # If the root is bracketed, return the side closer to zero
        if f_low * f_high > 0:
            if abs(f_low) < abs(f_high): return low
            else: return high
            
        a = low; b = high
        fa = f_low; fb = f_high
        
        root_est = b
        
        for _ in range(Solver.MAX_ITER):
            # Ensure 'b' is the better estimate (closer to 0)
            if abs(fa) < abs(fb):
                a, b = b, a
                fa, fb = fb, fa
            
            # Check for convergence
            if abs(fb - fa) < Solver.EPS_CONVERGE:
                root_est = b
                break
            else:
                # Try Inverse Quadratic Interpolation
                if abs(fa - fb) > 1e-14:
                    root_est = b - fb * (b - a) / (fb - fa)
                else:
                    root_est = b
            
            mid_val = 0.5 * (a + b)
            
            # Conditions to fall back to Bisection method
            cond1 = (root_est < mid_val and b < mid_val)
            cond2 = (root_est > mid_val and b > mid_val)
            
            if cond1 or cond2:
                root_est = mid_val
            
            # Ensure estimate stays within bounds
            min_val = a if a < b else b
            max_val = a if a > b else b
            
            if root_est < min_val or root_est > max_val:
                root_est = mid_val
                
            f_root = func(root_est)
            
            # Check if result is within tolerance
            if abs(f_root) < Solver.EPS_TOL or abs(b - a) < Solver.EPS_TOL:
                return root_est
            
            # Update the brackets based on the sign
            if fa * f_root < 0:
                b = root_est; fb = f_root
            else:
                a = root_est; fa = f_root
                
        return root_est


class MotionState:
    """
    Represents the kinematic state at a specific node in the trajectory.
    Includes position (p), velocity (v), acceleration (a), jerk (j),
    and the duration (dt) to maintain this jerk.
    """
    def __init__(self, dt=0.0, p=0.0, v=0.0, a=0.0, j=0.0):
        self.dt = float(dt) if dt > 0.0 else 0.0
        self.p = float(p)
        self.v = float(v)
        self.a = float(a)
        self.j = float(j)

    def copy(self):
        return MotionState(self.dt, self.p, self.v, self.a, self.j)

    def __repr__(self):
        return "State(dt=%.4f, P=%.3f, V=%.3f, A=%.3f, J=%.2f)" % \
               (self.dt, self.p, self.v, self.a, self.j)


class S7RTT:
    """
    S-Curve (7-Segment) Trajectory Planner.
    Generates jerk-limited motion profiles (S-curves) for position and velocity control.
    """
    
    EPS_TIME = 1e-9
    EPS_VAL  = 1e-6
    EPS_DIST = 1e-5

    def __init__(self):
        pass

    def _integrate_state(self, state, dt, j):
        """
        Calculates the next state after applying constant jerk 'j' for duration 'dt'.
        """
        p = state.p + state.v * dt + 0.5 * state.a * dt**2 + (1.0/6.0) * j * dt**3
        v = state.v + state.a * dt + 0.5 * j * dt**2
        a = state.a + j * dt
        return MotionState(0.0, p, v, a, 0.0)

    def _integrate_dist_only(self, v0, a0, state_list):
        """
        Lightweight integration used by the solver to calculate total displacement
        without the overhead of creating new MotionState objects.
        """
        dist = 0.0
        v = v0
        a = a0
        for s in state_list:
            t = s.dt
            j = s.j
            dist += v * t + 0.5 * a * t**2 + (1.0/6.0) * j * t**3
            v += a * t + 0.5 * j * t**2
            a += j * t
        return dist

    def _build_profile(self, v_start, a_start, v_target, a_max, j_max):
        """
        Generates the velocity change profile (acceleration or deceleration phase).
        Returns a list of segments (ramp up jerk, constant acc, ramp down jerk).
        Handles both Trapezoidal acceleration (reaching a_max) and Triangular acceleration.
        """
        shape_states = []
        
        # 1. Clamp start acceleration to limits for safety
        acc_clamped = a_start
        if acc_clamped > a_max: acc_clamped = a_max
        if acc_clamped < -a_max: acc_clamped = -a_max
        
        # 2. Calculate velocity reached if we immediately reduce acceleration to zero
        t_to_zero = abs(acc_clamped) / j_max
        j_to_zero = -j_max if acc_clamped > 0 else j_max
        dv_base = acc_clamped * t_to_zero + 0.5 * j_to_zero * t_to_zero**2
        v_base = v_start + dv_base
        
        # 3. Determine direction of the required profile
        direction = 1.0
        if v_target < v_base:
            direction = -1.0
            
        if direction > 0:
            j_up, j_down = j_max, -j_max
        else:
            j_up, j_down = -j_max, j_max
            
        # 4. Calculate constraints and thresholds
        v_req_total = v_target - v_start
        a_limit = a_max * direction
        
        # Time required to ramp from current acceleration to limit, and limit to zero
        t1_max = (a_limit - acc_clamped) / j_up
        t3_max = abs(a_limit) / j_max
        
        # Velocity change provided by a full trapezoidal profile without the constant acceleration phase
        dv_trapezoid = (acc_clamped * t1_max + 0.5 * j_up * t1_max**2) + \
                       (a_limit * t3_max + 0.5 * j_down * t3_max**2)
        
        # 5. Determine if we need a flat acceleration phase (Trapezoidal)
        needs_flat = False
        if direction > 0:
            if v_req_total > dv_trapezoid: needs_flat = True
        else:
            if v_req_total < dv_trapezoid: needs_flat = True
            
        if needs_flat:
            # Trapezoidal Profile (reaches max acceleration)
            v_missing = v_req_total - dv_trapezoid
            t_flat = v_missing / a_limit
            
            if t1_max > S7RTT.EPS_TIME: 
                shape_states.append(MotionState(t1_max, 0,0,0, j_up))
            if t_flat > S7RTT.EPS_TIME: 
                shape_states.append(MotionState(t_flat, 0,0,0, 0.0))
            if t3_max > S7RTT.EPS_TIME: 
                shape_states.append(MotionState(t3_max, 0,0,0, j_down))
        else:
            # Triangular Profile (does not reach max acceleration)
            term = v_req_total * j_up + 0.5 * acc_clamped**2
            if term < 0: term = 0.0
            a_peak_mag = math.sqrt(term)
            a_peak = a_peak_mag if direction > 0 else -a_peak_mag
                
            t1 = (a_peak - acc_clamped) / j_up
            t3 = (0.0 - a_peak) / j_down
            
            if t1 > S7RTT.EPS_TIME: 
                shape_states.append(MotionState(t1, 0,0,0, j_up))
            if t3 > S7RTT.EPS_TIME: 
                shape_states.append(MotionState(t3, 0,0,0, j_down))
            
        return shape_states

    def _compute_trajectory_dist(self, start_state, v_peak, target_v, a_max, j_max):
        """
        Simulates a full trajectory: Accelerate to v_peak, then Decelerate to target_v.
        Returns the total distance traveled and the shape segments.
        """
        # Phase 1: Start -> Peak Velocity
        acc_shape = self._build_profile(start_state.v, start_state.a, v_peak, a_max, j_max)
        
        # Integrate to find state at peak velocity
        v_mid = start_state.v
        a_mid = start_state.a
        for s in acc_shape:
            v_mid += a_mid * s.dt + 0.5 * s.j * s.dt**2
            a_mid += s.j * s.dt
            
        # Phase 2: Peak Velocity -> Target Velocity
        dec_shape = self._build_profile(v_mid, a_mid, target_v, a_max, j_max)
        
        # Calculate total distance
        d1 = self._integrate_dist_only(start_state.v, start_state.a, acc_shape)
        d2 = self._integrate_dist_only(v_mid, a_mid, dec_shape)
        
        return (d1 + d2), acc_shape, dec_shape

    def plan(self, start_state, target_p, target_v, v_max, a_max, j_max):
        """
        Planning Method: Position Control.
        Calculates a trajectory to move from start_state to target_p with final velocity target_v.
        """
        if v_max <= 0 or a_max <= 0 or j_max <= 0:
            return []

        current_state = start_state.copy()
        final_trajectory = []

        # --- 1. Acceleration Recovery ---
        # If current acceleration exceeds limits, ramp it down to safe limits first.
        if current_state.a > a_max + S7RTT.EPS_VAL:
            t_recover = (current_state.a - a_max) / j_max
            rec_state = current_state.copy()
            rec_state.dt = t_recover
            rec_state.j = -j_max
            final_trajectory.append(rec_state)
            
            current_state = self._integrate_state(current_state, t_recover, -j_max)
            current_state.a = a_max
            
        elif current_state.a < -a_max - S7RTT.EPS_VAL:
            t_recover = (-a_max - current_state.a) / j_max
            rec_state = current_state.copy()
            rec_state.dt = t_recover
            rec_state.j = j_max
            final_trajectory.append(rec_state)
            
            current_state = self._integrate_state(current_state, t_recover, j_max)
            current_state.a = -a_max

        dist_req = target_p - current_state.p

        # --- 2. Boundary Check ---
        # Calculate the distance traveled if we hit +V_max and -V_max respectively.
        d_upper, acc_up, dec_up = self._compute_trajectory_dist(
            current_state, v_max, target_v, a_max, j_max)
            
        d_lower, acc_lo, dec_lo = self._compute_trajectory_dist(
            current_state, -v_max, target_v, a_max, j_max)

        shape_list = []

        # --- 3. Logic Branching ---
        if dist_req > d_upper + S7RTT.EPS_DIST:
            # Distance is large: Accelerate to V_max, Cruise, then Decelerate.
            gap = dist_req - d_upper
            t_cruise = gap / v_max
            
            shape_list.extend(acc_up)
            if t_cruise > S7RTT.EPS_TIME: 
                shape_list.append(MotionState(t_cruise, 0,0,0, 0.0))
            shape_list.extend(dec_up)
            
        elif dist_req < d_lower - S7RTT.EPS_DIST:
            # Distance is large (negative): Accelerate to -V_max, Cruise, then Decelerate.
            gap = dist_req - d_lower
            t_cruise = gap / (-v_max)
            
            shape_list.extend(acc_lo)
            if t_cruise > S7RTT.EPS_TIME: 
                shape_list.append(MotionState(t_cruise, 0,0,0, 0.0))
            shape_list.extend(dec_lo)
            
        else:
            # Distance is within bounds: No cruising phase at V_max is needed.
            # Solve for the specific peak velocity that yields the exact distance.
            def get_error(v_p):
                d, _, _ = self._compute_trajectory_dist(
                    current_state, v_p, target_v, a_max, j_max)
                return d - dist_req
            
            max_abs_v = max(v_max, abs(current_state.v), abs(target_v)) * 2.0
            best_v = Solver.solve_monotonic_brent(get_error, -max_abs_v, max_abs_v)
            
            #best_v = Solver.solve_monotonic_brent(get_error, -v_max, v_max)
            #if abs(best_v - current_state.v) < 0.1:
            #    best_v = current_state.v
            
            _, acc_fin, dec_fin = self._compute_trajectory_dist(
                current_state, best_v, target_v, a_max, j_max)
            
            shape_list.extend(acc_fin)
            shape_list.extend(dec_fin)

        # --- 4. Stitching ---
        # Convert relative shape segments into absolute trajectory points.
        for shape in shape_list:
            node = current_state.copy()
            node.dt = shape.dt
            node.j = shape.j
            final_trajectory.append(node)
            
            current_state = self._integrate_state(current_state, shape.dt, shape.j)

        return final_trajectory

    def plan_velocity(self, start_state, target_v, v_max, a_max, j_max):
        """
        Planning Method: Velocity Control.
        Reaches target_v as fast as possible without considering position.
        """
        if v_max <= 0 or a_max <= 0 or j_max <= 0:
            return []

        target_v = max(-v_max, min(v_max, target_v))
        current_state = start_state.copy()
        final_trajectory = []

        # --- 1. Acceleration Recovery ---
        if current_state.a > a_max + S7RTT.EPS_VAL:
            t_recover = (current_state.a - a_max) / j_max
            rec_state = current_state.copy()
            rec_state.dt = t_recover
            rec_state.j = -j_max
            final_trajectory.append(rec_state)
            
            current_state = self._integrate_state(current_state, t_recover, -j_max)
            current_state.a = a_max
            
        elif current_state.a < -a_max - S7RTT.EPS_VAL:
            t_recover = (-a_max - current_state.a) / j_max
            rec_state = current_state.copy()
            rec_state.dt = t_recover
            rec_state.j = j_max
            final_trajectory.append(rec_state)
            
            current_state = self._integrate_state(current_state, t_recover, j_max)
            current_state.a = -a_max

        # --- 2. Build Profile ---
        shape_list = self._build_profile(
            current_state.v, 
            current_state.a, 
            target_v, 
            a_max, 
            j_max
        )

        # --- 3. Stitching ---
        for shape in shape_list:
            node = current_state.copy()
            node.dt = shape.dt
            node.j = shape.j
            final_trajectory.append(node)
            
            current_state = self._integrate_state(current_state, shape.dt, shape.j)

        return final_trajectory

    def at_time(self, trajectory, dt):
        """
        Samples the trajectory at a specific time 'dt' from the start.
        If dt exceeds trajectory duration, extrapolates with constant velocity.
        """
        if not trajectory:
            return MotionState()

        # Handle time before start (clamp to 0)
        if dt <= 1e-12:
            res = trajectory[0].copy()
            res.dt = 0 
            res.j = 0
            return res

        t_remaining = dt
        
        # Iterate segments to find the active one
        for node in trajectory:
            if t_remaining <= node.dt:
                t = t_remaining
                j = node.j
                
                p0, v0, a0 = node.p, node.v, node.a
                
                p = p0 + v0 * t + 0.5 * a0 * t**2 + (1.0/6.0) * j * t**3
                v = v0 + a0 * t + 0.5 * j * t**2
                a = a0 + j * t
                
                return MotionState(0.0, p, v, a, j)
            
            else:
                t_remaining -= node.dt

        # Handle Overshoot: Extrapolate from the last known state
        last_node = trajectory[-1]
        t = last_node.dt
        j = last_node.j
        
        p_end = last_node.p + last_node.v * t + 0.5 * last_node.a * t**2 + (1.0/6.0) * j * t**3
        v_end = last_node.v + last_node.a * t + 0.5 * j * t**2
        
        # Assume constant velocity after trajectory ends
        p_final = p_end + v_end * t_remaining
        
        return MotionState(0.0, p_final, v_end, 0.0, 0.0)
