# S7RTT.py
import math

# ==========================================
# 1. Mathematical Solver
# ==========================================

class Solver:
    """
    Utility class for numerical root finding.
    """
    
    # Constants for convergence criteria
    EPS_TOL = 1e-8
    EPS_CONVERGE = 1e-12
    MAX_ITER = 50

    @staticmethod
    def solve_monotonic_brent(func, low, high):
        """
        Brent's Method for finding a root (func(x) = 0) within [low, high].
        
        Args:
            func: Callable function taking a float and returning a float.
            low: Lower bound of interval.
            high: Upper bound of interval.
            
        Returns:
            float: The estimated root.
        """
        f_low = func(low)
        f_high = func(high)
        
        # If signs match, we cannot bracket the root strictly. 
        # Return the side closer to zero.
        if f_low * f_high > 0:
            if abs(f_low) < abs(f_high):
                return low
            else:
                return high
            
        a = low
        b = high
        fa = f_low
        fb = f_high
        
        root_est = b
        
        for _ in range(Solver.MAX_ITER):
            # Ensure 'b' is always the better guess
            if abs(fa) < abs(fb):
                # Swap a and b
                temp = a; a = b; b = temp
                temp = fa; fa = fb; fb = temp
            
            # Check convergence
            if abs(fb - fa) < Solver.EPS_CONVERGE:
                root_est = b
                break
            else:
                # Inverse quadratic interpolation / Secant step
                # Guard against division by zero
                if abs(fa - fb) > 1e-14:
                    root_est = b - fb * (b - a) / (fb - fa)
                else:
                    root_est = b
            
            mid_val = 0.5 * (a + b)
            
            # Bisection fallback conditions (standard Brent logic)
            cond1 = (root_est < mid_val and b < mid_val)
            cond2 = (root_est > mid_val and b > mid_val)
            
            if cond1 or cond2:
                root_est = mid_val
            
            # Safety: clamp strictly between a and b to prevent divergence
            min_val = a if a < b else b
            max_val = a if a > b else b
            
            if root_est < min_val or root_est > max_val:
                root_est = mid_val
                
            f_root = func(root_est)
            
            # Check if result is within tolerance
            if abs(f_root) < Solver.EPS_TOL or abs(b - a) < Solver.EPS_TOL:
                return root_est
            
            # Tighten the bracket
            if fa * f_root < 0:
                b = root_est
                fb = f_root
            else:
                a = root_est
                fa = f_root
                
        return root_est

# ==========================================
# 2. Data Structures
# ==========================================

class MotionState:
    """
    Represents kinematic state: P, V, A, J.
    """
    def __init__(self, p=0.0, v=0.0, a=0.0, j=0.0):
        self.p = float(p)
        self.v = float(v)
        self.a = float(a)
        self.j = float(j)

    def copy(self):
        return MotionState(self.p, self.v, self.a, self.j)

    def __repr__(self):
        return "State(P=%.3f, V=%.3f, A=%.3f)" % (self.p, self.v, self.a)

class Segment:
    """
    Represents a constant jerk segment (Jerk, Duration).
    """
    def __init__(self, jerk, duration):
        self.jerk = float(jerk)
        # Ensure time is never negative
        self.duration = duration if duration > 0.0 else 0.0

    def __repr__(self):
        return "Seg(J=%.2f, T=%.4f)" % (self.jerk, self.duration)

# ==========================================
# 3. Planner Class
# ==========================================

class S7RTT:
    """
    S-Curve Trajectory Planner (7-Segment Profile).
    Handles Asymmetric Jerk limits and non-zero initial acceleration.
    """
    
    # Global epsilons for decision making
    EPS_TIME = 1e-9
    EPS_VAL  = 1e-6
    EPS_DIST = 1e-5

    def __init__(self):
        pass

    def integrate_path(self, start_state, segments):
        """
        Forward integration of motion equations.
        Returns the final state after applying all segments.
        """
        curr = start_state.copy()
        for seg in segments:
            t = seg.duration
            j = seg.jerk
            
            # Standard kinematics integration (Taylor series)
            curr.p += curr.v * t + 0.5 * curr.a * t**2 + (1.0/6.0) * j * t**3
            curr.v += curr.a * t + 0.5 * j * t**2
            curr.a += j * t
            curr.j = j
        return curr

    def _integrate_dist(self, v0, a0, segments):
        """
        Helper: Calculate distance delta given initial v, a and segments.
        Avoids creating full State objects to save overhead.
        """
        dist = 0.0
        v = v0
        a = a0
        for seg in segments:
            t = seg.duration
            j = seg.jerk
            dist += v * t + 0.5 * a * t**2 + (1.0/6.0) * j * t**3
            v += a * t + 0.5 * j * t**2
            a += j * t
        return dist

    def _build_profile(self, v_start, a_start, v_target, a_max, j_max):
        """
        Constructs the velocity change profile (Acc/Decel segments)
        to go from (v_start, a_start) to (v_target, 0.0).
        """
        segments = []
        
        # 1. Clamp A_start to prevent numerical instability if input is bad
        acc_clamped = a_start
        if acc_clamped > a_max: acc_clamped = a_max
        if acc_clamped < -a_max: acc_clamped = -a_max
        
        # 2. Calculate "Base Velocity"
        # This is the velocity reached if we reduce acceleration to zero immediately.
        # It acts as a heuristic to decide direction.
        t_to_zero = abs(acc_clamped) / j_max
        j_to_zero = -j_max if acc_clamped > 0 else j_max
        
        dv_base = acc_clamped * t_to_zero + 0.5 * j_to_zero * t_to_zero**2
        v_base = v_start + dv_base
        
        # 3. Decide Direction (1.0 for Accel, -1.0 for Decel)
        direction = 1.0
        if v_target < v_base:
            direction = -1.0
            
        # Setup Jerk parameters based on direction
        if direction > 0:
            j_up = j_max
            j_down = -j_max
        else:
            j_up = -j_max
            j_down = j_max
            
        # 4. Solve for A_peak
        v_req_total = v_target - v_start
        a_limit = a_max * direction
        
        # Calculate capacity of a Trapezoidal profile (Acc -> Flat -> Dec)
        # T1: acc_clamped -> a_limit
        t1_max = (a_limit - acc_clamped) / j_up
        # T3: a_limit -> 0
        t3_max = abs(a_limit) / j_max
        
        dv_trapezoid = (acc_clamped * t1_max + 0.5 * j_up * t1_max**2) + \
                       (a_limit * t3_max + 0.5 * j_down * t3_max**2)
        
        # Check if we need a flat top (constant acceleration phase)
        needs_flat = False
        if direction > 0:
            if v_req_total > dv_trapezoid: needs_flat = True
        else:
            if v_req_total < dv_trapezoid: needs_flat = True
            
        if needs_flat:
            # --- Trapezoidal Profile ---
            v_missing = v_req_total - dv_trapezoid
            t_flat = v_missing / a_limit
            
            if t1_max > S7RTT.EPS_TIME: 
                segments.append(Segment(j_up, t1_max))
            if t_flat > S7RTT.EPS_TIME: 
                segments.append(Segment(0.0, t_flat))
            if t3_max > S7RTT.EPS_TIME: 
                segments.append(Segment(j_down, t3_max))
            
        else:
            # --- Triangular Profile ---
            # We cannot reach a_max. Solve for peak acceleration.
            # Robust calculation of a_peak magnitude
            term = v_req_total * j_up + 0.5 * acc_clamped**2
            
            if term < 0: 
                term = 0.0
                
            a_peak_mag = math.sqrt(term)
            a_peak = a_peak_mag if direction > 0 else -a_peak_mag
                
            t1 = (a_peak - acc_clamped) / j_up
            t3 = (0.0 - a_peak) / j_down
            
            if t1 > S7RTT.EPS_TIME: 
                segments.append(Segment(j_up, t1))
            if t3 > S7RTT.EPS_TIME: 
                segments.append(Segment(j_down, t3))
            
        return segments

    def _compute_trajectory_dist(self, start_state, v_peak, target_v, a_max, j_max):
        """
        Calculates total distance and generates segments for:
        Start -> v_peak -> Target
        """
        # Phase 1: Accelerate/Decelerate from start to v_peak
        acc_segs = self._build_profile(start_state.v, start_state.a, v_peak, a_max, j_max)
        
        # Simulate state at the middle point (v_peak achieved)
        # We only need the intermediate State to start the second half
        state_mid = self.integrate_path(start_state, acc_segs)
        
        # Phase 2: Accelerate/Decelerate from v_peak to target_v
        dec_segs = self._build_profile(state_mid.v, state_mid.a, target_v, a_max, j_max)
        
        # Calculate distances
        d1 = self._integrate_dist(start_state.v, start_state.a, acc_segs)
        d2 = self._integrate_dist(state_mid.v, state_mid.a, dec_segs)
        
        return (d1 + d2), acc_segs, dec_segs

    def plan(self, start_state, target_p, target_v, v_max, a_max, j_max):
        """
        Main planning interface.
        Returns a list of Segment objects.
        """
        # Input validation guard
        if v_max <= 0 or a_max <= 0 or j_max <= 0:
            return [] # Or raise Error

        # 'pre_segments' stores segments required to recover from invalid initial states
        pre_segments = [] 
        current_state = start_state.copy()

        # ========================================================
        # 1. Initial State Sanitization
        # ========================================================
        # If initial acceleration violates limits, reduce it first.
        if current_state.a > a_max + S7RTT.EPS_VAL:
            t_recover = (current_state.a - a_max) / j_max
            seg = Segment(-j_max, t_recover)
            pre_segments.append(seg)
            current_state = self.integrate_path(current_state, [seg])
            # Force clamp to avoid float noise
            current_state.a = a_max 
            
        elif current_state.a < -a_max - S7RTT.EPS_VAL:
            t_recover = (-a_max - current_state.a) / j_max
            seg = Segment(j_max, t_recover)
            pre_segments.append(seg)
            current_state = self.integrate_path(current_state, [seg])
            current_state.a = -a_max 

        dist_req = target_p - current_state.p

        # ========================================================
        # 2. Compute Physical Boundaries
        # ========================================================
        
        # Check Max Velocity (Positive limit)
        d_upper, acc_segs_up, dec_segs_up = self._compute_trajectory_dist(
            current_state, v_max, target_v, a_max, j_max)
            
        # Check Min Velocity (Negative limit)
        d_lower, acc_segs_lo, dec_segs_lo = self._compute_trajectory_dist(
            current_state, -v_max, target_v, a_max, j_max)

        plan_segments = []

        # ========================================================
        # 3. Three-Branch Decision Logic
        # ========================================================
        
        if dist_req > d_upper + S7RTT.EPS_DIST:
            # Case 1: Positive Saturation (Need to cruise at +v_max)
            gap = dist_req - d_upper
            t_cruise = gap / v_max
            
            plan_segments.extend(acc_segs_up)
            if t_cruise > S7RTT.EPS_TIME: 
                plan_segments.append(Segment(0.0, t_cruise))
            plan_segments.extend(dec_segs_up)
            
        elif dist_req < d_lower - S7RTT.EPS_DIST:
            # Case 2: Negative Saturation (Need to cruise at -v_max)
            gap = dist_req - d_lower
            t_cruise = gap / (-v_max)
            
            plan_segments.extend(acc_segs_lo)
            if t_cruise > S7RTT.EPS_TIME: 
                plan_segments.append(Segment(0.0, t_cruise))
            plan_segments.extend(dec_segs_lo)
            
        else:
            # Case 3: Peak Search (Trajectory fits within limits)
            # Define error function for the solver
            def get_error(v_p):
                d, _, _ = self._compute_trajectory_dist(
                    current_state, v_p, target_v, a_max, j_max)
                return d - dist_req
            
            # Expand search range (2.0x) to robustly handle initial state overshoot
            max_abs_v = max(v_max, abs(current_state.v), abs(target_v)) * 2.0
            
            best_v = Solver.solve_monotonic_brent(get_error, -max_abs_v, max_abs_v)
            
            _, acc_segs_fin, dec_segs_fin = self._compute_trajectory_dist(
                current_state, best_v, target_v, a_max, j_max)
            
            plan_segments.extend(acc_segs_fin)
            plan_segments.extend(dec_segs_fin)

        # Return combined list
        return pre_segments + plan_segments