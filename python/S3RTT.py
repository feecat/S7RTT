# ==============================================================================
# File Name:    S3RTT.py
# Author:       feecat
# Version:      V1.0
# Description:  Trapezoidal Velocity Profile Generator
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
# Description:
#   This module generates trajectory data using a Trapezoidal Velocity Profile 
#   (T-Curve) rather than an S-Curve (Sigmoid/Jerk-limited profile).
#   
#   The T-Curve approach significantly reduces computational complexity and 
#   CPU consumption, making it highly suitable for resource-constrained 
#   embedded devices.
#
# Note:
#   This Python version is intended for algorithm verification and testing.
#   For the production-ready embedded implementation, please refer to the 
#   C language header file 'S3RTT.h'.
# ==============================================================================
import tkinter as tk
from tkinter import ttk
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

S3_EPS = 1e-5

class MotionState:
    def __init__(self, dt=0.0, p=0.0, v=0.0, a=0.0):
        self.dt = dt
        self.p = p
        self.v = v
        self.a = a

class Path:
    def __init__(self):
        self.nodes = [] # List of MotionState

    def push(self, dt, p, v, a):
        if dt > S3_EPS:
            self.nodes.append(MotionState(dt, p, v, a))
            
    def total_time(self):
        return sum(n.dt for n in self.nodes)

def s3_sign(x):
    if x >= 0: return 1.0
    return -1.0

def s3_plan(start_state, target_p, target_v, v_max, a_max):
    path = Path()
    
    # Extract parameters and convert to float
    vi, pi = float(start_state.v), float(start_state.p)
    pf, vf = float(target_p), float(target_v)
    acc, v_cap = abs(float(a_max)), abs(float(v_max))
    
    # Clamp the target velocity within the allowable maximum velocity limits
    vf = max(-v_cap, min(v_cap, vf)) 

    # If acceleration is negligible, return a constant velocity path immediately
    if acc < S3_EPS:
        path.push(0.0, pi, vi, 0.0)
        return path

    dist = pf - pi
    # Determine the direction of motion (1 for positive, -1 for negative)
    s = s3_sign(dist) if abs(dist) > S3_EPS else 1.0

    # ==========================================
    # Phase 1: Kinematic Violation Handling
    # ==========================================
    
    # Check 1: Wrong Direction or Overshoot
    # Detect if the current velocity is moving away from the target or if inertia makes stopping impossible
    is_wrong_dir = (vi * s < -S3_EPS)
    is_overshoot = False
    
    # Only check for overshoot if we are moving towards the target
    if not is_wrong_dir and (vi * s > S3_EPS):
        # Case A: Moving towards target, but target velocity requires moving backwards (must stop and reverse)
        if vf * s < -S3_EPS:
            is_overshoot = True
        # Case B: Current kinetic energy is too high to slow down to 'vf' within distance 'dist'
        # Formula check: v_initial^2 > v_final^2 + 2 * a * distance
        elif vi**2 > vf**2 + 2.0 * acc * abs(dist) + S3_EPS:
            is_overshoot = True

    if is_wrong_dir or is_overshoot:
        # Strategy: Perform a full stop, then plan recursively from the stopped state.
        # Calculate the distance required to come to a complete stop
        d_stop = (vi**2) / (2.0 * acc) * s3_sign(vi)
        stop_pos = pi + d_stop
        
        # Add the braking segment to zero velocity
        path.push(abs(vi)/acc, pi, vi, -s3_sign(vi)*acc)
        # Recursively plan from the stop position to the original target
        path.nodes.extend(s3_plan(MotionState(0, stop_pos, 0, 0), pf, vf, v_max, a_max).nodes)
        return path

    # Check 2: Insufficient Run-up Distance
    # If we need to accelerate to a high 'vf' but don't have enough distance to reach it
    if abs(vf) > abs(vi) and (vf * s > 0):
        max_reachable_sq = vi**2 + 2.0 * acc * abs(dist)
        if vf**2 > max_reachable_sq + S3_EPS:
            # Strategy: Move backward first to gain runway (Run-up).
            # Calculate the turnaround point needed to generate enough acceleration distance
            p_turn = pf - s * ((vf**2)/(2.0*acc))
            
            # Plan segment to the turnaround point
            path.nodes.extend(s3_plan(MotionState(0, pi, vi, 0), p_turn, 0, v_max, a_max).nodes)
            # Plan segment from turnaround point to the target
            path.nodes.extend(s3_plan(MotionState(0, p_turn, 0, 0), pf, vf, v_max, a_max).nodes)
            return path

    # ==========================================
    # Phase 2: Standard Profile Generation
    # ==========================================
    
    # Calculate the theoretical peak velocity squared for a triangular profile
    # Derived from: 2 * a * d = (v_peak^2 - v_i^2) + (v_peak^2 - v_f^2)
    vp_sq = (2.0 * acc * abs(dist) + vi**2 + vf**2) / 2.0
    vp = math.sqrt(max(0.0, vp_sq))
    
    # Apply velocity constraints (Transition from Triangular to Trapezoidal profile)
    v_peak = min(vp, v_cap) * s

    # 1. Acceleration Phase: Move from current velocity to peak velocity
    t1 = abs(v_peak - vi) / acc
    path.push(t1, pi, vi, s3_sign(v_peak - vi) * acc)
    
    # Update current position after acceleration
    p_curr = pi + (vi + v_peak) * t1 * 0.5
    
    # 2. Cruising Phase: Constant velocity (exists only if profile is Trapezoidal)
    # Calculate distance required for the final deceleration phase
    t3 = abs(vf - v_peak) / acc
    d3 = (v_peak + vf) * t3 * 0.5
    
    # Calculate remaining distance available for cruising
    d_cruise = (pf - p_curr) - d3
    
    # If there is meaningful distance left, insert a cruise segment
    if d_cruise * s > S3_EPS:
        path.push(abs(d_cruise)/abs(v_peak), p_curr, v_peak, 0.0)
        p_curr += d_cruise

    # 3. Deceleration/Adjustment Phase: Move from peak velocity to target velocity
    path.push(t3, p_curr, v_peak, s3_sign(vf - v_peak) * acc)

    return path

def s3_plan_velocity(start_state, target_v, v_max, a_max):
    path = Path()
    
    vi = start_state.v
    pi = start_state.p
    acc = abs(a_max)
    v_cap = abs(v_max)
    
    vf = max(-v_cap, min(v_cap, target_v))
    if acc < S3_EPS or abs(vf - vi) < S3_EPS:
        path.push(0.0, pi, vi, 0.0)
        return path
    duration = abs(vf - vi) / acc
    sign = 1.0 if vf > vi else -1.0
    path.push(duration, pi, vi, sign * acc)
    
    return path

def s3_at_time(path, t):
    if not path.nodes: return MotionState()
    if t < 0: t = 0
    
    elapsed = 0.0
    for node in path.nodes:
        if t < elapsed + node.dt:
            dt = t - elapsed
            res = MotionState()
            res.a = node.a
            res.v = node.v + node.a * dt
            res.p = node.p + node.v * dt + 0.5 * node.a * dt * dt
            return res
        elapsed += node.dt
        
    # Extrapolate last
    last = path.nodes[-1]
    dt_seg = last.dt
    v_end = last.v + last.a * dt_seg
    p_end = last.p + last.v * dt_seg + 0.5 * last.a * dt_seg * dt_seg
    
    dt_ex = t - elapsed
    res = MotionState()
    res.a = 0.0
    res.v = v_end
    res.p = p_end + v_end * dt_ex
    return res

# ==========================================
# UI & Plotting
# ==========================================

class App:
    def __init__(self, root):
        self.root = root
        self.root.title("S3RTT C-Port Simulation (Overshoot Logic)")
        self.root.geometry("1200x800")
        
        # Controls Frame
        control_frame = ttk.LabelFrame(root, text="Parameters")
        control_frame.pack(fill="x", padx=10, pady=5)
        
        self.entries = {}
        params = [
            ("Start Pos", "0.0"),
            ("Start Vel", "0.0"),
            ("Target Pos", "1000"),
            ("Target Vel", "0.0"),
            ("Max Vel", "1000.0"),
            ("Max Acc", "10000.0")
        ]
        
        for i, (label, val) in enumerate(params):
            ttk.Label(control_frame, text=label).pack(side="left", padx=5)
            e = ttk.Entry(control_frame, width=8)
            e.insert(0, val)
            e.pack(side="left", padx=5)
            self.entries[label] = e
            
        btn = ttk.Button(control_frame, text="Calculate & Plot", command=self.calculate)
        btn.pack(side="left", padx=20)
        
        # Info Label
        self.info_lbl = ttk.Label(root, text="Ready", foreground="blue")
        self.info_lbl.pack(pady=5)
        
        # Plot Area
        self.fig, self.axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        
    def get_float(self, name):
        try:
            return float(self.entries[name].get())
        except ValueError:
            return 0.0

    def calculate(self):
        start_p = self.get_float("Start Pos")
        start_v = self.get_float("Start Vel")
        target_p = self.get_float("Target Pos")
        target_v = self.get_float("Target Vel")
        v_max = self.get_float("Max Vel")
        a_max = self.get_float("Max Acc")
        
        start_s = MotionState(0, start_p, start_v, 0)
        
        # Run Algorithm
        path = s3_plan(start_s, target_p, target_v, v_max, a_max)
        
        # Sampling
        total_t = path.total_time()
        sim_duration = total_t * 1.1 + 0.0001
        
        times = []
        pos = []
        vel = []
        acc = []
        
        steps = 1000
        dt = sim_duration / steps
        
        final_s = MotionState()
        
        for i in range(steps + 1):
            t = i * dt
            s = s3_at_time(path, t)
            times.append(t)
            pos.append(s.p)
            vel.append(s.v)
            acc.append(s.a)
            if i == steps: final_s = s
            
        # Update Info
        err = abs(final_s.p - target_p)
        self.info_lbl.config(text=f"Total Time: {total_t:.4f}s | Segments: {len(path.nodes)} | Final P: {final_s.p:.4f} (Err: {err:.4f}) | Final V: {final_s.v:.4f}")
        
        # Plotting
        for ax in self.axs: ax.clear()
        
        self.axs[0].plot(times, pos, label="Position", color="blue")
        self.axs[0].axhline(y=target_p, color="green", linestyle="--", label="Target P")
        self.axs[0].set_ylabel("Position")
        self.axs[0].grid(True)
        self.axs[0].legend()
        
        self.axs[1].plot(times, vel, label="Velocity", color="red")
        self.axs[1].axhline(y=target_v, color="orange", linestyle="--", label="Target V")
        self.axs[1].set_ylabel("Velocity")
        self.axs[1].grid(True)
        
        self.axs[2].plot(times, acc, label="Acceleration", color="purple")
        self.axs[2].set_ylabel("Acceleration")
        self.axs[2].set_xlabel("Time (s)")
        self.axs[2].grid(True)
        
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()
