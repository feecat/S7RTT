import tkinter as tk
from tkinter import ttk, messagebox
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

try:
    from ruckig import Ruckig, InputParameter, Trajectory, Result
    RUCKIG_AVAILABLE = True
except ImportError:
    RUCKIG_AVAILABLE = False
    print("Warning: 'ruckig' library not found. Please install via 'pip install ruckig'")

try:
    from S7RTT import S7RTT, MotionState
    S7RTT_AVAILABLE = True
except ImportError:
    S7RTT_AVAILABLE = False
    class MotionState:
        def __init__(self, p, v, a, j): pass
    class S7RTT:
        def plan(self, *args): raise NotImplementedError("S7RTT not found")

class CompareGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Planner Comparison: S7RTT vs Ruckig")
        self.root.geometry("1400x950")
        
        if S7RTT_AVAILABLE:
            self.planner_s7 = S7RTT()
        
        self._init_ui()
        self._init_plot()

    def _init_ui(self):
        control_frame = ttk.LabelFrame(self.root, text="Motion Parameters")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)

        def create_input(parent, label, row, col, default_val):
            ttk.Label(parent, text=label).grid(row=row, column=col, padx=5, pady=5, sticky="e")
            var = tk.DoubleVar(value=default_val)
            entry = ttk.Entry(parent, textvariable=var, width=10)
            entry.grid(row=row, column=col+1, padx=5, pady=5)
            return var

        ttk.Label(control_frame, text="[Constraints]", foreground="blue").grid(row=0, column=0, sticky="w", padx=5)
        self.v_max = create_input(control_frame, "Max Vel:", 0, 1, 10.0)
        self.a_max = create_input(control_frame, "Max Acc:", 0, 3, 10.0)
        self.j_max = create_input(control_frame, "Max Jerk:", 0, 5, 10.0)

        ttk.Label(control_frame, text="[Start State]", foreground="green").grid(row=1, column=0, sticky="w", padx=5)
        self.p_start = create_input(control_frame, "Start Pos:", 1, 1, 0.0)
        self.v_start = create_input(control_frame, "Start Vel:", 1, 3, 0.0)
        self.a_start = create_input(control_frame, "Start Acc:", 1, 5, 0.0)

        ttk.Label(control_frame, text="[Target State]", foreground="red").grid(row=2, column=0, sticky="w", padx=5)
        self.p_target = create_input(control_frame, "Target Pos:", 2, 1, 100.0)
        self.v_target = create_input(control_frame, "Target Vel:", 2, 3, 0.0)
        ttk.Label(control_frame, text="(End Acc = 0)").grid(row=2, column=5, sticky="w", padx=5)

        btn_frame = ttk.Frame(control_frame)
        btn_frame.grid(row=0, column=7, rowspan=3, padx=30, sticky="ns")
        
        btn_calc = ttk.Button(btn_frame, text="COMPARE & PLOT", command=self._on_calculate)
        btn_calc.pack(fill=tk.X, pady=5, ipady=5)

        stats_frame = ttk.Frame(self.root)
        stats_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=2)
        
        self.lbl_s7 = ttk.Label(stats_frame, text="S7RTT: N/A", font=("Arial", 12, "bold"), foreground="#1f77b4")
        self.lbl_s7.pack(side=tk.LEFT, padx=50)
        
        self.lbl_rk = ttk.Label(stats_frame, text="Ruckig: N/A", font=("Arial", 12, "bold"), foreground="#d62728")
        self.lbl_rk.pack(side=tk.RIGHT, padx=50)

    def _init_plot(self):
        self.fig, self.axs = plt.subplots(3, 2, figsize=(12, 8), sharex='col')
        self.fig.subplots_adjust(hspace=0.2, wspace=0.2, left=0.05, right=0.98, top=0.92, bottom=0.05)
        
        self.titles = ['Position', 'Velocity', 'Acceleration']
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c'] # Blue, Orange, Green
        
        self.axs[0, 0].set_title("S7RTT Planner", fontsize=14, fontweight='bold')
        self.axs[0, 1].set_title("Ruckig Planner", fontsize=14, fontweight='bold')

        self.lines = {}
        
        for col in range(2):
            for row in range(3):
                ax = self.axs[row, col]
                ax.set_ylabel(self.titles[row])
                ax.grid(True, linestyle=':', alpha=0.6)
                if row == 2:
                    ax.set_xlabel("Time (s)")
                
                line, = ax.plot([], [], color=self.colors[row], lw=2)
                self.lines[(row, col)] = line

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _calc_s7rtt(self, dt=0.01):
        if not S7RTT_AVAILABLE:
            return None, None, None, None, "Library Missing"
        
        try:
            start = MotionState(
                dt=0.0, 
                p=self.p_start.get(), 
                v=self.v_start.get(), 
                a=self.a_start.get(), 
                j=0.0
            )
            
            tgt_p = self.p_target.get()
            tgt_v = self.v_target.get()
            
            trajectory = self.planner_s7.plan(
                start, tgt_p, tgt_v, 
                self.v_max.get(), self.a_max.get(), self.j_max.get()
            )
            
            if not trajectory:
                print("S7RTT Plan Failed: Empty trajectory")
                return None, None, 0.0

            total_duration = sum(node.dt for node in trajectory)
            
            times = []
            pos = []
            vel = []
            acc = []
            
            t_curr = 0.0
            while t_curr <= total_duration + 1e-9:
                state = self.planner_s7.at_time(trajectory, t_curr)
                
                times.append(t_curr)
                pos.append(state.p)
                vel.append(state.v)
                acc.append(state.a)
                
                t_curr += dt
            
            return times, [pos, vel, acc], total_duration
            
        except Exception as e:
            print(f"S7RTT Error: {e}")
            import traceback
            traceback.print_exc()
            return None, None, 0.0

    def _calc_ruckig(self, dt=0.01):
        if not RUCKIG_AVAILABLE:
            return None, None, "Library Missing"
            
        try:
            otg = Ruckig(1) # 1 DoF
            inp = InputParameter(1)
            
            inp.current_position = [self.p_start.get()]
            inp.current_velocity = [self.v_start.get()]
            inp.current_acceleration = [self.a_start.get()]
            
            inp.target_position = [self.p_target.get()]
            inp.target_velocity = [self.v_target.get()]
            inp.target_acceleration = [0.0]
            
            inp.max_velocity = [self.v_max.get()]
            inp.max_acceleration = [self.a_max.get()]
            inp.max_jerk = [self.j_max.get()]
            
            traj = Trajectory(1)
            
            result = otg.calculate(inp, traj)
            
            if result.value < 0:
                raise ValueError(f"Ruckig calculation failed with error code: {result}")
            
            duration = traj.duration
            t_vals = np.arange(0, duration + dt, dt)
            if len(t_vals) > 0 and t_vals[-1] < duration:
                t_vals = np.append(t_vals, duration)
            elif len(t_vals) == 0:
                t_vals = np.array([0.0])
                
            pos, vel, acc = [], [], []
            for t in t_vals:
                t_clamped = min(t, duration)
                new_state = traj.at_time(t_clamped)
                pos.append(new_state[0]) # Position
                vel.append(new_state[1]) # Velocity
                acc.append(new_state[2]) # Acceleration
                
            return t_vals, [pos, vel, acc], duration

        except Exception as e:
            print(f"Ruckig Error: {e}")
            return None, None, 0.0

    def _on_calculate(self):
        s7_t, s7_data, s7_dur = self._calc_s7rtt()
        
        rk_t, rk_data, rk_dur = self._calc_ruckig()
        
        if s7_data:
            self.lbl_s7.config(text=f"S7RTT Duration: {s7_dur:.4f} s")
        else:
            self.lbl_s7.config(text="S7RTT: Error or N/A")
            
        if rk_data:
            self.lbl_rk.config(text=f"Ruckig Duration: {rk_dur:.4f} s")
        else:
            self.lbl_rk.config(text="Ruckig: Error or N/A")
            
        self._update_plot_column(0, s7_t, s7_data)
        self._update_plot_column(1, rk_t, rk_data)
        
        self.canvas.draw()

    def _update_plot_column(self, col_idx, t_data, y_datas):
        if t_data is None or y_datas is None:
            for row in range(3):
                self.lines[(row, col_idx)].set_data([], [])
                self.axs[row, col_idx].relim()
                self.axs[row, col_idx].autoscale_view()
            return

        for row in range(3):
            line = self.lines[(row, col_idx)]
            ax = self.axs[row, col_idx]
            
            line.set_data(t_data, y_datas[row])
            
            ax.relim()
            ax.autoscale_view()

if __name__ == "__main__":
    root = tk.Tk()
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except:
        pass
        
    app = CompareGUI(root)
    root.mainloop()