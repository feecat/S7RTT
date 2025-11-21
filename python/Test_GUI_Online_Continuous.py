# ==============================================================================
# File Name:    Test_GUI_Online_Continuous.py
# Description:  Online S-Curve Test with Continuous Replanning (Stress Test)
# ==============================================================================

import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import collections
import time

from S7RTT import S7RTT, MotionState

class OnlineTestGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("S7RTT Continuous Replanning (1ms Cycle)")
        self.root.geometry("1200x900")
        
        self.SIM_DT = 0.001  
        self.PLOT_WINDOW = 5.0 
        self.is_running = False
        
        self.last_plot_time = 0.0
        self.plot_interval = 0.05
        
        self.planner = S7RTT()
        
        self.current_state = MotionState(0, 0, 0, 0, 0)
        
        self.max_points = int(self.PLOT_WINDOW / self.SIM_DT)
        self.history_t = collections.deque(maxlen=self.max_points)
        self.history_p = collections.deque(maxlen=self.max_points)
        self.history_v = collections.deque(maxlen=self.max_points)
        self.history_a = collections.deque(maxlen=self.max_points)
        self.history_j = collections.deque(maxlen=self.max_points)
        
        self.global_time = 0.0

        self._init_ui()
        self._init_plot()

    def _init_ui(self):
        control_panel = ttk.Frame(self.root, width=250)
        control_panel.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        btn_frame = ttk.LabelFrame(control_panel, text="Simulation Control")
        btn_frame.pack(fill=tk.X, pady=5)
        
        self.btn_start = ttk.Button(btn_frame, text="START Loop", command=self._start_simulation)
        self.btn_start.pack(fill=tk.X, padx=5, pady=5)
        
        self.btn_stop = ttk.Button(btn_frame, text="STOP", command=self._stop_simulation, state=tk.DISABLED)
        self.btn_stop.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(btn_frame, text="Reset State", command=self._reset_state).pack(fill=tk.X, padx=5, pady=5)

        target_frame = ttk.LabelFrame(control_panel, text="Dynamic Target")
        target_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(target_frame, text="Target Position:").pack(anchor="w", padx=5)
        self.var_target_p = tk.DoubleVar(value=100.0)
        scale_p = tk.Scale(target_frame, variable=self.var_target_p, 
                           from_=-200, to=200, orient=tk.HORIZONTAL,
                           resolution=0.1)
        scale_p.pack(fill=tk.X, padx=5)
        
        ttk.Label(target_frame, text="Target Velocity:").pack(anchor="w", padx=5)
        self.var_target_v = tk.DoubleVar(value=0.0)
        scale_v = tk.Scale(target_frame, variable=self.var_target_v, 
                           from_=-20, to=20, orient=tk.HORIZONTAL,
                           resolution=0.1)
        scale_v.pack(fill=tk.X, padx=5)

        const_frame = ttk.LabelFrame(control_panel, text="Constraints")
        const_frame.pack(fill=tk.X, pady=10)
        
        self.v_max = self._create_input(const_frame, "Max Vel:", 1000.0)
        self.a_max = self._create_input(const_frame, "Max Acc:", 5000.0)
        self.j_max = self._create_input(const_frame, "Max Jerk:", 20000.0)

        stats_frame = ttk.LabelFrame(control_panel, text="Current State")
        stats_frame.pack(fill=tk.X, pady=10)
        self.lbl_cur_p = ttk.Label(stats_frame, text="P: 0.00")
        self.lbl_cur_p.pack(anchor="w", padx=5)
        self.lbl_cur_v = ttk.Label(stats_frame, text="V: 0.00")
        self.lbl_cur_v.pack(anchor="w", padx=5)
        self.lbl_cur_a = ttk.Label(stats_frame, text="A: 0.00")
        self.lbl_cur_a.pack(anchor="w", padx=5)

    def _create_input(self, parent, label, default_val):
        frame = ttk.Frame(parent)
        frame.pack(fill=tk.X, padx=5, pady=2)
        ttk.Label(frame, text=label, width=10).pack(side=tk.LEFT)
        var = tk.DoubleVar(value=default_val)
        ttk.Entry(frame, textvariable=var, width=10).pack(side=tk.RIGHT)
        return var

    def _init_plot(self):
        self.fig, self.axs = plt.subplots(4, 1, figsize=(8, 8), sharex=True)
        self.fig.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.05, hspace=0.2)
        titles = ['Position', 'Velocity', 'Acceleration', 'Jerk']
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        self.lines = []
        for i, ax in enumerate(self.axs):
            ax.set_ylabel(titles[i])
            ax.grid(True, linestyle=':', alpha=0.6)
            line, = ax.plot([], [], color=colors[i], lw=2)
            self.lines.append(line)
        plot_frame = ttk.Frame(self.root)
        plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _start_simulation(self):
        if self.is_running: return
        self.is_running = True
        self.btn_start.config(state=tk.DISABLED)
        self.btn_stop.config(state=tk.NORMAL)
        self._run_cycle()

    def _stop_simulation(self):
        self.is_running = False
        self.btn_start.config(state=tk.NORMAL)
        self.btn_stop.config(state=tk.DISABLED)

    def _reset_state(self):
        self.current_state = MotionState(0, 0, 0, 0, 0)
        self.history_t.clear()
        self.history_p.clear()
        self.history_v.clear()
        self.history_a.clear()
        self.history_j.clear()
        self.global_time = 0.0
        
        self._update_plot()
        self._update_labels()

    def _run_cycle(self):
        if not self.is_running:
            return

        loop_start_time = time.time()

        try:
            tgt_p = self.var_target_p.get()
            tgt_v = self.var_target_v.get()
            v_max = self.v_max.get()
            a_max = self.a_max.get()
            j_max = self.j_max.get()
            
            traj = self.planner.plan(
                self.current_state, 
                tgt_p, 
                tgt_v, 
                v_max, a_max, j_max
            )
            
            if traj:
                next_state = self.planner.at_time(traj, self.SIM_DT)
            else:
                next_state = self.current_state.copy()
                next_state.v = 0
                next_state.a = 0
                next_state.j = 0

            self.current_state = next_state
            self.global_time += self.SIM_DT

            self.history_t.append(self.global_time)
            self.history_p.append(self.current_state.p)
            self.history_v.append(self.current_state.v)
            self.history_a.append(self.current_state.a)
            self.history_j.append(self.current_state.j)

        except Exception as e:
            print(f"Error in cycle: {e}")
            self._stop_simulation()

        finally:
            current_real_time = time.time()
            
            if current_real_time - self.last_plot_time > self.plot_interval:
                if len(self.history_t) > 0:
                    self._update_plot()
                    self._update_labels()
                self.last_plot_time = current_real_time
            
            elapsed = current_real_time - loop_start_time
            
            wait_ms = int((self.SIM_DT - elapsed) * 1000)
            
            if wait_ms < 1: 
                wait_ms = 1
            
            if self.is_running:
                self.root.after(wait_ms, self._run_cycle)

    def _update_labels(self):
        s = self.current_state
        self.lbl_cur_p.config(text=f"P: {s.p:.3f}")
        self.lbl_cur_v.config(text=f"V: {s.v:.3f}")
        self.lbl_cur_a.config(text=f"A: {s.a:.3f}")

    def _update_plot(self):
        if not self.history_t: return
        
        t_data = list(self.history_t)
        
        step = 1
        if len(t_data) > 1000:
            step = int(len(t_data) / 500)
        
        t_plot = t_data[::step]
        p_plot = list(self.history_p)[::step]
        v_plot = list(self.history_v)[::step]
        a_plot = list(self.history_a)[::step]
        j_plot = list(self.history_j)[::step]

        data_sources = [p_plot, v_plot, a_plot, j_plot]
        
        min_t, max_t = t_plot[0], t_plot[-1]
        
        for i, line in enumerate(self.lines):
            line.set_data(t_plot, data_sources[i])
            self.axs[i].set_xlim(min_t, max_t + 0.05)
            
            y_data = data_sources[i]
            if y_data:
                mn, mx = min(y_data), max(y_data)
                span = mx - mn
                if span < 0.1: span = 0.1
                self.axs[i].set_ylim(mn - span*0.1, mx + span*0.1)
                
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = OnlineTestGUI(root)
    root.mainloop()
