# Test_GUI.py
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from S7RTT import S7RTT, MotionState

class TestGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("S7RTT Planner GUI")
        self.root.geometry("1200x900")
        self.planner = S7RTT()
        self._init_ui()
        self._init_plot()

    def _init_ui(self):
        control_frame = ttk.LabelFrame(self.root, text="Motion Parameters")
        control_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)

        def create_input(parent, label, row, col, default_val):
            ttk.Label(parent, text=label).grid(row=row, column=col, padx=5, pady=5, sticky="e")
            var = tk.DoubleVar(value=default_val)
            entry = ttk.Entry(parent, textvariable=var, width=12)
            entry.grid(row=row, column=col+1, padx=5, pady=5)
            return var

        # Constraints
        ttk.Label(control_frame, text="[Constraints]", foreground="blue").grid(row=0, column=0, sticky="w", padx=5)
        self.v_max = create_input(control_frame, "Max Vel:", 0, 1, 10.0)
        self.a_max = create_input(control_frame, "Max Acc:", 0, 3, 10.0)
        self.j_max = create_input(control_frame, "Max Jerk:", 0, 5, 10.0)

        # Start State
        ttk.Label(control_frame, text="[Start State]", foreground="green").grid(row=1, column=0, sticky="w", padx=5)
        self.p_start = create_input(control_frame, "Start Pos:", 1, 1, 0.0)
        self.v_start = create_input(control_frame, "Start Vel:", 1, 3, 0.0)
        self.a_start = create_input(control_frame, "Start Acc:", 1, 5, 0.0)

        # Target State
        ttk.Label(control_frame, text="[Target State]", foreground="red").grid(row=2, column=0, sticky="w", padx=5)
        self.p_target = create_input(control_frame, "Target Pos:", 2, 1, 100.0)
        self.v_target = create_input(control_frame, "Target Vel:", 2, 3, 0.0)

        # Buttons
        btn_frame = ttk.Frame(control_frame)
        btn_frame.grid(row=0, column=7, rowspan=3, padx=20, sticky="ns")
        ttk.Button(btn_frame, text="GENERATE", command=self._on_calculate).pack(fill=tk.X, pady=5)
        ttk.Button(btn_frame, text="RESET", command=self._on_reset).pack(fill=tk.X, pady=5)
        
        # Status Label
        self.lbl_stats = ttk.Label(control_frame, text="Ready", font=("Arial", 10, "bold"))
        self.lbl_stats.grid(row=3, column=0, columnspan=8, pady=5)

    def _init_plot(self):
        self.fig, self.axs = plt.subplots(4, 1, figsize=(10, 8), sharex=True)
        self.fig.subplots_adjust(hspace=0.15, left=0.08, right=0.95, top=0.95, bottom=0.05)
        self.titles = ['Position', 'Velocity', 'Acceleration', 'Jerk']
        self.colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
        self.lines = []
        for i, ax in enumerate(self.axs):
            ax.set_ylabel(self.titles[i])
            ax.grid(True, linestyle=':', alpha=0.6)
            line, = ax.plot([], [], color=self.colors[i], lw=2)
            self.lines.append(line)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _generate_plot_data(self, trajectory, sampling_dt=0.01):
        """
        Uses at_time to sample the trajectory for plotting.
        """
        if not trajectory:
            return [], [], [], [], []

        total_duration = sum(node.dt for node in trajectory)
        
        plot_duration = total_duration + 1e-5
        
        times = []
        pos = []
        vel = []
        acc = []
        jerk = []

        t_curr = 0.0
        while t_curr <= plot_duration:
            state = self.planner.at_time(trajectory, t_curr)
            
            times.append(t_curr)
            pos.append(state.p)
            vel.append(state.v)
            acc.append(state.a)
            jerk.append(state.j)
            
            t_curr += sampling_dt

        return times, pos, vel, acc, jerk

    def _on_calculate(self):
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
            
            trajectory = self.planner.plan(
                start, tgt_p, tgt_v, 
                self.v_max.get(), self.a_max.get(), self.j_max.get()
            )
            
            ts, ps, vs, as_, js = self._generate_plot_data(trajectory)
            
            data = [ps, vs, as_, js]
            for i, line in enumerate(self.lines):
                line.set_data(ts, data[i])
                self.axs[i].relim()
                self.axs[i].autoscale_view()
            self.canvas.draw()
            
            if trajectory:
                real_total_time = sum(n.dt for n in trajectory)
                final_state = self.planner.at_time(trajectory, real_total_time)
                peak_v = max(abs(min(vs)), max(vs)) if vs else 0.0
                
                self.lbl_stats.config(
                    text=f"Time: {real_total_time:.3f}s | End P: {final_state.p:.2f} | End V: {final_state.v:.2f} | Peak V: {peak_v:.2f}"
                )
            else:
                 self.lbl_stats.config(text="No movement required or invalid params.")
                 
        except Exception as e:
            import traceback
            traceback.print_exc()
            messagebox.showerror("Error", str(e))

    def _on_reset(self):
        for line in self.lines: line.set_data([], [])
        for ax in self.axs: ax.relim(); ax.autoscale_view()
        self.canvas.draw()
        self.lbl_stats.config(text="Reset.")

if __name__ == "__main__":
    root = tk.Tk()
    app = TestGUI(root)
    root.mainloop()
