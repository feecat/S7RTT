# ==============================================================================
# File Name:    Test_GUI_Online.py
# Description:  Online S-Curve Trajectory Generator Test with Real-time Plotting
# ==============================================================================

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import collections
import time

# 导入你提供的库
from S7RTT import S7RTT, MotionState

class OnlineTestGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("S7RTT Online Real-time Planner (Optimized)")
        self.root.geometry("1200x900")
        
        # --- 核心参数 ---
        self.SIM_DT = 0.02  # 模拟控制周期 20ms (0.02s)
        self.PLOT_WINDOW = 5.0 
        self.is_running = False
        
        # --- 规划器与状态 ---
        self.planner = S7RTT()
        
        # 当前机器人的状态
        self.current_state = MotionState(0, 0, 0, 0, 0)
        
        # --- 【新增】缓存变量用于校验 ---
        # 用于存储上一次周期的输入，判断是否需要重新规划
        self.last_inputs = None 
        # 当前正在执行的轨迹对象
        self.active_trajectory = None 
        # 当前轨迹已经执行的时间
        self.traj_exec_time = 0.0
        
        # 历史数据用于绘图
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
        # 左侧控制面板
        control_panel = ttk.Frame(self.root, width=250)
        control_panel.pack(side=tk.LEFT, fill=tk.Y, padx=10, pady=10)
        
        # 1. 控制按钮
        btn_frame = ttk.LabelFrame(control_panel, text="Simulation Control")
        btn_frame.pack(fill=tk.X, pady=5)
        
        self.btn_start = ttk.Button(btn_frame, text="START Online Loop", command=self._start_simulation)
        self.btn_start.pack(fill=tk.X, padx=5, pady=5)
        
        self.btn_stop = ttk.Button(btn_frame, text="STOP", command=self._stop_simulation, state=tk.DISABLED)
        self.btn_stop.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(btn_frame, text="Reset State to Zero", command=self._reset_state).pack(fill=tk.X, padx=5, pady=5)

        # 2. 动态目标
        target_frame = ttk.LabelFrame(control_panel, text="Dynamic Target")
        target_frame.pack(fill=tk.X, pady=10)
        
        ttk.Label(target_frame, text="Target Position:").pack(anchor="w", padx=5)
        self.var_target_p = tk.DoubleVar(value=100.0)
        scale_p = tk.Scale(target_frame, variable=self.var_target_p, from_=-200, to=200, orient=tk.HORIZONTAL)
        scale_p.pack(fill=tk.X, padx=5)
        
        ttk.Label(target_frame, text="Target Velocity:").pack(anchor="w", padx=5)
        self.var_target_v = tk.DoubleVar(value=0.0)
        scale_v = tk.Scale(target_frame, variable=self.var_target_v, from_=-20, to=20, orient=tk.HORIZONTAL)
        scale_v.pack(fill=tk.X, padx=5)

        # 3. 约束参数
        const_frame = ttk.LabelFrame(control_panel, text="Constraints")
        const_frame.pack(fill=tk.X, pady=10)
        
        self.v_max = self._create_input(const_frame, "Max Vel:", 1000.0)
        self.a_max = self._create_input(const_frame, "Max Acc:", 10000.0)
        self.j_max = self._create_input(const_frame, "Max Jerk:", 10000.0)

        # 4. 状态显示
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
        self.fig.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.05, hspace=0.2)
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
        # 启动时重置一下缓存逻辑，保证第一次必定规划
        self.last_inputs = None 
        self._run_cycle()

    def _stop_simulation(self):
        self.is_running = False
        self.btn_start.config(state=tk.NORMAL)
        self.btn_stop.config(state=tk.DISABLED)

    def _reset_state(self):
        self.current_state = MotionState(0, 0, 0, 0, 0)
        # 清除绘图数据
        self.history_t.clear()
        self.history_p.clear()
        self.history_v.clear()
        self.history_a.clear()
        self.history_j.clear()
        self.global_time = 0.0
        # --- 【新增】清除规划缓存 ---
        self.last_inputs = None
        self.active_trajectory = None
        self.traj_exec_time = 0.0
        
        self._update_plot()
        self._update_labels()

    def _run_cycle(self):
        """
        带输入校验优化的实时控制循环
        """
        if not self.is_running:
            return

        start_time = time.time()

        try:
            # 1. 获取当前的输入参数
            tgt_p = self.var_target_p.get()
            tgt_v = self.var_target_v.get()
            v_max = self.v_max.get()
            a_max = self.a_max.get()
            j_max = self.j_max.get()
            
            # 将所有影响规划的参数打包成元组
            current_inputs = (tgt_p, tgt_v, v_max, a_max, j_max)

            # 2. 【优化校验】检查输入是否变化
            # 如果 last_inputs 为空（第一次运行）或者 输入元组不相等
            inputs_changed = (self.last_inputs is None) or (current_inputs != self.last_inputs)

            if inputs_changed:
                # --- 变化了：重新规划 ---
                # print(f"Target/Constraints Changed at {self.global_time:.2f}, Re-planning...") # 调试用
                
                # 调用规划器生成新的轨迹对象
                # 注意：起点永远是当前这一刻的机器人状态 self.current_state
                self.active_trajectory = self.planner.plan(
                    self.current_state, 
                    tgt_p, 
                    tgt_v, 
                    v_max, a_max, j_max
                )
                
                # 重置该轨迹的执行时间计时器
                self.traj_exec_time = 0.0
                # 更新缓存的输入
                self.last_inputs = current_inputs
            
            else:
                # --- 没变化：不进行计算，仅延续之前的轨迹 ---
                pass

            # 3. 【步进与保护】计算下一时刻状态
            # 不管是否重新规划，我们都需要向前推进时间
            
            if not self.active_trajectory:
                # 异常情况或已到达完全静止：保持原地
                # 简单的物理阻尼模拟，防止数值飘移
                next_state = self.current_state.copy()
                if abs(next_state.v) < 0.001:
                    next_state.v = 0.0; next_state.a = 0.0; next_state.j = 0.0
                else:
                    next_state.p += next_state.v * self.SIM_DT
            else:
                # 累计轨迹执行时间
                self.traj_exec_time += self.SIM_DT
                
                # 【关键点】
                # 之前的代码是 at_time(traj, SIM_DT)，因为那是每次都算新轨迹，新轨迹t=0就是当前位置。
                # 现在我们复用轨迹，所以必须查询 "从轨迹开始到现在经过的时间" 的状态。
                next_state = self.planner.at_time(self.active_trajectory, self.traj_exec_time)

            # 4. 更新系统状态
            self.current_state = next_state
            self.global_time += self.SIM_DT

            # 5. 记录数据
            self.history_t.append(self.global_time)
            self.history_p.append(self.current_state.p)
            self.history_v.append(self.current_state.v)
            self.history_a.append(self.current_state.a)
            self.history_j.append(self.current_state.j)

            # 6. 更新文本
            self._update_labels()

        except Exception as e:
            print(f"Error: {e}")

        finally:
            # 7. 安排下一次循环
            elapsed = time.time() - start_time
            wait_ms = int((self.SIM_DT - elapsed) * 1000)
            if wait_ms < 1: wait_ms = 1
            
            # 绘图降频 (每0.1秒画一次)
            if len(self.history_t) > 0 and int(self.global_time / self.SIM_DT) % 5 == 0:
                self._update_plot()
            
            if self.is_running:
                self.root.after(wait_ms, self._run_cycle)

    def _update_labels(self):
        s = self.current_state
        # 增加显示当前是否正在重新计算的提示（可选）
        status_txt = "Hold" if (self.last_inputs and self.var_target_p.get() == self.last_inputs[0]) else "Plan"
        
        self.lbl_cur_p.config(text=f"P: {s.p:.2f} (Tgt: {self.var_target_p.get():.0f})")
        self.lbl_cur_v.config(text=f"V: {s.v:.2f}")
        self.lbl_cur_a.config(text=f"A: {s.a:.2f}")

    def _update_plot(self):
        if not self.history_t:
            return
        t_data = list(self.history_t)
        data_sources = [
            list(self.history_p), list(self.history_v),
            list(self.history_a), list(self.history_j)
        ]
        min_t, max_t = t_data[0], t_data[-1]
        if max_t - min_t < self.PLOT_WINDOW: min_t = max(0, max_t - self.PLOT_WINDOW)
            
        for i, line in enumerate(self.lines):
            line.set_data(t_data, data_sources[i])
            self.axs[i].set_xlim(min_t, max_t + 0.1)
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