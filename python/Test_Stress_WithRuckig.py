import random
import sys
import traceback
import time
import math


SIM_TOLERANCE = 1e-3 

REPORT_RUCKIG_ISSUES = False 

try:
    from S7RTT import S7RTT, MotionState
    S7RTT_AVAILABLE = True
except ImportError:
    print("Error: Cannot find S7RTT.py.")
    sys.exit(1)

try:
    from ruckig import Ruckig, InputParameter, Trajectory
    RUCKIG_AVAILABLE = True
except ImportError:
    print("Error: Ruckig not install. try 'pip install ruckig'")
    sys.exit(1)

def simulate_s7_end_state(start_s: MotionState, nodes):
    p, v, a = start_s.p, start_s.v, start_s.a
    
    for node in nodes:
        dt = node.dt
        j = node.j
        if dt <= 0: continue
        
        dt2 = dt * dt
        dt3 = dt2 * dt
        
        p += v * dt + 0.5 * a * dt2 + (1.0/6.0) * j * dt3
        v += a * dt + 0.5 * j * dt2
        a += j * dt
        
    return p, v, a

def run_stress_test_with_sim(num_tests=10000):
    print(f"Start test...")
    print(f"SIM_TOLERANCE: {SIM_TOLERANCE}")
    print("-" * 80)

    planner_s7 = S7RTT()
    
    cnt_total_run = 0
    
    cnt_s7_fail_plan  = 0
    cnt_s7_fail_sim   = 0
    
    cnt_rk_fail_plan  = 0
    cnt_rk_fail_sim   = 0
    
    cnt_s7_faster     = 0
    cnt_mismatches    = 0

    GLOBAL_J_LIMIT = 1000000.0
    GLOBAL_A_LIMIT = 100000.0
    GLOBAL_V_LIMIT = 10000.0

    for i in range(1, num_tests + 1):
        cnt_total_run += 1
        
        j_max = random.uniform(0.1, GLOBAL_J_LIMIT)
        
        a_upper_bound = min(GLOBAL_A_LIMIT, j_max * 10.0)
        a_max = random.uniform(0.1, a_upper_bound)
        
        v_upper_bound = min(GLOBAL_V_LIMIT, a_max * 10.0)
        v_max = random.uniform(0.1, v_upper_bound)
        
        p_limit = 5000.0
        p_start = random.uniform(-p_limit, p_limit)
        p_target = random.uniform(-p_limit, p_limit)
        
        v_start = random.uniform(-v_max, v_max)
        v_target = random.uniform(-v_max, v_max)
        
        a_start = random.uniform(-a_max, a_max)
        a_target = 0.0

        s7_ok = False
        s7_dur = 0.0
        s7_error_info = ""
        s7_sim_end_p = 0.0
        
        try:
            start_obj = MotionState(0.0, p_start, v_start, a_start, 0.0)
            traj = planner_s7.plan(start_obj, p_target, v_target, v_max, a_max, j_max)
            
            if traj and len(traj) > 0:
                s7_dur = sum(node.dt for node in traj)
                s7_sim_end_p, _, _ = simulate_s7_end_state(start_obj, traj)
                
                if abs(s7_sim_end_p - p_target) > SIM_TOLERANCE:
                    s7_ok = False
                    cnt_s7_fail_sim += 1
                    s7_error_info = f"Sim Accuracy Failed. Err={abs(s7_sim_end_p - p_target):.4f}"
                else:
                    s7_ok = True
            else:
                if abs(p_start - p_target) < 1e-3 and \
                   abs(v_start - v_target) < 1e-3:
                    s7_dur = 0.0
                    s7_ok = True
                else:
                    s7_ok = False
                    cnt_s7_fail_plan += 1
                    s7_error_info = "Result is empty list []"
        except Exception:
            s7_ok = False
            cnt_s7_fail_plan += 1
            s7_error_info = traceback.format_exc()

        rk_ok = False
        rk_dur = 0.0
        rk_sim_end_p = 0.0
        
        try:
            otg = Ruckig(1)
            inp = InputParameter(1)
            inp.current_position = [p_start]; inp.current_velocity = [v_start]; inp.current_acceleration = [a_start]
            inp.target_position = [p_target]; inp.target_velocity = [v_target]; inp.target_acceleration = [a_target]
            inp.max_velocity = [v_max]; inp.max_acceleration = [a_max]; inp.max_jerk = [j_max]
            
            res_traj = Trajectory(1)
            result = otg.calculate(inp, res_traj)
            
            if result.value >= 0:
                rk_dur = res_traj.duration
                final_state = res_traj.at_time(rk_dur)
                rk_sim_end_p = final_state[0][0] 
                
                if abs(rk_sim_end_p - p_target) > SIM_TOLERANCE:
                    rk_ok = False
                    cnt_rk_fail_sim += 1 
                else:
                    rk_ok = True
            else:
                rk_ok = False
                cnt_rk_fail_plan += 1
        except:
            rk_ok = False
            cnt_rk_fail_plan += 1

        should_print = False
        reason = ""

        if not s7_ok and "Sim Accuracy" in s7_error_info:
            should_print = True
            reason = f"S7 Accuracy Fail (Err={abs(s7_sim_end_p - p_target):.3f})"
            if not rk_ok:
                reason += " [Ruckig also failed]"

        elif not rk_ok:
            if s7_ok and REPORT_RUCKIG_ISSUES:
                should_print = True
                reason = "S7 Success, RK Failed/Inaccurate"

        elif not s7_ok:
            should_print = True
            reason = "S7 Plan Fail" 

        else:
            if s7_dur < rk_dur - 1e-6: 
                cnt_s7_faster += 1
            else:
                diff = abs(s7_dur - rk_dur)
                tol = 0.002 + 0.005 * max(s7_dur, rk_dur)
                
                if diff > tol:
                    should_print = True
                    reason = f"Duration Mismatch (S7 Slower by {diff:.3f}s)"

        if should_print:
            cnt_mismatches += 1
            print("\n" + "="*80)
            print(f"ISSUE FOUND at Test #{i}")
            print("="*80)
            print(f"Reason: {reason}")
            
            if not s7_ok and "Sim Accuracy" not in s7_error_info:
                if len(s7_error_info) > 5: 
                    print("-" * 40)
                    print("S7RTT TRACEBACK:")
                    print(s7_error_info.strip())
                    print("-" * 40)

            print("\nPARAMETERS:")
            print(f"  v_max   = {v_max:.3f}, a_max={a_max:.3f}, j_max={j_max:.3f}")
            print(f"  start   = (p={p_start:.3f}, v={v_start:.3f}, a={a_start:.3f})")
            print(f"  target  = (p={p_target:.3f}, v={v_target:.3f}, a=0.000)")
            
            print("-" * 40)
            print("RESULTS (Verified):")
            
            s7_res_str = f"{s7_dur:.3f}s" if s7_ok else f"FAIL" 
            rk_res_str = f"{rk_dur:.3f}s" if rk_ok else f"FAIL"
            
            if s7_dur > 0: s7_res_str += f" (EndErr={s7_sim_end_p-p_target:.1e})"
            if rk_dur > 0: rk_res_str += f" (EndErr={rk_sim_end_p-p_target:.1e})"
            
            print(f"  S7RTT  : {s7_res_str}")
            print(f"  Ruckig : {rk_res_str}")
            print("="*80 + "\n")

        if i % 500 == 0:
            print(f"Processing... {i}/{num_tests} cases. (Mismatches: {cnt_mismatches})")

    print("\n" + "#" * 80)
    print("FINAL STATISTICS (With Simulation Verification)")
    print("#" * 80)
    print(f"Total Tests         : {cnt_total_run}")
    print("-" * 40)
    print(f"S7RTT Failures:")
    print(f"  - Plan Crash/Empty: {cnt_s7_fail_plan}")
    print(f"  - Accuracy check  : {cnt_s7_fail_sim}")
    print("-" * 40)
    print(f"Ruckig Failures:")
    print(f"  - Plan Error      : {cnt_rk_fail_plan}")
    print(f"  - Accuracy check  : {cnt_rk_fail_sim}")
    print("-" * 40)
    print(f"Comparisons:")
    print(f"  - S7 Faster       : {cnt_s7_faster} (Optimizations)")
    print(f"  - Mismatches      : {cnt_mismatches} (S7 Slower/Failed)")
    print("#" * 80)

if __name__ == "__main__":
    try:
        run_stress_test_with_sim(10000) 
    except KeyboardInterrupt:
        print("\nTesting stopped by user.")