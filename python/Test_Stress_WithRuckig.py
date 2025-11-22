import random
import sys
import traceback
import time
import math

REPORT_RUCKIG_ISSUES = False
REPORT_S7RTT_ISSUES = True
REPORT_RUCKIG_SLOWS = False
REPORT_S7RTT_SLOWS = False

SIM_TOLERANCE = 1e-3

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

def simulate_s7_end_state(start_s, nodes):
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
    
    cnt_s7_fail_plan = 0
    cnt_s7_fail_sim = 0
    cnt_rk_fail_plan = 0
    cnt_rk_fail_sim = 0
    
    cnt_s7_faster = 0
    cnt_rk_faster = 0
    cnt_draw = 0
    
    cnt_printed_issues = 0

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
        s7_error_reason = ""
        s7_end_p = 0.0
        
        try:
            start_obj = MotionState(0.0, p_start, v_start, a_start, 0.0)
            traj = planner_s7.plan(start_obj, p_target, v_target, v_max, a_max, j_max)
            
            if traj and len(traj) > 0:
                s7_dur = sum(node.dt for node in traj)
                s7_end_p, _, _ = simulate_s7_end_state(start_obj, traj)
                
                if abs(s7_end_p - p_target) > SIM_TOLERANCE:
                    s7_ok = False
                    cnt_s7_fail_sim += 1
                    s7_error_reason = f"Accuracy Fail (Diff={abs(s7_end_p - p_target):.4f})"
                else:
                    s7_ok = True
            else:
                if abs(p_start - p_target) < 1e-3 and abs(v_start - v_target) < 1e-3:
                    s7_dur = 0.0
                    s7_ok = True
                else:
                    s7_ok = False
                    cnt_s7_fail_plan += 1
                    s7_error_reason = "Empty Trajectory"
        except Exception:
            s7_ok = False
            cnt_s7_fail_plan += 1
            s7_error_reason = traceback.format_exc()

        rk_ok = False
        rk_dur = 0.0
        rk_error_reason = ""
        rk_end_p = 0.0
        
        try:
            otg = Ruckig(1)
            inp = InputParameter(1)
            inp.current_position = [p_start]
            inp.current_velocity = [v_start]
            inp.current_acceleration = [a_start]
            inp.target_position = [p_target]
            inp.target_velocity = [v_target]
            inp.target_acceleration = [a_target]
            inp.max_velocity = [v_max]
            inp.max_acceleration = [a_max]
            inp.max_jerk = [j_max]
            
            res_traj = Trajectory(1)
            result = otg.calculate(inp, res_traj)
            
            if result.value >= 0:
                rk_dur = res_traj.duration
                final_state = res_traj.at_time(rk_dur)
                rk_end_p = final_state[0][0] 
                
                if abs(rk_end_p - p_target) > SIM_TOLERANCE:
                    rk_ok = False
                    cnt_rk_fail_sim += 1
                    rk_error_reason = f"Accuracy Fail (Diff={abs(rk_end_p - p_target):.4f})" 
                else:
                    rk_ok = True
            else:
                rk_ok = False
                cnt_rk_fail_plan += 1
                rk_error_reason = f"Plan Fail Code: {result.value}"
        except Exception:
            rk_ok = False
            cnt_rk_fail_plan += 1
            rk_error_reason = traceback.format_exc()

        winner_label = "DRAW"
        if s7_ok and rk_ok:
            diff = s7_dur - rk_dur
            if diff < -1e-6:
                cnt_s7_faster += 1
                winner_label = "S7RTT"
            elif diff > 1e-6:
                cnt_rk_faster += 1
                winner_label = "RUCKIG"
            else:
                cnt_draw += 1

        should_print = False
        report_reason = ""

        if not s7_ok:
            if REPORT_S7RTT_ISSUES:
                should_print = True
                report_reason = f"S7RTT Issue: {s7_error_reason}"
        elif not rk_ok:
            if REPORT_RUCKIG_ISSUES:
                should_print = True
                report_reason = f"Ruckig Issue: {rk_error_reason}"
        else:
            if winner_label == "S7RTT" and REPORT_RUCKIG_SLOWS:
                should_print = True
                report_reason = f"Ruckig Slower (+{s7_dur - rk_dur:.4f}s)"
            elif winner_label == "RUCKIG" and REPORT_S7RTT_SLOWS:
                should_print = True
                report_reason = f"S7RTT Slower (+{rk_dur - s7_dur:.4f}s)"

        if should_print:
            cnt_printed_issues += 1
            print("\n" + "="*80)
            print(f"REPORT at Test #{i} | Reason: {report_reason}")
            print("="*80)
            print(f"PARAMETERS:")
            print(f"  v_max={v_max:.3f}, a_max={a_max:.3f}, j_max={j_max:.3f}")
            print(f"  start=(p={p_start:.3f}, v={v_start:.3f}, a={a_start:.3f})")
            print(f"  target=(p={p_target:.3f}, v={v_target:.3f}, a={a_target:.3f})")
            print("-" * 40)
            print("RESULTS:")
            
            s7_str = f"OK ({s7_dur:.4f}s)" if s7_ok else "FAIL"
            rk_str = f"OK ({rk_dur:.4f}s)" if rk_ok else "FAIL"
            
            print(f"  S7RTT  : {s7_str}")
            if not s7_ok: print(f"           Info: {s7_error_reason}")
            
            print(f"  Ruckig : {rk_str}")
            if not rk_ok: print(f"           Info: {rk_error_reason}")
            
            print("="*80)

        if i % 1000 == 0:
            print(f"Progress: {i}/{num_tests} done.")

    print("\n" + "#" * 80)
    print("FINAL STATISTICS")
    print("#" * 80)
    print(f"Total Tests: {cnt_total_run}")
    print("-" * 40)
    print(f"{'CATEGORY':<20} | {'S7RTT':<10} | {'RUCKIG':<10}")
    print("-" * 40)
    print(f"{'Plan Failures':<20} | {cnt_s7_fail_plan:<10} | {cnt_rk_fail_plan:<10}")
    print(f"{'Sim Acc Failures':<20} | {cnt_s7_fail_sim:<10} | {cnt_rk_fail_sim:<10}")
    print("-" * 40)
    print(f"{'Faster Count':<20} | {cnt_s7_faster:<10} | {cnt_rk_faster:<10}")
    print(f"{'Draws':<20} | {cnt_draw:<10} | {cnt_draw:<10}")
    print("#" * 80)

if __name__ == "__main__":
    try:
        run_stress_test_with_sim(10000) 
    except KeyboardInterrupt:
        print("\nTesting stopped by user.")