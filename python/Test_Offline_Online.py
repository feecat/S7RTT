import sys
import random
import math

# Attempt to import the S7RTT module
try:
    from S7RTT import S7RTT, MotionState
except ImportError:
    print("Error: S7RTT.py not found")
    sys.exit(1)

def has_cruise_segment(nodes, epsilon=1e-5):
    """
    Checks if the trajectory contains a cruise segment.
    A cruise segment is defined by a=0, j=0, and a positive duration.
    """
    for n in nodes:
        # Strict check: acceleration and jerk must be near zero
        if abs(n.a) < 1e-2 and abs(n.j) < 1e-2 and n.dt > epsilon:
            return True
    return False

def get_total_duration(nodes):
    """Calculates the total duration of the trajectory nodes."""
    return sum(n.dt for n in nodes)

def run_stress_test():
    # 1. Parameter Settings
    V_MAX = 1000.0
    A_MAX = 10000.0
    J_MAX = 100000.0
    
    ITERATIONS = 10000
    SIM_DT = 0.001  # 1ms simulation step
    
    # Time tolerance threshold in seconds.
    # It is considered a bug if the online replanning remaining time differs
    # from the offline plan's remaining time by more than this value.
    TIME_TOLERANCE = 0.002  # Allow 2ms prediction error
    
    planner = S7RTT()
    
    print(f"=== S7RTT Stress Test v3.0 (Smart Jitter Filtering) ===")
    print(f"Params: V={V_MAX}, A={A_MAX}, J={J_MAX}")
    print(f"Criteria: 1. No strategy mutation (Cruise check)")
    print(f"          2. Consistent remaining time (Diff < {TIME_TOLERANCE*1000}ms)")
    print("-" * 60)

    bug_count = 0
    
    for i in range(ITERATIONS):
        # 2. Random Generation (Physically feasible values)
        p_start = random.uniform(-500, 500)
        v_start = random.uniform(-V_MAX * 0.9, V_MAX * 0.9)
        a_start = random.uniform(-A_MAX * 0.9, A_MAX * 0.9)
        
        p_target = random.uniform(-500, 500)
        v_target = random.uniform(-V_MAX * 0.9, V_MAX * 0.9)
        
        start_state = MotionState(0, p_start, v_start, a_start, 0)
        
        # 3. Offline Planning (Baseline)
        offline_traj = planner.plan(start_state, p_target, v_target, V_MAX, A_MAX, J_MAX)
        
        if not offline_traj: continue

        total_time_offline = get_total_duration(offline_traj)
        offline_has_cruise = has_cruise_segment(offline_traj)
        
        # Display progress
        if i % 10 == 0:
            sys.stdout.write(f"\rProgress: {i}/{ITERATIONS} | Bugs Found: {bug_count}")
            sys.stdout.flush()

        # 4. Online Simulation Verification
        # For efficiency, sample 20 checkpoints rather than every tick
        check_points = [t * total_time_offline / 20.0 for t in range(1, 20)]
        
        for t_curr in check_points:
            # Skip checking if we are too close to the end
            if t_curr >= total_time_offline - SIM_DT * 5: 
                continue

            # A. Get expected state from offline trajectory
            current_sim_state = planner.at_time(offline_traj, t_curr)
            
            # B. Online Re-planning
            online_traj = planner.plan(current_sim_state, p_target, v_target, V_MAX, A_MAX, J_MAX)
            
            if not online_traj:
                print(f"\n\n[Iter {i}] Online planning returned None! (t={t_curr:.4f})")
                bug_count += 1
                break
            
            # --- Smart Verification Logic ---
            
            # 1. Strategy Consistency Check (Solver Failure Detection)
            # If offline plan had no cruise segment (Bang-Bang), but online plan
            # introduces one, the time-optimal solver might be failing.
            online_has_cruise = has_cruise_segment(online_traj)
            
            solver_failed = False
            # Exception: If velocity is already near max, a cruise phase is physically mandatory.
            # We focus on cases where it should be sprinting but suddenly decides to cruise.
            if (not offline_has_cruise) and online_has_cruise:
                if abs(current_sim_state.v) < V_MAX * 0.99 and abs(v_target) < V_MAX * 0.99:
                    solver_failed = True

            # 2. Remaining Time Consistency Check
            rem_time_offline = total_time_offline - t_curr
            rem_time_online = get_total_duration(online_traj)
            
            time_diff = abs(rem_time_offline - rem_time_online)

            # --- Anomaly Detection ---
            
            # Only flag as a bug if the strategy collapses or time deviation is significant.
            # Minor micro-jitters in acceleration are ignored.
            if solver_failed or time_diff > TIME_TOLERANCE:
                
                print(f"\n\n=== Real Discrepancy Detected (Iter {i}) ===")
                print(f"Time t: {t_curr:.4f} / {total_time_offline:.4f} s")
                print(f"Input State: {current_sim_state}")
                print(f"Target: P={p_target:.2f}, V={v_target:.2f}")
                
                if solver_failed:
                    print(">>> [Severity] Solver Strategy Collapse (Crash)")
                    print("    Offline was Bang-Bang (No Cruise), Online degenerated into Cruise mode.")
                    print("    This causes unexpected zero acceleration.")
                
                if time_diff > TIME_TOLERANCE:
                    print(f">>> [Severity] Large Time Prediction Deviation")
                    print(f"    Offline Rem Time: {rem_time_offline:.5f} s")
                    print(f"    Online Pred Time: {rem_time_online:.5f} s")
                    print(f"    Diff: {time_diff*1000:.2f} ms (Threshold {TIME_TOLERANCE*1000} ms)")

                bug_count += 1
                
                # Print reproduction code
                print(f"--- Debug Info ---")
                print(f"start = MotionState(0, {current_sim_state.p}, {current_sim_state.v}, {current_sim_state.a}, {current_sim_state.j})")
                print(f"target_p = {p_target}; target_v = {v_target}")
                print("-" * 30)
                break # Break inner loop

    print(f"\n\nTest finished. Real logic discrepancies found: {bug_count}")

if __name__ == "__main__":
    run_stress_test()