# Test_100.py
import random
import time
from S7RTT import S7RTT, MotionState

def run_stress_test():
    print("============================================================")
    print("                 S7RTT STRESS TEST (100 Runs)               ")
    print("============================================================")
    
    planner = S7RTT()
    
    # Constraints
    v_lim = 10.0
    a_lim = 10.0
    j_lim = 10.0
    
    fail_count = 0
    pass_count = 0
    
    # Random seed for reproducibility (optional)
    # random.seed(42) 
    
    for i in range(1, 10001):
        # 1. Generate random initial state (within +/- 2x the limits)
        # Larger range for P
        p0 = random.uniform(-500, 500) 
        v0 = random.uniform(-2.0 * v_lim, 2.0 * v_lim)
        a0 = random.uniform(-2.0 * a_lim, 2.0 * a_lim)
        
        start_state = MotionState(p0, v0, a0, 0)
        
        # 2. Generate random target state (within limits)
        pt = random.uniform(-500, 500)
        vt = random.uniform(-v_lim, v_lim)
        # Target acceleration defaults to 0; the plan interface does not provide a target_acc parameter
        
        # 3. Run planning
        try:
            segments = planner.plan(start_state, pt, vt, v_lim, a_lim, j_lim)
            
            # 4. Verify results
            # Use S7RTT's built-in integration tool to verify the final state
            final_state = planner.integrate_path(start_state, segments)
            
            # Check errors
            p_err = abs(final_state.p - pt)
            v_err = abs(final_state.v - vt)
            a_err = abs(final_state.a) # Target acceleration should be 0
            
            # Allowed tolerance
            tol = 1e-3
            
            is_success = (p_err < tol) and (v_err < tol) and (a_err < tol)
            
            if is_success:
                pass_count += 1
                print(f"[RUN {i:03d}] PASS | Time:{sum(s.duration for s in segments):.2f}s")
            else:
                fail_count += 1
                print(f"[RUN {i:03d}] FAIL !!!")
                print(f"    Start : P={p0:.2f}, V={v0:.2f}, A={a0:.2f}")
                print(f"    Target: P={pt:.2f}, V={vt:.2f}")
                print(f"    Actual: P={final_state.p:.4f}, V={final_state.v:.4f}, A={final_state.a:.4f}")
                print(f"    Error : dP={p_err:.4f}, dV={v_err:.4f}, dA={a_err:.4f}")
                print("-" * 40)
                
        except Exception as e:
            fail_count += 1
            print(f"[RUN {i:03d}] EXCEPTION: {str(e)}")
            print(f"    Start : P={p0:.2f}, V={v0:.2f}, A={a0:.2f}")

    print("============================================================")
    print(f"SUMMARY: Total {pass_count + fail_count}, Pass {pass_count}, Fail {fail_count}")
    print("============================================================")

if __name__ == "__main__":
    run_stress_test()
