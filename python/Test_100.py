# Test_100.py
import random
from S7RTT import S7RTT, MotionState

def run_stress_test():
    print("============================================================")
    print("                 S7RTT STRESS TEST (100 Runs)             ")
    print("============================================================")
    
    planner = S7RTT()
    
    # Constraints
    v_lim = 10.0
    a_lim = 10.0
    j_lim = 10.0
    
    fail_count = 0
    pass_count = 0
    
    # random.seed(42) 
    
    for i in range(1, 101):
        # 1. Generate random initial state
        # Note: MotionState now requires dt. For start state, dt=0.
        p0 = random.uniform(-500, 500) 
        v0 = random.uniform(-2.0 * v_lim, 2.0 * v_lim)
        a0 = random.uniform(-2.0 * a_lim, 2.0 * a_lim)
        
        # Use keyword arguments for clarity with the new structure
        start_state = MotionState(dt=0.0, p=p0, v=v0, a=a0, j=0.0)
        
        # 2. Generate random target state
        pt = random.uniform(-500, 500)
        vt = random.uniform(-v_lim, v_lim)
        
        # 3. Run planning
        try:
            # 'trajectory' is now a list of MotionState objects
            trajectory = planner.plan(start_state, pt, vt, v_lim, a_lim, j_lim)
            
            if not trajectory:
                # If for some reason the plan is empty (e.g. invalid limits), count as fail
                raise ValueError("Planner returned empty trajectory.")

            # 4. Verify results
            # Calculate total duration from the MotionState nodes
            total_duration = sum(node.dt for node in trajectory)
            
            # Use the propagator to find the state at the EXACT end of the trajectory
            final_state = planner.at_time(trajectory, total_duration)
            
            # Check errors
            p_err = abs(final_state.p - pt)
            v_err = abs(final_state.v - vt)
            
            # Target acceleration and jerk should be 0 at the end of a standard S-Curve
            # (unless interrupted, but here we expect full completion)
            a_err = abs(final_state.a) 
            
            # Allowed tolerance
            tol = 1e-3
            
            is_success = (p_err < tol) and (v_err < tol) and (a_err < tol)
            
            if is_success:
                pass_count += 1
                print(f"[RUN {i:04d}] PASS | Time:{total_duration:.2f}s")
            else:
                fail_count += 1
                print(f"[RUN {i:04d}] FAIL !!!")
                print(f"    Start : P={p0:.2f}, V={v0:.2f}, A={a0:.2f}")
                print(f"    Target: P={pt:.2f}, V={vt:.2f}")
                print(f"    Actual: P={final_state.p:.4f}, V={final_state.v:.4f}, A={final_state.a:.4f}")
                print(f"    Error : dP={p_err:.4f}, dV={v_err:.4f}, dA={a_err:.4f}")
                print("-" * 40)
                
        except Exception as e:
            fail_count += 1
            print(f"[RUN {i:04d}] EXCEPTION: {str(e)}")
            print(f"    Start : P={p0:.2f}, V={v0:.2f}, A={a0:.2f}")

    print("============================================================")
    print(f"SUMMARY: Total {pass_count + fail_count}, Pass {pass_count}, Fail {fail_count}")
    print("============================================================")

if __name__ == "__main__":
    run_stress_test()
