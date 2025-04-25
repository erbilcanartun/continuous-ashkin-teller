import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpmath import mp, fabs, mpf
from renormalization_mp import get_initial_coefficients, rg_step


def identify_phase(J1, J2, M, n_max=10, n_rg_steps=20, b=2, d=2, history_length=3, dps=30):
    """
    Identify phase based on coefficient patterns with early stopping when convergence is detected.
    Uses mpmath for high precision calculations.

    Args:
        J1, J2, M: Model parameters
        n_max: Maximum Fourier mode
        n_rg_steps: Maximum number of RG steps
        b: Length rescaling factor
        d: Dimension
        history_length: Number of past steps to track for convergence
        dps: Decimal places of precision for mpmath

    Returns:
        String indicating the identified phase
    """
    # Set mpmath decimal precision
    mp.dps = dps

    # Use mpmath precision for the threshold
    e = mp.mpf('1e-3')  # Threshold for both phase identification and convergence

    # Convert inputs to mpmath precision
    J1, J2, M = mp.mpf(str(J1)), mp.mpf(str(J2)), mp.mpf(str(M))

    # Get initial coefficients (already in mpmath format)
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)

    # Define coefficient indices to track
    coefficient_indices = [
        (n_max+1, n_max),    # coeff_10
        (n_max, n_max+1),    # coeff_01
        (n_max+1, n_max+1),  # coeff_11
        (n_max+2, n_max),    # coeff_20
        (n_max, n_max+2),    # coeff_02
        (n_max+2, n_max+2)   # coeff_22
    ]
    
    # Keep track of convergence status for each coefficient
    stable_iterations = {idx: 0 for idx in coefficient_indices}
    
    # Perform RG flow with convergence check
    prev_lambda_nm = None
    
    for i in range(n_rg_steps):
        lambda_nm = rg_step(lambda_nm, b, d)

        # Only start checking conditions after 6 iterations
        if i >= 5:  # This ensures the loop runs at least 6 times (i = 0, 1, 2, 3, 4, 5)
            # Check if the disordered phase is reached
            all_near_zero = True
            for n_idx, m_idx in coefficient_indices:
                coeff_value = lambda_nm[n_idx, m_idx].real
                if mp.fabs(coeff_value) >= e:
                    all_near_zero = False
                    break
                    
            if all_near_zero:
                return "D_Phase"  # Early detection of disordered phase
                
            # Check for convergence if we have a previous state to compare with
            if prev_lambda_nm is not None:
                all_converged = True
                
                for n_idx, m_idx in coefficient_indices:
                    current_val = lambda_nm[n_idx, m_idx].real
                    prev_val = prev_lambda_nm[n_idx, m_idx].real
                    
                    # Calculate difference between iterations using mpmath
                    diff = mp.fabs(current_val - prev_val)
                    
                    if diff < e:
                        stable_iterations[(n_idx, m_idx)] += 1
                    else:
                        stable_iterations[(n_idx, m_idx)] = 0
                        all_converged = False
                
                # Break if all coefficients have been stable for enough iterations
                min_stable_count = min(stable_iterations.values())
                if all_converged and min_stable_count >= history_length-1:
                    #print(f"All converged for J={J1}, M={M} at step {i}")
                    break
                
        # Store current state for next iteration
        prev_lambda_nm = lambda_nm.copy()

    #print(f"RG finalized at step {i}")
    
    # Extract normalized coefficients in mpmath format
    coeff_10 = lambda_nm[n_max+1, n_max].real
    coeff_01 = lambda_nm[n_max, n_max+1].real
    coeff_11 = lambda_nm[n_max+1, n_max+1].real
    coeff_20 = lambda_nm[n_max+2, n_max].real
    coeff_02 = lambda_nm[n_max, n_max+2].real
    coeff_22 = lambda_nm[n_max+2, n_max+2].real

    # Phase conditions using mpmath's fabs
    if (mp.fabs(coeff_10) < e and  # (1,0) = 0
        mp.fabs(coeff_01) < e and  # (0,1) = 0
        mp.fabs(coeff_11) < e and  # (1,1) = 0
        mp.fabs(coeff_20) < e and  # (2,0) = 0
        mp.fabs(coeff_02) < e and  # (0,2) = 0
        mp.fabs(coeff_22) < e):    # (2,2) = 0
        return "D_Phase"  # Disordered
    elif (coeff_10 > e and  # (1,0) > 0
          coeff_01 > e and  # (0,1) > 0
          coeff_11 > e and  # (1,1) > 0
          coeff_20 > e and  # (2,0) > 0
          coeff_02 > e and  # (0,2) > 0
          coeff_22 > e):    # (2,2) > 0
        return "A_Phase"
    elif (mp.fabs(coeff_10) < e and  # (1,0) = 0
          mp.fabs(coeff_01) < e and  # (0,1) = 0
          coeff_11 > e and           # (1,1) > 0
          coeff_20 > e and           # (2,0) > 0
          coeff_02 > e and           # (0,2) > 0
          coeff_22 > e):             # (2,2) > 0
        return "B_Phase"
    elif (coeff_10 < -e and  # (1,0) < 0
          coeff_01 < -e and  # (0,1) < 0
          coeff_11 > e and   # (1,1) > 0
          coeff_20 > e and   # (2,0) > 0
          coeff_02 > e and   # (0,2) > 0
          coeff_22 > e):     # (2,2) > 0
        return "C_Phase"
    elif (mp.fabs(coeff_10) < e and  # (1,0) = 0
          mp.fabs(coeff_01) < e and  # (0,1) = 0
          coeff_11 < -e and          # (1,1) < 0
          coeff_20 > e and           # (2,0) > 0
          coeff_02 > e and           # (0,2) > 0
          coeff_22 > e):             # (2,2) > 0
        return "E_Phase"
    else:
        return "X_Phase"

def identify_phase_all_steps(J1, J2, M, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Identify phase based on coefficient patterns
    """
    # Use mpmath precision for the threshold
    e = mp.mpf('1e-3')

    # Convert inputs to mpmath precision
    J1, J2, M = mp.mpf(str(J1)), mp.mpf(str(J2)), mp.mpf(str(M))

    # Get initial coefficients (already in mpmath format)
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)

    # Perform RG flow
    for _ in range(n_rg_steps):
        lambda_nm = rg_step(lambda_nm, b, d)

    # Extract normalized coefficients in mpmath format (no conversion to float)
    coeff_10 = lambda_nm[n_max+1, n_max].real
    coeff_01 = lambda_nm[n_max, n_max+1].real
    coeff_11 = lambda_nm[n_max+1, n_max+1].real
    coeff_20 = lambda_nm[n_max+2, n_max].real
    coeff_02 = lambda_nm[n_max, n_max+2].real
    coeff_22 = lambda_nm[n_max+2, n_max+2].real

    # Phase conditions using mpmath's fabs
    if (mp.fabs(coeff_10) < e and  # (1,0) = 0
        mp.fabs(coeff_01) < e and  # (0,1) = 0
        mp.fabs(coeff_11) < e and  # (1,1) = 0
        mp.fabs(coeff_20) < e and  # (2,0) = 0
        mp.fabs(coeff_02) < e and  # (0,2) = 0
        mp.fabs(coeff_22) < e):    # (2,2) = 0
        return "D_Phase"  # Disordered

    elif (coeff_10 > e and  # (1,0) > 0
          coeff_01 > e and  # (0,1) > 0
          coeff_11 > e and  # (1,1) > 0
          coeff_20 > e and  # (2,0) > 0
          coeff_02 > e and  # (0,2) > 0
          coeff_22 > e):    # (2,2) > 0
        return "A_Phase"

    elif (mp.fabs(coeff_10) < e and  # (1,0) = 0
          mp.fabs(coeff_01) < e and  # (0,1) = 0
          coeff_11 < -e and          # (1,1) < 0
          coeff_20 > e and           # (2,0) > 0
          coeff_02 > e and           # (0,2) > 0
          coeff_22 > e):             # (2,2) > 0
        return "B_Phase"

    elif (coeff_10 < -e and  # (1,0) < 0
          coeff_01 < -e and  # (0,1) < 0
          coeff_11 > e and   # (1,1) > 0
          coeff_20 > e and   # (2,0) > 0
          coeff_02 > e and   # (0,2) > 0
          coeff_22 > e):     # (2,2) > 0
        return "C_Phase"

    else:
        return "X_Phase"

def generate_phase_diagram(J_values, M_values, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Generate phase diagram by scanning J, M parameter space for J = J1 = J2

    The marker size is dynamically adjusted based on the grid dimensions.
    """
    D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase = [],[],[],[],[],[]

    for i, J in enumerate(tqdm(J_values)):
        #J1 = J2 = J
        for j, M in enumerate(M_values):
            # Convert parameters to mpmath precision for identify_phase
            J_mp = mp.mpf(str(J))
            M_mp = mp.mpf(str(M))

            #phase = identify_phase(J_mp, J_mp, M_mp, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d)
            phase = identify_phase(J_mp, J_mp, M_mp, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d,
                                   convergence_threshold=1e-3, history_length=5, dps=30)

            # Store original float values for plotting
            if phase == "D_Phase":
                D_Phase.append([J, M])
            elif phase == "A_Phase":
                A_Phase.append([J, M])
            elif phase == "B_Phase":
                B_Phase.append([J, M])
            elif phase == "C_Phase":
                C_Phase.append([J, M])
            elif phase == "X_Phase":
                X_Phase.append([J, M])
            else:
                U_Phase.append([J, M])

    # Calculate appropriate marker size for plotting
    grid_size = len(J_values) * len(M_values)
    base_marker_size = 25  # Base size for a 10x10 grid
    reference_grid_size = 100
    ms = base_marker_size * np.sqrt(reference_grid_size / grid_size)
    ms = max(1, min(ms, 12))  # Limit between 1 and 12

    # Plot the results
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6), dpi=100)
    fig.set_facecolor("white")
    plt.rc(group="font", family="Arial", weight="bold", size=10)
    plt.rc(group="lines", linewidth=1)
    plt.rc(group="axes", linewidth=2)

    cdic = {"D_Phase":"grey",
            "A_Phase":"blue",
            "B_Phase":"red",
            "C_Phase":"pink",
            "X_Phase":"black",
            "U_Phase":"yellow"}

    if D_Phase: ax.plot(np.array(D_Phase)[:,0], np.array(D_Phase)[:,1], ls="", marker="s", mfc=cdic["D_Phase"], mec=cdic["D_Phase"], ms=ms, alpha=1)
    if A_Phase: ax.plot(np.array(A_Phase)[:,0], np.array(A_Phase)[:,1], ls="", marker="s", mfc=cdic["A_Phase"], mec=cdic["A_Phase"], ms=ms, alpha=1)
    if B_Phase: ax.plot(np.array(B_Phase)[:,0], np.array(B_Phase)[:,1], ls="", marker="s", mfc=cdic["B_Phase"], mec=cdic["B_Phase"], ms=ms, alpha=1)
    if C_Phase: ax.plot(np.array(C_Phase)[:,0], np.array(C_Phase)[:,1], ls="", marker="s", mfc=cdic["C_Phase"], mec=cdic["C_Phase"], ms=ms, alpha=1)
    if X_Phase: ax.plot(np.array(X_Phase)[:,0], np.array(X_Phase)[:,1], ls="", marker="s", mfc=cdic["X_Phase"], mec=cdic["X_Phase"], ms=ms, alpha=1)
    if U_Phase: ax.plot(np.array(U_Phase)[:,0], np.array(U_Phase)[:,1], ls="", marker="s", mfc=cdic["U_Phase"], mec=cdic["U_Phase"], ms=ms, alpha=1)

    ax.set_xlabel("J")
    ax.set_ylabel("M")
    ax.tick_params(axis="both", direction="in", width=2, length=4)
    fig.tight_layout()
    plt.show()

    return D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase