import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from renormalization_np import get_initial_coefficients, rg_step


def identify_phase(J1, J2, M, n_max=10, n_rg_steps=20, b=2, d=2,
                   convergence_threshold=1e-3, history_length=5):
    """
    Identify phase based on coefficient patterns with early stopping when convergence is detected.
    
    Args:
        J1, J2, M: Model parameters
        n_max: Maximum Fourier mode
        n_rg_steps: Maximum number of RG steps
        b: Length rescaling factor
        d: Dimension
        convergence_threshold: Threshold for convergence detection
        history_length: Number of past steps to track for convergence
    
    Returns:
        String indicating the identified phase
    """
    e = 1e-3  # Threshold for phase identification
    
    # Initialize coefficient array
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)
    
    # Initialize history arrays to track coefficient evolution
    history = []
    
    # Perform RG flow with convergence check
    for i in range(n_rg_steps):
        lambda_nm = rg_step(lambda_nm, b, d)
        
        # Extract key coefficients
        coeffs = {
            'coeff_10': lambda_nm[n_max+1, n_max].real,
            'coeff_01': lambda_nm[n_max, n_max+1].real,
            'coeff_11': lambda_nm[n_max+1, n_max+1].real,
            'coeff_20': lambda_nm[n_max+2, n_max].real,
            'coeff_02': lambda_nm[n_max, n_max+2].real,
            'coeff_22': lambda_nm[n_max+2, n_max+2].real
        }
        
        # Add current coefficients to history
        history.append(coeffs)
        
        # Check for convergence if we have enough history
        if len(history) >= history_length:
            # Check if all coefficients are close to zero (D_Phase)
            all_near_zero = all(
                abs(history[-1][key]) < e 
                for key in ['coeff_10', 'coeff_01', 'coeff_11', 'coeff_20', 'coeff_02', 'coeff_22']
            )
            
            if all_near_zero:
                #print(f"k = {i}")
                return "D_Phase"  # Early detection of disordered phase
            
            # Check for convergence in non-zero coefficients
            converged = True
            for key in coeffs:
                # Calculate maximum variation over last history_length steps
                variation = max(abs(history[-j][key] - history[-j-1][key]) 
                               for j in range(1, history_length))
                
                if variation > convergence_threshold:
                    converged = False
                    break
            
            if converged:
                # We've converged, so determine the phase based on final coefficients
                break
    #print(f"k = {i}")
    # Extract coefficients from the final step
    coeff_10 = lambda_nm[n_max+1, n_max].real
    coeff_01 = lambda_nm[n_max, n_max+1].real
    coeff_11 = lambda_nm[n_max+1, n_max+1].real
    coeff_20 = lambda_nm[n_max+2, n_max].real
    coeff_02 = lambda_nm[n_max, n_max+2].real
    coeff_22 = lambda_nm[n_max+2, n_max+2].real
    
    # Phase conditions
    if (abs(coeff_10) < e and abs(coeff_01) < e and  # (1,0) = 0, (0,1) = 0
        abs(coeff_11) < e and abs(coeff_20) < e and  # (1,1) = 0, (2,0) = 0
        abs(coeff_02) < e and abs(coeff_22) < e):    # (0,2) = 0, (2,2) = 0
        return "D_Phase"  # Disordered
    elif (coeff_10 > e and coeff_01 > e and  # (1,0) > 0, (0,1) > 0
          coeff_11 > e and coeff_20 > e and  # (1,1) > 0, (2,0) > 0
          coeff_02 > e and coeff_22 > e):    # (0,2) > 0, (2,2) > 0
        return "A_Phase"
    elif (abs(coeff_10) < e and abs(coeff_01) < e and  # (1,0) = 0, (0,1) = 0
          coeff_11 < -e and coeff_20 > e and         # (1,1) < 0, (2,0) > 0
          coeff_02 > e and coeff_22 > e):            # (0,2) > 0, (2,2) > 0
        return "B_Phase"
    elif (coeff_10 < -e and coeff_01 < -e and  # (1,0) < 0, (0,1) < 0
          coeff_11 > e and coeff_20 > e and   # (1,1) > 0, (2,0) > 0
          coeff_02 > e and coeff_22 > e):     # (0,2) > 0, (2,2) > 0
        return "C_Phase"
    else:
        return "X_Phase"

def identify_phase_0(J1, J2, M, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Identify phase based on coefficient patterns
    """
    e = 1e-3
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)
           
    # Perform RG flow
    for _ in range(n_rg_steps):
        lambda_nm = rg_step(lambda_nm, b, d)
   
    # Extract normalized coefficients
    coeff_10 = lambda_nm[n_max+1, n_max].real
    coeff_01 = lambda_nm[n_max, n_max+1].real
    coeff_11 = lambda_nm[n_max+1, n_max+1].real
    coeff_20 = lambda_nm[n_max+2, n_max].real
    coeff_02 = lambda_nm[n_max, n_max+2].real
    coeff_22 = lambda_nm[n_max+2, n_max+2].real
   
    # Phase conditions
    if (abs(coeff_10) < e and # (1,0) = 0
        abs(coeff_01) < e and # (0,1) = 0
        abs(coeff_11) < e and # (1,1) = 0
        abs(coeff_20) < e and # (2,0) = 0
        abs(coeff_02) < e and # (0,2) = 0
        abs(coeff_22) < e):   # (2,2) = 0
        return "D_Phase" # Disordered

    elif (coeff_10 > e and # (1,0) > 0
          coeff_01 > e and # (0,1) > 0
          coeff_11 > e and # (1,1) > 0
          coeff_20 > e and # (2,0) > 0
          coeff_02 > e and # (0,2) > 0
          coeff_22 > e):   # (2,2) > 0
        return "A_Phase" # 

    elif (abs(coeff_10) < e and # (1,0) = 0
          abs(coeff_01) < e and # (0,1) = 0
          coeff_11 < -e and     # (1,1) < 0
          coeff_20 > e and      # (2,0) > 0
          coeff_02 > e and      # (0,2) > 0
          coeff_22 > e):        # (2,2) > 0
        return "B_Phase" # 

    elif (coeff_10 < -e and # (1,0) < 0
          coeff_01 < -e and # (0,1) < 0
          coeff_11 > e and  # (1,1) > 0
          coeff_20 > e and  # (2,0) > 0
          coeff_02 > e and  # (0,2) > 0
          coeff_22 > e):    # (2,2) > 0
        return "C_Phase" # 

    else:
        return "X_Phase"

def generate_phase_diagram(J_values, M_values, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Generate phase diagram by scanning J, M parameter space for J = J1 = J2
    
    The marker size is dynamically adjusted based on the grid dimensions.
    """
    D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase = [],[],[],[],[],[]

    for i, J in enumerate(tqdm(J_values)):
        J1 = J2 = J
        for j, M in enumerate(M_values):
            #phase = identify_phase0(J1, J2, M, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d)
            phase = identify_phase(J1, J2, M, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d,
                                   convergence_threshold=1e-3, history_length=5)

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

    # Calculate appropriate marker size based on grid dimensions
    # We want marker size to be inversely proportional to the square root of grid size
    grid_size = len(J_values) * len(M_values)
    base_marker_size = 25  # Base size for a 10x10 grid
    reference_grid_size = 100
    ms = base_marker_size * np.sqrt(reference_grid_size / grid_size)
    
    # Ensure marker size is within reasonable bounds
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

def generate_phase_diagram_inverse(inv_J_values, inv_M_values, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Generate phase diagram by scanning 1/J, 1/M parameter space
    
    The marker size is dynamically adjusted based on the grid dimensions.
    """
    D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase = [],[],[],[],[],[]

    for i, inv_J in enumerate(tqdm(inv_J_values)):
        for j, inv_M in enumerate(inv_M_values):
            J = 1 / inv_J if inv_J != 0 else float('inf')
            J1 = J2 = J
            M = 1 / inv_M if inv_M != 0 else float('inf')
            phase = identify_phase(J1, J2, M, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d)

            if phase == "D_Phase":
                D_Phase.append([inv_J, inv_M])
            elif phase == "A_Phase":
                A_Phase.append([inv_J, inv_M])
            elif phase == "B_Phase":
                B_Phase.append([inv_J, inv_M])
            elif phase == "C_Phase":
                C_Phase.append([inv_J, inv_M])
            elif phase == "X_Phase":
                X_Phase.append([inv_J, inv_M])
            else:
                U_Phase.append([inv_J, inv_M])

    # Calculate appropriate marker size based on grid dimensions
    # We want marker size to be inversely proportional to the square root of grid size
    grid_size = len(inv_J_values) * len(inv_M_values)
    base_marker_size = 25  # Base size for a 10x10 grid
    reference_grid_size = 100
    ms = base_marker_size * np.sqrt(reference_grid_size / grid_size)
    
    # Ensure marker size is within reasonable bounds
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
    
    ax.set_xlabel("1/J")
    ax.set_ylabel("1/M")
    ax.tick_params(axis="both", direction="in", width=2, length=4)
    fig.tight_layout()
    plt.show()
                
    return D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase