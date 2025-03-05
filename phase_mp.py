import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from mpmath import mp, fabs, mpf
from renormalization_mp import get_initial_coefficients, rg_step


def identify_phase(J1, J2, M, n_max=10, n_rg_steps=20, b=2, d=2, convergence_threshold=1e-3, history_length=5, dps=30):
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
    # Set mpmath decimal precision
    mp.dps = dps
    
    # Use mpmath precision for the threshold
    e = mp.mpf('1e-3')
    convergence_threshold = mp.mpf(str(convergence_threshold))
    
    # Convert inputs to mpmath precision
    J1, J2, M = mp.mpf(str(J1)), mp.mpf(str(J2)), mp.mpf(str(M))
    
    # Get initial coefficients (already in mpmath format)
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)
    
    # Initialize history list to track coefficient evolution
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
            all_near_zero = True
            for key in ['coeff_10', 'coeff_01', 'coeff_11', 'coeff_20', 'coeff_02', 'coeff_22']:
                if mp.fabs(history[-1][key]) >= e:
                    all_near_zero = False
                    break
            
            if all_near_zero:
                return "D_Phase"  # Early detection of disordered phase
            
            # Check for convergence in non-zero coefficients
            converged = True
            for key in coeffs:
                # Calculate maximum variation over last history_length steps
                max_variation = mp.mpf('0')
                for j in range(1, history_length):
                    variation = mp.fabs(history[-j][key] - history[-j-1][key])
                    if variation > max_variation:
                        max_variation = variation
                
                if max_variation > convergence_threshold:
                    converged = False
                    break
            
            if converged:
                # We've converged, so determine the phase based on final coefficients
                break
    
    # Extract normalized coefficients in mpmath format (no conversion to float)
    coeff_10 = lambda_nm[n_max+1, n_max].real
    coeff_01 = lambda_nm[n_max, n_max+1].real
    coeff_11 = lambda_nm[n_max+1, n_max+1].real
    coeff_20 = lambda_nm[n_max+2, n_max].real
    coeff_02 = lambda_nm[n_max, n_max+2].real
    coeff_22 = lambda_nm[n_max+2, n_max+2].real
    
    # Phase conditions
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
    elif (mp.fabs(coeff_10) < e and # (1,0) = 0
          mp.fabs(coeff_01) < e and # (0,1) = 0
          coeff_11 < -e and         # (1,1) < 0
          coeff_20 > e and          # (2,0) > 0
          coeff_02 > e and          # (0,2) > 0
          coeff_22 > e):            # (2,2) > 0
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

def identify_phase_0(J1, J2, M, n_max=10, n_rg_steps=20, b=2, d=2):
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

    elif (mp.fabs(coeff_10) < e and # (1,0) = 0
          mp.fabs(coeff_01) < e and # (0,1) = 0
          coeff_11 < -e and         # (1,1) < 0
          coeff_20 > e and          # (2,0) > 0
          coeff_02 > e and          # (0,2) > 0
          coeff_22 > e):            # (2,2) > 0
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

def generate_phase_diagram_inverse(inv_J_values, inv_M_values, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Generate phase diagram by scanning 1/J, 1/M parameter space
    
    The marker size is dynamically adjusted based on the grid dimensions.
    """
    D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase = [],[],[],[],[],[]

    for i, inv_J in enumerate(tqdm(inv_J_values)):
        for j, inv_M in enumerate(inv_M_values):
            # Handle division by zero with mpmath
            if inv_J != 0:
                J_mp = mp.mpf('1') / mp.mpf(str(inv_J))
            else:
                # Set to a very large value
                J_mp = mp.mpf('1e100')
                
            if inv_M != 0:
                M_mp = mp.mpf('1') / mp.mpf(str(inv_M))
            else:
                # Set to a very large value
                M_mp = mp.mpf('1e100')
            
            # All calculations with mpmath precision
            phase = identify_phase(J_mp, J_mp, M_mp, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d)

            # Store original float values for plotting
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

    # Calculate appropriate marker size for plotting
    grid_size = len(inv_J_values) * len(inv_M_values)
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
    
    ax.set_xlabel("1/J")
    ax.set_ylabel("1/M")
    ax.tick_params(axis="both", direction="in", width=2, length=4)
    fig.tight_layout()
    plt.show()
                
    return D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase