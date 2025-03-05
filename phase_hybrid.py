import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
from phase_np import identify_phase as identify_phase_np
from phase_mp import identify_phase as identify_phase_mp


def identify_phase_hybrid(J1, J2, M, n_max_np=20, n_max_mp=10, n_rg_steps=30, b=2, d=2, dps=30):
    """
    Hybrid approach to phase identification that uses numpy for speed
    and falls back to mpmath for precision when needed.
    
    Args:
        J1, J2, M: Model parameters
        n_max: Maximum Fourier mode
        n_rg_steps: Number of RG steps
        b: Length rescaling factor
        d: Dimension
    
    Returns:
        String indicating the identified phase
    """
    # First try with the faster numpy implementation
    phase = identify_phase_np(J1, J2, M, n_max_np, n_rg_steps, b, d)
    
    # If the numpy implementation returns "X_Phase" or fails to converge,
    # fall back to the more precise mpmath implementation
    if phase == "X_Phase":
        print(f"Falling back to mpmath for J1={J1}, J2={J2}, M={M}")
        
        # Convert parameters to mpmath precision
        J1_mp = mp.mpf(str(J1))
        J2_mp = mp.mpf(str(J2))
        M_mp = mp.mpf(str(M))
        
        # Use the mpmath implementation for increased precision
        phase = identify_phase_mp(J1_mp, J2_mp, M_mp, n_max_mp, n_rg_steps, b, d, dps=dps)
    
    return phase

def generate_phase_diagram_hybrid(J_values, M_values, n_max_np=20, n_max_mp=10, n_rg_steps=20, b=2, d=2, dps=30, path=None):
    """
    Generate phase diagram by scanning J, M parameter space for J = J1 = J2
    using a hybrid approach (numpy for speed, mpmath for precision when needed).
    
    Args:
        J_values: Array of J values
        M_values: Array of M values
        n_max: Maximum Fourier mode
        n_rg_steps: Number of RG steps
        b: Length rescaling factor
        d: Dimension
    
    Returns:
        Tuple of phase points (D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase)
    """
    D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase = [],[],[],[],[],[]
    
    # Calculate the total number of points for progress tracking
    total_points = len(J_values) * len(M_values)
    processed_points = 0

    for i, J in enumerate(J_values):
        J1 = J2 = J
        for j, M in enumerate(M_values):
            # Update progress counter
            processed_points += 1
            if processed_points % 100 == 0:
                print(f"Progress: {processed_points}/{total_points} points processed ({processed_points/total_points*100:.1f}%)")
            
            # Use the hybrid phase identification function
            phase = identify_phase_hybrid(J1, J2, M, n_max_np=n_max_np, n_max_mp=n_max_mp, n_rg_steps=n_rg_steps, b=b, d=d, dps=dps)

            # Store the result in the appropriate list
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
    if path:
        fig.savefig(path)
    
    # Report statistics
    total = len(D_Phase) + len(A_Phase) + len(B_Phase) + len(C_Phase) + len(X_Phase) + len(U_Phase)
    print(f"Phase distribution:")
    print(f"D_Phase: {len(D_Phase)}/{total} ({len(D_Phase)/total*100:.1f}%)")
    print(f"A_Phase: {len(A_Phase)}/{total} ({len(A_Phase)/total*100:.1f}%)")
    print(f"B_Phase: {len(B_Phase)}/{total} ({len(B_Phase)/total*100:.1f}%)")
    print(f"C_Phase: {len(C_Phase)}/{total} ({len(C_Phase)/total*100:.1f}%)")
    print(f"X_Phase: {len(X_Phase)}/{total} ({len(X_Phase)/total*100:.1f}%)")
    print(f"U_Phase: {len(U_Phase)}/{total} ({len(U_Phase)/total*100:.1f}%)")
                
    return D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase

def generate_phase_diagram_inverse_hybrid(inv_J_values, inv_M_values, n_max=10, n_rg_steps=20, b=2, d=2):
    """
    Generate phase diagram by scanning 1/J, 1/M parameter space
    using a hybrid approach (numpy for speed, mpmath for precision when needed).
    
    Args:
        inv_J_values: Array of 1/J values
        inv_M_values: Array of 1/M values
        n_max: Maximum Fourier mode
        n_rg_steps: Number of RG steps
        b: Length rescaling factor
        d: Dimension
    
    Returns:
        Tuple of phase points (D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase)
    """
    D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase = [],[],[],[],[],[]
    
    # Calculate the total number of points for progress tracking
    total_points = len(inv_J_values) * len(inv_M_values)
    processed_points = 0

    for i, inv_J in enumerate(inv_J_values):
        for j, inv_M in enumerate(inv_M_values):
            # Update progress counter
            processed_points += 1
            if processed_points % 100 == 0:
                print(f"Progress: {processed_points}/{total_points} points processed ({processed_points/total_points*100:.1f}%)")
            
            # Convert from inverse parameters to actual parameters
            if inv_J != 0:
                J = 1.0 / inv_J
            else:
                J = float('inf')
                
            if inv_M != 0:
                M = 1.0 / inv_M
            else:
                M = float('inf')
            
            J1 = J2 = J
            
            # Use the hybrid phase identification function
            phase = identify_phase_hybrid(J1, J2, M, n_max=n_max, n_rg_steps=n_rg_steps, b=b, d=d)

            # Store the result in the appropriate list
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