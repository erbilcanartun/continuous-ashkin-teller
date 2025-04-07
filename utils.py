import pickle
import numpy as np
import matplotlib.pyplot as plt

def save_phase_diagram(filename, J_values, M_values, D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase):
    """
    Save phase diagram data to a pickle file, including parameter values.

    Args:
        filename (str): Name of the file to save the data to
        J_values (array): Array of J parameter values
        M_values (array): Array of M parameter values
        D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase: Lists of points for each phase
    """
    # Create a dictionary to hold all phase data
    phase_data = {
        'J_values': J_values,
        'M_values': M_values,
        'D_Phase': D_Phase,
        'A_Phase': A_Phase,
        'B_Phase': B_Phase,
        'C_Phase': C_Phase,
        'E_Phase': E_Phase,
        'X_Phase': X_Phase,
        'U_Phase': U_Phase
    }

    # Ensure the filename has a .pkl extension
    if not filename.endswith('.pkl'):
        filename += '.pkl'

    # Save the data to file
    with open(filename, 'wb') as f:
        pickle.dump(phase_data, f)

    print(f"Phase diagram data saved to {filename}")


def load_phase_diagram(filename):
    """
    Load phase diagram data from a pickle file.

    Args:
        filename (str): Name of the file to load the data from

    Returns:
        tuple: (J_values, M_values, D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase)
    """
    # Ensure the filename has a .pkl extension
    if not filename.endswith('.pkl'):
        filename += '.pkl'

    # Load the data from file
    with open(filename, 'rb') as f:
        phase_data = pickle.load(f)

    # Extract parameter values
    J_values = phase_data['J_values']
    M_values = phase_data['M_values']

    # Extract individual phase lists
    D_Phase = phase_data['D_Phase']
    A_Phase = phase_data['A_Phase']
    B_Phase = phase_data['B_Phase']
    C_Phase = phase_data['C_Phase']
    E_Phase = phase_data['E_Phase']
    X_Phase = phase_data['X_Phase']
    U_Phase = phase_data['U_Phase']

    print(f"Phase diagram data loaded from {filename}")

    return J_values, M_values, D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase


def plot_phase_diagram(J_values, M_values, D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase):
    """
    Plot the phase diagram using saved data.

    Args:
        J_values (array): Array of J parameter values
        M_values (array): Array of M parameter values
        D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase: Lists of points for each phase
    """
    # Calculate appropriate marker size
    grid_size = len(J_values) * len(M_values)
    base_marker_size = 25  # Base size for a 10x10 grid
    reference_grid_size = 100
    ms = base_marker_size * np.sqrt(reference_grid_size / grid_size)
    ms = max(1, min(ms, 12))  # Limit between 1 and 12

    # Plot the results
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 8), dpi=100)
    fig.set_facecolor("white")
    plt.rc(group="font", family="Arial", weight="bold", size=10)
    plt.rc(group="lines", linewidth=1)
    plt.rc(group="axes", linewidth=2)

    cdic = {"D_Phase":"grey",
            "A_Phase":"lightblue", # F
            "B_Phase":"green",     # A
            "C_Phase":"red",       # Fsigma
            "E_Phase":"orange",    # Asigma
            "X_Phase":"black",
            "U_Phase":"yellow"}

    # Add a legend with colored squares for each phase
    handles = []
    labels = []

    # Plot each phase and add to legend if the phase exists
    if D_Phase:
        ax.plot(np.array(D_Phase)[:,1], np.array(D_Phase)[:,0], ls="", marker="s", mfc=cdic["D_Phase"], mec=cdic["D_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["D_Phase"], markersize=10))
        labels.append("D_Phase")
    if A_Phase:
        ax.plot(np.array(A_Phase)[:,1], np.array(A_Phase)[:,0], ls="", marker="s", mfc=cdic["A_Phase"], mec=cdic["A_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["A_Phase"], markersize=10))
        labels.append("A_Phase")
    if B_Phase:
        ax.plot(np.array(B_Phase)[:,1], np.array(B_Phase)[:,0], ls="", marker="s", mfc=cdic["B_Phase"], mec=cdic["B_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["B_Phase"], markersize=10))
        labels.append("B_Phase")
    if C_Phase:
        ax.plot(np.array(C_Phase)[:,1], np.array(C_Phase)[:,0], ls="", marker="s", mfc=cdic["C_Phase"], mec=cdic["C_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["C_Phase"], markersize=10))
        labels.append("C_Phase")
    if E_Phase:
        ax.plot(np.array(E_Phase)[:,1], np.array(E_Phase)[:,0], ls="", marker="s", mfc=cdic["E_Phase"], mec=cdic["E_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["E_Phase"], markersize=10))
        labels.append("E_Phase")
    if X_Phase:
        ax.plot(np.array(X_Phase)[:,1], np.array(X_Phase)[:,0], ls="", marker="s", mfc=cdic["X_Phase"], mec=cdic["X_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["X_Phase"], markersize=10))
        labels.append("X_Phase")
    if U_Phase:
        ax.plot(np.array(U_Phase)[:,1], np.array(U_Phase)[:,0], ls="", marker="s", mfc=cdic["U_Phase"], mec=cdic["U_Phase"], ms=ms, alpha=1)
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=cdic["U_Phase"], markersize=10))
        labels.append("U_Phase")

    # Add the legend
    ax.legend(handles, labels, loc='best')

    # Set labels and styling
    ax.set_xlabel("Four-spin Interaction M", fontsize=14)
    ax.set_ylabel("Two-spin Interaction J", fontsize=14)
    ax.tick_params(axis="both", direction="in", width=2, length=4)

    # Set axis limits to show full parameter range
    ax.set_xlim(min(M_values), max(M_values))
    ax.set_ylim(min(J_values), max(J_values))

    # Add a title with parameter information
    plt.title(f"Phase Diagram (M vs. J)\n{len(J_values)}×{len(M_values)} grid", fontsize=12)

    fig.tight_layout()
    return fig