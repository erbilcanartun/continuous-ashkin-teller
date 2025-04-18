import matplotlib.pyplot as plt
from joblib import Parallel, delayed
from utils import plot_phase_diagram
from phase_hybrid import identify_phase_hybrid


def generate_phase_diagram_joblib(J_values, M_values, n_max_np=20, n_max_mp=10, 
                                n_rg_steps=20, b=2, d=2, history_length=3, dps=30, path=None, n_jobs=1, only_mp=False):
    """
    Generate phase diagram by scanning J, M parameter space using parallel processing with joblib.
    
    Args:
        J_values: Array of J values
        M_values: Array of M values
        n_max_np: Maximum Fourier mode for numpy
        n_max_mp: Maximum Fourier mode for mpmath
        n_rg_steps: Number of RG steps
        b: Length rescaling factor
        d: Dimension
        dps: Decimal places for mpmath
        path: Path to save the figure
        n_jobs: Number of parallel jobs (default -1 uses all cores)

    Returns:
        Tuple of phase points (D_Phase, A_Phase, B_Phase, C_Phase, X_Phase, U_Phase)
    """

    # Function to process one point in parameter space
    def process_one_point(J, M):
        J1 = J2 = J
        phase = identify_phase_hybrid(J1, J2, M, n_max_np=n_max_np,
                                      n_max_mp=n_max_mp, n_rg_steps=n_rg_steps,
                                      b=b, d=d, history_length=history_length,
                                      dps=dps, only_mp=only_mp)
        return J, M, phase

    # Create list of all parameter combinations
    total_points = len(J_values) * len(M_values)
    print(f"Starting calculation for {total_points} points with {n_jobs} jobs")

    # Process all points in parallel
    results = Parallel(n_jobs=n_jobs, verbose=10)(
        delayed(process_one_point)(J, M)
        for J in J_values
        for M in M_values
    )

    # Organize results into phase groups
    D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase = [], [], [], [], [], [], []

    for J, M, phase in results:
        if phase == "D_Phase":
            D_Phase.append([J, M])
        elif phase == "A_Phase":
            A_Phase.append([J, M])
        elif phase == "B_Phase":
            B_Phase.append([J, M])
        elif phase == "C_Phase":
            C_Phase.append([J, M])
        elif phase == "E_Phase":
            E_Phase.append([J, M])
        elif phase == "X_Phase":
            X_Phase.append([J, M])
        else:
            U_Phase.append([J, M])

    # Plot the results
    fig = plot_phase_diagram(J_values, M_values, D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase)
    if path:
        plt.savefig(path, dpi=300, bbox_inches='tight')

    return D_Phase, A_Phase, B_Phase, C_Phase, E_Phase, X_Phase, U_Phase