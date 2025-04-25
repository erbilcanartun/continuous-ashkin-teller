import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from tqdm import tqdm

def get_initial_coefficients(J1, J2, M, n_max=20):
    """
    Calculate initial double Fourier coefficients for generalized Ashkin-Teller model.

    Args:
        J1: First 2-spin coupling constant
        J2: Second 2-spin coupling constant
        M: 4-spin coupling constant
        n_max: Maximum Fourier mode

    Returns:
        2D array of coefficients λ[n,m] from -n_max to n_max for both n,m
    """
    # Create grid for fast integration
    n_points = 100  # Number of points for integration
    theta = np.linspace(0, 2*np.pi, n_points)
    phi = np.linspace(0, 2*np.pi, n_points)
    THETA, PHI = np.meshgrid(theta, phi)
    dtheta = theta[1] - theta[0]
    dphi = phi[1] - phi[0]

    # Calculate potential on grid
    u = np.exp(J1*np.cos(THETA) + J2*np.cos(PHI) + M*np.cos(THETA)*np.cos(PHI))

    # Initialize coefficient array
    size = 2*n_max + 1
    lambda_nm = np.zeros((size, size), dtype=complex)

    # Calculate coefficients using fast 2D integration
    for n in range(-n_max, n_max + 1):
        for m in range(-n_max, n_max + 1):
            n_idx = n + n_max
            m_idx = m + n_max

            # Fourier integrand
            integrand = u * np.exp(-1j*(n*THETA + m*PHI))

            # Fast 2D integration using numpy sum
            lambda_nm[n_idx, m_idx] = np.sum(integrand) * dtheta * dphi / (4*np.pi**2)

    # Normalize
    lambda_nm = lambda_nm / np.max(np.abs(lambda_nm))

    return lambda_nm

def bond_moving(lambda_A, lambda_B):
    """
    Perform bond-moving step for Ashkin-Teller model

    Args:
        lambda_A: First 2D array of Fourier coefficients
        lambda_B: Second 2D array of Fourier coefficients

    Returns:
        Bond-moved coefficients
    """
    n_max = (lambda_A.shape[0] - 1) // 2
    size = 2*n_max + 1
    lambda_tilde = np.zeros((size, size), dtype=complex)

    for k1 in range(-n_max, n_max + 1):
        for k2 in range(-n_max, n_max + 1):
            k1_idx = k1 + n_max
            k2_idx = k2 + n_max
            sum_terms = 0.0 + 0.0j  # Initialize as complex number

            # Convolution in both dimensions
            for n1 in range(-n_max, n_max + 1):
                for n2 in range(-n_max, n_max + 1):
                    m1 = k1 - n1
                    m2 = k2 - n2
                    if abs(m1) <= n_max and abs(m2) <= n_max:
                        n1_idx = n1 + n_max
                        n2_idx = n2 + n_max
                        m1_idx = m1 + n_max
                        m2_idx = m2 + n_max
                        sum_terms += (lambda_A[n1_idx, n2_idx] *
                                      lambda_B[m1_idx, m2_idx])

            lambda_tilde[k1_idx, k2_idx] = sum_terms

    # Normalize
    max_val = np.max(np.abs(lambda_tilde))
    lambda_tilde = lambda_tilde / max_val

    return lambda_tilde

def decimation(lambda_12, lambda_23):
    """
    Perform decimation step for Ashkin-Teller model

    Args:
        lambda_12: Fourier coefficients of first interaction
        lambda_23: Fourier coefficients of second interaction

    Returns:
        Decimated coefficients
    """
    n_max = (lambda_12.shape[0] - 1) // 2
    size = 2*n_max + 1
    lambda_prime = np.zeros((size, size), dtype=complex)

    for n in range(-n_max, n_max + 1):
        for m in range(-n_max, n_max + 1):
            n_idx = n + n_max
            m_idx = m + n_max
            # For AT model, the middle spin integration gives same condition
            # for both theta and phi
            lambda_prime[n_idx, m_idx] = (lambda_12[n_idx, m_idx] *
                                         lambda_23[n_idx, m_idx])

    # Normalize
    max_val = np.max(np.abs(lambda_prime))
    lambda_prime = lambda_prime / max_val

    return lambda_prime

def rg_step(lambda_nm, b=2, d=2):
    """
    Perform one complete RG step for Ashkin-Teller model

    Args:
        lambda_nm: Input 2D coefficients
        b: Length rescaling factor
        d: Dimension

    Returns:
        Renormalized coefficients
    """
    # Calculate number of bonds to move/combine
    m = b**(d-1)

    # Bond moving
    lambda_tilde = lambda_nm.copy()  # Create a copy to avoid modifying the original
    for _ in range(m-1):
        lambda_tilde = bond_moving(lambda_tilde, lambda_nm)

    # Decimation
    lambda_prime = lambda_tilde.copy()  # Create a copy
    for _ in range(b-1):
        lambda_prime = decimation(lambda_prime, lambda_tilde)

    return lambda_prime

def track_rg_flow(J1, J2, M, b=2, d=2, n_max=20, n_steps=20):
    """
    Track RG flow for generalized Ashkin-Teller model

    Args:
        J1: First 2-spin coupling
        J2: Second 2-spin coupling
        M: 4-spin coupling
        b: Length rescaling factor
        d: Dimension
        n_max: Maximum Fourier mode
        n_steps: Number of RG steps

    Returns:
        flow_history: Array of shape (n_steps+1, 2*n_max+1, 2*n_max+1)
    """
    # Store coefficient history
    size = 2*n_max + 1
    flow_history = np.zeros((n_steps+1, size, size), dtype=complex)

    # Get initial coefficients
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)

    # Store initial coefficients
    flow_history[0] = lambda_nm

    # Perform RG steps
    for i in tqdm(range(n_steps)):
        lambda_nm = rg_step(lambda_nm, b, d)
        flow_history[i+1] = lambda_nm

    return flow_history

def display_coefficient_matrix(flow_history, step, decimals=4):
    """
    Display the coefficient matrix at a specific RG step with desired precision.

    Args:
        flow_history: Array containing coefficient matrices at each RG step
        step: The RG step to display (0 for initial, negative indices work for counting from the end)
        decimals: Number of decimal places to display
    """
    # Get the matrix at the specified step
    matrix = flow_history[step]
    n_max = (matrix.shape[0] - 1) // 2

    # Create a format string with alignment and decimal precision
    format_str = f"{{:{decimals+8}.{decimals}f}}"

    # Create row and column labels
    col_labels = [f"{m-n_max}" for m in range(matrix.shape[1])]
    row_labels = [f"{n-n_max}" for n in range(matrix.shape[0])]

    # Print column headers
    print(" " * 6, end="")
    for label in col_labels:
        print(f"{label:>{decimals+8}}", end="")
    print()

    # Print each row with row label
    for i, row_label in enumerate(row_labels):
        print(f"{row_label:<6}", end="")
        for j in range(matrix.shape[1]):
            # Convert to float and format
            value = float(matrix[i, j].real)
            print(format_str.format(value), end="")
        print()
    print()

def display_coefficient_matrix_colormap(matrix, decimals=4, component='real',
                               title=None, size_inches=(10, 8)):
    """
    Display the coefficient matrix using a heatmap plot with specified decimal precision.

    Args:
        matrix: 2D array of Fourier coefficients
        decimals: Number of decimal places to display (default: 4)
        component: Which component to display ('real', 'imag', or 'abs')
        title: Optional title for the plot
        size_inches: Size of the figure in inches (width, height)
    """
    # Make sure we're working with a 2D matrix
    if len(matrix.shape) > 2:
        raise ValueError("Input should be a 2D matrix, not higher dimensional")

    max_s = (matrix.shape[0] - 1) // 2

    # Extract the specified component with desired precision
    display_matrix = np.zeros((matrix.shape[0], matrix.shape[1]), dtype=float)

    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix[i, j]

            # Handle different types of values (complex, float)
            if np.iscomplexobj(val):  # Complex number
                if component == 'real':
                    display_matrix[i, j] = val.real
                elif component == 'imag':
                    display_matrix[i, j] = val.imag
                else:  # 'abs'
                    display_matrix[i, j] = abs(val)
            else:  # Already a float or similar
                display_matrix[i, j] = float(val)

    # Create figure
    plt.figure(figsize=size_inches)
    plt.imshow(display_matrix, cmap='viridis')
    plt.colorbar(format=f'%.{decimals}f')

    # Set labels and ticks
    plt.xlabel('m index')
    plt.ylabel('n index')

    # Create custom tick labels that show the actual Fourier indices
    tick_positions = np.arange(0, 2*max_s+1, max(1, (2*max_s+1)//10))  # Show up to 10 ticks
    tick_labels = [str(pos - max_s) for pos in tick_positions]

    plt.xticks(tick_positions, tick_labels)
    plt.yticks(tick_positions, tick_labels)

    if title:
        plt.title(title)
    else:
        plt.title(f'Coefficient Matrix ({component} component)')

    plt.tight_layout()
    plt.show()

def print_coefficient_matrix(matrix, decimals=4, component='real', max_size=5):
    """
    Print the coefficient matrix values with specified decimal precision.

    Args:
        matrix: 2D array of Fourier coefficients
        decimals: Number of decimal places to display (default: 4)
        component: Which component to display ('real', 'imag', or 'abs')
        max_size: Maximum size to display (central portion of the matrix)
    """
    # Make sure we're working with a 2D matrix
    if len(matrix.shape) > 2:
        raise ValueError("Input should be a 2D matrix, not higher dimensional")

    max_s = (matrix.shape[0] - 1) // 2

    # Determine the indices to display (central portion)
    display_size = min(max_size, 2*max_s+1)
    start_idx = max_s - display_size // 2
    end_idx = start_idx + display_size

    # Create format string
    fmt = f'{{:.{decimals}f}}'

    # Print the matrix indices
    print(f"\nCoefficient Matrix ({component} component), decimals={decimals}:")

    # Print column headers
    print(f"{'n\\m':<6}", end="")
    for m in range(start_idx, end_idx):
        print(f"{m-max_s:^{decimals+4}}", end="")
    print()

    # Print separator
    print("-" * (6 + (decimals+4) * display_size))

    # Print rows
    for n in range(start_idx, end_idx):
        print(f"{n-max_s:<6}", end="")
        for m in range(start_idx, end_idx):
            val = matrix[n, m]

            # Handle different types of values (complex, float)
            if np.iscomplexobj(val):  # Complex number
                if component == 'real':
                    val_to_display = val.real
                elif component == 'imag':
                    val_to_display = val.imag
                else:  # 'abs'
                    val_to_display = abs(val)
            else:  # Already a float or similar
                val_to_display = float(val)

            print(fmt.format(val_to_display), end="  ")
        print()

    # If we're showing a subset, indicate that
    if display_size < 2*max_s+1:
        print(f"\nNote: Showing central {display_size}×{display_size} portion of {matrix.shape[0]}×{matrix.shape[1]} matrix")

def reconstruct_potential(lambda_nm):
    """
    Reconstruct AT potential from Fourier coefficients

    Args:
        lambda_nm: 2D array of Fourier coefficients (complex)

    Returns:
        theta, phi: Mesh grids of angles
        u: Potential values (complex)
    """
    n_max = (lambda_nm.shape[0] - 1) // 2

    # Create angle grids
    theta = np.linspace(0, 2*np.pi, 100)
    phi = np.linspace(0, 2*np.pi, 100)
    THETA, PHI = np.meshgrid(theta, phi)

    # Initialize potential
    u = np.zeros_like(THETA, dtype=complex)

    # Sum over all modes
    for n in range(-n_max, n_max + 1):
        for m in range(-n_max, n_max + 1):
            n_idx = n + n_max
            m_idx = m + n_max
            u += lambda_nm[n_idx, m_idx] * np.exp(1j * (n*THETA + m*PHI))

    return THETA, PHI, u

def plot_evolution(flow_history, J1, J2, M, steps_to_show=None):
    """
    Visualize evolution of generalized AT model under RG

    Args:
        flow_history: Array of coefficients history (complex)
        J1: First 2-spin coupling
        J2: Second 2-spin coupling
        M: 4-spin coupling
        steps_to_show: List of RG steps to show
    """
    if steps_to_show is None:
        steps_to_show = [0, 1, 5, -1]

    n_steps = len(steps_to_show)
    fig = plt.figure(figsize=(16, 4*n_steps))

    for i, step in enumerate(steps_to_show):
        # Coefficient heatmap - using absolute values
        ax1 = plt.subplot(n_steps, 2, 2*i + 1)
        # Take magnitude of complex coefficients
        coeff_magnitude = np.abs(flow_history[step])
        im = ax1.imshow(coeff_magnitude, 
                       norm=LogNorm(vmin=max(1e-10, coeff_magnitude.min())),
                       cmap='viridis')
        plt.colorbar(im, ax=ax1)
        ax1.set_title(f'Coefficient Magnitudes at Step {step}')
        ax1.set_xlabel('m index')
        ax1.set_ylabel('n index')

        # Potential surface - using absolute values
        ax2 = plt.subplot(n_steps, 2, 2*i + 2, projection='3d')
        THETA, PHI, u = reconstruct_potential(flow_history[step])
        # Take magnitude of complex potential
        u = np.abs(u)
        surf = ax2.plot_surface(THETA, PHI, u, cmap='viridis')
        plt.colorbar(surf, ax=ax2)
        ax2.set_title(f'Potential Magnitude at Step {step}')
        ax2.set_xlabel('θ')
        ax2.set_ylabel('φ')
        ax2.set_zlabel('|u(θ,φ)|')

    plt.suptitle(f'Ashkin-Teller Model Evolution (J1={J1}, J2={J2}, M={M})')
    plt.tight_layout()
    return fig

def plot_coefficient_evolution(flow_history, J1, J2, M):
    """
    Plot evolution of selected coefficients over RG steps

    Args:
        flow_history: Array of coefficients history (complex)
        J1: First 2-spin coupling
        J2: Second 2-spin coupling
        M: 4-spin coupling
    """
    n_steps = flow_history.shape[0]
    n_max = (flow_history.shape[1] - 1) // 2

    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))

    # Select coefficients to track
    # Track (n,m) pairs: (0,0), (1,0), (0,1), (1,1), (2,0), (0,2), (2,2)
    coeff_indices = [(0,0), (1,0), (0,1), (1,1), (2,0), (0,2), (2,2)]

    # Plot magnitude evolution
    for n, m in coeff_indices:
        n_idx = n + n_max
        m_idx = m + n_max
        #values = np.abs([flow_history[step, n_idx, m_idx] for step in range(n_steps)])
        values = [flow_history[step, n_idx, m_idx].real for step in range(n_steps)]
        ax.plot(range(n_steps), values, 'o-', label=f'({n},{m})')

    # Force x-axis to use only integers
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('RG Step')
    ax.set_ylabel('|λ(n,m)|')
    ax.grid(True)
    ax.legend(title="(n,m)")
    ax.set_title('Evolution of Coefficient Magnitudes')

    plt.suptitle(f'Coefficient Evolution under RG (J1={J1}, J2={J2}, M={M})')
    plt.tight_layout()
    return fig

def coefficient_sink(flow_history, coeff_indices=[(0,0), (1,0), (0,1), (1,1), (2,0), (0,2), (2,2)], rg_step=-1):
    """
    Print coefficient values at the specified RG step (default: last step)

    Args:
        flow_history: Array of coefficients history
        coeff_indices: List of (n,m) pairs to display
        rg_step: Which RG step to examine
    """
    n_max = (flow_history.shape[1] - 1) // 2
    for n, m in coeff_indices:
        n_idx = n + n_max
        m_idx = m + n_max
        print(f'({n},{m}) = {round(flow_history[rg_step, n_idx, m_idx].real, 5)}')