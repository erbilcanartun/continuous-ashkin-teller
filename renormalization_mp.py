import numpy as np
from mpmath import mp, sin, cos, pi, exp, matrix, conj
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from tqdm import tqdm


def get_initial_coefficients(J1, J2, M, n_max=20):
    """
    Calculate initial double Fourier coefficients for generalized Ashkin-Teller model
    using high precision mpmath calculations.
    
    Args:
        J1: First 2-spin coupling constant (will be converted to mpmath)
        J2: Second 2-spin coupling constant (will be converted to mpmath)
        M: 4-spin coupling constant (will be converted to mpmath)
        n_max: Maximum Fourier mode
    
    Returns:
        2D array of coefficients λ[n,m] from -n_max to n_max for both n,m
    """
    # Convert input parameters to mpmath precision
    J1, J2, M = mp.mpf(J1), mp.mpf(J2), mp.mpf(M)
    
    # Create grid for integration with mpmath precision
    n_points = 100  # Number of points for integration
    dtheta = mp.mpf('2') * mp.pi / mp.mpf(n_points)
    dphi = mp.mpf('2') * mp.pi / mp.mpf(n_points)
    
    # Initialize coefficient array with mpmath complex numbers
    size = 2*n_max + 1
    lambda_nm = np.zeros((size, size), dtype='object')
    
    # Calculate coefficients using numerical integration
    for n in range(-n_max, n_max + 1):
        n_idx = n + n_max
        for m in range(-n_max, n_max + 1):
            m_idx = m + n_max
            
            # Initialize sum for numerical integration
            integral_sum = mp.mpc('0')
            
            # Perform 2D numerical integration
            for i in range(n_points):
                theta = mp.mpf(i) * dtheta
                for j in range(n_points):
                    phi = mp.mpf(j) * dphi
                    
                    # Calculate potential at this point
                    u_val = mp.exp(J1*mp.cos(theta) + J2*mp.cos(phi) + M*mp.cos(theta)*mp.cos(phi))
                    
                    # Fourier integrand
                    integrand = u_val * mp.exp(mp.mpc('0', '-1') * (mp.mpf(n)*theta + mp.mpf(m)*phi))
                    
                    # Add to integral sum
                    integral_sum += integrand
            
            # Finish integration with proper scaling
            lambda_nm[n_idx, m_idx] = integral_sum * dtheta * dphi / (mp.mpf('4') * mp.pi * mp.pi)
    
    # Normalize by maximum absolute value
    max_val = mp.mpf('0')
    for i in range(size):
        for j in range(size):
            abs_val = mp.fabs(lambda_nm[i, j])
            if abs_val > max_val:
                max_val = abs_val
    
    for i in range(size):
        for j in range(size):
            lambda_nm[i, j] = lambda_nm[i, j] / max_val
    
    return lambda_nm

def bond_moving(lambda_A, lambda_B):
    """
    Perform bond-moving step for Ashkin-Teller model with mpmath precision
    
    Args:
        lambda_A: First 2D array of Fourier coefficients (mpmath complex)
        lambda_B: Second 2D array of Fourier coefficients (mpmath complex)
    
    Returns:
        Bond-moved coefficients (mpmath complex)
    """
    n_max = (lambda_A.shape[0] - 1) // 2
    size = 2*n_max + 1
    lambda_tilde = np.zeros((size, size), dtype='object')
    
    for k1 in range(-n_max, n_max + 1):
        k1_idx = k1 + n_max
        for k2 in range(-n_max, n_max + 1):
            k2_idx = k2 + n_max
            sum_terms = mp.mpc('0')  # Initialize as mpmath complex number
            
            # Convolution in both dimensions
            for n1 in range(-n_max, n_max + 1):
                n1_idx = n1 + n_max
                for n2 in range(-n_max, n_max + 1):
                    n2_idx = n2 + n_max
                    m1 = k1 - n1
                    m2 = k2 - n2
                    if abs(m1) <= n_max and abs(m2) <= n_max:
                        m1_idx = m1 + n_max
                        m2_idx = m2 + n_max
                        sum_terms += (lambda_A[n1_idx, n2_idx] *
                                      lambda_B[m1_idx, m2_idx])
            
            lambda_tilde[k1_idx, k2_idx] = sum_terms
    
    # Normalize by maximum absolute value
    max_val = mp.mpf('0')
    for i in range(size):
        for j in range(size):
            abs_val = mp.fabs(lambda_tilde[i, j])
            if abs_val > max_val:
                max_val = abs_val
    
    for i in range(size):
        for j in range(size):
            lambda_tilde[i, j] = lambda_tilde[i, j] / max_val
    
    return lambda_tilde

def decimation(lambda_12, lambda_23):
    """
    Perform decimation step for Ashkin-Teller model with mpmath precision
    
    Args:
        lambda_12: Fourier coefficients of first interaction (mpmath complex)
        lambda_23: Fourier coefficients of second interaction (mpmath complex)
    
    Returns:
        Decimated coefficients (mpmath complex)
    """
    n_max = (lambda_12.shape[0] - 1) // 2
    size = 2*n_max + 1
    lambda_prime = np.zeros((size, size), dtype='object')
    
    for n in range(-n_max, n_max + 1):
        n_idx = n + n_max
        for m in range(-n_max, n_max + 1):
            m_idx = m + n_max
            # For AT model, the middle spin integration gives same condition
            # for both theta and phi
            lambda_prime[n_idx, m_idx] = (lambda_12[n_idx, m_idx] *
                                         lambda_23[n_idx, m_idx])
    
    # Normalize by maximum absolute value
    max_val = mp.mpf('0')
    for i in range(size):
        for j in range(size):
            abs_val = mp.fabs(lambda_prime[i, j])
            if abs_val > max_val:
                max_val = abs_val
    
    for i in range(size):
        for j in range(size):
            lambda_prime[i, j] = lambda_prime[i, j] / max_val
    
    return lambda_prime

def rg_step(lambda_nm, b=2, d=2):
    """
    Perform one complete RG step for Ashkin-Teller model with mpmath precision
    
    Args:
        lambda_nm: Input 2D coefficients (mpmath complex)
        b: Length rescaling factor
        d: Dimension
    
    Returns:
        Renormalized coefficients (mpmath complex)
    """
    
    # Calculate number of bonds to move/combine
    m = b**(d-1)
    
    # Bond moving (create a deep copy to avoid modifying the original)
    lambda_tilde = np.zeros_like(lambda_nm, dtype='object')
    for i in range(lambda_nm.shape[0]):
        for j in range(lambda_nm.shape[1]):
            lambda_tilde[i, j] = lambda_nm[i, j]
            
    for _ in range(m-1):
        lambda_tilde = bond_moving(lambda_tilde, lambda_nm)
    
    # Decimation (create a deep copy)
    lambda_prime = np.zeros_like(lambda_tilde, dtype='object')
    for i in range(lambda_tilde.shape[0]):
        for j in range(lambda_tilde.shape[1]):
            lambda_prime[i, j] = lambda_tilde[i, j]
            
    for _ in range(b-1):
        lambda_prime = decimation(lambda_prime, lambda_tilde)
    
    return lambda_prime

def track_rg_flow(J1, J2, M, b=2, d=2, n_max=20, n_steps=20):
    """
    Track RG flow for generalized Ashkin-Teller model with mpmath precision
    
    Args:
        J1: First 2-spin coupling
        J2: Second 2-spin coupling
        M: 4-spin coupling
        b: Length rescaling factor
        d: Dimension
        n_max: Maximum Fourier mode
        n_steps: Number of RG steps
    
    Returns:
        flow_history: Array of shape (n_steps+1, 2*n_max+1, 2*n_max+1) with mpmath values
    """
    # Store coefficient history
    size = 2*n_max + 1
    flow_history = np.zeros((n_steps+1, size, size), dtype='object')
    free_energy_constants = np.zeros(n_steps+1, dtype='object')
    
    # Get initial coefficients with mpmath precision
    print("Calculating initial coefficients...")
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)
    
    # Store initial coefficients
    for i in range(size):
        for j in range(size):
            flow_history[0, i, j] = lambda_nm[i, j]
    
    # Perform RG steps
    for i in tqdm(range(n_steps), desc="Performing RG steps"):
        lambda_nm = rg_step(lambda_nm, b, d)
        
        # Store the results
        for j in range(size):
            for k in range(size):
                flow_history[i+1, j, k] = lambda_nm[j, k]
    
    return flow_history

def track_rg_flow_free_energy(J1, J2, M, b=2, d=2, n_max=20, n_steps=20):
    """
    Track RG flow for generalized Ashkin-Teller model with mpmath precision
    
    Args:
        J1: First 2-spin coupling
        J2: Second 2-spin coupling
        M: 4-spin coupling
        b: Length rescaling factor
        d: Dimension
        n_max: Maximum Fourier mode
        n_steps: Number of RG steps
    
    Returns:
        tuple: (flow_history, free_energy_constants)
            flow_history: Array of shape (n_steps+1, 2*n_max+1, 2*n_max+1) with mpmath values
            free_energy_constants: Array of free energy constants at each RG step
    """
    # Store coefficient history
    size = 2*n_max + 1
    flow_history = np.zeros((n_steps+1, size, size), dtype='object')
    free_energy_constants = np.zeros(n_steps+1, dtype='object')
    
    # Convert parameters to mpmath precision
    b, d = mp.mpf(b), mp.mpf(d)
    
    # Get initial coefficients with mpmath precision
    print("Calculating initial coefficients...")
    lambda_nm = get_initial_coefficients(J1, J2, M, n_max)
    
    # Store initial coefficients
    for i in range(size):
        for j in range(size):
            flow_history[0, i, j] = lambda_nm[i, j]
    
    # Initial free energy constant is 0
    free_energy_constants[0] = mp.mpf('0')
    
    # Perform RG steps
    for i in tqdm(range(n_steps), desc="Performing RG steps"):
        # Before RG transformation, store the maximum coefficient value
        max_val_before = mp.mpf('0')
        for j in range(size):
            for k in range(size):
                abs_val = mp.fabs(lambda_nm[j, k])
                if abs_val > max_val_before:
                    max_val_before = abs_val
        
        # Perform RG step
        lambda_nm = rg_step(lambda_nm, b, d)
        
        # After RG transformation, find the maximum coefficient value
        max_val_after = mp.mpf('0')
        for j in range(size):
            for k in range(size):
                abs_val = mp.fabs(lambda_nm[j, k])
                if abs_val > max_val_after:
                    max_val_after = abs_val
        
        # Calculate free energy constant for this step
        # G^(n) = ln(λ_max) when normalization is done by dividing by max value
        # This represents the additive constant in the Hamiltonian
        if max_val_before > mp.mpf('1e-50'):  # Avoid log of very small numbers
            free_energy_constant = mp.log(max_val_before)
        else:
            free_energy_constant = mp.mpf('0')
        
        # Store the free energy constant, scaled by the appropriate power of b
        # to account for the rescaling of the system
        free_energy_constants[i+1] = free_energy_constant / mp.power(b, i * d)
        
        # Store the renormalized coefficients
        for j in range(size):
            for k in range(size):
                flow_history[i+1, j, k] = lambda_nm[j, k]
    
    return flow_history, free_energy_constants

def reconstruct_potential(lambda_nm):
    """
    Reconstruct AT potential from Fourier coefficients with mpmath precision
    
    Args:
        lambda_nm: 2D array of Fourier coefficients (mpmath complex)
    
    Returns:
        theta, phi: Mesh grids of angles (numpy arrays)
        u: Potential values (numpy array of complex values)
    """
    n_max = (lambda_nm.shape[0] - 1) // 2
    
    # Create angle grids (using numpy for visualization purposes)
    theta = np.linspace(0, float(2*mp.pi), 100)
    phi = np.linspace(0, float(2*mp.pi), 100)
    THETA, PHI = np.meshgrid(theta, phi)
    
    # Initialize potential with numpy complex for visualization
    u = np.zeros_like(THETA, dtype=complex)
    
    # Sum over all modes, converting mpmath complex to numpy complex for visualization
    for n in range(-n_max, n_max + 1):
        n_idx = n + n_max
        for m in range(-n_max, n_max + 1):
            m_idx = m + n_max
            # Convert mpmath complex to numpy complex
            coeff = complex(lambda_nm[n_idx, m_idx])
            u += coeff * np.exp(1j * (n*THETA + m*PHI))
    
    return THETA, PHI, u

def plot_evolution(flow_history, J1, J2, M, steps_to_show=None):
    """
    Visualize evolution of generalized AT model under RG
    
    Args:
        flow_history: Array of coefficients history (mpmath complex)
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
        
        # Convert mpmath values to numpy for visualization
        size = flow_history.shape[1]
        coeff_magnitude = np.zeros((size, size))
        for j in range(size):
            for k in range(size):
                coeff_magnitude[j, k] = float(mp.fabs(flow_history[step, j, k]))
        
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
        flow_history: Array of coefficients history (mpmath complex)
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
        # Convert mpmath values to numpy for plotting
        values = np.zeros(n_steps)
        for step in range(n_steps):
            values[step] = float(mp.fabs(flow_history[step, n_idx, m_idx]))
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
        flow_history: Array of coefficients history (mpmath complex)
        coeff_indices: List of (n,m) pairs to display
        rg_step: Which RG step to examine
    """
    n_max = (flow_history.shape[1] - 1) // 2
    for n, m in coeff_indices:
        n_idx = n + n_max
        m_idx = m + n_max
        # Use mpmath's nstr for controlled precision output
        print(f'({n},{m}) = {mp.nstr(flow_history[rg_step, n_idx, m_idx].real, n=6)}')


def calculate_thermodynamic_quantities(J1_values, J2, M, b=2, d=2, n_max=20, n_steps=20):
    """
    Calculate free energy, internal energy, and specific heat for a range of J1 values.
    
    Args:
        J1_values: Array of J1 values (first coupling constant)
        J2: Second coupling constant
        M: 4-spin coupling constant
        b: Length rescaling factor
        d: Dimension
        n_max: Maximum Fourier mode
        n_steps: Number of RG steps
    
    Returns:
        tuple: (J1_values, free_energies, internal_energies, specific_heats)
    """
    # Convert parameters to mpmath
    b, d = mp.mpf(b), mp.mpf(d)
    J2, M = mp.mpf(J2), mp.mpf(M)
    
    # Initialize arrays for thermodynamic quantities
    num_points = len(J1_values)
    free_energies = np.zeros(num_points, dtype='object')
    internal_energies = np.zeros(num_points, dtype='object')
    specific_heats = np.zeros(num_points, dtype='object')
    
    # Calculate thermodynamic quantities for each J1 value
    for i, J1 in enumerate(tqdm(J1_values, desc="Calculating thermodynamics")):
        # Convert J1 to mpmath
        J1_mp = mp.mpf(J1)
        
        # Run RG flow and get free energy constants
        _, free_energy_constants = track_rg_flow(J1_mp, J2, M, b, d, n_max, n_steps)
        
        # Calculate free energy per bond by summing the constants
        free_energy = mp.mpf('0')
        for n in range(len(free_energy_constants)):
            free_energy += free_energy_constants[n]
        
        free_energies[i] = free_energy
        
        # Calculate derivatives for internal energy and specific heat
        # For internal energy: dF/dJ1
        # For specific heat: d²F/dJ1²
        
        # We'll use finite difference for derivatives
        if i > 0 and i < num_points - 1:
            # Central difference for internal energy
            dJ1 = mp.mpf(J1_values[i+1] - J1_values[i-1])
            dF = free_energies[i+1] - free_energies[i-1]
            internal_energies[i] = -dF / dJ1  # Negative sign for U = -dF/dβ (β=1/kT)
            
            # Central difference for specific heat
            if i > 1 and i < num_points - 2:
                d2F = (free_energies[i+2] - 2*free_energies[i] + free_energies[i-2])
                d2J1 = mp.mpf(J1_values[i+2] - J1_values[i-2])
                specific_heats[i] = J1_mp * J1_mp * d2F / (d2J1 * d2J1)  # C = kT² * d²F/dT²
    
    # Handle endpoints for derivatives (forward/backward differences)
    if num_points > 1:
        # For i=0 (first point)
        dJ1 = mp.mpf(J1_values[1] - J1_values[0])
        dF = free_energies[1] - free_energies[0]
        internal_energies[0] = -dF / dJ1
        
        # For i=n-1 (last point)
        dJ1 = mp.mpf(J1_values[-1] - J1_values[-2])
        dF = free_energies[-1] - free_energies[-2]
        internal_energies[-1] = -dF / dJ1
        
        # Specific heat at endpoints
        if num_points > 3:
            # First point
            d2F = (free_energies[2] - 2*free_energies[1] + free_energies[0])
            d2J1 = mp.mpf(J1_values[2] - J1_values[0])
            specific_heats[0] = mp.mpf(J1_values[0])**2 * d2F / (d2J1 * d2J1)
            
            # Last point
            d2F = (free_energies[-1] - 2*free_energies[-2] + free_energies[-3])
            d2J1 = mp.mpf(J1_values[-1] - J1_values[-3])
            specific_heats[-1] = mp.mpf(J1_values[-1])**2 * d2F / (d2J1 * d2J1)
    
    return J1_values, free_energies, internal_energies, specific_heats

def plot_thermodynamic_quantities(J1_values, free_energies, internal_energies, specific_heats, J2, M):
    """
    Plot thermodynamic quantities as functions of J1.
    
    Args:
        J1_values: Array of J1 values
        free_energies: Array of free energy values
        internal_energies: Array of internal energy values
        specific_heats: Array of specific heat values
        J2: Second coupling constant
        M: 4-spin coupling constant
    
    Returns:
        matplotlib figure
    """
    # Convert mpmath values to float for plotting
    free_energies_float = np.array([float(val) for val in free_energies])
    internal_energies_float = np.array([float(val) for val in internal_energies])
    specific_heats_float = np.array([float(val) for val in specific_heats])
    
    # Create figure with 3 subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
    
    # Plot free energy
    ax1.plot(J1_values, free_energies_float, 'o-', color='blue', linewidth=2)
    ax1.set_ylabel('Free Energy per Bond (F/kN)')
    ax1.set_title('Thermodynamic Quantities vs. J1')
    ax1.grid(True)
    
    # Plot internal energy
    ax2.plot(J1_values, internal_energies_float, 'o-', color='red', linewidth=2)
    ax2.set_ylabel('Internal Energy per Bond (U/kN)')
    ax2.grid(True)
    
    # Plot specific heat
    ax3.plot(J1_values, specific_heats_float, 'o-', color='green', linewidth=2)
    ax3.set_xlabel('Coupling Constant J1')
    ax3.set_ylabel('Specific Heat per Bond (C/kN)')
    ax3.grid(True)
    
    # Mark any potential phase transitions (look for peaks in specific heat)
    if len(specific_heats_float) > 3:
        peaks = []
        for i in range(1, len(specific_heats_float)-1):
            if (specific_heats_float[i] > specific_heats_float[i-1] and 
                specific_heats_float[i] > specific_heats_float[i+1]):
                peaks.append(i)
        
        for i in peaks:
            # Critical point at peak in specific heat
            ax1.axvline(x=J1_values[i], color='k', linestyle='--', alpha=0.5)
            ax2.axvline(x=J1_values[i], color='k', linestyle='--', alpha=0.5)
            ax3.axvline(x=J1_values[i], color='k', linestyle='--', alpha=0.5)
            ax3.plot(J1_values[i], specific_heats_float[i], 'ko', markersize=8)
    
    plt.tight_layout()
    plt.suptitle(f'Ashkin-Teller Model Thermodynamics (J2={J2}, M={M})', y=0.98)
    plt.subplots_adjust(top=0.92)
    
    return fig


def analyze_ashkin_teller_phase_transition():
    """
    Analyze the phase transition in the Ashkin-Teller model by tracking
    thermodynamic quantities as a function of temperature.
    """
    # Set parameters
    J2 = 0.5     # Fixed second coupling
    M = 0.2      # Fixed 4-spin coupling
    
    # We'll vary J1 from low to high (equivalent to varying temperature)
    J1_values = np.linspace(0.2, 1.5, 20)
    
    # System parameters
    b = 2        # Length rescaling factor
    d = 2        # Spatial dimension
    n_max = 10   # Maximum Fourier mode (reduced for faster computation)
    n_steps = 10 # Number of RG steps (reduced for faster computation)
    
    # Calculate thermodynamic quantities
    J1_values, free_energies, internal_energies, specific_heats = calculate_thermodynamic_quantities(
        J1_values, J2, M, b, d, n_max, n_steps
    )
    
    # Plot results
    fig = plot_thermodynamic_quantities(J1_values, free_energies, internal_energies, specific_heats, J2, M)
    plt.savefig('ashkin_teller_thermodynamics.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Find and print the phase transition point
    if len(specific_heats) > 3:
        max_index = np.argmax([float(val) for val in specific_heats])
        critical_J1 = J1_values[max_index]
        print(f"Critical point estimate: J1_c ≈ {critical_J1}")
        print(f"Specific heat at critical point: C/kN ≈ {float(specific_heats[max_index])}")
    
    # Example of examining a specific point on the phase diagram
    critical_index = len(J1_values) // 2  # Just an example point
    J1_critical = J1_values[critical_index]
    print(f"\nExamining RG flow at J1 = {J1_critical}:")
    
    # Run the RG flow and plot the evolution of coefficients
    flow_history, _ = track_rg_flow(J1_critical, J2, M, b, d, n_max, n_steps)
    fig_coeff = plot_coefficient_evolution(flow_history, J1_critical, J2, M)
    plt.savefig('ashkin_teller_coefficient_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Visualize the potential at different RG steps
    fig_evo = plot_evolution(flow_history, J1_critical, J2, M, steps_to_show=[0, 2, 5, -1])
    plt.savefig('ashkin_teller_rg_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print the coefficients at the final RG step
    print("\nFinal coefficients at the RG sink:")
    coefficient_sink(flow_history)
    
    return {
        'J1_values': J1_values,
        'free_energies': free_energies,
        'internal_energies': internal_energies,
        'specific_heats': specific_heats,
        'flow_history': flow_history
    }