    import numpy as np
from mpmath import mp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator

from tqdm import tqdm

def get_initial_coefficients(J, M, n_max=20):
    """
    Calculate initial double Fourier coefficients for Ashkin-Teller model using fast numpy integration
    
    Args:
        J: 2-spin coupling constant
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
    u = np.exp(J*(np.cos(THETA) + np.cos(PHI)) + M*np.cos(THETA)*np.cos(PHI))
    
    # Initialize coefficient array
    size = 2*n_max + 1
    lambda_nm = np.zeros((size, size), dtype=complex)
    
    # Calculate coefficients using fast 2D integration
    for n in tqdm(range(-n_max, n_max + 1)):
        for m in range(-n_max, n_max + 1):
            n_idx = n + n_max
            m_idx = m + n_max
            
            # Fourier integrand
            integrand = u * np.exp(-1j*(n*THETA + m*PHI))
            
            # Fast 2D integration using numpy sum
            lambda_nm[n_idx, m_idx] = np.sum(integrand) * dtheta * dphi / (4*np.pi**2)
    
    # Normalize
    lambda_nm = lambda_nm / np.max(np.abs(lambda_nm))

    # Convert array to mpmath
    lambda_nm_mp = np.zeros_like(lambda_nm, dtype='object')
    for i in range(size):
        for j in range(size):
            lambda_nm_mp[i,j] = mp.mpc(lambda_nm[i,j])
    
    return lambda_nm_mp

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
    lambda_tilde = np.zeros((size, size), dtype='object')

    for k1 in range(-n_max, n_max + 1):
        for k2 in range(-n_max, n_max + 1):
            k1_idx = k1 + n_max
            k2_idx = k2 + n_max
            sum_terms = mp.mpf('0')
            
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
    max_val = max(abs(x) for row in lambda_tilde for x in row)
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
    lambda_prime = np.zeros((size, size), dtype='object')
    
    for n in range(-n_max, n_max + 1):
        for m in range(-n_max, n_max + 1):
            n_idx = n + n_max
            m_idx = m + n_max
            # For AT model, the middle spin integration gives same condition
            # for both theta and phi
            lambda_prime[n_idx, m_idx] = (lambda_12[n_idx, m_idx] * 
                                        lambda_23[n_idx, m_idx])
    
    # Normalize
    max_val = max(abs(x) for row in lambda_prime for x in row)
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
    lambda_tilde = lambda_nm
    for _ in range(m-1):
        lambda_tilde = bond_moving(lambda_tilde, lambda_nm)
        
    # Decimation
    lambda_prime = lambda_tilde
    for _ in range(b-1):
        lambda_prime = decimation(lambda_prime, lambda_tilde)
    
    return lambda_prime

def track_rg_flow(J, M, b=2, d=2, n_max=20, n_steps=20):
    """
    Track RG flow for Ashkin-Teller model
    
    Args:
        J: Initial 2-spin coupling
        M: Initial 4-spin coupling
        b: Length rescaling factor
        d: Dimension
        n_max: Maximum Fourier mode
        n_steps: Number of RG steps
        
    Returns:
        flow_history: Array of shape (n_steps+1, 2*n_max+1, 2*n_max+1)
    """
    # Store coefficient history - now using complex type
    size = 2*n_max + 1
    flow_history = np.zeros((n_steps+1, size, size), dtype=complex)
    
    # Get initial coefficients
    lambda_nm = get_initial_coefficients(J, M, n_max)
    # Store complex values directly
    flow_history[0] = lambda_nm
    
    # Perform RG steps
    for i in tqdm(range(n_steps)):
        lambda_nm = rg_step(lambda_nm, b, d)
        # Convert mpmath numbers to complex
        flow_history[i+1] = np.array([[complex(x) for x in row] for row in lambda_nm])
    
    return flow_history

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

def plot_evolution(flow_history, J, M, steps_to_show=None):
    """
    Visualize evolution of AT model under RG
    
    Args:
        flow_history: Array of coefficients history (complex)
        J: Initial 2-spin coupling
        M: Initial 4-spin coupling
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
    
    plt.suptitle(f'Ashkin-Teller Model Evolution (J={J}, M={M})')
    plt.tight_layout()
    return fig

def plot_coefficient_evolution(flow_history, J, M):
    """
    Plot evolution of selected coefficients over RG steps
    
    Args:
        flow_history: Array of coefficients history (complex)
        J: Initial 2-spin coupling
        M: Initial 4-spin coupling
    """
    n_steps = flow_history.shape[0]
    n_max = (flow_history.shape[1] - 1) // 2
    
    # Create figure with two subplots
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    #fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Select coefficients to track
    # Track (n,m) pairs: (0,0), (1,0), (0,1), (1,1), (2,0), (0,2), (2,2)
    coeff_indices = [(0,0), (1,0), (0,1), (1,1), (2,0), (0,2), (2,2)]
    
    # Plot magnitude evolution
    for n, m in coeff_indices:
        n_idx = n + n_max
        m_idx = m + n_max
        values = [flow_history[step, n_idx, m_idx] for step in range(n_steps)]
        ax.plot(range(n_steps), values, 'o-', label=f'(n,m)=({n},{m})')

    # Force x-axis to use only integers
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel('RG Step')
    ax.set_ylabel('λ(n,m)')
    #ax.set_yscale('log')
    ax.grid(True)
    ax.legend()
    ax.set_title('Evolution of Coefficient Magnitudes')
    
    plt.suptitle(f'Coefficient Evolution under RG (J={J}, M={M})')
    plt.tight_layout()
    return fig

def coefficient_sink(flow_history, coeff_indices = [(0,0), (1,0), (0,1), (1,1), (2,0), (0,2), (2,2)], rg_step = -1):
    n_max = (flow_history.shape[1] - 1) // 2
    for n, m in coeff_indices:
        n_idx = n + n_max
        m_idx = m + n_max
        print(f'({n},{m}) = {round(flow_history[rg_step, n_idx, m_idx].real, 5)}')