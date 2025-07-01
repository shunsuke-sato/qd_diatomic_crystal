import numpy as np
import matplotlib.pyplot as plt

def initialize_spins(N):
    """Initialize spins randomly to +1 or -1"""
    return np.random.choice([-1, 1], size=N)

def energy(spins, J):
    """Compute the total energy of the 1D spin configuration with periodic boundary"""
    return -J * np.sum(spins * np.roll(spins, 1))

def metropolis_step(spins, J, beta):
    """Perform one Metropolis update step"""
    N = len(spins)
    for _ in range(N):
        i = np.random.randint(N)
        # Periodic boundary neighbors
        left = spins[i - 1]
        right = spins[(i + 1) % N]
        dE = 2 * J * spins[i] * (left + right)
        # Metropolis criterion
        if dE < 0 or np.random.rand() < np.exp(-beta * dE):
            spins[i] *= -1
    return spins

def simulate(N, J, T, steps):
    beta = 1.0 / T
    spins = initialize_spins(N)
    energies = []

    for step in range(steps):
        spins = metropolis_step(spins, J, beta)
        if step % 10 == 0:
            E = energy(spins, J)
            energies.append(E / N)

    return spins, energies

def calc_correlation(spins):
    N = len(spins)
    corr = np.zeros(N)
    dist = np.zeros(N)

    for i in range(N):
        dist[i] = i
        corr[i] = (-1)**i*np.sum(spins * np.roll(spins, i))/N

    return dist, corr

# Parameters
N = 200           # Number of spins
J = -1.0          # Antiferromagnetic interaction (J < 0)
T = 0.2           # Temperature
steps = 1000      # Monte Carlo steps

final_spins, energy_trace = simulate(N, J, T, steps)
dist, correlation = calc_correlation(final_spins)

# Plot the energy over time
plt.plot(dist, correlation)
plt.xlabel('Length')
plt.ylabel('Correlation')

plt.show()
