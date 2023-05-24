import csv

import numpy as np


class GrossPitaevskiiProblem:
    """Class that generates a Gross-Pitaevskii problem instance,
    in harmonic oscillator units, for a 3D harmonic trap.
    """

    def __init__(
        self,
        particle_number,
        grid_size,
        grid_step,
        scattering_length,
        sigma,
        time_step,
        iterations,
        interacting_system=True,
    ):
        # Put default values in the declaration above
        self.particle_number = particle_number
        self.grid_size = grid_size
        self.grid_step = grid_step
        self.scattering_length = scattering_length
        self.sigma = sigma
        self.time_step = time_step
        self.iterations = iterations
        self.interacting_system = interacting_system

    def evolution(self):
        if self.interacting_system:
            interaction = self.scattering_length * self.particle_number
        else:
            interaction = 0

        r_vector = np.arange(0., self.grid_size, self.grid_step)
        psi = self.ansatz()
        mu = np.zeros(len(r_vector))
        psi0 = psi
        for i in range(self.iterations):
            normalization = self.simpson_integral(psi0 ** 2, self.grid_step)
            psi0 = psi0 / np.sqrt(normalization)
            is_normal = self.simpson_integral(psi0 ** 2, self.grid_step)
            if abs(is_normal) > 1.1:
                raise ValueError(f"Wave function is not normalized, Normalization: {is_normal}")
            ddpsi = self.second_derivative(psi0, self.grid_step)
            mu = self._calculate_mu(r_vector=r_vector, 
                                    psi=psi0, 
                                    dpsi=ddpsi, 
                                    interaction=interaction)
            psi = psi0 - self.time_step * mu * psi0
            psi0 = psi
        normalization = self.simpson_integral(psi0 ** 2, self.grid_step)
        psi0 = psi0 / np.sqrt(normalization)
        return psi0

    def integrate_problem(self):
        if self.interacting_system:
            interaction = self.scattering_length * self.particle_number
        else:
            interaction = 0
        r_vector = np.arange(0., self.grid_size, self.grid_step)
        with open("results.csv", "w") as f:
            csv_writer = csv.writer(f)
            # Initialize the wave function
            csv_writer.writerow(["Energy", "Kinetic Energy", "Harmonic Potential Energy", "Interaction Energy", "Radii"])
            psi = self.ansatz()

            mu = np.zeros(len(r_vector))
            # energy = np.zeros(self.iterations)
            # kinetic_energy = np.zeros(self.iterations)
            # harmonic_potential_energy = np.zeros(self.iterations)
            # interaction_energy = np.zeros(self.iterations)
            # radii = np.zeros(self.iterations)
            # kinetic_term = np.zeros(len(r_vector))
            # harmonic_potential_term = np.zeros(len(r_vector))
            # interaction_energy_term = np.zeros(len(r_vector))
            # radii_term = np.zeros(len(r_vector))

            for i in range(self.iterations):
                normalization = self.simpson_integral(psi ** 2, self.grid_step)
                print(f"Normalization: {normalization}")
                psi = psi / np.sqrt(normalization)
                ddpsi = self.second_derivative(psi, self.grid_step)
                mu = self._calculate_mu(r_vector=r_vector, psi=psi, dpsi=ddpsi, interaction=interaction)
                # for i, x in enumerate(r_vector):
                #     harmonic_potential_term[i] = 0.5 * psi[i] ** 2 * x ** 2
                #     if i == 0:
                #         interaction_energy_term[i] = 0
                #         continue
                #     interaction_energy_term[i] = interaction / 2 * psi[i] ** 4 / x ** 2

                # kinetic_energy[i] = self.simpson_integral(
                #     -0.5 * psi * ddpsi, self.grid_step
                # )

                # harmonic_potential_energy[i] = self.simpson_integral(
                #     harmonic_potential_term, self.grid_step
                # )
                # interaction_energy[i] = self.simpson_integral(
                #     interaction_energy_term, self.grid_step
                # )
                # energy[i] = (
                #     kinetic_energy[i] 
                #     + harmonic_potential_energy[i]
                #     + interaction_energy[i]
                # )
                
                # radii[i] = self.simpson_integral(radii_term, self.grid_step)
                # csv_writer.writerow([energy[i],
                #                     kinetic_energy[i],
                #                     harmonic_potential_energy[i],
                #                     interaction_energy[i],
                #                     radii[i]]
                #                     )
                psi = psi - self.time_step * mu * psi
        return 0

    def _calculate_mu(self, r_vector, psi, dpsi, interaction):
        """Method to calculate the chemical potential

        Args:
            r_vector (array like): values of the r-vector
            psi (array like): values of the wave function
            dpsi (array like): values of the second derivative of the wave function

        Returns:
            float: chemical potential
        """
        mu = np.zeros(len(r_vector))
        for i, x in enumerate(r_vector):
            if i == 0:
                mu[i] = 0
                continue
            mu[i] = (
                - 0.5 * dpsi[i] / psi[i]
                + 0.5 * x ** 2 * psi[i]
                + 0.5 * interaction * psi[i] ** 2
            )
        return mu

    @staticmethod
    def second_derivative(f, h):
        derivative = np.zeros(len(f))
        for i, _ in enumerate(f):
            if i == 0:
                derivative[i] = 0
            elif i == len(f) - 1:
                derivative[i] = (f[i - 1] - 2 * f[i]) / (h ** 2)
            else:
                derivative[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / (h ** 2)
        return derivative

    @staticmethod
    def simpson_integral(f, h):
        """Method to calculate integrals using Simpson's rule

        Args:
            f (array like): values of the function to integrate, separated by h
            h (float): step size

        Returns:
            float: integral value
        """
        n = len(f) - 1
        integral = sum(2 * f[i] if i % 2 == 0 else 4 * f[i] for i in range(1, n))
        integral += f[0] + f[n]
        return integral * h / 3

    def ansatz(self):
        r = np.arange(0, self.grid_size, self.grid_step)
        psi = np.zeros(len(r))
        cvar = 2 * np.sqrt(self.sigma) ** 3 / np.sqrt(np.sqrt(np.pi))
        psi = cvar*r*np.exp(-0.5*self.sigma**2*r**2)
        return psi

    def __str__(self) -> str:
        return f"{self.particle_number} Bosons in a spherical trap \n r-grid in {self.grid_size} points, r-step {self.grid_step} \n A0={self.scattering_length}, sigma={self.sigma} \n time={self.time_step}, number-iter={self.iterations}"


if __name__ == "__main__":
    # Create a problem instance
    problem = GrossPitaevskiiProblem(100, 100, 0.1, 1, 1, 0.1, 100)
    print(problem)
