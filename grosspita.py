import numpy as np


class GrossPitaevskiiProblem:
    def __init__(self,particle_number, grid_size, grid_step, scattering_length, sigma, time_step, iterations):
        # Put default values in the declaration above
        self.particle_number = particle_number
        self.grid_size = grid_size
        self.grid_step = grid_step
        self.scattering_length = scattering_length
        self.sigma = sigma
        self.time_step = time_step
        self.iterations = iterations

    def integrate_problem(self):
        pass
        
    def ansatz(self):
        for i in range(self.grid_size):
            x = self.grid_step*i
            # self.psi[i] = np.exp(-((x - self.grid[self.grid_size//2])**2)/(2*self.sigma**2))
        return 
    
    def __str__(self) -> str:
        return f"{self.particle_number} Bosons in a spherical trap \n r-grid in {self.grid_size} points, r-step {self.grid_step} \n A0={self.scattering_length}, sigma={self.sigma} \n time={self.time_step}, number-iter={self.iterations}"



if __name__ == "__main__":
    # Create a problem instance
    problem = GrossPitaevskiiProblem(100, 100, 0.1, 1, 1, 0.1, 100)
    print(problem)