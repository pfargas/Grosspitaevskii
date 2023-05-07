import numpy as np


class GrossPitaevskiiProblem:
    """Class that generates a Gross-Pitaevskii problem instance, in harmonic oscillator units, for a 3D harmonic trap.
    """
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
        interaction = self.scattering_length*self.particle_number
        return 0
        
    @staticmethod
    def RungeKutta4(f, h, y0):
        """Method to solve a differential equation using the Runge-Kutta 4th order method

        Args:
            f (function): function to solve
            h (float): step size
            y0 (float): initial value

        Returns:
            float: solution
        """        
        k1 = f(y0)
        k2 = f(y0 + 0.5*h*k1)
        k3 = f(y0 + 0.5*h*k2)
        k4 = f(y0 + h*k3)
        return y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    @staticmethod
    def second_derivative(f,h):
        derivative=np.zeros(len(f))
        for i in len(f):
            if i == 0:
                derivative[i] = (f[i+1] - 2*f[i])/(h**2)
            elif i == len(f)-1:
                derivative[i] = (f[i-1] - 2*f[i])/(h**2)
            else:
                derivative[i] = (f[i+1] - 2*f[i] + f[i-1])/(h**2)
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
        n = len(f)-1
        integral = sum(2*f[i] if i % 2 == 0 else 4*f[i] for i in range(1,n))
        integral+= f[0] + f[n]
        return integral * h / 3

    def ansatz(self):
        psi=np.zeros(self.grid_size)
        r=np.linspace(0,self.grid_size,self.grid_step)
        cvar=2*np.sqrt(self.sigma)**3/np.sqrt(np.sqrt(np.pi))
        for i in range(self.grid_size):
            x = self.grid_step*i
            psi[i] = cvar*r[i]*np.exp(-0.5*self.sigma**2*r[i]**2)
        return psi
    
    def __str__(self) -> str:
        return f"{self.particle_number} Bosons in a spherical trap \n r-grid in {self.grid_size} points, r-step {self.grid_step} \n A0={self.scattering_length}, sigma={self.sigma} \n time={self.time_step}, number-iter={self.iterations}"



if __name__ == "__main__":
    # Create a problem instance
    problem = GrossPitaevskiiProblem(100, 100, 0.1, 1, 1, 0.1, 100)
    print(problem)