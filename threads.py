import threading
from grosspita import GrossPitaevskiiProblem
import csv
number_of_particles_list = [100,1000,10000,100000,1000000]
section_a = []
threads = []
for i,number_of_particles in enumerate(number_of_particles_list):
    problem = GrossPitaevskiiProblem(
                                particle_number=number_of_particles,
                                grid_size=10, 
                                grid_step=0.02, 
                                scattering_length=0.00433, 
                                sigma=0.5, 
                                time_step = 0.0001, 
                                iterations=20000, 
                                thomas_fermi=False,
                                interacting_system=True
                                )
    section_a.append(problem)
    threads.append(threading.Thread(target=section_a[i].evolution))


for i in range(len(number_of_particles_list)):
    threads[i].start()
    
for i in range(len(number_of_particles_list)):
    threads[i].join()
    
with open("prova.csv","w") as file:
    writer = csv.writer(file)
    writer.writerow(["number_of_particles","energy"])
    for i,number_of_particles in enumerate(number_of_particles_list):
        writer.writerow([number_of_particles,section_a[i].energy])