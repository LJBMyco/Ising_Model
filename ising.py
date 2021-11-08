"""
Ising Model using Glauber or Kawasaki dynamics
"""

import sys
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

np.random.seed(12345)

class IsingModel(object):

    def __init__(self, method, shape, temp, sweeps):

        self.method = method
        self.shape = shape
        self.temp = temp
        self.sweeps = sweeps
        self.create_lattice()

    """Create an NxN lattice with randomly assigned 0s and 1s"""
    def create_lattice(self):

        # Create NxN lattice
        self.lattice = np.zeros((self.shape, self.shape))
        for i in range(self.shape):
            for j in range(self.shape):
                # Randomly assign +1 or -1  to each lattice site
                random_number = np.random.random()
                if random_number >= 0.5:
                    self.lattice[i][j] = 1.0
                else:
                    self.lattice[i][j] = -1.0


    """Impliment periodic boundary conditions """
    def pbc(self, i):

        if (i > self.shape-1) or (i<0):
            i = np.mod(i, self.shape)
            return i
        else:
            return i

########## Calculate properties of lattice ##########

    """Calculate total magnetisation """
    def total_mag(self):

        self.mag = np.abs(np.sum(self.lattice))
        return self.mag

    """Calculate total energy"""
    def total_energy(self):

        e_tot = 0.0
        for i in range(self.shape):
            for j in range(self.shape):
                #Energy at each site is -J sum nearest neighbours (J=1 in this sim)
                e_ij = 0.0
                e_ij += self.lattice[self.pbc(i+1)][j]
                e_ij += self.lattice[i][self.pbc(j-1)]
                e_ij = -1*e_ij*self.lattice[i][j]
                e_tot += e_ij
        return e_tot

    """Calculate Heat Capacity"""
    def heat_capacity(self, e, e_sqr):
        return (e_sqr-e**2.0)/((self.temp**2.0)*(self.shape*self.shape))

    """Calculate Suseptability"""
    def chi(self, m, m_sqr):
        return  (m_sqr-m**2.0)/((self.temp*self.shape)**2.0)

########## Impliment Glauber dynamics #########

    """Use Glauber dynamics to pick which lattice site to change"""
    def glauber(self):

        i = np.random.randint(self.shape)
        j = np.random.randint(self.shape)
        self.glauber_energy(i,j)

    """Work out the energy difference caused by the change"""
    def glauber_energy(self, i, j):

        e_new = 0.0
        e_new += self.lattice[self.pbc(i-1)][j]
        e_new += self.lattice[self.pbc(i+1)][j]
        e_new += self.lattice[i][self.pbc(j+1)]
        e_new += self.lattice[i][self.pbc(j-1)]
        # Swap sign so energy change is just 2x the new energy
        self.delta_e  = 2*e_new*self.lattice[i][j]
        self.glauber_metropolis_test(i, j)

    """Run the Metropolis Test to determine if the spin should be flipped"""
    def glauber_metropolis_test(self, i, j):

        #If energy state if preferable make the flip
        if self.delta_e <= 0:
            self.lattice[i][j] *= -1.0

        #If not flip is accepted with prob exp(-deltaE/k_B T)
        else:
            random_number = np.random.random()
            probability = math.exp(-1.0*(self.delta_e/self.temp))
            if random_number <= probability:
                self.lattice[i][j] *= -1.0


########## Impliment Kawasaki dynamics #########

    """Use Kawasaki dynamics to pick two random lattice sites to swap """
    def kawasaki(self):

        i1 = np.random.randint(self.shape)
        i2 = np.random.randint(self.shape)
        j1 = np.random.randint(self.shape)
        j2 = np.random.randint(self.shape)

        #Keep doing random selection if the sites match
        while (i1==i2) and (j1==j2):
            i1 = np.random.randint(self.shape)
            i2 = np.random.randint(self.shape)
            j1 = np.random.randint(self.shape)
            j2 = np.random.randint(self.shape)

        #If the sites have the opposite spin swap them
        if self.lattice[i1][j1] != self.lattice[i2][j2]:
            self.kawasaki_energy(i1, i2, j1, j2)

    """Work out the energy difference caused by the change"""
    def kawasaki_energy(self, i1, i2, j1, j2):

        #Calculate energy at both sites
        E1 = 0.0
        E1 += self.lattice[self.pbc(i1+1)][j1]
        E1 += self.lattice[self.pbc(i1-1)][j1]
        E1 += self.lattice[i1][self.pbc(j1+1)]
        E1 += self.lattice[i1][self.pbc(j1-1)]
        E1 = 2*E1*self.lattice[i1][j1]

        E2 = 0.0
        E2 += self.lattice[self.pbc(i2+1)][j2]
        E2 += self.lattice[self.pbc(i2-1)][j2]
        E2 += self.lattice[i2][self.pbc(j2+1)]
        E2 += self.lattice[i2][self.pbc(j2-1)]
        E2 = 2*E2*self.lattice[i2][j2]

        #Calculate delta_e
        self.delta_e = E1 + E2

        #Account for double counting nearest neighbours
        if (i1 == i2) and ((j1 == np.mod(j2+1, self.shape)) or (j1 == np.mod(j2-1, self.shape))):
             self.delta_e -= 2.0
        elif (j1 == j2) and ((i1 == np.mod(i2+1, self.shape)) or (i1 == np.mod(i2-1, self.shape))):
            self.delta_e -= 2.0

        self.kawasaki_metropolis_test(i1, i2, j1, j2)

    """Run the Metropolis Test to determine if the spin should be flipped"""
    def kawasaki_metropolis_test(self, i1, i2, j1, j2):

        #If energy state if preferable make the flip
        if self.delta_e <= 0:
            self.lattice[i1][j1] *= -1.0
            self.lattice[i2][j2] *= -1.0

        #If not flip is accepted with prob exp(-deltaE/k_B T)
        else:
            random_number = np.random.random()
            probability = math.exp(-1.0*(self.delta_e/self.temp))
            if random_number <= probability:
                self.lattice[i1][j1] *= -1.0
                self.lattice[i2][j2] *= -1.0

########## Calculate Errors ##########

    """Calculate standard error on the mean (used for E and Mag)"""
    def standard_error(self, values, values_sqr, list):
        return math.sqrt((values_sqr-values**2.0)/len(list))

    """Use bootstrap method to Calculate errors (used for Heat Capacity and Chi)"""
    def bootstrap_error(self, input, input_sqr, mode):

        #List for sample data poins
        values = []
        values_sqr = []

        #Randomly sample 101 times 
        for n in range(101):
            #Select indices
            samples = [np.random.randint(len(input)) for i in range(len(input))]
            #Select correspoinding data
            input_samples = [input[sample] for sample in samples]
            input_sqr_samples = [input_sqr[sample] for sample in samples]

            #Find means
            sample_mean = np.mean(input_samples)
            sample_sqr_mean = np.mean(input_sqr_samples)
            
            #Find values
            if mode == 'Heat Capacity':
                value = self.heat_capacity(sample_mean, sample_sqr_mean)
            elif mode == 'Chi':
                value = self.chi(sample_mean, sample_sqr_mean)

            values.append(value)
            values_sqr.append(value**2.0)
            
            
        #Find error
        mean = np.mean(values)
        mean_sqr = np.mean(values_sqr)

        error = math.sqrt(mean_sqr-mean**2.0)

        return error

########## Data Collection  ##########

    def data_collection(self):

        # Set up data collection lists
        cv = []
        chi = []
        cv_error = []
        chi_error = []
        e_mean = []
        mag_mean = []
        e_error = []
        mag_error = []

        #Create list of temps
        temp_range = np.arange(1.0, 3.25, 0.25)

        for T in reversed(temp_range):
            self.temp = T
            print(f'\n Temperature = {T} \n')
            #Only need to create an inital lattice, subsequent temps can use last sweep from previous temp
            if T == temp_range[-1]:
                self.create_lattice()

            #collect data for Errors
            e_values = []
            e_sqr_values = []
            mag_values = []
            mag_sqr_values = []

            #run for set number of sweeps
            for i in range(self.sweeps):
                #in each seep need to complete a test for every site in lattice
                for n in range(self.shape*self.shape):
                    if self.method == 'Glauber':
                        self.glauber()
                    elif self.method == 'Kawasaki':
                        self.kawasaki()
                if i ==0:
                    e = self.total_energy()
                    e_values.append(e)
                    e_sqr_values.append(e**2.0)
                    mag = self.total_mag()
                    mag_values.append(mag)
                    mag_sqr_values.append(mag**2.0)

                #Leave for 100 sweeps to equilibiriate
                if i >100:
                    #Every 10 sweeps collect data
                    if i%10==0:
                        e = self.total_energy()
                        e_values.append(e)
                        e_sqr_values.append(e**2.0)
                        mag = self.total_mag()
                        mag_values.append(mag)
                        mag_sqr_values.append(mag**2.0)

                        print(f'Sweep: {i}, Energy = {e}, Mag = {mag}')

            #At the end of all sweeps calculate means
            energy_mean = np.mean(e_values)
            energy_sqr_mean = np.mean(e_sqr_values)

            m_mean = np.mean(mag_values)
            m_sqr_mean = np.mean(mag_sqr_values)

            # Append to lists
            e_mean.append(energy_mean)
            mag_mean.append(m_mean)

            #Calculate Cv and Chi
            cv.append(self.heat_capacity(energy_mean, energy_sqr_mean))
            chi.append(self.chi(m_mean, m_sqr_mean))

            #Calculate errors
            e_error.append(self.standard_error(energy_mean, energy_sqr_mean, e_values))
            mag_error.append(self.standard_error(m_mean, m_sqr_mean, mag_values))
            cv_error.append(self.bootstrap_error(e_values, e_sqr_values, mode='Heat Capacity'))
            chi_error.append(self.bootstrap_error(mag_values, mag_sqr_values, mode='Chi'))

        #Write all data to files
        self.write_files(e_mean, e_error, mag_mean, mag_error, cv, cv_error, chi, chi_error, temp_range)

    """Write all data to .dat files """
    def write_files(self, e, e_error, mag, mag_error, cv, cv_error, chi, chi_error, T):

        energy_file = open('energy_output.dat', 'w')
        mag_file = open('mag_output.dat', 'w')
        cv_file = open('cv_output.dat', 'w')
        chi_file= open('chi_output.dat', 'w')

        for i in range(len(T)):
            energy_file.write(str(T[i]) + ',')
            energy_file.write(str(np.flip(e)[i]) + ',')
            energy_file.write(str(np.flip(e_error)[i]) + '\n')

            mag_file.write(str(T[i]) + ',')
            mag_file.write(str(np.flip(mag)[i]) + ',')
            mag_file.write(str(np.flip(mag_error)[i]) + '\n')

            cv_file.write(str(T[i]) + ',')
            cv_file.write(str(np.flip(cv)[i]) + ',')
            cv_file.write(str(np.flip(cv_error)[i]) + '\n')

            chi_file.write(str(T[i]) + ',')
            chi_file.write(str(np.flip(chi)[i]) + ',')
            chi_file.write(str(np.flip(chi_error)[i]) + '\n')

        energy_file.close()
        mag_file.close()
        cv_file.close()
        chi_file.close()

########## Run animation #########

    def update(self):

        for i in range(self.shape*self.shape):
            if self.method == 'Glauber':
                self.glauber()
            elif self.method == 'Kawasaki':
                self.kawasaki()

    def animate(self, i):
        self.update()
        self.mat.set_data(self.lattice)
        return self.mat,

    def run_animation(self):
        fig, ax = plt.subplots()
        self.mat = ax.imshow(self.lattice, cmap = 'seismic')
        ani = FuncAnimation(fig, self.animate, interval= 1, blit = False)

        plt.show()

if __name__ == "__main__":

    if len(sys.argv) != 6:
        print("Incorrect Number of Arguments Presented.")
        print("Usage: " + sys.argv[0] + "Method, Lattice Size, Temperature, Sweeps, Data/Animate")
        quit()
    elif sys.argv[1] not in ['Glauber', 'Kawasaki']:
        print("Please enter either Glauber or Kawasaki")
        quit()
    elif sys.argv[5] not in ['Data', 'Animate']:
        print('Please enter either Data or animate')
        quit()
    else:
        method = sys.argv[1]
        shape = int(sys.argv[2])
        temp = float(sys.argv[3])
        sweeps = int(sys.argv[4])

    ising_model = IsingModel(method, shape, temp, sweeps)
    if sys.argv[5] == 'Data':
        ising_model.data_collection()
    if sys.argv[5] == 'Animate':
        ising_model.run_animation()
