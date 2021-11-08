import matplotlib.pyplot as plt

def Energy():

    Temp = []
    Energy = []
    Error = []
    EnergyFile = open('energy_output.dat', 'r')
    Lines = EnergyFile.readlines()
    for i in Lines:
        Tokens = i.split(',')
        Temp.append(float(Tokens[0]))
        Energy.append(float(Tokens[1]))
        Error.append(float(Tokens[2]))

    plt.figure(1)
    plt.errorbar(Temp, Energy, yerr=Error, ecolor='r', markersize=2, capsize=3, elinewidth=0.5, barsabove=True)
    plt.savefig('Temp_Energy.png')

    EnergyFile.close()

def Mag():

    Temp = []
    Mag = []
    Error = []
    MagFile = open('mag_output.dat', 'r')
    Lines = MagFile.readlines()
    for i in Lines:
        Tokens = i.split(',')
        Temp.append(float(Tokens[0]))
        Mag.append(float(Tokens[1]))
        Error.append(float(Tokens[2]))

    plt.figure(2)
    plt.errorbar(Temp, Mag, yerr=Error, ecolor='r', markersize=2, capsize=3, elinewidth=0.5, barsabove=True)
    plt.savefig('Temp_Mag.png')

    MagFile.close()

def C_v():

    Temp = []
    C_v = []
    Error = []
    C_vFile = open('cv_output.dat', 'r')
    Lines = C_vFile.readlines()
    for i in Lines:
        Tokens = i.split(',')
        Temp.append(float(Tokens[0]))
        C_v.append(float(Tokens[1]))
        Error.append(float(Tokens[2]))

    plt.figure(3)
    plt.errorbar(Temp, C_v, yerr=Error, ecolor='r', markersize=2, capsize=3, elinewidth=0.5, barsabove=True)
    plt.savefig('Temp_C_v.png')

    C_vFile.close()

def Chi():

    Temp = []
    Chi = []
    Error = []
    ChiFile = open('chi_output.dat', 'r')
    Lines = ChiFile.readlines()
    for i in Lines:
        Tokens = i.split(',')
        Temp.append(float(Tokens[0]))
        Chi.append(float(Tokens[1]))
        Error.append(float(Tokens[2]))

    plt.figure(4)
    plt.errorbar(Temp, Chi, yerr=Error, ecolor='r', markersize=2, capsize=3, elinewidth=0.5, barsabove=True)
    plt.savefig('Temp_Chi.png')

    ChiFile.close()

def main():

    Energy()
    Mag()
    C_v()
    Chi()

main()
