import numpy as np
import matplotlib.pyplot as plt
import os

def get_data(filename, header_rows=1, **kwargs):
    path_to_file = os.path.realpath(filename)
    data = np.genfromtxt(path_to_file, delimiter=",", skip_header=header_rows, **kwargs)
    if header_rows > 0:
        f = open(path_to_file, "r")
        params_str = f.readline()
        params = get_header_data(params_str)
        f.close()
        print params
        return data, params
    else:
        return data

def plot_timecourse(data):
    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax1.set_ylabel('x')
    ax2 = fig.add_subplot(212, sharex=ax1)
    ax2.set_ylabel('y')
    ax2.set_xlabel('t')
    maxtime = np.max(data[:,4])
    ax1.scatter(data[:,4], data[:,0], marker='.', c=(data[:,4]), cmap="Blues", lw=0)
    ax1.hold(True)
    ax1.scatter(data[:,4], data[:,5], marker='.', c=(data[:,4]), cmap="RdPu", lw=0)
    ax2.scatter(data[:,4], data[:,1], marker='.', c=(data[:,4]), cmap="YlOrRd", lw=0)
    ax2.hold(True)
    ax2.scatter(data[:,4], data[:,6], marker='.', c=(data[:,4]), cmap="RdPu", lw=0)
    plt.savefig("tc.png")

def plot_phaseplane(data):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    maxtime = np.max(data[:,4])
    ax.scatter(data[:,0], data[:,1], marker='.', c=(data[:,4]), cmap="Purples", lw=0)
    ax.hold(True)
    ax.scatter(data[:,5], data[:,6], marker='.', c=(data[:,4]), cmap="Greens", lw=0)
    plt.savefig("pp.png")

if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()
    #change after properly including header in data files
    for file in args.input_files:
        data = get_data(file, header_rows=0)
        plot_timecourse(data)
        plot_phaseplane(data)
