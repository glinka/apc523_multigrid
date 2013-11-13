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
    plt.savefig("xytc.png")
    ax1.hold(False)
    ax2.hold(False)
    ax1.plot(data[:,4], data[:,7], c='b')
    ax2.plot(data[:,4], data[:,8], c='g')
    ax1.set_ylabel('energy')
    ax2.set_ylabel('angular momentum')
    ax2.set_xlabel('t')
    plt.savefig("eamtc.png")

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

def plot_sln_evo(data):
    fig = plt.figure();
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    ax.hold(True)
    n = data.shape[0]
    m = data.shape[1]
    x_vals = np.linspace(0, 1, n)
    labels = ["20 iters", "100 iters", "1000 iters", "moar iters", "analytical sln"]
    for i in range(m):
        ax.plot(x_vals, data[:,i], label=labels[i])
    ax.legend(loc=2)
    plt.savefig("jacobi_evo.png")

def plot_sln(data):
    fig = plt.figure();
    ax = fig.add_subplot(111)
    ax.set_xlabel('x')
    ax.set_ylabel('u')
    n = data.shape[0]
    x_vals = np.linspace(0, 1, n)
    ax.plot(x_vals, data)
    plt.savefig("multigrid.png")



if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input_files', nargs='+')
    args = parser.parse_args()
    #change after properly including header in data files
    for file in args.input_files:
        data = get_data(file, header_rows=0)
        plot_sln_evo(data)

