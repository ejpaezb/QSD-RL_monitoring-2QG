import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.special import erf

colormap = mpl.cm.get_cmap('rainbow')

def result_parser(filename):
    
    file = open(filename,'r')
    
    expectation_value = []
    for idx, line in enumerate(file): 
        if idx > 0:
            expectation_value.append(float(line.split()[1]) + 1.0j * float(line.split()[2]))
    
    return expectation_value

# def plot_pulse_sequence(sequence, gate_time, title, ylabel, xlabel, ylim=None):
    
#     # Plotting the pulse sequence    
#     fig = plt.figure()

#     pulse_width = gate_time/len(sequence)
#     for idx in range(len(sequence)):
#         plt.bar(pulse_width * idx, np.real(sequence[idx]), pulse_width, align='edge', color=colormap(np.real(idx/len(sequence))))        

#     if ylim != None:
#         plt.ylim(ylim)
        
#     plt.ylabel(ylabel, fontsize=20)
#     plt.xlabel(xlabel, fontsize=20)
#     plt.title(title, fontsize=20)

#     plt.xticks(fontsize=15)
#     plt.yticks(fontsize=15)

#     plt.show()

#     return fig

def plot_pulse_sequence(sequences, labels, linestyles, t_g, title, ylabel, xlabel, ylim=None, legend=False):
    
    # Plotting the pulse sequence    
    fig = plt.figure()

    _sequences = []

    t = np.linspace(0.0, t_g, 1000)
    for sequence in sequences:
        seq = []
        for time in t:
            seq.append(pulse_const(sequence, t_g, time))
        _sequences.append(np.array(seq))
            
    
    for idx, val in enumerate(_sequences):
        plt.plot(t, val, label=labels[idx], linestyle=linestyles[idx])

    if ylim != None:
        plt.ylim(ylim)
        
    plt.ylabel(ylabel, fontsize=20)
    plt.xlabel(xlabel, fontsize=20)
    plt.title(title, fontsize=20)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    if legend:
        plt.legend()
    
    # plt.show()

    return fig

def pulse_const(pulse, t_g, t):

    pulse_steps = len(pulse)
    
    idx = int(t / t_g * pulse_steps)
    if (idx == pulse_steps):
        idx -= 1
    
    return pulse[idx]
    
def pulse_erf(pulse, t_g, t):
    
    pulse_steps = len(pulse)

    dt = t_g / pulse_steps
    t = t - dt / 2
    
    idx = int(t / t_g * pulse_steps)
    if (idx == pulse_steps):
        idx -= 1
        
    t_i = idx * dt
    t_i_1 = (idx + 1) * dt
        
    if (idx == pulse_steps - 1):
        return pulse[idx]
    else:
        return (pulse[idx] + pulse[idx + 1]) / 2 + (pulse[idx + 1] - pulse[idx]) / 2 * erf(5.0 / dt * (t - (t_i + t_i_1) / 2))

def filter_pulses(pulse_list, threshold):
    for idx, val in enumerate(pulse_list):
        if np.max(np.abs(np.real(val))) < threshold:
            print(idx, np.max(np.abs(np.real(val))))

def save_pulse(pulse_list, filename):
    _pulse_list = np.real(np.sign(pulse_list)) * np.abs(pulse_list)
    np.save(filename, np.real(_pulse_list))