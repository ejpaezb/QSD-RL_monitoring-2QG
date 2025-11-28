import pickle
import numpy as np

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy.special import erf
from scipy.signal import savgol_filter
from scipy.fft import rfft, rfftfreq, irfft
from scipy.fft import fft, fftfreq, ifft

from pathlib import Path

colormap = mpl.cm.get_cmap('rainbow')

def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

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

def plot_pulse_sequence(sequences, labels, linestyles, colors, t_g, title, ylabel, xlabel, ylim=None, legend=False, legend_loc=0, mask=None):
    
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
        if mask == None:
            plt.plot(t, val, label=labels[idx], linestyle=linestyles[idx], color=colors[idx])    
        elif mask[idx]:
            plt.plot(t, val, label=labels[idx], linestyle=linestyles[idx], color=colors[idx])

    if ylim != None:
        plt.ylim(ylim)
        
    plt.ylabel(ylabel, fontsize=20)
    plt.xlabel(xlabel, fontsize=20)
    plt.title(title, fontsize=20)

    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)

    if legend:
        plt.legend(loc=legend_loc)
    
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

def sort_f_vs_d(input_file, output_file):

    fidelities = []
    indexes = []

    for idx, line in enumerate(input_file):
        if line.split()[0] == 'idx:':
            fidelities.append(line.split())
            indexes.append(int(line.split()[1]))

    fidelities.sort(key=lambda x: int(x[1]), reverse=False)
    indexes.sort(reverse=False)

    for idx in range(0, 99):
        if idx in indexes:
            output_file.write(' '.join(fidelities[indexes.index(idx)]) + '\n')
        else:
            output_file.write('idx: {} 0.0 (took: 0)\n'.format(idx))

def extract_f_vs_d_robustness(results_path):

    f_list = []
    d_list = []

    input_file = open(results_path, 'r')
    for idx, line in enumerate(input_file):
        if line.split()[0] == 'Offset:':
            f_list.append(np.float(line.split()[3]))
            d_list.append(np.float(line.split()[1]))

    return f_list, d_list

def config_generator(experiment_cfg, mode_cfg, qsd_simulation_cfg, output_path):
 
    gate_time = experiment_cfg['gate_time']
    simulation_mode = experiment_cfg['simulation_mode']

    cfg_list = [experiment_cfg, mode_cfg, qsd_simulation_cfg]
    cfg_path_list = ['experiment', 'mode-{}'.format(simulation_mode), 'qsd_simulation']

    for cfg in zip(cfg_list, cfg_path_list):

        input_cfg = open('./configs/default/{}.cfg'.format(cfg[1]), 'r')

        # index
        new_cfg = input_cfg.read() # .replace('?gate_time?', '{}'.format(gate_time))

        # Parse the config file
        for key, val in cfg[0].items():
            new_cfg = new_cfg.replace('?{}?'.format(key), '{}'.format(val))
        
        Path('{}/{}'.format(output_path, gate_time)).mkdir(parents=True, exist_ok=True)
        output_cfg = open('{}/{}/{}.cfg'.format(output_path, gate_time, cfg[1]), 'w')
        output_cfg.write(new_cfg)

def pulse_fft(t_g, pulse, sample_rate, n_repeats=1, filter=False):
    """
    t_g [us]
    pulse [Mrad/s]
    sample_rate [MHz]
    """
       
    _pulse = []
    for _ in range(n_repeats):
        _pulse.append(pulse)
    
    _pulse = np.concatenate(_pulse)
    _pulse = _pulse / (2.0 * np.pi) # _pulse [MHz]
    
    _t_g = t_g * n_repeats

    # Number of samples
    n_samples = sample_rate * _t_g 

    t = np.linspace(0.0, _t_g, n_samples, endpoint=False)

    sequence = []
    for time in t:
        sequence.append(pulse_const(_pulse, _t_g, time))
    
    if filter:
        sequence = savgol_filter(sequence, 11, 7)

    # Note the extra 'r' at the front
    yf = fft(sequence)
    xf = fftfreq(n_samples, 1 / sample_rate)

    # plt.plot(t, np.array(sequence) * 1000, '-', label="Global optimization (3,4)")
    # plt.plot(t, np.array(sequence) * 1000, '-', label="State of the art (3,4)")

    # plt.xlabel("$t$ ($\mu$s)", fontsize=20)
    # plt.ylabel("$\Omega$ (kHz)", fontsize=20)
    # plt.legend()
    
    # plt.savefig("pulse_15_200us(stoa).pdf", bbox_inches="tight")
    # plt.savefig("pulse_15_200us(go).pdf", bbox_inches="tight")
    # plt.show()
    
    # plt.plot(xf * 1e3, np.abs(yf), 'o-', label="Global optimization (3,4)")
    # plt.plot(xf * 1e3, np.abs(yf), 'o-', label="State of the art (3,4)")

    # plt.ylabel('Power', fontsize=20)
    # plt.xlabel('Frequency (kHz)', fontsize=20)
    # plt.legend()
    
    # plt.ylim([0, 5])
    # plt.xlim([20, 300])
    # plt.xlim([10, 150])

    # plt.savefig("pulse_15_200us_freq(go).pdf", bbox_inches="tight")
    # plt.savefig("pulse_15_200us_freq(stoa).pdf", bbox_inches="tight")
    # plt.show()
    
    
    # new_sig = irfft(yf)
    # plt.plot(t, new_sig * 1000)
    # plt.show()
    
    return xf, yf

def mrad2mhz(mrad):
    return mrad / (2.0 * np.pi)

def mrad2khz(mrad):
    return mrad / (2.0 * np.pi) * 1000

def khz2mrad(khz):
    return khz * (2.0 * np.pi) / 1000

def detuning_calc(mu1, mu2, r):
    
    det = (r * mu2 - mu1) / (r - 1.0)
    
    return det

def sort_f_vs_t(gate_times, res_path, pulse_list_path, verbose=False):

    flist = []
    idx_list = []
    max_rabi_list = []

    flist_2 = []
    idx_list_2 = []
    max_rabi_list_2 = []
    
    for gate_time in gate_times:

        if verbose:
            print('gate_time:', gate_time)
        
        input_file = open('{}/simulation_log_{}_us.txt'.format(res_path, gate_time), 'r')

        fidelities = []
        indexes = []

        for idx, line in enumerate(input_file):
            if line.split()[0] == 'idx:':
                fidelities.append(line.split())
                indexes.append(int(line.split()[1]))

        fidelities.sort(key=lambda x: float(x[2][1:-1]), reverse=True)
        # fidelities.sort(key=lambda x: float(x[2]), reverse=True)

        # Load pre-generated pulses
        pulse_list = np.load(pulse_list_path.format(gate_time))

        try:
            flist.append(np.float(fidelities[0][2][1:-1]))
            # flist.append(np.float(fidelities[0][2]))
            idx_list.append(int(np.float(fidelities[0][1])))
            max_rabi_list.append(np.max(np.abs(pulse_list[int(fidelities[0][1])])))
        except:
            flist.append(0.0)
            idx_list.append(0)
            max_rabi_list.append(0.0)

        flist_bests = []
        flist_bests_filtered = []
        for cnt in range(5): # Print the top pulses

            try:
                flist_bests.append(fidelities[cnt])
                
                if verbose:
                    max_rabi = np.max(np.abs(pulse_list[int(fidelities[cnt][1])]))
                    print(' '.join(fidelities[cnt]), 'Max Rabi: {} ({} kHz)'.format(max_rabi, max_rabi / (2.0 * np.pi) * 1000))#, pulse_list[int(fidelities[cnt][1])])
            except:
                if verbose:
                    print('list index out of range')
                
        # process the flist_bests
        if len(flist_bests) > 0:
            
            max_fidelity = np.float(flist_bests[0][2][1:-1])
            # max_fidelity = np.float(flist_bests[0][2])
            for cnt in range(len(flist_bests)):
                if (max_fidelity - np.float(flist_bests[cnt][2][1:-1])) / max_fidelity * 100 < 0.1:
                # if (max_fidelity - np.float(flist_bests[cnt][2])) / max_fidelity * 100 < 0.1:
                    flist_bests_filtered.append(flist_bests[cnt])

            flist_bests_filtered.sort(key=lambda x: np.max(np.abs(pulse_list[int(x[1])])), reverse=False)

            if verbose:
                print(len(flist_bests_filtered))
                print(flist_bests_filtered[0][2][1:-1])
                # print(flist_bests_filtered[0][2])

            flist_2.append(np.float(flist_bests_filtered[0][2][1:-1]))
            # flist_2.append(np.float(flist_bests_filtered[0][2]))
            idx_list_2.append(int(np.float(flist_bests_filtered[0][1])))
            max_rabi_list_2.append(np.max(np.abs(pulse_list[int(flist_bests_filtered[0][1])])))

        else:
            flist_2.append(0.0)
            idx_list_2.append(0)
            max_rabi_list_2.append(0.0)

    return flist_2, idx_list_2, max_rabi_list_2

def extract_error_bars(gate_times, ntraj, traj_path, verbose=False):

    f_list = []
    f_list_std = []

    for gate_time in gate_times:
        
        f_list_traj = []
        for idx in range(0, ntraj, 1):

            input_file_1 = open('{}/{}/filename_1_worker_{}'.format(traj_path, gate_time, idx), 'r')
            input_file_2 = open('{}/{}/filename_2_worker_{}'.format(traj_path, gate_time, idx), 'r')        
            
            last_line_1 = input_file_1.read().splitlines()[-1]
            last_line_2 = input_file_2.read().splitlines()[-1]
                    
            fidelity_1 = np.float(last_line_1.split()[1])
            fidelity_2 = np.float(last_line_2.split()[1])

            f_list_traj.append(np.max([fidelity_1, fidelity_2]))

        if verbose:
            print('gate_time =', gate_time)
            print('mean(f_list_traj) =', np.mean(np.array(f_list_traj)))
            print('std(f_list_traj) =', np.std(np.array(f_list_traj)))
            print('')

        f_list.append(np.mean(np.array(f_list_traj)))
        f_list_std.append(np.std(np.array(f_list_traj)))

    return f_list, f_list_std