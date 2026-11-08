### Simulation (C++)

# QSD-RL monitoring for two-qubit gates

C++ quantum-trajectory simulators for trapped-ion gates, with a reinforcement-learning loop that designs the laser amplitude while a gate is running. Dynamics are integrated with the Quantum State Diffusion (QSD) library. A Python agent, embedded through pybind11, supplies the pulse.

The reinforcement-learning target is a three-ion Mølmer–Sørensen gate. The policy chooses each segment of a piecewise-constant Rabi pulse. The simulator returns motional quadratures during the gate and a processor fidelity at the end. Proximal policy optimization (PPO-Clip) is the agent selected in the shipped configuration; DDPG is also implemented.

## Layout

| Path | Contents |
| --- | --- |
| `include/qsd/` | QSD trajectory library (state, operators, stochastic integration) |
| `include/cnpy/` | NumPy `.npy` reader used by the C++ side |
| `src/ion-trap/` | Gate simulator. Builds `IonTrap-QSD` and `IonTrap-RL` |
| `src/ion-trap/ion-trap-rl/` | PPO, DDPG, and the Gymnasium environment called from C++ |
| `src/ion-trap/experiment/` | Ion-trap Hamiltonian, Lindblad operators, pulse, and fidelity |
| `src/cooling-dynamics/` | EIT and Doppler cooling. Builds `IonTrap-COOLING` |
| `src/dynamical-decoupling/` | Shelving fidelities. Builds `IonTrap-DD` |
| `src/SQUID/` | Circuit-QED example. Builds `CV-QSD` |
| `src/qsd-playground/` | Small QSD driver. Builds `QSD-Playground` |
| `data/ion-trap/configs/` | Boost INFO configs for gates, cooling, shelving, and RL |
| `utilities/src/` | Python helpers for pulses, Wigner functions, and trap Hamiltonians |

There is no top-level `CMakeLists.txt`. Each program is configured from its own directory under `src/`.

## Dependencies

C++ (C++17 for the ion-trap targets):

- CMake 3.6 or newer
- A C++17 compiler
- [Boost](https://www.boost.org/) (`filesystem`, `serialization`, `asio`, `property_tree`)
- [Eigen](https://eigen.tuxfamily.org/) 3.3 or newer
- [pybind11](https://github.com/pybind/pybind11) and Python 3 development headers, for `IonTrap-RL`
- NumPy C headers (`find_package(Python3 COMPONENTS Interpreter Development NumPy)`)

Python, for the agent imported by `IonTrap-RL`:

- `numpy`
- `torch`
- `gymnasium`

`src/ion-trap/CMakeLists.txt` appends `-L /opt/homebrew/lib` to the link line. On Linux, point that flag at the directory that holds the Boost libraries, or drop it if the linker already finds them.

Configs under `data/ion-trap/configs/` load trap frequencies, Lamb–Dicke parameters, detunings, and pulses from `.npy` files such as `data/ion-trap/trap-parameters/` and `data/ion-trap/pulse/`. Those arrays are not in this repository. A run needs them on the paths written in `trap.cfg` and `solutions.cfg`, or the manual lists in those files.

## Build

From the directory of the program you want:

```bash
cd src/ion-trap
cmake -S . -B build
cmake --build build --target IonTrap-RL
cmake --build build --target IonTrap-QSD
```

The same pattern applies to the other trees:

```bash
cd src/cooling-dynamics   && cmake -S . -B build && cmake --build build --target IonTrap-COOLING
cd src/dynamical-decoupling && cmake -S . -B build && cmake --build build --target IonTrap-DD
cd src/SQUID              && cmake -S . -B build && cmake --build build --target CV-QSD
cd src/qsd-playground     && cmake -S . -B build && cmake --build build --target QSD-Playground
```

## Run

Every ion-trap executable takes one argument: a directory of Boost INFO configs, with a trailing slash. The binary concatenates that prefix with names such as `experiment.cfg`.

```bash
cd src/ion-trap/build
./IonTrap-RL ../../../data/ion-trap/configs/default_rl/
./IonTrap-QSD ../../../data/ion-trap/configs/default_3/
```

`IonTrap-RL` imports `ppo.py` and `ion_trap_env.py` from `dirname(cwd) + py_dir`. With the shipped `mode-6.cfg` value `py_dir = "/ion-trap-rl/"`, the working directory must be a child of `src/ion-trap` (a `build/` directory there satisfies that). Launching the binary from another cwd will fail to import the agent.

Logs go to the directory in `experiment.log_dir` (for the RL config, `./cont_mon/simulation_log.txt`, relative to the working directory). Training also writes `rewards.txt`, `fidelities.txt`, and `actor.h5` under `mode-6.save_dir` (`./results/` in the shipped file).

### IonTrap-QSD modes

`experiment.simulation_mode` in `experiment.cfg` selects the run:

| Mode | Program |
| --- | --- |
| 0 | Single trajectory, fidelity from `solutions.cfg` |
| 1 | Fidelity versus detuning |
| 3 | Robustness (detuning, pulse scale, or gate time) |
| 5 | Playground |

### IonTrap-RL modes

`IonTrap-RL` only accepts `experiment.simulation_mode = 6`. Behavior is then set by `mode-6.mode` in `mode-6.cfg`:

| `mode-6.mode` | Behavior |
| --- | --- |
| 0 | Train. Sample actions from the Gaussian policy. Save `actor.h5` at the end, or when fidelity reaches `termination_condition` (0.9993). |
| 1 | Test. Load `save_dir/actor.h5` and act with the Gaussian mean. The shipped config is mode 1, and this repository does not contain a trained `actor.h5`. |
| 2 | Debug. One trajectory, print motional quadratures, do not call the agent. |

`mode-6.rl_agent.select` is `ppo` or `ddpg`. PPO hyperparameters in the shipped `default_rl` config: actor learning rate `1e-4`, critic learning rate `4e-4`, clip `ε = 0.185`, GAE `γ = 0.95` and `λ = 0.8`, 128 epochs, batch size 128, rollout length 500.

Other executables follow the same config-directory convention. `IonTrap-COOLING` uses mode 7 (EIT when `mode-7.processor` is 0, Doppler when it is 1) or mode 5. `IonTrap-DD` uses mode 8 (shelving) or mode 5. `CV-QSD` reads a single config file rather than a directory; see the usage string printed when it is started with no arguments.

## What the PPO loop does

One episode is one quantum trajectory of the gate in `default_rl`: three ions, pair `(0, 1)`, nine amplitude segments, gate duration `0.1 × 1000 = 100`.

1. The agent resets its observation to the zero vector.
2. At `t = 0`, `take_action` samples a Rabi amplitude and the callback installs it as segment 0.
3. On each later segment boundary the callback passes the joint quadratures `(x0, x1, x2, p0, p1, p2)` and the segment index to `observe_and_learn`, then asks for the next amplitude.
4. After the trajectory, processor fidelity is the larger of the `rho_ideal_1` and `rho_ideal_2` expectations. That call is the terminal step.
5. Every 500 stored transitions, the agent computes generalized advantage estimates and applies the PPO-Clip update. The actor loss is the clipped surrogate plus an entropy bonus. The critic loss is the squared error against the GAE return.

The per-step reward is the negative distance between a weighted motional modulus and `π / (4k)`, where `k` is the segment index. The terminal step adds a fidelity-dependent term defined in `IonTrapEnv.step`.

## Configuration files

A gate config directory typically contains:

| File | Role |
| --- | --- |
| `experiment.cfg` | Log path, trajectory count, seed, output operators, `simulation_mode` |
| `interaction.cfg` | Ion number, gate pair, Hamiltonian terms, Lindblad channels, pulse shape, initial state, time step |
| `trap.cfg` | Motional frequencies and Lamb–Dicke parameters |
| `solutions.cfg` | Detuning and the pulse used when the agent is not driving the amplitude |
| `mode-N.cfg` | Parameters for simulation mode `N` |

`interaction.pulse.AM.profile = 0` is piecewise constant. `steps` is the number of amplitude segments the agent fills.

## Credits

Ion-trap and circuit-QED drivers: Seyed Shakib Vedaie and Eduardo J. Páez (2019). Cooling and dynamical-decoupling drivers: Eduardo J. Páez (2023).

The QSD sources in `include/qsd/` are the stochastic-simulation library of Todd Brun and Rüdiger Schack. Their headers carry the GNU General Public License, version 2 or later. This repository has no separate top-level license file.

The PPO implementation in `src/ion-trap/ion-trap-rl/ppo.py` is adapted from [pg-is-all-you-need](https://github.com/MrSyee/pg-is-all-you-need), following Schulman et al., “Proximal Policy Optimization Algorithms,” arXiv:1707.06347 (2017). Generalized advantage estimation follows Schulman et al., arXiv:1506.02438 (2016).
