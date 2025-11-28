### Manning's 11-segment pulse
```md
hamiltonian   "FIVE_ION_KERR"
lindblads     "FULL"
stark_status  1
ion_pair      {0, 1}
time          {0.1, 10, 100}
delta_idx     64
seed          74298

ntraj 1  -> Result: 0.991716 (took: 123s)
ntraj 5  -> Result: 0.987278 (took: 578s)
ntraj 30 -> Result: 0.986543 (took: 3763s)
ntraj 50 -> Result: 0.987396 (took: 6281s)
```
### DE's 15-segment pulse [0.0, 2.0]
```md
hamiltonian   "FIVE_ION_KERR"
lindblads     "FULL"
stark_status  1
ion_pair      {0, 1}
time          {0.1, 10, 100}
delta_idx     64
seed          74298

ntraj 1  -> Result: 0.986269 (took: 120s)
ntraj 5  -> Result: 0.988027 (took: 591s)
ntraj 30 -> Result: 0.9858   (took: 3796s)
ntraj 50 -> 
```
### DE's 15-segment pulse [-2.0, 0.0]
```md
hamiltonian   "FIVE_ION_KERR"
lindblads     "FULL"
stark_status  1
ion_pair      {0, 1}
time          {0.1, 10, 100}
delta_idx     64
seed          74298

ntraj 1  -> Result: 0.993177 (took: 122s)
ntraj 5  -> Result: 0.984657 (took: 583s)
ntraj 30 -> 
ntraj 50 -> 
```
