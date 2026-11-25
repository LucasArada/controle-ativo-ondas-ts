# Active control of Tollmien-Schlichting waves (IFFC)

MATLAB code of the undergraduate thesis *Controle Ativo de Ondas de Tollmien-Schlichting em Camada Limite* (Mechatronics Engineering, Escola Politécnica da USP, 2026), by Lucas Arada Alves and Luiz Henrique Piffer Marques, advised by Prof. Marlon Sproesser Mathias.

A compressible Navier-Stokes solver simulates a two-dimensional flat plate boundary layer at Re = 1000 and Ma = 0.5. An oscillator at the wall generates a T-S wave, a pressure sensor upstream reads it, and an actuator downstream cancels it at a target sensor. The inverse feed-forward controller (IFFC) is identified from two simulation runs.

The flow solver (`nseq`, `initDerivatives`, `generateMesh`, `blasiusFlow` and the base of `runSolver` and `boundary`) was developed by Prof. Marlon Sproesser Mathias (Mathias, 2021). The oscillator, actuator, controller, calibration and post-processing were added in this work.

Tag `v1.0-monografia` is the version used for the results of the thesis.

## Requirements

MATLAB R2018b or later, no toolboxes.

## Quick start

All scripts run from the repository folder. Choose the case at the top of `runSolver.m`:

```matlab
run_number = 3;                     % 0, 1, 2 or 3
lx = 250;                           % Physical domain length: 250 or 500
actuator_profile = 'dipole';        % 'uniform' or 'dipole'
control_mode = 'IFFC_buffer_comp';  % Run 3: 'IFFC_phasor', 'IFFC_buffer' or 'IFFC_buffer_comp'
```

Full sequence:

1. `run_number = 1`, run `runSolver.m`: open loop, oscillator on, actuator off.
2. `run_number = 2`, run `runSolver.m`: actuator test tone, oscillator off. Repeat for each actuator profile.
3. Run `iffcCalibration.m` (same `lx` and `actuator_profile`): identifies the transfer functions and saves the controller parameters.
4. `run_number = 3`, run `runSolver.m`: closed loop with the chosen control mode.
5. Run `analyzeResults.m` for the attenuation and `generateGraphics.m` for the figures.

To go straight to step 4, `runSolver.m` uses the parameters identified in the thesis (`data/controller-params-<lx>-dipole.mat`) when no calibration file exists. Steps 1 and 2 are still needed for `analyzeResults.m` and `generateGraphics.m`.

Run 0 regenerates a base flow from the Blasius profile. The base flows used in the thesis are in `data/`.

## Files

| File | Role | Reads | Writes |
|---|---|---|---|
| `runSolver.m` | Case setup, RK4 time loop, sponge zones and virtual sensors | `data/baseflow-<lx>.mat`, controller parameters (run 3) | `nt_*.mat`, `history.mat` |
| `nseq.m` | Right-hand side of the compressible Navier-Stokes equations | | |
| `boundary.m` | Boundary conditions, oscillator and actuator | | |
| `controller.m` | Actuator command at each time step | | |
| `initDerivatives.m` | Compact finite differences and 10th order filter | | |
| `generateMesh.m` | Mesh and buffer zones | | |
| `blasiusFlow.m` | Blasius profile (run 0) | | |
| `iffcCalibration.m` | Phasors, transfer functions and IFFC synthesis | `history.mat` of runs 1 and 2 | `controller_params_<profile>.mat` |
| `analyzeResults.m` | Attenuation at the target sensor and reference contamination | `history.mat` of runs 1 and 3 | |
| `generateGraphics.m` | Sensor signals, wall pressure and u' along the plate | runs 1, 2 and 3 | `figures/lx<lx>/*.png` |

Each run writes to `results/lx<lx>/<case>/`:

- `nt_XXXXXXXXXX.mat`: fields `U`, `V`, `R`, `E`, mesh `X`, `Y` and time `t`, every `nsave` steps;
- `history.mat`: wall pressure perturbation at the reference and target sensors (`hist_p_ref`, `hist_p_target`) and actuator command (`hist_u`) at every step.

A run to t = 2000 takes 200 000 steps and about 0.9 GB of snapshots with `nsave = 100`. `generateGraphics.m` needs the snapshots of the last two periods; the other scripts only need `history.mat`.

## Thesis cases

Main results use `lx = 250`. The cases with the dipole actuator ended at t = 1500 (`tf = 1500` in `runSolver.m` and `N_periods = 5` in `analyzeResults.m`). The extended domain cases use `lx = 500`.

## Reference

Mathias, M. S. *Computational study of the hydrodynamic stability of gaps and cavities in a subsonic compressible boundary layer*. PhD thesis, Escola de Engenharia de São Carlos, Universidade de São Paulo, 2021. doi:10.11606/T.18.2021.tde-29032022-085318
