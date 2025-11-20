---
sidebar_position: 3
---
# Self-Consist Field

> Calculation type: `calculation = 'scf'`

This document provides a **standard self-consistent field (SCF)** input file template for `pw.x` in Quantum ESPRESSO.  
It includes both essential parameters and optional advanced parameters (commented out with explanations).

The SCF calculation is typically performed **after structure optimization** (`relax` or `vc-relax`) to obtain the self-consistent charge density and total energy for the optimized structure.

---

## 1. Standard SCF Template

```pw
&control
    calculation = 'scf',            ! self-consistent field calculation
    prefix      = 'MY_SYS',         ! name tag used for saving data
    outdir      = './tmp',          ! scratch directory for wavefunctions, etc.
    pseudo_dir  = './pseudo',       ! directory for pseudopotentials
    tstress     = .true.,           ! print stress tensor
    tprnfor     = .true.,           ! print atomic forces (optional for scf)
    verbosity   = 'high',           ! more detailed output

    ! disk_io    = 'none',          ! minimizes disk I/O; use only if enough RAM
    ! restart_mode = 'from_scratch',! start from scratch (default)
    ! forc_conv_thr = 1.0d-3,       ! threshold on forces (for consistency check)
    ! etot_conv_thr = 1.0d-5,       ! convergence threshold on total energy (Ry)
/

&system
    ibrav   = 0,                    ! free-form cell (use CELL_PARAMETERS)
    nat     = 4,                    ! number of atoms
    ntyp    = 2,                    ! number of atomic species

    ecutwfc = 60.0,                 ! kinetic energy cutoff for wavefunctions (Ry)
    ecutrho = 480.0,                ! cutoff for charge density (typically 4–8× ecutwfc)

    occupations = 'fixed',
    ! occupations = 'smearing',     ! smearing required for metals
    ! smearing    = 'mp',           ! Methfessel–Paxton scheme
    ! degauss     = 0.02,           ! smearing width (Ry)

    ! nbnd     = 40,                ! manually set number of bands (optional)

    ! input_dft = 'PBE',            ! explicitly define the XC functional (overrides pseudopotential header)
    ! vdw_corr = 'dft-d3',          ! include DFT-D3 dispersion correction
    ! dftd3_version = 3,            ! specify D3 method variant
    ! nspin = 2,                    ! enable spin-polarized calculation
    ! starting_magnetization(1) = 0.5, ! initial magnetic moment for species 1
    ! lda_plus_u = .true.,          ! enable DFT+U correction
    ! Hubbard_U(1) = 4.0,           ! on-site Coulomb interaction in eV (converted to Ry internally)
    ! ecutfock = 60.0,              ! Fock-exchange cutoff (for hybrid functionals)
    ! nqx1 = 2, nqx2 = 2, nqx3 = 2, ! q-grid sampling for Fock exchange
    ! assume_isolated = 'mp',       ! Makov–Payne correction for charged systems
    ! force_symmorphic = .false.,   ! disable symmetry enforcement (if necessary for distorted systems)
/

&electrons
    conv_thr    = 1.0d-8,           ! convergence threshold on total energy (Ry)
    mixing_beta = 0.3,              ! mixing factor; lower = more stable SCF
    mixing_mode = 'plain',          ! mixing algorithm (plain/Broyden/local-TF)
    diagonalization = 'david',      ! Davidson diagonalization (robust choice)
    electron_maxstep = 200,         ! maximum SCF iterations

    ! mixing_ndim = 8,              ! number of steps kept in mixing history
    ! conv_thr_init = 1.0d-5,       ! relaxed initial convergence criterion (used with adaptive convergence)
    ! adaptive_thr = .true.,        ! enable adaptive convergence thresholds
    ! startingwfc = 'file',         ! read starting wavefunctions from previous run
    ! startingpot = 'file',         ! read starting potential from previous run
    ! real_space_q = .true.,        ! use real-space FFT for hybrid functionals
/

ATOMIC_SPECIES
Si  28.0855  Si.pbe-n-kjpaw_psl.1.0.0.UPF
O   15.9994  O.pbe-n-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
  5.431   0.000   0.000
  0.000   5.431   0.000
  0.000   0.000   5.431

ATOMIC_POSITIONS angstrom
Si  0.000000  0.000000  0.000000
Si  2.715500  2.715500  2.715500
O   1.357750  1.357750  1.357750
O   4.073250  4.073250  4.073250

K_POINTS automatic
  8 8 8  0 0 0         ! Monkhorst–Pack grid for SCF; denser than relaxation step
```

⸻

## 2. Recommended Workflow and Notes

1. Purpose:
    - The SCF step provides a converged charge density (`prefix.save/charge-density.dat`).
    - This density is used for subsequent bands, nscf, and dos calculations.
2. Key Parameters:
    - `ecutwfc`, `ecutrho`: must be converged with respect to total energy and stress.
    - `K_POINTS`: use a denser grid than for relaxation.
    - `conv_thr`: tighten to 1.0d-8 or smaller for precise band/DOS analysis.
    - `smearing`: required for metals; for insulators, set `occupations = 'fixed'`.
    - `vdw_corr`: optional for weakly bonded solids or molecular crystals.
    - `nspin` and `starting_magnetization`: needed for magnetic systems.

3. Spin-Polarized SCF Example Snippet:

    ```3.
    &system
        nspin = 2,                            ! spin-polarized
        starting_magnetization(1) = 0.5,      ! initial guess for atom type 1
        starting_magnetization(2) = 0.0,
    /
    ```

4. DFT+U Example Snippet:

    ```4.
    &system
        lda_plus_u = .true.,                  ! enable DFT+U
        Hubbard_U(1) = 5.0,                   ! U parameter for atom type 1 (in eV)
    /
    ```

5. Hybrid Functional Example Snippet:

    ```5
    &system
        input_dft = 'HSE',                    ! hybrid functional
        ecutfock = 60.0,                      ! cutoff for Fock exchange
        nqx1 = 2, nqx2 = 2, nqx3 = 2,         ! q-grid for Fock term
    /
    ```

6. Restarting SCF:

    ```6
    &control
        restart_mode = 'restart',             ! restart from previous run
    /
    &electrons
        startingwfc = 'file',
        startingpot = 'file',
    /
    ```

---

## 3. Compact SCF Template for Automation

A minimal template suitable for scripting or batch generation.

```template
&control
    calculation = 'scf',
    prefix = 'PREFIX',
    outdir = './tmp',
    pseudo_dir = './pseudo',
/

&system
    ibrav = 0,
    nat = NAT,
    ntyp = NTYP,
    ecutwfc = ECUTWFC,
    ecutrho = ECUTRHO,
    occupations = 'smearing',
    smearing = 'mp',
    degauss = 0.02,
/

&electrons
    conv_thr = 1.0d-8,
    mixing_beta = 0.3,
    diagonalization = 'david',
/

ATOMIC_SPECIES
  ...   (auto-generated)

CELL_PARAMETERS angstrom
  ...   (auto-generated)

ATOMIC_POSITIONS angstrom
  ...   (auto-generated)

K_POINTS automatic
  KX KY KZ  0 0 0
```

---

## Summary

| Step | Purpose | Key Input Keywords | Notes |
| ---- | ---- | ---- | ---- |
| `relax` / `vc-relax` | Optimize atomic structure (and optionally lattice) | `ion_dynamics`, `cell_dynamics`, `tstress` | Must precede SCF |
| `scf` | Compute self-consistent charge density and total energy | `ecutwfc`, `ecutrho`, `K_POINTS`, `conv_thr` | Used for bands/DOS calculations |
| `nscf` | Non-self-consistent band/DOS sampling | `calculation='nscf'` | Uses charge density from SCF |
| `bands` / `dos` | Post-processing | `bands.x`, `dos.x`, `projwfc.x` | Read from SCF results |

This SCF template is designed for high compatibility across solid-state, molecular, and surface systems.
You may adapt it for hybrid DFT, DFT+U, spin-polarized, or dispersion-corrected calculations as needed.
