# Geometry Optimization

> Calculation type: `relax` / `vc-relax`

This document provides two standard Quantum ESPRESSO input templates for **geometry optimization** using `pw.x`.  
Both include commented advanced parameters with explanations for flexible and educational use.

---

## 1. Fixed Cell Relaxation (`calculation = 'relax'`)

Used when optimizing **atomic positions only**, keeping the lattice parameters fixed (e.g., experimental cell or already optimized structure).

```pw
&control
    calculation = 'relax',          ! ionic relaxation at fixed cell
    prefix      = 'SYSTEM_NAME',
    outdir      = './tmp',          ! scratch directory for wavefunctions, charge density, etc.
    pseudo_dir  = './pseudo',       ! directory containing UPF pseudopotentials
    nstep       = 200,              ! maximum number of ionic steps

    tstress     = .true.,           ! print stress tensor
    tprnfor     = .true.,           ! print forces

    ! disk_io   = 'none',           ! reduces disk I/O; requires enough memory
    ! verbosity = 'high',           ! verbose output, good for debugging
    ! restart_mode = 'restart',     ! restart from a previous run
/

&system
    ibrav   = 0,                    ! use CELL_PARAMETERS for full flexibility
    nat     = 4,                    ! number of atoms
    ntyp    = 2,                    ! number of atomic species

    ecutwfc = 60.0,                 ! kinetic energy cutoff for wavefunctions (Ry)
    ecutrho = 480.0,                ! cutoff for charge density (typically 4–8× ecutwfc)

    occupations = 'fixed',          ! for semiconductor or insulator
    ! occupations = 'smearing',     ! useful for metals or small-gap systems
    ! smearing    = 'mp',           ! Methfessel–Paxton smearing
    ! degauss     = 0.02,           ! smearing width (Ry)

    ! nbnd     = 40,                ! manually set number of bands (optional)

    ! input_dft = 'PBE',            ! specify XC functional explicitly
    ! vdw_corr = 'grimme-d3',       ! enable DFT-D3 dispersion correction
    ! london_s6 = 0.75d0,           ! scaling factor for vdW correction
    ! nspin    = 2,                 ! spin-polarized calculation
    ! starting_magnetization(1) = 0.5, ! initial magnetization for atom type 1
    ! lda_plus_u = .true.,          ! enable DFT+U correction
    ! Hubbard_U(1) = 4.0,           ! U parameter in eV for atom type 1
/

&electrons
    diagonalization = 'david',      ! Davidson diagonalization (default)
    mixing_mode = 'plain',          ! simple mixing mode
    mixing_beta = 0.3,              ! mixing factor; smaller = more stable
    
    ! conv_thr    = 1.0d-8,         ! convergence threshold for SCF (Ry), default presetted 1e-6
    ! electron_maxstep = 200,       ! max SCF iterations per ionic step
    ! mixing_ndim = 8,              ! history size for Broyden mixing
    ! startingpot = 'file',         ! read potential from file
    ! startingwfc = 'file',         ! read wavefunctions from file
/

&ions
    ! ion_dynamics = 'bfgs',        ! BFGS algorithm for ionic relaxation
    ! ion_dynamics = 'damp',        ! damped dynamics (alternative)
    ! trust_radius_max = 0.8d0,     ! maximum ionic step (Bohr)
    ! trust_radius_min = 1.0d-3,    ! minimum ionic step (Bohr)
    ! bfgs_ndim = 3,                ! number of BFGS history steps
/

&cell
    ! cell_dynamics = 'bfgs'        ! BFGS quasi-newton algorithm (default), ion_dynamics must be 'bfgs' too
    ! press = 0.0
    ! press_conv_thr = 0.5
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
  6 6 6  0 0 0         ! Monkhorst–Pack grid; test convergence
```

- Optimizes **only atomic positions**.
- `tstress` and `tprnfor` should be `.true.` to print stress and forces.
- Useful advanced options: DFT+U, vdW corrections, spin-polarization, and restart control.

---

## 2. Variable Cell Relaxation (`calculation = 'vc-relax'`)

Used to optimize both atomic positions and lattice parameters.
Essential for equilibrium structures and materials without known lattice constants.

```vc-relax
&control
    calculation = 'vc-relax',       ! optimize both ions and cell
    prefix      = 'MY_SYS',
    outdir      = './tmp',
    pseudo_dir  = './pseudo',
    nstep       = 200,

    tstress     = .true.,           ! required for variable-cell relaxation
    tprnfor     = .true.,

    ! disk_io   = 'none',           ! reduce I/O
    ! etot_conv_thr = 1.0D-4,       ! energy convergence between ionic steps (Ry)
    ! forc_conv_thr = 1.0D-3,       ! force convergence (Ry/Bohr)
/

&system
    ibrav   = 0,
    nat     = 296,
    ntyp    = 6,

    ecutwfc = 50.0,
    ecutrho = 400.0,

    occupations = 'fixed',
    ! occupations = 'smearing',
    ! smearing    = 'mp',
    ! degauss     = 0.02,

    vdw_corr    = 'dft-d3',         ! include dispersion correction
    dftd3_version = 3,              ! version of DFT-D3 correction

    ! input_dft = 'PBEsol',         ! choose alternative GGA (PBEsol for solids)
    ! ecutfock  = 60.0,             ! Fock exchange cutoff (for hybrid DFT)
    ! nqx1 = 2, nqx2 = 2, nqx3 = 2, ! q-grid for hybrid calculations
    ! nspin    = 2,                 ! spin-polarized vc-relax
    ! starting_magnetization(1) = 0.5,
/

&electrons
    conv_thr    = 1.0d-8,
    mixing_beta = 0.3,
    mixing_mode = 'plain',
    diagonalization = 'david',

    ! electron_maxstep = 200,
    ! startingpot = 'file',
    ! startingwfc = 'file',
/

&ions
    ion_dynamics = 'bfgs',
    ! ion_dynamics = 'damp',
    ! trust_radius_max  = 0.8d0,
    ! trust_radius_min  = 1.0d-3,
    ! bfgs_ndim         = 3,
/

&cell
    cell_dynamics = 'bfgs',         ! BFGS for cell optimization
    press         = 0.0,            ! target pressure (kbar)
    press_conv_thr = 0.5,           ! pressure convergence threshold (kbar)

    ! cell_dofree = 'all',          ! choose cell degrees of freedom:
    !                               !  'all'     : full relaxation (default)
    !                               !  'volume'  : fixed shape, volume only
    !                               !  'shape'   : shape at fixed volume
    !                               !  '2Dxy'    : in-plane only (for slabs)
    !                               !  'xyz'     : custom constraints
/

ATOMIC_SPECIES
C   12.011  C.pbe-n-kjpaw_psl.1.0.0.UPF
N   14.007  N.pbe-n-kjpaw_psl.1.0.0.UPF
H   1.0079  H.pbe-kjpaw_psl.1.0.0.UPF
Pb  207.2   Pb.pbe-dn-kjpaw_psl.1.0.0.UPF
I   126.90  I.pbe-n-kjpaw_psl.1.0.0.UPF
Cs  132.91  Cs.pbe-spn-kjpaw_psl.1.0.0.UPF

CELL_PARAMETERS angstrom
  ... (3×3 cell matrix)

ATOMIC_POSITIONS angstrom
  ... (atomic coordinates)

K_POINTS automatic
  3 3 3  0 0 0       ! coarse grid for large supercells
```

- Optimizes **both atomic positions and lattice vectors**.
- `tstress = .true.` is **mandatory** for `vc-relax`.
- press sets the **target external pressure** (in kbar).
- `cell_dofree` can restrict allowed cell degrees of freedom:
  - `all` — full optimization (default)
  - `volume` — isotropic scaling only
  - `shape` — fixed volume, shape change only
  - `2Dxy` — only in-plane cell relaxation
- Use `DFT-D3` or `PBEsol` for solids to improve accuracy.

---

## 3. Compact Reusable Skeleton Template

A simplified version for automation (e.g., scripting with ASE or SeekPath).  

```extra
&control
    calculation = 'vc-relax',       ! or 'relax'
    prefix      = 'PREFIX',
    outdir      = './tmp',
    pseudo_dir  = './pseudo',
    nstep       = 300,
    tstress     = .true.,
    tprnfor     = .true.,
/

&system
    ibrav   = 0,
    nat     = NAT,
    ntyp    = NTYP,
    ecutwfc = ECUTWFC,
    ecutrho = ECUTRHO,
    occupations = 'smearing',
    smearing    = 'mp',
    degauss     = 0.02,
    ! vdw_corr  = 'dft-d3',
    ! dftd3_version = 3,
    ! nspin    = 2,
/

&electrons
    conv_thr    = 1.0d-8,
    mixing_beta = 0.3,
    mixing_mode = 'plain',
    diagonalization = 'david',
/

&ions
    ion_dynamics = 'bfgs',
/

&cell
    cell_dynamics = 'bfgs',
    press         = 0.0,
    press_conv_thr = 0.5,
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

### Summary

- `relax`: atomic positions only.
- `vc-relax`: atomic positions + cell shape/volume.
- Typical relaxation workflow:
  - `vc-relax` (fully relaxed structure)
  - `scf` (self-consistent total energy)
  - `bands` / `nscf` / `dos` as post-processing.

These templates are designed to be directly usable and extendable for research or high-throughput workflows.

## 4. Take Result

```bash
awk '/Begin final coordinates/{flag=1} flag; /End final coordinates/{flag=0}' vc-relax.out > final.coords
```
