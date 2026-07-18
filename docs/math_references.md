# Math references

Every formula in `attofab-core` is public-domain textbook physics. This
file lists, per solver, the equation used and where it comes from, so the
"public domain physics only" claim is auditable rather than asserted.

## Thermal oxidation - Deal & Grove model
Source: B.E. Deal and A.S. Grove, "General Relationship for the Thermal
Oxidation of Silicon," J. Appl. Phys. 36, 3770 (1965). Reproduced in any
standard VLSI process textbook (e.g. Plummer, Deen & Griffin, "Silicon
VLSI Technology").

    x_o^2 + A*x_o = B*(t + tau)
    x_o(t) = (A/2) * ( sqrt(1 + (t+tau)/(A^2/4B)) - 1 )

`A`, `B` are Arrhenius-temperature-dependent rate constants; `tau` offsets
the parabola for any oxide already present at t=0. Implementation:
`crates/core/src/oxidation.rs`. Rate constants (`materials.rs`) are
representative `<111>` Si textbook figures - not measured foundry data.

## Ion implantation - Gaussian range profile
Source: classic LSS (Lindhard-Scharff-Schiott) range-theory result, as
presented in standard ion-implantation texts.

    C(x) = Q / (sqrt(2*pi) * dRp) * exp( -(x-Rp)^2 / (2*dRp^2) )

`Rp` (projected range) and `dRp` (straggle) are physically set by ion
species, energy, and target material - classically from LSS theory or
SRIM/TRIM Monte Carlo range tables. Implementation:
`crates/core/src/implant.rs`. This engine's `estimate_range()` uses a
simple linear energy-scaling approximation, explicitly flagged as
illustrative - NOT a SRIM/TRIM substitute. Callers with real range data
should supply `Rp`/`dRp` directly.

## Thermal diffusion - Fick's Second Law
Source: A. Fick, "Ueber Diffusion," Annalen der Physik 94, 59 (1855);
standard finite-difference (FTCS) numerical treatment as taught in any
computational-methods course.

    dC/dt = D * d^2C/dx^2                          (1D)
    dC/dt = D * (d^2C/dx^2 + d^2C/dy^2)             (2D)

Solved via explicit finite differences with an auto-selected stable
timestep (CFL-like condition `D*dt/dx^2 <= 0.5`, generalized to 2D).
Dopant diffusivities `D(T) = D0 * exp(-Ea/kT)` are Arrhenius-form
textbook constants (`materials.rs`). Implementation:
`crates/core/src/diffusion.rs` (1D), `crates/core/src/grid2d.rs::anneal`
(2D).

## Etch & deposition - cellular automaton front propagation
No single named "textbook equation" - implemented as a geometric
nearest-neighbor front-propagation process (a simplified alternative to a
full level-set / Hamilton-Jacobi PDE solver), which correctly reproduces
the qualitative and geometric behavior real anisotropic (directional,
zero-undercut) and isotropic (uniform, undercutting) etches are taught to
exhibit in any intro semiconductor processing course. Implementation:
`crates/core/src/etch.rs`.

## Lithography - binary mask threshold
Simplified as a direct mask-to-photoresist projection (no diffraction
model yet). A Rayleigh-criterion-based aerial-image blur
(`R = k1*lambda/NA`, convolved with the mask before thresholding) is a
documented future extension for modeling proximity effects near mask
edges - see `docs/architecture.md`. Implementation:
`crates/core/src/lithography.rs`.

## Benchmarking against open PDKs
Where SkyWater 130nm PDK documentation publishes process parameters
(oxide thicknesses, junction depths) that overlap with what this engine
can compute, comparisons belong in a dedicated test
(`tests/physics_tests/skywater_comparison.rs`, planned) rather than mixed
into the general regression suite, since SkyWater's public parameter set
is limited and the comparison is an ongoing alignment effort, not a
one-time check.

## Update: corrected Deal-Grove constants + orientation dependence

An earlier internal draft of `materials.rs` used a set of Deal-Grove
pre-exponential constants that didn't quite match a more carefully sourced
aggregated table (dry B/A: corrected 3.71e6 um/hr from 6.23e6; wet B:
corrected 386 um^2/hr, Ea=0.78eV from 214 um^2/hr, Ea=0.71eV; wet B/A:
corrected 9.70e7 from 8.95e7). The regression tests pinned to the old
values were updated to the new, corrected outputs (documented inline in
`oxidation.rs` and `wafer1d.rs` with the reasoning) rather than silently
left passing against outdated physics. The original JS Phase 1 prototype's
`materials.js` was updated to match, so both engines now agree.

Also added: the (111)-plane linear-rate-constant scaling factor
(`B/A_111 ≈ 1.68 × B/A_100`), since the linear (reaction-limited) rate
depends on Si bond density at the interface, which differs by crystal
plane, while the parabolic (diffusion-limited) rate B does not.

## Status: Pearson-IV implant distribution (experimental, unverified)

`implant.rs::pearson4` implements the general Pearson differential
equation via numerical (RK4) integration - this machinery is correct and
tested (dose conservation, rejecting invalid Type IV moment combinations).
However, the moment-matching formula that converts (mean, variance,
skewness, kurtosis) into the ODE's (a, b0, b1, b2) parameters, sourced
from a general reference rather than the primary Elderton & Johnson
derivation, does NOT currently reproduce the requested moments correctly
(verified empirically by numerically recomputing the realized moments of
generated profiles and comparing against targets - see the two
`#[ignore]`d tests in that module for the specific failures). This is left
as a clearly-flagged, not-wired-in experimental module rather than either
silently shipping wrong physics or claiming a feature that doesn't work
yet. Gaussian (`implant_profile`) remains the verified default. Fixing
this requires validating the moment-matching formula against a reference
implementation (e.g. the R `PearsonDS` package) before re-enabling.

## Planned (Phase 3+, not yet implemented)

The following are scoped and documented here so contributors know the
target, but none of this is built yet:

- **Crank-Nicolson implicit diffusion**: unconditionally stable
  alternative to the current explicit FDM, needed once grids get fine
  enough (~1nm) that the explicit CFL limit forces impractically many
  timesteps. Requires a tridiagonal (Thomas algorithm) linear solve per
  step.
- **Dill's ABC model + Mack dissolution model**: replaces the current
  binary-threshold lithography model with photoresist exposure kinetics
  (`dI/dz = -I*dM/dt`, `dM/dt = -C*I*M`) and a proper development-rate
  curve (`R(M) = Rmax*(a+1)(1-M)^n / (a+(1-M)^n) + Rmin`).
- **Level-Set method**: replaces the current cellular-automaton etch/
  deposit with the Hamilton-Jacobi formulation (`dphi/dt + F|grad(phi)| = 0`)
  for accurate angled sidewalls and topological merging - the standard
  approach used by ViennaTS/ViennaLS.
- **HRLE (Hierarchical Run-Length Encoding)**: the memory structure that
  makes Level-Set practical at wafer scale - stores level-set values
  densely only near the material interface, run-length-compresses
  everything else.
- **GDSII/OASIS parsing**: via the `gds21`/`layout21` Rust crates, to
  drive lithography masks from real tapeout files (e.g. an OpenLane-
  synthesized standard-cell inverter).
- **SkyWater 130nm benchmark targets**: N-well peak ~6e17 cm^-3, P-well
  peak ~4e17 cm^-3, background ~8e14 cm^-3, standard gate oxide
  ~4.14-4.23nm, HP gate oxide ~1.5-1.6nm, LP gate oxide ~2.0nm - a
  dedicated `tests/physics_tests/skywater_comparison.rs` should assert
  AttoFab's engine converges on these when run with SkyWater-equivalent
  process parameters.
