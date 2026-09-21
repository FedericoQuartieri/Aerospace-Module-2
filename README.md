# A two-temperature model for reacting hypersonic flows in OpenFOAM-13

An implementation of Park's two-temperature model in OpenFOAM-13, coupled with
the thermochemistry library [Mutation++](https://github.com/mutationpp/Mutationpp),
and its verification against the zero-dimensional (heat-bath) test cases of

> V. Casseau, R. C. Palharini, T. J. Scanlon, R. E. Brown,
> *A Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows,
> Part One: Zero-Dimensional Analysis*, Aerospace **3** (2016) 34.

The paper is available at
[applications/test/nonEqTTv/aerospace-03-00034-1.pdf](applications/test/nonEqTTv/aerospace-03-00034-1.pdf).

Tommaso Marchesini, Federico Quartieri, Daniele Salvi —
Politecnico di Milano, Department of Aerospace Science and Technology.

This README describes the code on the `revised` branch.

## Overview

The gas has two temperatures: `T` (translational-rotational) and `Tve`
(vibro-electronic). The solver advances the conserved variables (partial
densities, sensible energy `e`, vibro-electronic energy `eve`) and recovers the
two temperatures through Mutation++, which also provides the thermodynamic
properties.

- **`highEnthalpyThermo` thermophysical model**
  (`libhighEnthalpyThermophysicalModels.so`): the interface to Mutation++. It
  initialises `e` and `eve` from `(p, T, Tve)`, evaluates the source terms and,
  after each step, recovers `T`, `Tve` and `psi` from the energies.
- **Source terms** (`mutationSources.H`): Landau-Teller V-T exchange with the
  relaxation times of eqs. 9-17 of the paper, chemistry at Park's temperature
  `T^0.7 Tv^0.3`, and preferential or non-preferential chemistry-vibration
  coupling. The solver and the zero-dimensional programs use the same
  functions.
- **`shockThermo` solver** (`libshockThermo.so`): a `foamRun` module derived from
  `shockFluid` that solves the species, sensible-energy and `eve` equations.
- **Zero-dimensional programs** (`Test-N2`, `Test-N2N`, `Test-N2O2`,
  `Test-air5`): integrate the heat baths of the paper using Mutation++ alone.
- **`solverHeatBath` case**: the same heat baths solved with `shockThermo` on a
  single cell, to check the integration of the model into the solver.

## Results

All zero-dimensional cases of Section 3 of the paper are reproduced. For
reacting nitrogen, the Park configuration is reproduced; the CVDV model is not
implemented.

| figure | case | zero-dimensional program | solver |
|---|---|---|---|
| 3a, 3b | N2, heating and cooling | `Test-N2` | yes |
| 4 | N2 from 30000 K, with and without electronic energy | `Test-N2` | yes |
| 5 | N2 + N, non-reacting | `Test-N2N` | yes |
| 6 | N2 + O2, with and without V-V exchange | `Test-N2O2` | no |
| 7, 8 | N2 + N, reacting | `Test-N2N` | yes (fig. 7) |
| 9 | five-species air | `Test-air5` | no |

- The solver and the zero-dimensional programs agree to within 0.3 %,
  chemistry included, except for nitrogen at 30000 K (fig. 4), where the
  difference reaches 1.2 % in the first 100 ns.
- Against the curves of the paper, extracted from the vector figures of the
  pdf, final temperatures agree to within 0.5 %, the composition of the
  nitrogen mixtures to within 1.7 % and the composition of air to within 15 %.
- The largest remaining differences are in fig. 3b and fig. 6; they are
  discussed in the report.

## Requirements

- OpenFOAM-13 (OpenFOAM Foundation)
- Mutation++, built with CMake
- Python 3 with `numpy` and `matplotlib`, for the comparisons and plots
- TeX Live with `latexmk`, only to build the report

The easiest option is the Docker image in [docker/](docker/), which already
contains OpenFOAM-13 and Mutation++ at a pinned revision.

## Building

### With Docker

```sh
./docker/build.sh                                        # once (~15 min: builds Mutation++)
./docker/run.sh ./Allwmake                               # libraries and solver
./docker/run.sh applications/test/nonEqTTv/Allwmake      # zero-dimensional programs
./docker/run.sh                                          # shell in the container
```

The project is mounted at `/project`, and the compiled libraries are kept in
the `aero-m2-home` volume, so nothing is rebuilt at each start. See
[docker/README.md](docker/README.md) for details.

### Without Docker

```sh
source /opt/openfoam13/etc/bashrc       # OpenFOAM first
source etc/bashrc                       # then the project: POLIMI_* and MPP_* variables

(cd thirdParty && ./makeMutationpp)     # clones and builds Mutation++ in thirdParty/Mutationpp
./Allwmake
applications/test/nonEqTTv/Allwmake
```

The order of the two `source` commands matters: `etc/bashrc` relies on the
OpenFOAM variables. `makeMutationpp` fetches the latest Mutation++ revision,
whereas the Docker image uses `e8edf4f3a7f2bc22b4c62a23060dbe86967f1b3e`.

Libraries and executables are installed in `$FOAM_USER_LIBBIN` and
`$FOAM_USER_APPBIN`. `etc/bashrc` also copies `etc/codeTemplates/` to
`~/.OpenFOAM/13/`, which `dynamicCode` needs to accept the `rrho` base
thermophysical model.

## Reproducing the verification

Run the following commands in the container (`./docker/run.sh`) or in a shell
with the environment loaded.

**Zero-dimensional programs**, all figures and model variants, with the
comparisons:

```sh
applications/test/nonEqTTv/Allrun
```

Each case takes about one second. The results are written to
`applications/test/nonEqTTv/output/`: one csv per case, plus the plots
`<figure>.png` (temperatures) and `<figure>-n.png` (number densities).

**Solver on a single cell**, one figure at a time (`fig3a`, `fig3b`,
`fig4-noEl`, `fig4-el`, `fig5`, `fig7`):

```sh
applications/test/nonEqTTv/solverHeatBath/Allrun fig3a
applications/test/nonEqTTv/solverHeatBath/Allrun fig7 0.5   # second argument: Park exponent
python3 applications/test/nonEqTTv/compare.py fig3a
```

`Allrun` writes the initial conditions of the figure to `conditions`, runs
`blockMesh` and `foamRun`, and converts the probes to
`output/<figure>-solver.csv`. `compare.py` compares the paper, the
zero-dimensional program and the solver, and prints the maximum differences.
Fig. 7 (one million steps of 1 ns) takes about 26 minutes.

**A single case by hand**, for example fig. 7 with non-preferential C-V
coupling:

```sh
cd applications/test/nonEqTTv
MPP_DATA_DIRECTORY=$PWD/mutation-data-noElectronic \
    Test-N2N 30000 1000 1e-3 N2_Park 0.7 output/fig7-nonPref.csv paper nonPreferential
```

Run without arguments, each program prints its usage. `MPP_DATA_DIRECTORY`
selects the data: `mutation-data` contains the electronic levels of Table A2 of
the paper, `mutation-data-noElectronic` only the ground state.

## Model selection in the solver

The models are selected in the `highEnthalpyMutation` sub-dictionary of
`constant/physicalProperties`. The defaults are those of the paper, and the
`foamRun` log prints the models in use.

```
highEnthalpyMutation
{
    mixture                 air_5;
    stateModel              ChemNonEqTTv;
    thermodynamicDatabase   RRHO;
    mechanism               N2_Park;        // none: no chemistry
    relaxationTime          paper;          // eqs. 9-17; mutation: Mutation++ formula
    chemistryVibration      preferential;   // eq. 32; nonPreferential: eq. 31
    preferentialFactor      0.3;
    parkExponent            0.7;            // T_P = T^a Tv^(1-a), eq. 29
}
```

In the zero-dimensional programs the same choices are command-line arguments:
`Test-N2N` takes the Park exponent, the relaxation time and the C-V coupling,
`Test-N2O2` the relaxation time.

## Repository layout

```
Allwmake                          builds src/ and applications/modules/
etc/                              project environment (bashrc, config.sh/, codeTemplates/)
src/thermophysicalModels/
  multicomponentThermo/highEnthalpyMulticomponentThermo/
                                  highEnthalpyThermo model and mutationSources.H
  specie/thermo/rrho/             rrho, a copy of janafThermo used as the base thermo
applications/modules/
  shockFluid/                     OpenFOAM-13 density-based module, base of shockThermo
  shockThermo/                    two-temperature solver
applications/test/nonEqTTv/       verification against the zero-dimensional cases of the paper
  N2/ N2N/ N2O2/ air5/            one zero-dimensional program per mixture
  common/heatBath.H               energies, equilibrium and output times of the programs
  solverHeatBath/                 single-cell OpenFOAM case
  mutation-data/                  Mutation++ data with Tables A1 and A2 of the paper
  mutation-data-noElectronic/     the same, without excited electronic levels
  hystrath-data/                  QK mechanism for air, taken from hyStrath
  paper-data/                     curves extracted from the paper pdf
  output/                         csv files and plots
  compare.py                      paper / zero-dimensional / solver comparison
applications/test/mutation++/     preliminary tests, not used in the verification
applications/test/thermoMixturePark2T/
test/chemistry/                   one-temperature chemistry cases, not used in the verification
tutorials/shockThermo/shockTube/  1D shock tube, not yet valid (see Limitations)
thirdParty/                       build scripts for Mutation++ and GSL, Mutation++ Doxygen documentation
docker/                           OpenFOAM-13 + Mutation++ image
explainations/                    working notes and dependency diagrams
report/                           LaTeX report
```

## Documentation

- **Report**: [report/](report/). Build it with `make` inside `report/`;
  instructions in [report/README.md](report/README.md).
- **Working notes** (in Italian), in [explainations/](explainations/):
  - [milestone-1-VT-relaxation.md](explainations/milestone-1-VT-relaxation.md):
    thermophysical model, source terms and solver;
  - [milestone-2-heat-bath-paper.md](explainations/milestone-2-heat-bath-paper.md):
    figure-by-figure comparison, sources of every modelling choice,
    conservation checks;
  - [janaf-rrho-a-cosa-serve.md](explainations/janaf-rrho-a-cosa-serve.md):
    what the `rrho` base thermo does and what Mutation++ does instead;
  - [dipendenze/](explainations/dipendenze/): build, class and execution
    diagrams.
- **Paper curves**: [paper-data/README.md](applications/test/nonEqTTv/paper-data/README.md).
- **Data from hyStrath**: [hystrath-data/README.md](applications/test/nonEqTTv/hystrath-data/README.md).

## Limitations

- The verification covers only uniform heat baths, with no flow, on a single
  cell.
- The solver advances `eve` without the transport term, so the model cannot
  yet be used in a case with flow, including the `shockTube` tutorial.
- `Cp`, `Cv`, viscosity and thermal conductivity come from the `rrho` base
  thermo; they are computed at start-up and never updated. They play no role
  in the heat baths; in a flow they will have to be computed by Mutation++ at
  the two temperatures.
- Above 20000 K the `rrho` base thermo prints out-of-range warnings. They are
  harmless, because `T`, `Tve` and `psi` are computed by Mutation++.

The extensions needed to run on a real mesh are described in Sections 4 and 5
of the report.

## Third-party data

- The QK mechanism for air (`hystrath-data/`) comes from the
  [hyStrath](https://github.com/vincentcasseau/hyStrath) repository, commit
  `984e3000a5f8`, licensed under the GNU GPL v3.
- The reference curves in `paper-data/` are extracted from the figures of
  Casseau et al. (2016).
