# myReichardtWallFunction

A custom OpenFOAM wall function for turbulent kinematic viscosity (`nut`) based on Reichardt's law of the wall. Compatible with OpenFOAM 10.

## Description

This boundary condition replaces the standard `nutWallFunction` on no-slip walls. It computes `nut` by iteratively solving for the friction velocity `u_τ` using **Reichardt's composite velocity profile**:

$$u^+ = \frac{1}{\kappa} \ln(1 + \kappa\, y^+) + C\!\left(1 - e^{-y^+/11} - \frac{y^+}{11}\,e^{-y^+/3}\right)$$

where:

| Symbol | Meaning |
|--------|---------|
| `u+` | non-dimensional velocity `U / u_τ` |
| `y+` | non-dimensional wall distance `y u_τ / ν` |
| `κ`  | von Kármán constant (inherited from base class, default 0.41) |
| `C`  | Reichardt's constant (default 7.8) |

Unlike the standard log-law, Reichardt's formula is **continuous from the wall** through the viscous sublayer, buffer layer, and log-law region. The friction velocity is found by Newton-Raphson iteration (up to 10 iterations, convergence tolerance 1 %).

Once `u_τ` is known, the wall value of `nut` is back-calculated as:

$$\nu_t = \frac{u_\tau^2}{|\partial U / \partial n|} - \nu$$

## Files

```
myReichardtWallFunctionFvPatchScalarField.H   # Class declaration
myReichardtWallFunctionFvPatchScalarField.C   # Implementation
Make/files                                    # wmake build list
Make/options                                  # Compiler/linker flags
```

## Compilation

1. Source your OpenFOAM environment (if not already done):

   ```bash
   source /opt/openfoam10/etc/bashrc
   # or wherever your OpenFOAM 10 installation lives
   ```

2. From the library source directory, run `wmake`:

   ```bash
   cd $WM_PROJECT_USER_DIR/src/myReichardtWallFunction
   wmake libso
   ```

   The shared library is installed to `$FOAM_USER_LIBBIN/libmyReichardtWallFunction.so`.

To clean build artefacts:

```bash
wclean
```

## Usage

### 1. Load the library in `system/controlDict`

```c++
libs
(
    "libmyReichardtWallFunction.so"
);
```

### 2. Apply the boundary condition in `0/nut`

```c++
<patchName>
{
    type            myReichardtWallFunction;
    value           uniform 0;
}
```

The optional parameter `C` can be used to override Reichardt's constant:

```c++
<patchName>
{
    type            myReichardtWallFunction;
    C               7.8;        // default; adjust if needed
    value           uniform 0;
}
```

### 3. Turbulence model compatibility

The wall function works with any RANS model that exposes a `nut` field through the `momentumTransportModel` interface (e.g. `kEpsilon`, `kOmegaSST`). No changes to `k` or `epsilon`/`omega` boundary conditions are required beyond the usual wall function settings.

## Reference

Reichardt, H. (1951). *Vollständige Darstellung der turbulenten Geschwindigkeitsverteilung in glatten Leitungen*. ZAMM, 31(7), 208–219.
