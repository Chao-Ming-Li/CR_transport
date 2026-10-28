# Python simulation analysis

These project guidelines are based on `output/spallation_check.ipynb` and
`output/adiabatic_loss_analysis.ipynb`, together with the owner's explicit
instructions. Apply them when creating or modifying simulation analyses.
Match nearby code and preserve existing notebook cells and results when making
focused changes.

## Basic syntax and readable calculations

- Follow basic Python syntax and grammar. Prefer straightforward assignments,
  explicit intermediate variables, indexing, slices, `for` loops, and simple
  functions. Avoid high-level commands, compact method chains, clever syntax,
  and unnecessary abstractions when basic operations express the calculation
  clearly.
- Ordinary NumPy array arithmetic, broadcasting, and familiar operations such
  as `np.sqrt`, `np.exp`, `np.logspace`, `np.sum`, and `np.histogram` fit this
  style. Do not replace clear array equations with loops merely to avoid NumPy.
- If a high-level command is necessary, provide enough nearby comments to
  explain **why it is necessary and what it really does**. Describe its inputs,
  outputs, and the mathematical or physical operation it performs; explain
  non-obvious arguments, conventions, and units. Naming the routine alone is
  not an explanation.
- For example, explain that `expm(Q * T_yr)` computes the matrix exponential
  solving a constant-rate coupled fragmentation system, including successive
  reactions. State the row/column convention and whether the input and output
  vectors represent particle counts or fractions. Explain why elementwise
  `np.exp` cannot replace this operation.
- Similarly, when using numerical ODE solvers, interpolation, optimization, or
  chained pandas operations, explain their purpose and actual operation. Use
  simpler direct calculations where practical; do not introduce these tools
  solely to shorten the code.

## Notebook organization

- Use a straightforward sequence of cells: imports and run parameters, binary
  loading and column extraction, physical reference calculation, counts or
  spectra, and comparison plots or residuals.
- Keep filenames, initial energy bounds, particle counts, propagation times,
  densities, and conversion constants visible near the beginning. Use explicit
  filenames for runs being compared rather than an elaborate discovery layer.
- Use `numpy as np` and `matplotlib.pyplot as plt`. Add other dependencies only
  when the analysis needs them; avoid unused imports.
- Write scientific equations directly with readable intermediate variables,
  following names such as `E_start`, `p_start`, `E_t1e8`, `theory_t1e8`, and
  `data_E_t1e8`. Retain established physical symbols and distinguish simulation
  arrays from theoretical predictions.
- Use small functions for calculations that are reused or evaluated for
  different parameters, such as `predicted_fractions(n_h)`. Keep one-off
  analysis steps directly in cells rather than building a framework.

## Data, comparisons, and plots

- Read binary output directly with `np.fromfile`, specify the writer's dtype,
  reshape to the documented record width, and extract columns explicitly.
  Check the C++ writer before assigning column meanings or energy units.
- Keep isotope order and mass numbers explicit. Distinguish total nuclear
  energy, total energy per nucleon, and kinetic energy. Treat ID `-1` as the
  untracked-fragment sink rather than an ordinary isotope.
- Make theoretical expectations visible as equations or explicit rate-matrix
  construction. Use the same bins for simulated and predicted spectra.
- Use direct `plt.subplots`, `step`, `scatter`, and simple loops over axes for
  comparisons. Label cases, physical quantities, and units, and select log
  scales explicitly where appropriate.
- Compute residuals directly, for example `(observed - predicted) / predicted`,
  with zero predictions handled when they actually occur.
- Preserve correct normalization: counts in equal logarithmic bins are
  proportional to `E*dN/dE`; divide by the logarithmic bin width when reporting
  that quantity numerically. Do not copy incidental labeling or normalization
  mistakes as style.
- Follow the reference notebooks' focused validation style. Add checks needed
  to prevent concrete reading or calculation failures, rather than broad
  unsolicited diagnostics. Report what was changed, what was checked, and any
  missing inputs or dependencies.
