# Transporter–CYP Zonation Simulator

A minimal compartmental model of the rat hepatic lobule that combines two
independently zonated processes — an uptake transporter and a CYP-mediated
metabolic sink — to predict spatial concentration profiles for a parent drug
*and* its metabolite across the portal-vein-to-central-vein axis.

**Code for:**
> *Transporter and CYP Zonation Leave Separable Signatures on Parent Drug
> versus Metabolite: A Minimal Hepatic Lobule Model* — submitted to
> *Biopharmaceutics & Drug Disposition*.

## Why this exists

Hepatic zonation is usually modeled for a transporter or an enzyme in
isolation. Real drugs go through both: an uptake transporter gets them into
the hepatocyte, and a CYP enzyme clears them once there — and the two are not
necessarily zonated the same way, or even in the same direction. This model
asks: when transporter zonation and CYP zonation are set independently, what
does the *parent drug* actually show across the lobule, versus its
*metabolite* — and does that combination give a practical way to tell whether
a given drug-drug interaction (DDI) is happening at the transporter or at the
enzyme, without invasive regional sampling?

## What's in this repo

- **`model.py`** — the model itself: a 15-compartment sinusoidal transport
  model (PV → CV; compartment count set from Ruijter et al. 2004
  stereological rat-liver data, see the docstring) with saturable transporter
  influx, reversible transporter efflux, a saturable CYP metabolism sink, and
  a metabolite compartment. Transporter and CYP capacity are each drawn
  independently from one of three zonation shapes (pericentral, non-zonated/
  flat, periportal), giving nine combinations. Two independent competitive-
  inhibition perpetrators (a transporter inhibitor and a CYP inhibitor) can
  be simulated. Pure dependency-free Python plus the standard library; run
  directly with `python3 model.py` to print the baseline transporter×CYP
  zonation grid.
- **`make_figure1.py`** — generates the manuscript's Figure 1 (schematic,
  parent/metabolite CV:PV grids across all nine zonation combinations, and
  the transporter-inhibitor-vs-CYP-inhibitor DDI discrimination panels).
  Requires `numpy` and `matplotlib`. Run with `python3 make_figure1.py`;
  output is `figure1.png` (included in this repo).
- **`sensitivity_suite.py`** — the six robustness/sensitivity checks
  reported in the manuscript's Figure 2 and "Robustness checks" section:
  (1) metabolite export zonated instead of uniform, (2) rank-order
  robustness across fm/Km,cyp/kef/flow sweeps, (3) DDI-signature robustness
  across fm × zonation configurations, (4) coarse-graining to a 3-zone
  experimental readout, (5) decoupling CYP capacity from kef (bypassing the
  fm parameterization), (6) CV:PV boundary-width (n_edge) sensitivity. Run
  directly with `python3 sensitivity_suite.py`. Requires `numpy`.
- **`make_figure2_sensitivity.py`** — generates the manuscript's Figure 2
  from `sensitivity_suite.py`. Requires `numpy` and `matplotlib`. Run with
  `python3 make_figure2_sensitivity.py`; output is `figure2_sensitivity.png`
  (included in this repo).
- **`requirements.txt`** — `numpy`/`matplotlib` versions used to generate
  the manuscript's figures (only needed for the sensitivity suite and
  figure scripts; `model.py` itself has no dependencies).

## Headline findings

- **Parent drug tracks transporter zonation; metabolite tracks CYP
  zonation — largely independently.** A non-zonated transporter combined
  with a pericentral CYP produces a periportal-biased parent signal
  alongside a strongly pericentral-biased metabolite signal — opposite
  spatial directions from the same simulated animal.
- **Co-localized zonation can mask the parent signal.** When transporter
  and CYP zonation are both pericentral, increasing the metabolic fraction
  suppresses the parent's spatial contrast even though the transporter
  itself remains sharply zonated.
- **This model reaches a genuine, non-trivial steady state** under constant
  portal inflow, because CYP metabolism here is a real, irreversible sink
  rather than a purely reversible transporter-efflux process, which would
  instead relax to a flat equilibrium for every zonation hypothesis.
- **Transporter- and CYP-mediated DDIs leave distinguishable signatures in
  most (8/9), but not all, tested configurations**: a transporter inhibitor
  collapses total exposure while barely reshaping parent spatial contrast;
  a CYP inhibitor raises parent exposure and substantially reshapes its
  spatial contrast. The one failure case (fm = 0.2, pericentral/pericentral)
  has a mechanistic explanation in the manuscript (Robustness checks) — a
  near-miss, not a reversal.

## What's robust, and what isn't

Six sensitivity/robustness checks (`sensitivity_suite.py`, Figure 2)
distinguish which claims above are generic model behavior versus contingent
on the specific illustrative parameter set:

- **Robust (rank order/direction):** parent-tracks-transporter and
  metabolite-tracks-CYP hold across a >10-fold sweep of fm, Km,cyp, kef, and
  flow; survive decoupling CYP capacity from kef entirely; survive
  coarse-graining to a 3-zone experimental readout; survive redefining the
  CV:PV boundary width (n_edge = 1/3/5); survive re-running at N = 24
  compartments (the compartment count used in earlier drafts).
- **Not robust (absolute magnitude):** metabolite CV:PV shifts by up to
  ~270-fold depending on whether metabolite export is itself zonated — a
  process this minimal model holds uniform by construction. Absolute
  metabolite gradients should not be over-interpreted quantitatively even
  where their direction is reliable.

## Status

This is a hypothesis-generating model: transporter/CYP capacities, the three
zonation shape families, and the illustrative `fm` (fraction metabolized)
values are chosen to be physiologically plausible, not fitted to a matched
kinetic dataset for a specific drug — see the manuscript's Limitations
section. All numerical results are dt-convergence-checked (explicit Euler,
dt = 0.0025 min, verified against a 16-fold finer step at baseline and at
every tested sensitivity-sweep extreme). A follow-up study, replacing these
illustrative parameters with real rat atorvastatin/cyclosporine kinetics
from two independent literature datasets, has since been completed; it
found the qualitative parent-tracks-transporter/metabolite-tracks-CYP
signature survives at real parameter scale (robust to parameter,
compartment-count, and inhibition-mechanism uncertainty), and separately
identified genuine numerical stiffness as a consequence of using real
transporter kinetics rather than illustrative ones.

## License

MIT — see [LICENSE](LICENSE).
