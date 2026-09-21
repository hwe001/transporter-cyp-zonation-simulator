# Transporter–CYP Zonation Simulator

A minimal compartmental model of the rat hepatic lobule that combines two
independently zonated processes — an uptake transporter and a CYP-mediated
metabolic sink — to predict spatial concentration profiles for a parent drug
*and* its metabolite across the portal-vein-to-central-vein axis.

**Companion code for:**
> *Transporter and CYP Zonation Leave Separable Signatures on Parent Drug
> versus Metabolite: A Minimal Hepatic Lobule Model* — submitted to
> *Biopharmaceutics & Drug Disposition*.

This model extends the transporter-only model in the companion repository
[oatp1a4-zonation-simulator](https://github.com/hwe001/oatp1a4-zonation-simulator)
(itself the companion tool for a Perspective on Oatp1a4 zonation, submitted to
*CPT: Pharmacometrics & Systems Pharmacology*) by adding a second, saturable
CYP metabolism flux and a tracked metabolite compartment.

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

- **`model.py`** — the model itself: a 24-compartment sinusoidal transport
  model (PV → CV) with saturable transporter influx, reversible transporter
  efflux, a saturable CYP metabolism sink, and a metabolite compartment.
  Transporter and CYP capacity are each drawn independently from one of three
  zonation shapes (pericentral, non-zonated/flat, periportal), giving nine
  combinations. Two independent competitive-inhibition perpetrators (a
  transporter inhibitor and a CYP inhibitor) can be simulated. Pure
  dependency-free Python plus the standard library; run directly with
  `python3 model.py` to print the baseline transporter×CYP zonation grid.
- **`make_figure1.py`** — generates the manuscript's Figure 1 (schematic,
  parent/metabolite CV:PV grids across all nine zonation combinations, and
  the transporter-inhibitor-vs-CYP-inhibitor DDI discrimination panel).
  Requires `numpy` and `matplotlib`. Run with `python3 make_figure1.py`;
  output is `figure1.png` (included in this repo).

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
  portal inflow — unlike the transporter-only companion model, whose true
  equilibrium is exactly flat for every hypothesis, because CYP metabolism
  here is a real, irreversible sink.
- **Transporter- and CYP-mediated DDIs leave distinguishable signatures**:
  a transporter inhibitor collapses total exposure while barely reshaping
  parent spatial contrast; a CYP inhibitor raises parent exposure and
  substantially reshapes its spatial contrast — offering a non-invasive way
  to localize a DDI's mechanistic locus.

## Status

This is a hypothesis-generating model: transporter/CYP capacities, the three
zonation shape families, and the illustrative `fm` (fraction metabolized)
values are chosen to be physiologically plausible, not fitted to a matched
kinetic dataset for a specific drug — see the manuscript's Limitations
section. A follow-up model, fitting real rat atorvastatin/cyclosporine
kinetics, is in scoping.

## License

MIT — see [LICENSE](LICENSE).
