# How ilreco reconstructs showers

This section describes the physics algorithms inside the library, for a
calorimetry or nuclear-physics reader — what each stage assumes, what it
computes, and where its numbers come from. No code knowledge is needed;
the API pages cover usage.

The input is one event: a list of fired cells (column, row, energy) of a
square-cell calorimeter. The output is a list of shower candidates with
energy, position, a fit-quality value, and a fiducial classification.

The algorithms assume **electromagnetic showers** — electrons, positrons,
and photons — because everything except the clustering itself is computed
against a measured *transverse shower profile*. Nothing forbids running
them on other deposits (a hadronic cluster, a MIP), but energies and
positions are then corrected with a profile that does not describe those
showers, and the fit quality loses its meaning. Where this section says
"shower" it means an electromagnetic one; the two-shower separation is
written for the π⁰ → γγ case that motivates it.

## The pipeline

| Stage | What it does | Page |
|---|---|---|
| 1. Island clustering | groups touching cells into connected clusters | [Clustering & peaks](./clustering) |
| 2. Peak finding | locates local maxima = shower candidates within a cluster | [Clustering & peaks](./clustering) |
| 3. Energy sharing | divides a multi-peak cluster's cells between its peaks using the shower profile | [Overlapping showers](./two-gamma) |
| 4. Position & energy fit | fits each (sub)cluster as one shower: moments + containment correction + metric minimization | [χ² fit](./chi2-method) / [log-likelihood fit](./log-likelihood) |
| 5. Second-order separation | tests whether a *single-peak* cluster is really two merged showers | [Overlapping showers](./two-gamma) |

Two ingredients feed everything: the **[transverse shower profile](./profiles)**
(mean deposited-energy fraction and its variance versus distance from the
shower axis) and, for the log-likelihood method, the
**[tabulated probability distributions](./log-likelihood)** of those
fractions. Both are detector-specific; producing them for your calorimeter
is described under Development
([profiles](./generating-profiles),
[likelihood tables](./generating-ll-tables)).

## The two reconstruction metrics

The historic method compares each cell with the *mean* profile, weighted by
a parametrized variance — a **χ²**. The log-likelihood method replaces the
mean-and-variance summary with the full measured probability distribution
of the cell's energy fraction — it asks "how probable is this cell content
for a shower here?" rather than "how far is it from the average?". The
metric enters twice: in the position fit (stage 4) and in the split
decision (stage 5). Energy reconstruction of isolated showers is common to
both.

Measured differences (20×20 PbWO₄, 2.05 cm cells, photons, full response
simulation): the log-likelihood fit gives 9–13 % better single-shower
position resolution at 1–10 GeV with a much smaller systematic pull toward
cell centers, recognizes more merged two-shower clusters at the same
false-split rate, and divides the energy of a merged cluster between the
two showers more accurately. Details and plots: the
[log-likelihood page](./log-likelihood) and
[overlapping showers](./two-gamma).

## Historic origin

The clustering, sharing, and χ² algorithms are Ilya Larin's PrimEx/HyCal
code, preserved expression-for-expression (results are pinned bit-exact by
regression tables). The log-likelihood method and the reconstructed
second-order separation are the 2026 additions; both are optional and off
by default.
