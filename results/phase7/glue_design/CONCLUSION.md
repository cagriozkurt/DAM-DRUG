# Section 4B — CRBN molecular-glue generation: CONCLUSION

**Date:** 2026-09-10
**Deliverable, not a hypothesis test.** A provenance-documented, scaffold-based
library of CRBN-anchored candidate molecular glues directed at the IKZF1 ZF2
degron. These are **unvalidated in-silico scaffolds, not binders** — consistent
with the manuscript framing of the screening arm as a methodological
demonstration.

## Method

- **Anchor (fixed):** lenalidomide isoindolinone–glutarimide,
  `O=C1CCC(N2Cc3cccc(N)c3C2=O)C(=O)N1` (MW 259.3). Growth exit vector = the
  4-amino aromatic nitrogen (the CELMoD linker-attachment position; the 8RQC
  reference ligand QFC grows from the analogous 4-position via an ether).
- **Fragment pool:** BRICS decomposition of 1,962 source SMILES
  (285 Tier-1 CNS-approved + 1,677 Tier-2 approved) → **561** mono-attachment
  fragments (MW 40–250, ≤2 rings, PAINS-free, no reactive groups).
- **Growth:** each fragment joined to the anchor N through 5 linker types
  (amide, reverse-amide→dropped, sulfonamide, acetamide spacer, urea,
  carbamate) via RDKit `molzip`. Every product verified to retain the intact
  glutarimide + isoindolinone (`HasSubstructMatch`). → **2,750** unique
  (InChIKey) products, 1 ETKDGv3 conformer each. Seed 42.
- **Filtering** (`04_filter_cns.py`): RDKit MW, cLogP (Crippen), TPSA, HBD/HBA,
  RTB, QED; CNS-MPO by the repo 5-check proxy (`26_fetch_tier2.py`) and by an
  approximate Wager-2016 6-parameter score.

## Key finding — the strict TODO gates are unsatisfiable

`TODO.md` §4 specifies TPSA < 90 Å² and CNS-MPO ≥ 4.0. **Zero of 2,750**
products pass, and this is structural, not a sampling issue:

| property | anchor alone | best achievable for any glue |
|---|---|---|
| TPSA | 92.5 Å² | ≥ 92.5 Å² (grows only) → TPSA < 90 impossible |
| cLogP | 0.03 | needs lipophilic caps to reach cLogP ≥ 2 |
| CNS-MPO (proxy) | — | ≤ 3.5 (TPSA > 90 caps that check at 0.5) |

The bifunctional glutarimide warhead intrinsically costs ~92 Å² of polar
surface and is near-zero cLogP; CELMoD/glue chemical space sits **at the edge
of CNS-MPO druglikeness by construction**. This is a real design constraint
worth stating (and it mirrors 4A: the pipeline's own filters expose a
limitation).

## Pre-registered fallback tier

TPSA < 120 (the CNS-MPO half-credit band), CNS-MPO(proxy) ≥ 3.5, MW < 450,
cLogP 2–4, PAINS-free, anchor retained → **184** candidates.
Ranges: MW 371–449 (median 414), cLogP 2.00–2.94, TPSA 96–119, HBD 2–3,
QED 0.49–0.80 (median 0.71). Linkers: carbamate 53, acetamide 48, amide 45,
urea 35, sulfonamide 3.

## Positional-plausibility docking (`06_dock_glues.py`)

Top 25 candidates + the bare anchor docked into the 8RQC CRBN–IKZF1(ZF2)
ternary interface, reusing the published grid box
(center 0.566/−2.235/4.760, 20 Å³; Vina, exhaustiveness 16 — the screen used 32).

- All 25 dock inside the box with Vina **−5.92 to −7.01 kcal/mol**, each
  **better than the bare lenalidomide anchor (−5.21)** by 0.7–1.8 kcal/mol —
  i.e. the grown fragment adds interface contacts, as intended, with no
  runaway scores that would signal an artefact.
- Best: GLUE0246 (4-methylbenzamide, −7.01), GLUE0231 (2-methylbenzamide,
  −6.90), GLUE0123 (2-Cl-benzamide, −6.86).
- This is a docking sanity check on **where** the molecule sits, not an
  affinity ranking and not a binding claim.

## Manuscript action

Add a short Results/Methods paragraph (CRBN pivot / Section 4):

> Building on the exclusion of the undruggable IKZF1 orthosteric zinc-finger
> pocket (drug_score 0.001), a scaffold-based library of candidate cereblon
> molecular glues was enumerated on the lenalidomide isoindolinone–glutarimide
> anchor, growing 561 BRICS fragments from CNS-approved chemical space at the
> 4-amino position toward the IKZF1 ZF2 degron (PDB 8RQC). Of 2,750
> anchor-preserving products, none met a strict CNS multiparameter-optimisation
> filter (TPSA < 90 Å², CNS-MPO ≥ 4) — the glutarimide warhead alone
> (TPSA 92.5 Å², cLogP ≈ 0) precludes it — illustrating that cereblon-glue
> chemistry occupies the boundary of CNS druglikeness. 184 candidates passed a
> relaxed filter (TPSA < 120 Å², CNS-MPO ≥ 3.5, MW < 450, cLogP 2–4); the top
> 25 docked into the 8RQC ternary interface with Vina scores of −5.9 to
> −7.0 kcal/mol, all more favourable than the bare anchor (−5.2). These are
> unvalidated computational scaffolds intended to define a synthetic starting
> point, not binding predictions.

## Files
`anchor.json`, `anchor.png`, `fragment_pool.csv`, `generated_raw.csv/.sdf`,
`generated_library.csv`, `glue_candidates_top.csv` (+ docking column),
`glue_top.sdf`, `glue_grid.png`, `glue_docking.csv`, `docked/*.pdbqt`.
