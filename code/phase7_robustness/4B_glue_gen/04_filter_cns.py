"""
Section 4B.4 — CNS-druglikeness filtering of generated glues
==========================================================
Gates (TODO.md Section 4):
  MW < 450  AND  cLogP in [2, 4]  AND  TPSA < 90  AND
  CNS-MPO(repo 5-check proxy) >= 4.0  AND  PAINS-free  AND  anchor retained.

CNS-MPO of record = the repo's 5-check proxy (code/slurm/26_fetch_tier2.py:133),
reproduced exactly. Wager-2016 6-parameter MPO is reported alongside as
secondary (uses a rule-based pKa estimate — flagged approximate).

Input:  results/phase7/glue_design/generated_raw.csv
Output: results/phase7/glue_design/generated_library.csv   (all, with gate flags)
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem import FilterCatalog

RDLogger.DisableLog("rdApp.*")

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
GD = PROJECT / "results/phase7/glue_design"
ANCHOR_CORE = "O=C1CCC(N2Cc3ccccc3C2=O)C(=O)N1"

pains = FilterCatalog.FilterCatalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)


def mpo_proxy(mw, logp, tpsa, hbd, rtb):
    """Repo 5-check CNS-MPO proxy (26_fetch_tier2.py); each check 0 / 0.5 / 1."""
    c = 0.0
    c += 1 if mw <= 360 else (0.5 if mw <= 500 else 0)
    c += 1 if logp <= 3.0 else (0.5 if logp <= 5.0 else 0)
    c += 1 if 40 <= tpsa <= 90 else (0.5 if tpsa <= 120 else 0)
    c += 1 if hbd == 0 else (0.5 if hbd <= 3 else 0)
    c += 1 if rtb <= 8 else (0.5 if rtb <= 10 else 0)
    return c


def _hump(x, lo_good, hi_good, lo_zero, hi_zero):
    if lo_good <= x <= hi_good:
        return 1.0
    if x <= lo_zero or x >= hi_zero:
        return 0.0
    if x < lo_good:
        return (x - lo_zero) / (lo_good - lo_zero)
    return (hi_zero - x) / (hi_zero - hi_good)


def mpo_wager(mw, logp, tpsa, hbd, pka_basic):
    """Wager 2016 6-parameter desirability sum (0-6). pka_basic is approximate."""
    d_logp = 1.0 if logp <= 3 else (0.0 if logp >= 5 else (5 - logp) / 2)
    d_logd = 1.0 if logp - 0.5 <= 2 else (0.0 if logp - 0.5 >= 4 else (4 - (logp - 0.5)) / 2)
    d_mw = 1.0 if mw <= 360 else (0.0 if mw >= 500 else (500 - mw) / 140)
    d_tpsa = _hump(tpsa, 40, 90, 20, 120)
    d_hbd = 1.0 if hbd <= 0.5 else (0.0 if hbd >= 3.5 else (3.5 - hbd) / 3)
    d_pka = 1.0 if pka_basic <= 8 else (0.0 if pka_basic >= 10 else (10 - pka_basic) / 2)
    return d_logp + d_logd + d_mw + d_tpsa + d_hbd + d_pka


def est_basic_pka(mol):
    """Crude: highest if aliphatic tertiary/secondary amine present, else low."""
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;!$(N=*);!$(N-C=O);!$(N-S=O);!$([n])][C]")):
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H0;!$(N=*);!$(N-C=O);!$([n])]")):
            return 9.5
        return 9.0
    if mol.HasSubstructMatch(Chem.MolFromSmarts("[NX3;H2][c]")):
        return 4.5   # aniline
    return 5.0


def main():
    core = Chem.MolFromSmarts(ANCHOR_CORE)
    df = pd.read_csv(GD / "generated_raw.csv")
    rows = []
    for _, r in df.iterrows():
        m = Chem.MolFromSmiles(r["smiles"])
        if m is None:
            continue
        mw = Descriptors.MolWt(m)
        logp = Crippen.MolLogP(m)
        tpsa = rdMolDescriptors.CalcTPSA(m)
        hbd = rdMolDescriptors.CalcNumHBD(m)
        hba = rdMolDescriptors.CalcNumHBA(m)
        rtb = rdMolDescriptors.CalcNumRotatableBonds(m)
        arom = rdMolDescriptors.CalcNumAromaticRings(m)
        fsp3 = rdMolDescriptors.CalcFractionCSP3(m)
        qed = QED.qed(m)
        pka = est_basic_pka(m)
        mpo5 = mpo_proxy(mw, logp, tpsa, hbd, rtb)
        wager = mpo_wager(mw, logp, tpsa, hbd, pka)
        pains_hit = pains.HasMatch(m)
        anchor_ok = m.HasSubstructMatch(core)

        g_mw = mw < 450
        g_logp = 2.0 <= logp <= 4.0
        g_tpsa = tpsa < 90
        g_mpo = mpo5 >= 4.0
        g_pains = not pains_hit
        passes = all([g_mw, g_logp, g_tpsa, g_mpo, g_pains, anchor_ok])
        # Pre-registered fallback: the glutarimide+isoindolinone warhead alone
        # has TPSA ~92 and cLogP ~ -0.4, so TPSA<90 is unreachable and MPO>=4 is
        # anchor-limited. Fallback: TPSA to the MPO half-credit band (<120),
        # MPO >=3.5; MW<450 and cLogP 2-4 unchanged.
        g_tpsa_fb = tpsa < 120
        g_mpo_fb = mpo5 >= 3.5
        passes_fallback = all([g_mw, g_logp, g_tpsa_fb, g_mpo_fb, g_pains, anchor_ok])

        rows.append({
            **r.to_dict(), "mw": round(mw, 2), "clogp": round(logp, 2),
            "tpsa": round(tpsa, 1), "hbd": hbd, "hba": hba, "rtb": rtb,
            "aromatic_rings": arom, "fsp3": round(fsp3, 3), "qed": round(qed, 3),
            "pka_basic_est": pka, "cns_mpo_proxy": mpo5,
            "cns_mpo_wager_approx": round(wager, 2),
            "gate_mw_lt450": g_mw, "gate_clogp_2_4": g_logp,
            "gate_tpsa_lt90": g_tpsa, "gate_mpo_ge4": g_mpo,
            "gate_pains_free": g_pains, "gate_anchor": anchor_ok,
            "gate_tpsa_lt120_fb": g_tpsa_fb, "gate_mpo_ge3p5_fb": g_mpo_fb,
            "passes_all_gates": passes, "passes_fallback": passes_fallback,
        })

    out = pd.DataFrame(rows)
    out.to_csv(GD / "generated_library.csv", index=False)
    n_pass = int(out["passes_all_gates"].sum())
    n_fb = int(out["passes_fallback"].sum())
    print(f"{len(out)} molecules; {n_pass} pass strict TODO gates; {n_fb} pass fallback")
    print("\ngate pass rates:")
    for g in ["gate_mw_lt450", "gate_clogp_2_4", "gate_tpsa_lt90",
              "gate_tpsa_lt120_fb", "gate_mpo_ge4", "gate_mpo_ge3p5_fb", "gate_pains_free"]:
        print(f"  {g:22s} {out[g].mean():.1%}")
    tier = "passes_all_gates" if n_pass >= 10 else "passes_fallback"
    n_use = out[tier].sum()
    if n_use:
        top = out[out[tier]].sort_values(
            ["cns_mpo_proxy", "qed"], ascending=False).head(12)
        print(f"\ntop of '{tier}' ({n_use}):\n",
              top[["glue_id", "linker", "mw", "clogp", "tpsa", "hbd", "cns_mpo_proxy",
                   "cns_mpo_wager_approx", "qed"]].to_string(index=False))
    print(f"\nWrote {GD}/generated_library.csv")


if __name__ == "__main__":
    main()
