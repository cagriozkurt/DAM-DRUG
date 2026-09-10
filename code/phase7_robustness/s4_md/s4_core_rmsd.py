"""
WP3.4 — Ligand core-RMSD analysis of the T-REMD trajectories
==========================================================
For each ternary complex and each of the 8 replicas: fit on the receptor
backbone, compute ligand core-RMSD on the 10% lowest-RMSF ligand heavy atoms
over the final 20 ns. Advancement gate: mean core-RMSD < 3.5 A across ALL
replicas (per TODO.md Section 4).

Run inside a GROMACS/MDAnalysis-capable env on TRUBA:
  python s4_core_rmsd.py

Outputs -> results/phase7/glue_md/
  core_rmsd_by_replica.csv
  core_rmsd_summary.csv
  CONCLUSION.md
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
except ImportError:
    mda = None

PROJECT = Path(os.environ.get("DAM_DRUG_DIR", str(Path.cwd())))
MD = PROJECT / "results/phase7/glue_md"
TOPCSV = PROJECT / "results/phase7/glue_design/glue_candidates_top.csv"
N_REP = 8
LAST_NS = 20.0
GATE = 3.5


def analyse_replica(tpr, xtc):
    u = mda.Universe(str(tpr), str(xtc))
    total_ns = (u.trajectory.n_frames - 1) * u.trajectory.dt / 1000.0
    start_frame = int(max(0, (total_ns - LAST_NS) / total_ns) * (u.trajectory.n_frames - 1))

    rec = u.select_atoms("protein and backbone")
    lig = u.select_atoms("resname LIG or resname UNL or resname MOL")
    if len(lig) == 0:
        lig = u.select_atoms("not protein and not resname SOL NA CL HOH")

    ref = u.copy()
    align.AlignTraj(u, ref, select="protein and backbone", in_memory=True).run()

    # RMSF of ligand heavy atoms over the analysis window -> pick 10% lowest
    lig_heavy = lig.select_atoms("not name H*")
    coords = []
    for ts in u.trajectory[start_frame:]:
        coords.append(lig_heavy.positions.copy())
    coords = np.array(coords)
    mean_pos = coords.mean(axis=0)
    rmsf = np.sqrt(((coords - mean_pos) ** 2).sum(axis=2).mean(axis=0))
    k = max(3, int(0.10 * len(rmsf)))
    core_idx = np.argsort(rmsf)[:k]

    ref_core = coords[0][core_idx]
    dev = np.sqrt(((coords[:, core_idx, :] - ref_core) ** 2).sum(axis=2).mean(axis=1))
    return float(dev.mean()), float(dev.std()), total_ns, k, len(lig_heavy)


def main():
    if mda is None:
        print("MDAnalysis not available — install it in the analysis env")
        return
    top = pd.read_csv(TOPCSV)
    gids = top["glue_id"].head(3).tolist()

    rows = []
    for gid in gids:
        wd = MD / gid
        for i in range(N_REP):
            tpr = wd / f"rep{i}" / "prod.tpr"
            xtc = wd / f"rep{i}" / "prod.xtc"
            if not xtc.exists():
                rows.append({"glue_id": gid, "replica": i, "status": "no trajectory"})
                continue
            mean_r, sd_r, ns, ncore, nheavy = analyse_replica(tpr, xtc)
            rows.append({"glue_id": gid, "replica": i, "status": "ok",
                         "sim_ns": round(ns, 1), "n_core_atoms": ncore,
                         "n_lig_heavy": nheavy,
                         "core_rmsd_A": round(mean_r, 2), "core_rmsd_sd": round(sd_r, 2)})
    by_rep = pd.DataFrame(rows)
    by_rep.to_csv(MD / "core_rmsd_by_replica.csv", index=False)
    print(by_rep.to_string(index=False))

    summ = []
    for gid in gids:
        g = by_rep[(by_rep["glue_id"] == gid) & (by_rep["status"] == "ok")]
        if len(g) == 0:
            summ.append({"glue_id": gid, "status": "incomplete"})
            continue
        mx = g["core_rmsd_A"].max()
        summ.append({
            "glue_id": gid, "n_replicas_ok": len(g),
            "mean_core_rmsd_A": round(g["core_rmsd_A"].mean(), 2),
            "max_replica_core_rmsd_A": round(mx, 2),
            "advancement_gate_pass": bool(mx < GATE and len(g) == N_REP),
        })
    sdf = pd.DataFrame(summ)
    sdf.to_csv(MD / "core_rmsd_summary.csv", index=False)

    lines = ["# WP3 — glue ternary-complex T-REMD stability: CONCLUSION\n",
             f"Advancement gate: mean ligand core-RMSD < {GATE} Å over the final "
             f"{LAST_NS:.0f} ns in ALL {N_REP} replicas.\n"]
    for _, r in sdf.iterrows():
        if r.get("status") == "incomplete":
            lines.append(f"- **{r['glue_id']}**: incomplete — trajectories missing.")
        else:
            v = "PASS" if r["advancement_gate_pass"] else "FAIL"
            lines.append(f"- **{r['glue_id']}**: {v} "
                         f"(mean {r['mean_core_rmsd_A']} Å, worst replica "
                         f"{r['max_replica_core_rmsd_A']} Å, {r['n_replicas_ok']}/{N_REP} replicas)")
    lines.append("\nManuscript: glues that PASS advance as kinetically-stable "
                 "ternary-complex scaffolds; those that FAIL are reported as "
                 "not stable under explicit-solvent T-REMD.")
    (MD / "CONCLUSION.md").write_text("\n".join(lines) + "\n")
    print("\n" + "\n".join(lines))


if __name__ == "__main__":
    main()
