"""
Step 11: Boltz-2 protein-ligand structure predictions for top drug candidates.

Generates YAML input files for each drug-protein pair and runs Boltz-2 via
uvx (isolated Python 3.12 environment — boltz pins biopython==1.84 which
requires Python <3.13).

Boltz-2 is MIT-licensed and supports commercial use. It predicts protein-ligand
co-folded structures and binding affinity, making it suitable for the structural
validation panel in the p01 manuscript.

Usage (generate inputs + run all predictions):
    uv run python scripts/11_boltz_predict.py

Usage (generate inputs only, run manually):
    uv run python scripts/11_boltz_predict.py --prepare-only

Usage (run a single drug):
    uv run python scripts/11_boltz_predict.py --drug nintedanib

Output:
    results/boltz/inputs/<drug>_<target>.yaml      — Boltz-2 input files
    results/boltz/predictions/<drug>_<target>/     — Boltz-2 output per job
      boltz_results_<drug>_<target>/
        predictions/<drug>_<target>_model_0.cif   — predicted structure
        predictions/<drug>_<target>_model_0.json  — confidence scores

Hardware:
    Runs on CPU (M-series Mac): ~5–15 min per prediction.
    Set BOLTZ_ACCELERATOR=gpu to use MPS on Apple Silicon (experimental).
"""

import argparse
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "shared/scripts"))
from pipeline_utils import find_repo_root, load_config

DEFAULT_CONFIG = "shared/configs/igan.yaml"

# Drug-protein pairs from final_candidates.tsv (top 5)
PAIRS = [
    {
        "drug": "nintedanib",
        "target": "LYN",
        "uniprot": "P07948",
        "rank": 1,
        "smiles": "CN1CCN(CC1)CC(=O)N(C)C2=CC=C(C=C2)N=C(C3=CC=CC=C3)C4=C(NC5=C4C=CC(=C5)C(=O)OC)O",
        "rationale": "Rank 1 GNN candidate — LYN is a direct Kiryluk 2023 GWAS locus",
    },
    {
        "drug": "bosutinib",
        "target": "LYN",
        "uniprot": "P07948",
        "rank": 2,
        "smiles": "CN1CCN(CC1)CCCOC2=C(C=C3C(=C2)N=CC(=C3NC4=CC(=C(C=C4Cl)Cl)OC)C#N)OC",
        "rationale": "Rank 2 — LYN-cluster convergence",
    },
    {
        "drug": "masitinib",
        "target": "LYN",
        "uniprot": "P07948",
        "rank": 3,
        "smiles": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC(=CS4)C5=CN=CC=C5",
        "rationale": "Rank 3 — LYN-cluster convergence",
    },
    {
        "drug": "tranilast",
        "target": "TGFB1",
        "uniprot": "P01137",
        "rank": 4,
        "smiles": "COC1=C(C=C(C=C1)C=CC(=O)NC2=CC=CC=C2C(=O)O)OC",
        "rationale": "Rank 4 — TGFB1 anti-fibrotic cluster",
    },
    {
        "drug": "dasatinib",
        "target": "LYN",
        "uniprot": "P07948",
        "rank": 5,
        "smiles": "CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)CCO",
        "rationale": "Rank 5 (substituted for adaprev) — classic LYN inhibitor",
    },
]

_SEQ_CACHE: dict[str, str] = {}


def fetch_uniprot_seq(uniprot_id: str) -> str:
    """Fetch canonical sequence from UniProt REST API (cached)."""
    if uniprot_id in _SEQ_CACHE:
        return _SEQ_CACHE[uniprot_id]
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    with urllib.request.urlopen(url) as r:
        lines = r.read().decode().splitlines()
    seq = "".join(l for l in lines if not l.startswith(">"))
    _SEQ_CACHE[uniprot_id] = seq
    return seq


def generate_yaml(pair: dict, output_dir: Path) -> Path:
    """Write a Boltz-2 YAML input file for one drug-protein pair."""
    job_name = f"{pair['drug']}_{pair['target']}"
    seq = fetch_uniprot_seq(pair["uniprot"])

    content = f"""\
# Boltz-2 input: {pair['drug']} vs {pair['target']} ({pair['uniprot']})
# Rank {pair['rank']} candidate — {pair['rationale']}
# MIT license — commercial use permitted

sequences:
  - protein:
      id: A
      sequence: "{seq}"
  - ligand:
      id: B
      smiles: "{pair['smiles']}"

version: 1
"""
    out_path = output_dir / f"{job_name}.yaml"
    out_path.write_text(content)
    return out_path


def run_boltz(yaml_path: Path, out_dir: Path, accelerator: str = "cpu") -> bool:
    """Run boltz predict for a single input YAML via uvx."""
    cmd = [
        "uvx", "--python", "3.12", "--from", "boltz",
        "boltz", "predict", str(yaml_path),
        "--out_dir", str(out_dir),
        "--model", "boltz2",
        "--accelerator", accelerator,
        "--use_msa_server",
        "--diffusion_samples", "3",
        "--output_format", "mmcif",
        "--override",
    ]
    print(f"  Running: boltz predict {yaml_path.name} (accelerator={accelerator})")
    result = subprocess.run(cmd, capture_output=False, text=True)
    return result.returncode == 0


def parse_confidence(pred_dir: Path, job_name: str) -> dict | None:
    """Extract confidence metrics from boltz output JSON (model_0 representative)."""
    # Boltz writes: <out_dir>/boltz_results_<name>/predictions/<name>/confidence_<name>_model_0.json
    json_path = (pred_dir / f"boltz_results_{job_name}" / "predictions"
                 / job_name / f"confidence_{job_name}_model_0.json")
    if not json_path.exists():
        return None
    with open(json_path) as f:
        return json.load(f)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    parser.add_argument("--prepare-only", action="store_true",
                        help="Only write YAML inputs; do not run boltz")
    parser.add_argument("--drug", help="Run only this drug (e.g. nintedanib)")
    parser.add_argument("--accelerator", default="gpu",
                        choices=["cpu", "gpu"],
                        help="Accelerator: 'gpu' uses MPS on Apple Silicon (default). "
                             "Fall back to 'cpu' only if MPS causes errors.")
    args = parser.parse_args()

    _, repo_root = load_config(args.config)
    base = repo_root / "papers/p01-gnn-repurposing/results/boltz"
    input_dir = base / "inputs"
    pred_dir = base / "predictions"
    input_dir.mkdir(parents=True, exist_ok=True)
    pred_dir.mkdir(parents=True, exist_ok=True)

    pairs = PAIRS
    if args.drug:
        pairs = [p for p in PAIRS if p["drug"] == args.drug]
        if not pairs:
            print(f"Unknown drug '{args.drug}'. Available: {[p['drug'] for p in PAIRS]}")
            sys.exit(1)

    print(f"Generating {len(pairs)} Boltz-2 input files...\n")
    yaml_paths = []
    for pair in pairs:
        p = generate_yaml(pair, input_dir)
        yaml_paths.append((pair, p))
        print(f"  Wrote {p.name}")

    if args.prepare_only:
        print(f"\nYAML inputs ready in {input_dir}")
        print("\nTo run all predictions:")
        print(f"  uv run python scripts/11_boltz_predict.py --accelerator cpu")
        return

    print(f"\nRunning Boltz-2 predictions (accelerator={args.accelerator})...")
    print("First run downloads model weights (~3 GB) — this takes a few minutes.\n")

    results = []
    for pair, yaml_path in yaml_paths:
        job_name = f"{pair['drug']}_{pair['target']}"
        print(f"\n[{pair['rank']}/5] {pair['drug']} vs {pair['target']}")
        ok = run_boltz(yaml_path, pred_dir, args.accelerator)
        conf = parse_confidence(pred_dir, job_name) if ok else None
        results.append({"pair": pair, "ok": ok, "confidence": conf})
        if conf:
            lig_iptm = round(conf.get("ligand_iptm", conf.get("iptm", 0)), 4)
            ptm      = round(conf.get("ptm", 0), 4)
            plddt    = round(conf.get("complex_plddt", 0), 4)
            print(f"  ligand_ipTM={lig_iptm}  pTM={ptm}  pLDDT={plddt}")

    print("\n=== Summary ===")
    print(f"{'Drug':<12} {'Target':<8} {'Rank':<6} {'Status':<8} {'ligand_ipTM':<14} {'pTM':<8} {'pLDDT'}")
    print("-" * 65)
    for r in results:
        p = r["pair"]
        status = "OK" if r["ok"] else "FAILED"
        conf = r["confidence"] or {}
        lig_iptm = round(conf.get("ligand_iptm", conf.get("iptm", 0)), 4)
        ptm      = round(conf.get("ptm", 0), 4)
        plddt    = round(conf.get("complex_plddt", 0), 4)
        print(f"{p['drug']:<12} {p['target']:<8} {p['rank']:<6} {status:<8} {str(lig_iptm):<14} {str(ptm):<8} {plddt}")

    print(f"\nStructures in: {pred_dir}")
    print("Open .cif files in ChimeraX or PyMOL for visualisation.")


if __name__ == "__main__":
    main()
