"""Reproducible LogP evaluation for MolecuLens.

Runs the QSAR training/validation pipeline and writes results/logp_eval.json with
full provenance (date, library versions, sample sizes). Run from the project root:

    python scripts/evaluate_logp.py

Commit the generated JSON so that every metric quoted in the docs can be traced to a run.
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

_HERE = Path(__file__).resolve().parent
ROOT = _HERE if (_HERE / "chemistry_engine").exists() else _HERE.parent
sys.path.insert(0, str(ROOT))

from chemistry_engine.molecular_intelligence import MolecularIntelligence  # noqa: E402


def main() -> int:
    bundle = MolecularIntelligence.train_and_evaluate_qsar_logp_model()
    if not bundle.get("available"):
        print("RDKit and NumPy/scikit-learn are required:", bundle.get("message_ar"))
        return 1

    versions = {"python": platform.python_version()}
    for mod in ("rdkit", "numpy", "sklearn"):
        try:
            versions[mod] = __import__(mod).__version__
        except Exception:
            versions[mod] = None

    skip = {"pipeline", "test_set_table", "_numpy_ref_mu", "_numpy_ref_std", "_numpy_ref_beta", "_numpy_ref_y_mean", "_numpy_xtx_inv", "_hybrid_res_predictor"}
    result = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "versions": versions,
        "dataset_citation": bundle["dataset_citation"],
        "metrics": {k: v for k, v in bundle.items() if k not in skip},
        "test_set_table": bundle["test_set_table"],
    }
    out = ROOT / "results" / "logp_eval.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    m = result["metrics"]
    print(f"n_train={m['n_train']} n_test={m['n_test']}")
    print(f"holdout   R2={m['test_r2']}  RMSE={m['test_rmse']}  (Crippen R2={m['crippen_test_r2']} RMSE={m['crippen_test_rmse']})")
    print(f"5-fold CV R2={m['cv_r2']}  RMSE={m['cv_rmse']}")
    print(f"class-out R2={m['group_cv_r2']}  RMSE={m['group_cv_rmse']}  (Crippen all R2={m['crippen_all_r2']} RMSE={m['crippen_all_rmse']})")
    print(f"written: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
