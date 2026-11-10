import logging
import os
import sys
import pandas as pd
from config.settings import RAW_DATA_PATH, TARGET_INFO

logger = logging.getLogger(__name__)


def fetch_dpp4_data(output_path: str = RAW_DATA_PATH, limit: int = 1000) -> pd.DataFrame:
    target_id = TARGET_INFO["target_id"]
    try:
        from chembl_webresource_client.new_client import new_client
        activity = new_client.activity
        query = activity.filter(
            target_chembl_id=target_id,
            standard_type="IC50",
        )[:limit]
        records = list(query)
    except Exception as exc:
        msg = f"خطأ حرج: تعذر الاتصال بقاعدة بيانات ChEMBL لسحب البيانات الحقيقية ({exc}). تم إيقاف العملية لمنع استخدام بيانات اصطناعية."
        logger.error(msg)
        raise RuntimeError(msg) from exc

    if not records:
        msg = f"خطأ حرج: لم يتم استرجاع أي سجلات من ChEMBL للهدف {target_id}."
        logger.error(msg)
        raise RuntimeError(msg)

    df = pd.DataFrame(records)
    required_cols = {"molecule_chembl_id", "canonical_smiles", "standard_type", "standard_value", "standard_units"}
    if not required_cols.issubset(set(df.columns)) or df["canonical_smiles"].nunique() < 10:
        msg = "خطأ حرج: البيانات المسترجعة غير مكتملة أو لا تستوفي معايير الأصالة العلمية."
        logger.error(msg)
        raise RuntimeError(msg)

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)
    return df


if __name__ == "__main__":
    try:
        fetch_dpp4_data()
    except RuntimeError as err:
        print(str(err), file=sys.stderr)
        sys.exit(1)

