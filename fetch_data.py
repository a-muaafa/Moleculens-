import os
import sys
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def fetch_dpp4_data(limit=1000):
    output_dir = os.path.join(BASE_DIR, "data", "raw")
    output_path = os.path.join(output_dir, "dpp4_raw.csv")

    try:
        from chembl_webresource_client.new_client import new_client
        activity = new_client.activity
        query = activity.filter(
            target_chembl_id="CHEMBL284",
            standard_type="IC50"
        )[:limit]
        records = list(query)
    except Exception as exc:
        print(f"خطأ حرج: تعذر الاتصال بقاعدة بيانات ChEMBL لسحب البيانات الحقيقية ({exc}).", file=sys.stderr)
        print("تم إيقاف العملية لمنع استخدام أو توليد أي بيانات اصطناعية.", file=sys.stderr)
        sys.exit(1)

    if not records:
        print("خطأ حرج: لم يتم استرجاع أي سجلات من ChEMBL للهدف CHEMBL284.", file=sys.stderr)
        sys.exit(1)

    df = pd.DataFrame(records)
    required_cols = {"molecule_chembl_id", "canonical_smiles", "standard_type", "standard_value", "standard_units"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols or df["canonical_smiles"].nunique() < 10:
        print("خطأ حرج: البيانات المسترجعة غير مكتملة أو لا تستوفي معايير الأصالة العلمية.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(output_path, index=False)
    print(f"تم سحب وحفظ {len(df)} سجل حقيقي من ChEMBL في: {output_path}")


if __name__ == "__main__":
    fetch_dpp4_data()

