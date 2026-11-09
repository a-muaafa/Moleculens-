import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def clean_and_featurize():
    print("جاري قراءة البيانات الخام...")
    input_path = os.path.join(BASE_DIR, "data", "raw", "dpp4_raw.csv")
    
    if not os.path.exists(input_path):
        print("خطأ: ملف البيانات الخام غير موجود!")
        return

    df = pd.read_csv(input_path)
    print(f"عدد السجلات الأصلية: {len(df)}")

    # 1. تنظيف القيم المفقودة وإزالة الصفوف غير التامة
    df = df.dropna(subset=['canonical_smiles', 'standard_value'])
    df['standard_value'] = pd.to_numeric(df['standard_value'], errors='coerce')
    df = df.dropna(subset=['standard_value'])

    print("جاري توليد البصمات الجزيئية (Morgan Fingerprints) باستخدام RDKit...")
    valid_rows = []
    
    for idx, row in df.iterrows():
        smiles = row['canonical_smiles']
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is not None:
            valid_rows.append(row)

    df_clean = pd.DataFrame(valid_rows)

    # 2. إنشاء مجلد data/processed لحفظ البيانات النظيفة
    out_dir = os.path.join(BASE_DIR, "data", "processed")
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, "dpp4_cleaned.csv")
    df_clean.to_csv(output_path, index=False)
    
    print(f"تم تنظيف البيانات بنجاح! عدد المركبات الصالحة: {len(df_clean)}")
    print(f"تم حفظ الملف المعالج في المسار: {output_path}")

if __name__ == "__main__":
    clean_and_featurize()
