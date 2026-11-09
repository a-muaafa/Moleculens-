import os
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
import joblib

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def train_dpp4_model():
    print("جاري تحميل البيانات المعالجة لتوليد مصفوفة التدريب...")
    input_path = os.path.join(BASE_DIR, "data", "processed", "dpp4_cleaned.csv")
    
    if not os.path.exists(input_path):
        print("خطأ: الملف المعالج غير موجود!")
        return

    df = pd.read_csv(input_path)
    
    X = []
    y = []
    
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row['canonical_smiles']))
        if mol is not None:
            # توليد بصمات Morgan Fingerprints بحجم 1024 بت
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
            arr = np.zeros((1024,))
            DataStructs.ConvertToNumpyArray(fp, arr)
            X.append(arr)
            y.append(row['standard_value'])
            
    X = np.array(X)
    y = np.array(y)
    
    print(f"أبعاد مصفوفة المدخلات للنموذج: {X.shape}")
    
    # تقسيم البيانات (80% تدريب، 20% اختبار)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    
    print("جاري تدريب نموذج Random Forest...")
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X_train, y_train)
    
    # قياس الأداء
    y_pred = model.predict(X_test)
    r2 = r2_score(y_test, y_pred)
    print(f"تم التدريب بنجاح! دقة النموذج (R² Score): {r2:.3f}")
    
    # حفظ النموذج المدرب
    models_dir = os.path.join(BASE_DIR, "models")
    os.makedirs(models_dir, exist_ok=True)
    model_path = os.path.join(models_dir, "dpp4_rf_model.pkl")
    joblib.dump(model, model_path)
    print(f"تم حفظ النموذج في المسار: {model_path}")

if __name__ == "__main__":
    train_dpp4_model()
