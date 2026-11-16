# MolecuLens — وثيقة عقد واجهة البرمجة الموحد (Backend–Frontend API Contract v1.0.0)

تحدد هذه الوثيقة العقد البرمجي الموحد بين **الخادم العلمي (Backend: `core/api_contract.py` & `core/orchestrator.py`)** و**واجهة المستخدم (Frontend: `frontend/` أو أي عميل ويب مستقبلي)** لضمان التكامل السلس عند تسليم الواجهة أو تطويرها بشكل مستقل.

---

## 1. المبادئ المعمارية الأساسية (Core Architectural Principles)
1. **منع وضع منطق الحسابات العلمية داخل الفرونت إند:**
   - جميع الحسابات الكيميائية (الموازنة الجبرية، حفظ الذرات والشحنة، المتفاعل المحدد والفائض، المردود النظري، تحويل الصيغ الكيميائية وحساب الكتل المولية، حسابات `cLogP` ونموذج `Ridge QSAR`، ومراحل التجربة البصرية الأربع) تُحسب حصراً داخل الـ **Backend** وتُسلَّم جاهزة للعرض في حقل `data`.
2. **غلاف استجابة موحد (Unified Response Envelope):**
   - تُرجع جميع نقاط النهاية (`Endpoints`) في `MolecuLensAPI` نفس الهيكل القياسي دون استثناء:
```json
{
  "ok": true,
  "status_code": 200,
  "ui_state": "ready",
  "endpoint": "/api/v1/reactions/predict",
  "message_ar": "تم إجراء الحسابات العلمية واسترجاع البيانات الموثقة بنجاح.",
  "message_en": "Scientific calculations and verified data retrieved successfully.",
  "units": {
    "molar_mass": "g/mol",
    "amount": "mol",
    "mass": "g",
    "volume_stp": "L_STP",
    "temperature": "°C",
    "pressure": "atm",
    "enthalpy": "kJ/mol",
    "tpsa": "Å²",
    "concentration_molar": "M",
    "volume_ml": "mL"
  },
  "data": { ... },
  "error": null,
  "meta": {
    "api_version": "1.0.0",
    "timestamp_utc": "2026-11-16T17:30:00+00:00",
    "computed_in_backend": true,
    "cache_hit": false
  }
}
```

---

## 2. رموز حالة HTTP وحالات الواجهة (`status_code` & `ui_state`)

| رمز HTTP (`status_code`) | حالة النجاح (`ok`) | حالة الواجهة (`ui_state`) | المعنى والاستخدام العلمي |
|---|---|---|---|
| `200 OK` | `true` | `"ready"` | نجاح استرجاع المادة أو حساب التفاعل أو تقييم النموذج. |
| `400 Bad Request` | `false` | `"error"` | مدخلات ناقصة، فارغة، غير رقمية، تتجاوز الطول المسموح، أو تحتوي على وسوم HTML/JS مرفوضة أمنياً. |
| `404 Not Found` | `false` | `"empty"` | المادة المطلوبة غير موجودة في السجل المحلي أو المصادر الخارجية (مع إرجاع مقترحات `did_you_mean`). |
| `422 Unprocessable Entity` | `false` | `"empty"` أو `"error"` | مدخلات صحيحة شكلياً لكنها مرفوضة علمياً (مثل طلب حجم غازي `L_STP` لمادة صلبة/سائلة، أو غياب أدلة كافية لتوقع تفاعل دون اختلاق نواتج). |
| `503 Service Unavailable` | `false` | `"error"` | تعذر الوصول لخدمة خارجية عند طلب التحديث الشبكي الصريح. |

### رسائل حالات الواجهة القياسية (`UI_STATE_MESSAGES`):
- **`loading`**: `"جارٍ معالجة البيانات العلمية وحساب النتائج في الخادم..."`
- **`empty`**: `"لا توجد نتائج مطابقة للمدخلات الحالية. جرّب البحث بالاسم العلمي أو الصيغة الكيميائية."`
- **`error`**: `"تعذر إتمام العملية بسبب خطأ في المدخلات أو عدم كفاية الأدلة العلمية."`
- **`ready`**: `"تم إجراء الحسابات العلمية واسترجاع البيانات الموثقة بنجاح."`

---

## 3. توحيد أسماء الحقول والوحدات (Standardized Field Names & Units)

| اسم الحقل الموحد | النوع | الوحدة القياسية | الوصف |
|---|---|---|---|
| `symbol_or_formula` / `formula` | `str` | — | الصيغة الكيميائية أو الرمز القياسي (مثل `H2O`, `NaCl`, `Fe`). |
| `canonical_smiles` / `smiles` | `str` | — | التمثيل الخطي المعياري للجزيء (`Canonical SMILES`). |
| `molar_mass` | `float` | `g/mol` | الكتلة المولية المحسوبة من الأوزان الذرية القياسية (`IUPAC`). |
| `initial_mol`, `consumed_mol`, `remaining_mol` | `float` | `mol` | كميات المادة الابتدائية والمستهلكة والمتبقية بالمول. |
| `initial_mass_g`, `consumed_mass_g`, `remaining_mass_g` | `float` | `g` | الكتل الابتدائية والمستهلكة والمتبقية بالجرام. |
| `theoretical_yield_mol`, `theoretical_yield_g` | `float` | `mol`, `g` | المردود النظري للناتج بالمول والجرام. |
| `estimated_delta_h_kj` / `delta_h_kj` | `float \| null` | `kJ` | التغير في المحتوى الحراري الموثق أو المحسوب (`null` في الحالات التوضيحية). |
| `calculated_clogp`, `qsar_predicted_logp`, `experimental_logp` | `float \| null` | وحدة لوغاريتمية | قيم معامل التوزع الدهني-المائي ($\log P$). |

---

## 4. توثيق نقاط النهاية (Endpoints Specification)

### 4.1 استكشاف مادة علمية (`/api/v1/explore`)
- **الدالة:** `MolecuLensAPI.explore_entity(query, search_type="auto", allow_network=True, force_network=False)`
- **المدخلات:**
  - `query` (`str`): اسم المادة (عربي/إنجليزي)، الصيغة الكيميائية، أو `SMILES`.
  - `search_type` (`str`): `"auto" | "name" | "formula" | "smiles"`.
- **المخرجات (`data` عند `200 OK`):**
  - `entity`: السجل الموحد للمادة (`key`, `name_ar`, `name_en`, `symbol_or_formula`, `canonical_smiles`, `molar_mass`, `entity_type`).
  - `profile`: الملف العلمي الكامل شريطة التحقق البنيوي، خصائص ليبينسكي، وبطاقة سلامة ما قبل المختبر.

### 4.2 توقع التفاعل الكيميائي والمختبر البصري (`/api/v1/reactions/predict`)
- **الدالة:** `MolecuLensAPI.predict_reaction(reactants_input, domain=None, custom_conditions=None)`
- **المدخلات:**
  - `reactants_input`: نص حر مثل `"CH3COOH + K2CO3"` أو قائمة كائنات `[{"substance": {...}, "amount": 1.0, "unit": "mol"}]`.
- **المخرجات (`data` عند `200 OK`):**
  - `mode`: `"documented" | "rule_inferred"`.
  - `equation_unicode`, `equation_latex`, `net_ionic`: المعادلات الموزونة والصافية.
  - `atom_conservation`: تدقيق حفظ الذرات والكتلة.
  - `reactants_table`, `theoretical_yields`, `limiting_reagent`: الحسابات الستوكيومترية الكاملة.
  - `visual_experiment_stages`: المراحل البصرية الأربع (`stage1_preparation`, `stage2_active_vessel`, `stage3_mass_balance`, `stage4_final_inventory`).
- **حالة عدم كفاية الأدلة (`422 Unprocessable Entity`):**
  - `mode = "insufficient_evidence"`, `products = []` (لا يتم اختلاق أي نواتج غير موثقة).

### 4.3 تحليل التفاعلات الموثقة والستوكيومتري (`/api/v1/reactions/analyze`)
- **الدالة:** `MolecuLensAPI.analyze_reaction(reactant_inputs, domain="general_chemistry", custom_conditions=None)`
- **رموز الخطأ الخاصة:**
  - يرجع `422` مع `error.code = "INVALID_STP_PHASE"` عند محاولة استخدام وحدة `L_STP` لمادة صلبة أو سائلة أو محلول مائي.

### 4.4 تقييم معامل التوزع ونموذج QSAR (`/api/v1/molecular-intelligence/logp`)
- **الدالة:** `MolecuLensAPI.evaluate_logp(smiles)`
- **المخرجات (`data` عند `200 OK`):**
  - `deterministic_clogp`: حساب `Wildman-Crippen` بالتجزئة الذرية دون فواصل ثقة اصطناعية.
  - `qsar_model_evaluation`: تقييم نموذج `Ridge QSAR` على 120 مركباً (`90 train / 30 test`), والتحقق المتقاطع الخماسي (`5-Fold CV`)، والتحقق الخالي من التسرب بحجب الفئة الكيميائية (`Leak-Free LOCO`)، وجدول المقارنة المباشر `model_comparison_table` مع `Wildman-Crippen` (`R2`, `RMSE`, `MAE`).

### 4.5 التعديل البنيوي الافتراضي (`/api/v1/molecular-intelligence/modify`)
- **الدالة:** `MolecuLensAPI.modify_structure(smiles, modification_key)`
- **المخرجات (`data` عند `200 OK`):**
  - `is_virtual_modification = True`, `pharmacological_efficacy_claimed = False`.
  - `before`, `after`, `deltas`, `chemical_validity`, `site_info`, `scientific_caveats_ar`.

### 4.6 تقرير التدقيق العلمي المستقل (`/api/v1/reactions/audit`)
- **الدالة:** `MolecuLensAPI.scientific_audit()`
- **المخرجات (`data` عند `200 OK`):**
  - نتائج فحص الـ 87 تفاعلاً والـ 13 قاعدة استنتاج، مع التمييز الصريح بين التفاعلات المراجعة يدوياً من المصدر (`manually_verified_references = 34`) والتفاعلات المدققة برمجياً (`algorithmically_verified_only_references = 53`).
