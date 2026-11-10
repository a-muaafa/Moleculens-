# المصادر العلمية وتراخيص المكتبات والبيانات (Sources and Licenses)

يوثق هذا الملف جميع المكتبات البرمجية ومصادر البيانات العلمية المعتمدة في مشروع **MolecuLens** (فريق **BioCode**)، مع روابطها الرسمية، وتراخيصها، وشروط استخدامها، ومتطلبات الإسناد العلمي (Attribution Requirements).

> **تنبيه قانوني ومنهجي:** ترخيص الكود المصدري للمشروع (MIT License في ملف `LICENSE`) لا يلغي ولا يحل محل تراخيص وشروط استخدام قواعد البيانات العلمية الخارجية أو المكتبات البرمجية المدرجة أدناه.

---

## أولا: المكتبات البرمجية المفتوحة المصدر (Open-Source Software Libraries)

| اسم المكتبة | الرابط الرسمي | نوع الترخيص | شروط الاستخدام والإسناد |
|---|---|---|---|
| **RDKit** | https://www.rdkit.org | BSD 3-Clause License | يسمح بالاستخدام والتعديل والتوزيع للأغراض البحثية والتعليمية والتجارية مع الإبقاء على إشعار حقوق النشر الأصلي. |
| **Streamlit** | https://streamlit.io | Apache License 2.0 | يسمح بالاستخدام والتوزيع مع الاحتفاظ بإشعارات الترخيص وحقوق الملكية الفكرية لشركة Snowflake / Streamlit. |
| **Scikit-Learn** | https://scikit-learn.org | BSD 3-Clause License | يسمح بالاستخدام البحثي والتعليمي والتطويري مع الإبقاء على إشعار حقوق النشر والإسناد الأكاديمي عند النشر العلمي. |
| **Pandas** | https://pandas.pydata.org | BSD 3-Clause License | يسمح بالاستخدام والتوزيع مع الاحتفاظ بإشعار حقوق النشر. |
| **NumPy** | https://numpy.org | BSD 3-Clause License | يسمح بالاستخدام والتوزيع مع الاحتفاظ بإشعار حقوق النشر. |
| **Plotly** | https://plotly.com/python | MIT License | يسمح بالاستخدام والنسخ والتعديل مع تضمين إشعار الترخيص. |
| **Pytest** | https://docs.pytest.org | MIT License | يستخدم لاختبارات الوحدة الآلية والتكامل المستمر. |

---

## ثانيا: مصادر البيانات الكيميائية والحيوية (Scientific Data Sources & Attribution)

### 1. قاعدة بيانات ChEMBL (EMBL-EBI)
- **اسم المصدر:** ChEMBL Database (European Molecular Biology Laboratory - European Bioinformatics Institute)
- **الرابط الرسمي:** https://www.ebi.ac.uk/chembl/
- **الهدف البروتيني المستخدم:** `CHEMBL284` (Human Dipeptidyl peptidase 4)
- **الترخيص:** Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0)
- **شروط الاستخدام ومتطلبات الإسناد العلمي:**
  - يجب نسبة البيانات صراحة إلى قاعدة بيانات ChEMBL (EMBL-EBI) والاحتفاظ بمعرفات المركبات (`molecule_chembl_id`) ومعرفات التجارب (`assay_chembl_id`) ومعرفات الوثائق (`document_chembl_id`).
  - أي قاعدة بيانات مشتقة أو معدلة من بيانات ChEMBL يجب أن تلتزم بشرط الترخيص بالمثل (ShareAlike) وإسناد المصدر الأصلي:
    *Zdrazil, B. et al. The ChEMBL Database in 2023: a drug discovery platform spanning multiple bioactivity data types and time periods. Nucleic Acids Research (2024).*

### 2. قاعدة بيانات UniProt (Universal Protein Resource)
- **اسم المصدر:** UniProt Consortium (EMBL-EBI, SIB, PIR)
- **الرابط الرسمي:** https://www.uniprot.org/uniprotkb/P27487/entry
- **المعرف المرجعي المستخدم:** `P27487` (`DPP4_HUMAN` - Homo sapiens)
- **الترخيص:** Creative Commons Attribution 4.0 International (CC BY 4.0)
- **شروط الاستخدام ومتطلبات الإسناد العلمي:**
  - يسمح بمشاركة واستخدام البيانات الوظيفية والبنيوية للبروتين مع ذكر المصدر:
    *The UniProt Consortium, UniProt: the Universal Protein Knowledgebase, Nucleic Acids Research.*

### 3. قاعدة بيانات PubChem (NCBI / NLM / NIH)
- **اسم المصدر:** PubChem - National Center for Biotechnology Information (U.S. National Library of Medicine)
- **الرابط الرسمي:** https://pubchem.ncbi.nlm.nih.gov/
- **واجهة البرمجة المستخدمة:** PubChem PUG REST API
- **الترخيص وشروط الاستخدام:** Public Domain / NCBI Website and Data Usage Policies
- **متطلبات الإسناد والاستخدام التقني:**
  - الالتزام بحدود معدل الطلبات (Rate Limits) الخاصة بـ PUG REST API وعدم إرسال طلبات مكثفة تعيق الخدمة.
  - إسناد معرف المركب `PubChem CID` وذكر المصدر العلمي:
    *Kim, S. et al. PubChem 2023 update. Nucleic Acids Research (2023).*

