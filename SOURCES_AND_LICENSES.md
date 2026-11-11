# المصادر العلمية والتراخيص (Sources & Licenses) — MolecuLens

## 1. ترخيص الكود المصدري (Source Code License)
يخضع الكود المصدري لمنصة **MolecuLens** لترخيص **MIT License** المرفق في ملف [LICENSE](./LICENSE).

---

## 2. المصادر العلمية المعتمدة في المنصة (Verified Scientific Data Sources)
تعتمد طبقة `data_integrations/data_manager.py` والفهرس المحلي الموثق (`SmartLocalDatabase`) على المصادر العلمية المرجعية التالية، مع الفصل الصريح بين السجلات المحلية الموثقة بمراجع علمية (`local_curated_registry`) والسجلات المسترجعة مباشرة من واجهات البرمجة الخارجية (`External APIs`):

| المصدر العلمي (Source) | الدور في المنصة (Role in Platform) | الرابط المرجعي (Official URL) | الترخيص وشروط الإسناد (License) |
|---|---|---|---|
| **PubChem (NCBI)** | الأسماء، الصيغ الجزيئية، البنية الجزيئية (SMILES / InChIKey)، والخصائص الفيزيائية | [https://pubchem.ncbi.nlm.nih.gov](https://pubchem.ncbi.nlm.nih.gov) | Public Domain (US Government Work — NCBI/NLM) |
| **ChEMBL (EMBL-EBI)** | البيانات الدوائية، الهدف الجزيئي، والنشاط الحيوي التجريبي الموثق | [https://www.ebi.ac.uk/chembl](https://www.ebi.ac.uk/chembl) | Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0) |
| **ChEBI (EMBL-EBI)** | التصنيف الكيميائي الحيوي للجزيئات الصغيرة والمرادفات القياسية | [https://www.ebi.ac.uk/chebi](https://www.ebi.ac.uk/chebi) | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| **UniProt Knowledgebase** | البروتينات، التسلسلات الأمينية، والوظائف الحيوية | [https://www.uniprot.org](https://www.uniprot.org) | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| **RCSB Protein Data Bank (PDB)** | التراكيب الفراغية البلورية ثلاثية الأبعاد للبروتينات والمعقدات الدوائية | [https://www.rcsb.org](https://www.rcsb.org) | CC0 1.0 Universal (Public Domain Dedication) |
| **Rhea Database** | التفاعلات الكيميائية الحيوية والإنزيمية الموثقة (مع التحقق من حفظ الذرات) | [https://www.rhea-db.org](https://www.rhea-db.org) | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| **Open Reaction Database (ORD)** | استعلام التفاعلات العضوية الخارجية بمعرفات `ord-` الموثقة بعد التحقق من موازنتها | [https://open-reaction-database.org](https://open-reaction-database.org) | Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0) |
| **IUPAC & NIST-JANAF** | الكتل الذرية القياسية للعناصر الكيميائية وجداول الثوابت الديناميكية الحرارية | [https://iupac.org](https://iupac.org) / [https://webbook.nist.gov](https://webbook.nist.gov) | IUPAC Standard Atomic Weights & NIST Standard Reference Data |

---

## 3. المراجع المطبوعة والدستورية لسجل التفاعلات المحلي (Literature Citations for the 87 Reactions)
يُوثق كل تفاعل من التفاعلات الـ87 في `reaction_engine/reaction_registry.py` بمعرف مرجعي فريد (`source_id`) وتصنيف توثيقي صريح (`experimentally_documented` | `literature_referenced` | `theoretical_model`) استناداً إلى المراجع القياسية التالية:
1. **United States Pharmacopeia and National Formulary (USP 43–NF 38)** & **European Pharmacopoeia (Ph. Eur. 10th Ed.)**: معايرة الأحماض والقواعد القياسية واختبارات الهوية الدوائية.
2. **Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis (5th Ed., 1979)** & **Vogel's Textbook of Quantitative Chemical Analysis (5th Ed., 1989)**: تفاعلات الترسيب النوعي والوزني، وتفاعلات الكربونات والبيكربونات، وإحلال الفلزات.
3. **Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements (2nd Ed., Butterworth-Heinemann, 1997)**: تفاعلات الأكسدة والاختزال غير العضوية وتحضير الأكاسيد والهاليدات.
4. **Chase, M.W., NIST-JANAF Thermochemical Tables (4th Ed., J. Phys. Chem. Ref. Data, Monograph 9, 1998)**: التفاعلات الحرارية القياسية وتكوين الأكاسيد والأملاح.
5. **Ullmann's Encyclopedia of Industrial Chemistry (Wiley-VCH, 2012)**: العمليات الصناعية القياسية (عملية هابر-بوش، طريقة التلامس، كلسنة الحجر الجيري، التحليل الكهربائي للماء).
6. **Smith, M.B., March's Advanced Organic Chemistry (8th Ed., Wiley, 2020)** & **Vogel's Textbook of Practical Organic Chemistry (5th Ed., 1989)**: تفاعلات أسترة فيشر، الهدرجة الحفزية، الهلجنة، والنترتة العطرية.
7. **Nelson, D.L. & Cox, M.M., Lehninger Principles of Biochemistry (8th Ed., W.H. Freeman, 2021)**: المعادلات الستوكيومترية الكلية للتنفس الخلوي الهوائي، البناء الضوئي، والتخمر الكحولي.

---

## 4. المكتبات البرمجية مفتوحة المصدر (Open-Source Software Libraries)
- **RDKit Cheminformatics:** BSD 3-Clause License ([https://www.rdkit.org](https://www.rdkit.org))
- **Streamlit Framework:** Apache License 2.0 ([https://streamlit.io](https://streamlit.io))
- **NumPy / Pandas / Plotly / Requests / Pillow / Pytest:** BSD / MIT / Apache 2.0 Licenses
