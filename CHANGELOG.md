# سجل التغييرات (Changelog) — MolecuLens

## [2.0.0] - 2026-11-11 — إعادة بناء شاملة للمنصة (Platform Architecture Rebuild)
- إعادة بناء مشروع **MolecuLens** من الصفر كمنصة تعليمية واستكشافية علمية (`Scientific Learning & Exploration Platform`) قائمة على قسمين رئيسيين فقط: `Explore & Understand` و `Reaction Lab`.
- إزالة الواجهات والملفات القديمة المتمحورة حول مراحل DPP4 الخمس.
- إنشاء الهيكل المعماري الجديد المكون من 11 طبقة (`config`, `core`, `search`, `data_integrations`, `entity_profiles`, `chemistry_engine`, `reaction_engine`, `pharmacy`, `biology`, `education`, `frontend`).
- بناء محرك تصنيف الكيانات العلمية (`عنصر كيميائي`، `مركب كيميائي`، `مادة دوائية`، `بروتين أو جزيء حيوي`) وبناء ملفات علمية مخصصة لكل فئة.
- بناء `Data Manager` موحد للاتصال بمصادر `PubChem`, `ChEMBL`, `UniProt`, `RCSB PDB`, و `Rhea` مع نظام تخزين مؤقت (`Cache Store`).
- بناء محرك التفاعلات (`Reaction Engine`) لموازنة المعادلات جبريا، والتحقق من حفظ الذرات، وحساب المتفاعل المحدد والفائض والناتج النظري، ورفض اختلاق التفاعلات غير المدعومة.
- بناء حزمة اختبارات آلية شاملة (`24/24` اختبار ناجح في `tests/`).
