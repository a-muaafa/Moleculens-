"""
Drug-Drug Interaction (DDI) engine for MolecuLens.
Provides pairwise interaction checks, severity grading, pharmacological mechanisms,
clinical management recommendations, evidence citations, and mandatory safety disclaimers.
"""

from typing import Dict, Any, List, Optional, Tuple
import html


class DDIEngine:
    """Evaluates drug-drug interactions with scientific rigor and clinical safety disclaimers."""

    DISCLAIMER_AR = (
        "تنبيه سريري وأكاديمي هام: هذه الأداة مخصصة للأغراض التعليمية والاسترشادية فقط. "
        "إن عدم وجود تفاعل مسجل في قاعدة البيانات الحالية لا يعني بالضرورة الأمان التام "
        "أو غياب التداخل الدوائي السريري. لا تُستخدم هذه الأداة كبديل عن الاستشارة الصيدلانية "
        "أو الطبية المباشرة، ويجب دائماً مراجعة المراجع الدوائية المعتمدة لكل حالة سريرية."
    )

    DISCLAIMER_EN = (
        "Important Clinical & Academic Notice: This tool is designed strictly for educational and guidance purposes. "
        "The absence of a documented interaction in this database does NOT prove safety or the complete absence "
        "of a clinical interaction. This tool does not replace professional medical or clinical pharmacy consultation; "
        "always cross-reference official prescribing monographs for individualized patient care."
    )

    # Keyed by frozenset of canonical names
    _INTERACTIONS: Dict[frozenset, Dict[str, Any]] = {
        frozenset(["sitagliptin", "metformin"]): {
            "severity": "Minor",
            "severity_ar": "بسيط / مراقبة روتينية",
            "mechanism_en": "Synergistic glycemic control without direct pharmacokinetic competition; both drugs are cleared primarily via renal tubular excretion.",
            "mechanism_ar": "تأثير تآزري مفيد في خفض سكر الدم دون تداخل حركي مباشر، كلاهما يُطرح عبر الكلى، مما يتطلب مراقبة وظائف الكلى الدورية.",
            "clinical_effect_en": "Enhanced reduction in HbA1c; mild increase in gastrointestinal side effects common to metformin; low risk of hypoglycemia unless combined with insulin secretagogues.",
            "clinical_effect_ar": "تحسن ملحوظ في ضبط السكر التراكمي (HbA1c) مع احتمالية زيادة خفيفة في الأعراض الهضمية الشائعة للميتفورمين، وخطر منخفض لهبوط السكر.",
            "management_en": "Commonly prescribed as fixed-dose combination. Monitor renal function (eGFR) periodically. Dosage adjustment required if eGFR falls below 45 mL/min/1.73m2.",
            "management_ar": "يُصرفان معاً بشكل شائع كمركب دوائي مشترك. تجب مراقبة وظائف الكلى دورياً وتعديل الجرعة إذا انخفض معدل الترشيح الكبيبي عن 45 مل/دقيقة.",
            "evidence_level": "Established",
            "citations": [
                "FDA Prescribing Information: Janumet (Sitagliptin/Metformin HCl)",
                "Diabetes Care 2023; Standards of Medical Care in Diabetes",
            ],
        },
        frozenset(["aspirin", "ibuprofen"]): {
            "severity": "Major",
            "severity_ar": "خطير / رئيسي",
            "mechanism_en": "Ibuprofen competitively inhibits platelet cyclooxygenase-1 (COX-1), blocking the irreversible acetylation of Ser529 by low-dose aspirin and diminishing cardioprotection; additive inhibition of gastric mucosal COX-1 increases bleeding ulcer risk.",
            "mechanism_ar": "يثبط الإيبوبروفين إنزيم الأكسدة الحلقية-1 (COX-1) في الصفائح الدموية تنافسياً، مما يمنع الأسبرين من تثبيطه نهائياً ويقلل الحماية القلبية؛ كما يضاعف الخطر على بطانة المعدة.",
            "clinical_effect_en": "Attenuation of aspirin's antiplatelet cardioprotective effect; significantly heightened risk of gastrointestinal mucosal erosion, ulceration, and hemorrhage.",
            "clinical_effect_ar": "إضعاف التأثير المضاد للتخثر الواقي للقلب للأسبرين، وارتفاع كبير في مخاطر قرحة ونزيف الجهاز الهضمي.",
            "management_en": "Avoid concurrent use if possible. If ibuprofen is required, take immediate-release ibuprofen at least 30 minutes after or at least 8 hours before low-dose aspirin.",
            "management_ar": "يُفضل تجنب التناول المتزامن. إذا لزم الإيبوبروفين، يجب تناوله بعد 30 دقيقة على الأقل من الأسبرين أو قبل الأسبرين بـ 8 ساعات على الأقل.",
            "evidence_level": "Established",
            "citations": [
                "FDA Drug Safety Communication: Concomitant use of Ibuprofen and Aspirin (2006)",
                "Catella-Lawson F et al. N Engl J Med 2001; 345:1809-1817",
            ],
        },
        frozenset(["aspirin", "warfarin"]): {
            "severity": "Major",
            "severity_ar": "خطير / رئيسي",
            "mechanism_en": "Combined anticoagulant (vitamin K antagonism) and antiplatelet (irreversible COX-1 inhibition) actions produce profound synergistic impairment of hemostasis.",
            "mechanism_ar": "تأثير مضاعف يجمع بين تثبيط عوامل التخثر المعتمدة على فيتامين K وتثبيط الصفائح الدموية، مما يؤدي إلى تثبيط شديد لمنظومة التجلط.",
            "clinical_effect_en": "Substantially increased risk of severe, potentially life-threatening gastrointestinal or intracranial bleeding.",
            "clinical_effect_ar": "ارتفاع حاد في احتمالية حدوث نزيف هضمي أو دماغي شديد ومهدد للحياة.",
            "management_en": "Only combine when strictly indicated (e.g., mechanical heart valve plus acute coronary syndrome). Close INR monitoring and co-prescription of gastroprotective PPI recommended.",
            "management_ar": "لا يُجمع بينهما إلا بدواعٍ سريرية محددة وبإشراف استشاري دقيق مع مراقبة لصيقة لمؤشر INR وإضافة واقٍ للمعدة مثل مثبطات مضخة البروتون.",
            "evidence_level": "Established",
            "citations": [
                "Chest 2012; 141(2 Suppl):e152S-e184S (Antithrombotic Therapy Guidelines)",
                "Circulation 2018; Dual Antithrombotic Therapy Consensus",
            ],
        },
        frozenset(["aspirin", "paracetamol"]): {
            "severity": "Minor",
            "severity_ar": "طفيف / مراقبة",
            "mechanism_en": "Both inhibit prostaglandin synthesis but through different central and peripheral mechanisms. Concomitant therapeutic doses generally do not cause major pharmacokinetic interference.",
            "mechanism_ar": "كلاهما يؤثر في تثبيط البروستاغلاندينات عبر آليات مركزية ومحيطية مختلفة، ولا يسببان تداخلاً حركياً كبيراً عند الجرعات الاعتيادية.",
            "clinical_effect_en": "Additive analgesic effect. High chronic doses may increase cumulative renal toxicity risk.",
            "clinical_effect_ar": "تأثير مسكن مضاف، ولكن الاستخدام المزمن بجرعات عالية قد يرفع العبء التراكمي على الكلى.",
            "management_en": "Safe for short-term combination. Do not exceed maximum daily limits (4000 mg paracetamol / 4000 mg aspirin).",
            "management_ar": "آمن للاستخدام المؤقت قصير المدى. يجب عدم تجاوز الحد الأقصى اليومي المسموح به لكل منهما.",
            "evidence_level": "Established",
            "citations": [
                "British National Formulary (BNF 85): Analgesic Combinations",
            ],
        },
        frozenset(["paracetamol", "warfarin"]): {
            "severity": "Moderate",
            "severity_ar": "متوسط / يستلزم المتابعة",
            "mechanism_en": "Paracetamol metabolite (NAPQI) may inhibit enzymes in the vitamin K cycle (VKORC1), causing mild inhibition of clotting factor synthesis when taken regularly.",
            "mechanism_ar": "قد يثبط ناتج أيض الباراسيتامول الإنزيمات المسؤولة عن تدوير فيتامين K، مما يؤدي إلى تثبيط إضافي لإنتاج عوامل التجلط عند تناوله بانتظام.",
            "clinical_effect_en": "Regular high-dose paracetamol (>2 g/day for >3 days) may elevate INR and increase bleeding tendency.",
            "clinical_effect_ar": "تناول الباراسيتامول بجرعات عالية (>2 جم/يوم لأكثر من 3 أيام) قد يؤدي لارتفاع ملحوظ في قراءة INR وزيادة قابلية النزيف.",
            "management_en": "Occasional single doses have minimal impact. Monitor INR if paracetamol is used regularly for more than 3-4 days.",
            "management_ar": "الجرعات المفردة العارضة آمنة ولا تؤثر، ولكن عند الاستمرار لأكثر من 3 أيام يجب فحص INR وتعديل جرعة الوارفارين إن لزم.",
            "evidence_level": "Probable",
            "citations": [
                "Mahé I et al. Haematologica 2006; 91(12):1621-1627",
                "Ann Pharmacother 2007; 41(7):1220-1224",
            ],
        },
        frozenset(["caffeine", "aspirin"]): {
            "severity": "Minor",
            "severity_ar": "طفيف / تآزر إيجابي",
            "mechanism_en": "Caffeine accelerates gastric emptying and improves gastrointestinal absorption of aspirin; also exerts central adenosine receptor antagonism.",
            "mechanism_ar": "يعجل الكافيين إفراغ المعدة ويزيد معدل امتصاص الأسبرين، كما يحصر مستقبلات الأدينوزين المركزية.",
            "clinical_effect_en": "Enhanced analgesic onset and potency (often formulated together in OTC pain relievers). Higher doses may cause tachycardia or gastric distress.",
            "clinical_effect_ar": "تسريع وبدء مفعول تسكين الألم وزيادة فاعليته (يُدمجان عمداً في أدوية الصداع)، وقد يسبب زيادة خفيفة في نبضات القلب أو تهيج المعدة.",
            "management_en": "Frequently utilized clinically. Avoid in patients with acute peptic ulcers or uncontrolled tachyarrhythmias.",
            "management_ar": "مزيج مستخدم سريرياً بشكل واسع. يُنصح بالحذر لدى مرضى القرحة النشطة أو اضطرابات نبضات القلب.",
            "evidence_level": "Established",
            "citations": [
                "Cochrane Database Syst Rev 2014; (12):CD009281",
            ],
        },
        frozenset(["omeprazole", "amoxicillin"]): {
            "severity": "Minor",
            "severity_ar": "طفيف / تآزر علاجي مستهدف",
            "mechanism_en": "Omeprazole elevates gastric pH, preventing acid degradation of amoxicillin and improving its local antibacterial stability against Helicobacter pylori.",
            "mechanism_ar": "يرفع الأوميبرازول درجة حموضة المعدة، مما يحمي الأموكسيسيلين من التحلل الحمضي ويزيد من استقراره وفاعليته الموضعية ضد جرثومة المعدة.",
            "clinical_effect_en": "Synergistic eradication of Helicobacter pylori in standard triple-therapy regimens.",
            "clinical_effect_ar": "تعزيز القضاء على بكتيريا الملوية البوابية (جرثومة المعدة) ضمن بروتوكولات العلاج الثلاثي القياسية.",
            "management_en": "Standard approved combination for H. pylori eradication therapy. Follow recommended treatment course (typically 10-14 days).",
            "management_ar": "مزيج قياسي معتمد عالمياً لعلاج جرثومة المعدة. يتم الالتزام بالمدة المقررة للعلاج (غالباً 10 إلى 14 يوماً).",
            "evidence_level": "Established",
            "citations": [
                "ACG Clinical Guideline: Treatment of Helicobacter pylori Infection (Am J Gastroenterol 2017)",
            ],
        },
        frozenset(["digoxin", "amiodarone"]): {
            "severity": "Major",
            "severity_ar": "خطير / رئيسي",
            "mechanism_en": "Amiodarone inhibits P-glycoprotein (P-gp / ABCB1) efflux transporter and decreases renal and non-renal clearance of digoxin.",
            "mechanism_ar": "يثبط الأميودارون مضخة P-glycoprotein المسؤولة عن إخراج الديجوكسين، مما يقلل طرحه الكلوي والكبدي ويرفع تركيزه بالدم.",
            "clinical_effect_en": "Serum digoxin levels approximately double within 3-7 days, creating high risk of fatal digitalis toxicity (arrhythmias, heart block, visual disturbances).",
            "clinical_effect_ar": "مضاعفة تركيز الديجوكسين في مصل الدم خلال أيام، مما يسبب تسمماً خطيراً بالديجوكسين واضطرابات قلبية قاتلة.",
            "management_en": "Reduce digoxin maintenance dose by 50% upon initiating amiodarone. Measure serum digoxin concentrations after 1-2 weeks and monitor ECG.",
            "management_ar": "تخفيض جرعة الديجوكسين بنسبة 50% فور بدء الأميودارون، مع قياس مستواه في الدم دورياً وتخطيط القلب الكهربائي.",
            "evidence_level": "Established",
            "citations": [
                "FDA Prescribing Information: Lanoxin (Digoxin)",
                "Circulation 2010; Digoxin-Amiodarone Pharmacokinetic Interactions",
            ],
        },
    }

    _SYNONYMS = {
        "paracetamol": "paracetamol",
        "acetaminophen": "paracetamol",
        "panadol": "paracetamol",
        "tylenol": "paracetamol",
        "aspirin": "aspirin",
        "acetylsalicylic acid": "aspirin",
        "asa": "aspirin",
        "ibuprofen": "ibuprofen",
        "advil": "ibuprofen",
        "motrin": "ibuprofen",
        "sitagliptin": "sitagliptin",
        "januvia": "sitagliptin",
        "metformin": "metformin",
        "glucophage": "metformin",
        "caffeine": "caffeine",
        "amoxicillin": "amoxicillin",
        "amoxil": "amoxicillin",
        "omeprazole": "omeprazole",
        "prilosec": "omeprazole",
        "losec": "omeprazole",
        "warfarin": "warfarin",
        "coumadin": "warfarin",
        "digoxin": "digoxin",
        "lanoxin": "digoxin",
        "amiodarone": "amiodarone",
        "cordarone": "amiodarone",
    }

    @classmethod
    def normalize_drug_name(cls, name: str) -> str:
        """Normalizes drug name or brand name to standard canonical key."""
        clean = (name or "").strip().lower()
        return cls._SYNONYMS.get(clean, clean)

    @classmethod
    def get_supported_drugs(cls) -> List[str]:
        """Returns sorted list of distinct canonical drugs available in interaction database."""
        drugs = set()
        for pair in cls._INTERACTIONS.keys():
            drugs.update(pair)
        return sorted(list(drugs))

    @classmethod
    def check_interaction(cls, drug_a: str, drug_b: str) -> Dict[str, Any]:
        """
        Evaluates potential interaction between two drugs.
        Returns clinical interaction details or safe negative report with mandatory disclaimer.
        """
        norm_a = cls.normalize_drug_name(drug_a)
        norm_b = cls.normalize_drug_name(drug_b)

        if not norm_a or not norm_b:
            return {
                "success": False,
                "error": "Two valid drug names are required.",
                "disclaimer_ar": cls.DISCLAIMER_AR,
                "disclaimer_en": cls.DISCLAIMER_EN,
            }

        if norm_a == norm_b:
            return {
                "success": True,
                "is_same_drug": True,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "severity": "Duplicate Therapy",
                "severity_ar": "تكرار دوائي لنفس المادة",
                "mechanism_en": "Duplicate active ingredient therapy.",
                "mechanism_ar": "استخدام نفس المادة الفعالة أو مشتق مطابق مما يؤدي لمضاعفة الجرعة واحتمال السمية.",
                "clinical_effect_en": "Risk of unintentional overdose and enhanced adverse events.",
                "clinical_effect_ar": "خطر الجرعة الزائدة غير المقصودة ومضاعفة الآثار الجانبية السمية.",
                "management_en": "Verify that patient is not unknowingly taking multiple brands containing identical active substance.",
                "management_ar": "التأكد من عدم تناول المريض لعلامات تجارية مختلفة تحتوي على نفس المادة الفعالة.",
                "evidence_level": "Established",
                "citations": ["Clinical Pharmacology Guidelines: Duplicate Therapy Prevention"],
                "disclaimer_ar": cls.DISCLAIMER_AR,
                "disclaimer_en": cls.DISCLAIMER_EN,
            }

        key = frozenset([norm_a, norm_b])
        interaction = cls._INTERACTIONS.get(key)

        if interaction:
            return {
                "success": True,
                "has_interaction": True,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "canonical_a": norm_a,
                "canonical_b": norm_b,
                "severity": interaction["severity"],
                "severity_ar": interaction["severity_ar"],
                "mechanism_en": interaction["mechanism_en"],
                "mechanism_ar": interaction["mechanism_ar"],
                "clinical_effect_en": interaction["clinical_effect_en"],
                "clinical_effect_ar": interaction["clinical_effect_ar"],
                "management_en": interaction["management_en"],
                "management_ar": interaction["management_ar"],
                "evidence_level": interaction["evidence_level"],
                "citations": interaction["citations"],
                "disclaimer_ar": cls.DISCLAIMER_AR,
                "disclaimer_en": cls.DISCLAIMER_EN,
            }

        # Safe negative report with explicit disclaimer
        return {
            "success": True,
            "has_interaction": False,
            "drug_a": drug_a,
            "drug_b": drug_b,
            "canonical_a": norm_a,
            "canonical_b": norm_b,
            "severity": "No Interaction Documented",
            "severity_ar": "لا يوجد تداخل مسجل بقاعدة البيانات",
            "message_en": f"No documented clinical interaction between '{drug_a}' and '{drug_b}' in the local database.",
            "message_ar": f"لا يوجد تداخل دوائي مسجل بين '{drug_a}' و '{drug_b}' في قاعدة البيانات المحلية.",
            "evidence_level": "Not Found",
            "citations": [],
            "disclaimer_ar": cls.DISCLAIMER_AR,
            "disclaimer_en": cls.DISCLAIMER_EN,
        }

    @classmethod
    def screen_polypharmacy(cls, drug_list: List[str]) -> Dict[str, Any]:
        """Screens all pairwise combinations in a multi-drug regimen."""
        clean_list = [d.strip() for d in drug_list if d and d.strip()]
        unique_drugs = list(dict.fromkeys(clean_list))
        pairs_checked = []
        interactions_found = []

        for i in range(len(unique_drugs)):
            for j in range(i + 1, len(unique_drugs)):
                d1 = unique_drugs[i]
                d2 = unique_drugs[j]
                res = cls.check_interaction(d1, d2)
                pairs_checked.append((d1, d2))
                if res.get("has_interaction"):
                    interactions_found.append(res)

        severity_counts = {"Major": 0, "Moderate": 0, "Minor": 0}
        for inter in interactions_found:
            sev = inter.get("severity")
            if sev in severity_counts:
                severity_counts[sev] += 1

        return {
            "success": True,
            "total_drugs": len(unique_drugs),
            "total_pairs_evaluated": len(pairs_checked),
            "total_interactions_found": len(interactions_found),
            "severity_summary": severity_counts,
            "interactions": interactions_found,
            "disclaimer_ar": cls.DISCLAIMER_AR,
            "disclaimer_en": cls.DISCLAIMER_EN,
        }

