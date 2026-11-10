import difflib
import logging
import urllib.parse
import requests
from chemistry.molecular_structure import get_molecular_formula, validate_smiles
from search.local_dictionary import LOCAL_COMPOUND_DICTIONARY, normalize_arabic_text

logger = logging.getLogger(__name__)


def _entry_matches_field(entry: dict, norm_q: str, raw_q: str, search_type: str) -> tuple[bool, str]:
    q_lower = raw_q.strip().lower()

    if search_type in ("auto", "name_ar"):
        if normalize_arabic_text(entry["name_ar"]) == norm_q:
            return True, "name_ar"
        for syn in entry.get("synonyms_ar", []):
            if normalize_arabic_text(syn) == norm_q:
                return True, "synonym_ar"

    if search_type in ("auto", "name_en"):
        if entry["name_en"].lower() == q_lower:
            return True, "name_en"

    if search_type in ("auto", "iupac"):
        if entry.get("iupac_name", "").lower() == q_lower:
            return True, "iupac"

    if search_type in ("auto", "formula"):
        if entry.get("molecular_formula", "").lower() == q_lower:
            return True, "formula"

    if search_type in ("auto", "abbreviation"):
        for abbr in entry.get("abbreviations", []):
            if abbr.lower() == q_lower:
                return True, "abbreviation"

    if search_type in ("auto", "synonym"):
        for syn_en in entry.get("synonyms_en", []):
            if syn_en.lower() == q_lower:
                return True, "synonym_en"
        for syn_ar in entry.get("synonyms_ar", []):
            if normalize_arabic_text(syn_ar) == norm_q:
                return True, "synonym_ar"

    if search_type in ("auto", "brand"):
        for br in entry.get("brand_names", []):
            if br.lower() == q_lower or normalize_arabic_text(br) == norm_q:
                return True, "brand"

    if search_type in ("auto", "smiles"):
        mol, canon_smi, _ = validate_smiles(raw_q)
        if mol is not None and canon_smi == entry.get("canonical_smiles"):
            return True, "smiles"

    return False, ""


def _compute_fuzzy_score(entry: dict, raw_q: str, norm_q: str) -> float:
    q_lower = raw_q.strip().lower()
    terms_en = (
        [entry["name_en"], entry.get("iupac_name", ""), entry.get("molecular_formula", "")]
        + entry.get("synonyms_en", [])
        + entry.get("abbreviations", [])
        + entry.get("brand_names", [])
    )
    terms_ar = [entry["name_ar"]] + entry.get("synonyms_ar", []) + entry.get("brand_names", [])

    best = 0.0
    for t in terms_en:
        if not t:
            continue
        t_low = t.lower()
        if q_lower in t_low:
            best = max(best, 0.86)
        ratio = difflib.SequenceMatcher(None, q_lower, t_low).ratio()
        best = max(best, ratio)

    for t in terms_ar:
        if not t:
            continue
        t_norm = normalize_arabic_text(t)
        if norm_q and norm_q in t_norm:
            best = max(best, 0.86)
        ratio = difflib.SequenceMatcher(None, norm_q, t_norm).ratio()
        best = max(best, ratio)

    return best


def query_pubchem_online(query_str: str, timeout: float = 4.0) -> list[dict]:
    encoded = urllib.parse.quote(query_str.strip())
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded}/"
        "property/Title,MolecularFormula,IUPACName,CanonicalSMILES,IsomericSMILES/JSON"
    )
    resp = requests.get(url, timeout=timeout)
    if resp.status_code != 200:
        return []
    data = resp.json()
    props = data.get("PropertyTable", {}).get("Properties", [])
    results = []
    for p in props[:5]:
        raw_smi = p.get("IsomericSMILES") or p.get("CanonicalSMILES") or ""
        mol, canon_smi, _ = validate_smiles(raw_smi)
        if mol is None:
            continue
        cid = p.get("CID")
        title = p.get("Title") or query_str.strip()
        formula = p.get("MolecularFormula") or get_molecular_formula(mol)
        iupac = p.get("IUPACName") or title
        results.append(
            {
                "id": f"pubchem_{cid}",
                "name_ar": title,
                "name_en": title,
                "synonyms_ar": [],
                "synonyms_en": [title],
                "abbreviations": [],
                "brand_names": [],
                "iupac_name": iupac,
                "molecular_formula": formula,
                "pubchem_cid": cid,
                "chembl_id": None,
                "canonical_smiles": canon_smi,
                "info_ar": "تم جلب بيانات البنية الجزيئية والمعرفات من قاعدة بيانات PubChem والتحقق من صحتها عبر RDKit.",
                "source": f"قاعدة بيانات PubChem عبر الإنترنت (CID: {cid})",
            }
        )
    return results


def verify_compound_identity(entry: dict) -> dict:
    smi = entry.get("canonical_smiles", "")
    mol, canon_smi, err = validate_smiles(smi)
    if mol is None:
        return {
            "is_verified": False,
            "error": err or "تعذر التحقق من البنية الجزيئية للمركب.",
            "compound": entry,
        }
    calc_formula = get_molecular_formula(mol)
    verified_entry = dict(entry)
    verified_entry["canonical_smiles"] = canon_smi
    if not verified_entry.get("molecular_formula"):
        verified_entry["molecular_formula"] = calc_formula
    return {
        "is_verified": True,
        "error": None,
        "compound": verified_entry,
        "verification_steps": {
            "compound_identified": True,
            "name_and_synonyms_checked": bool(verified_entry.get("name_en") and verified_entry.get("name_ar")),
            "formula_verified": calc_formula == verified_entry.get("molecular_formula", calc_formula),
            "structure_2d_ready": True,
            "smiles_validated": True,
            "identifier_present": bool(verified_entry.get("pubchem_cid") or verified_entry.get("chembl_id")),
        },
    }


def search_compound(
    query: str,
    search_type: str = "auto",
    allow_network: bool = True,
) -> dict:
    if not query or not isinstance(query, str) or not query.strip():
        return {
            "found": False,
            "match_type": "empty",
            "is_ambiguous": False,
            "ambiguity_message": None,
            "did_you_mean": [],
            "candidates": [],
            "selected": None,
            "data_source": "القاموس المحلي",
            "connection_status": "جاهز (القاموس المحلي الموثق)",
            "message": "يرجى إدخال اسم المركب، صيغته الكيميائية، اختصاره، أو صيغة SMILES للبحث.",
        }

    raw_q = query.strip()
    norm_q = normalize_arabic_text(raw_q)

    exact_matches = []
    for item in LOCAL_COMPOUND_DICTIONARY:
        matched, field_used = _entry_matches_field(item, norm_q, raw_q, search_type)
        if matched:
            copy_item = dict(item)
            copy_item["matched_by"] = field_used
            exact_matches.append(copy_item)

    if raw_q.upper() == "O" and search_type == "auto":
        o2_entry = next((dict(x) for x in LOCAL_COMPOUND_DICTIONARY if x["id"] == "molecular_oxygen"), None)
        if o2_entry and all(m["id"] != "molecular_oxygen" for m in exact_matches):
            o2_entry["matched_by"] = "related_oxygen"
            exact_matches.append(o2_entry)
        return {
            "found": True,
            "match_type": "exact_ambiguous",
            "is_ambiguous": True,
            "ambiguity_message": (
                "المدخل (O) يحتمل أكثر من معنى علمي: الرمز O يمثل ذرة أكسجين منفردة في الصيغ الكيميائية، "
                "ويمثل جزيء الماء (H2O) في ترميز SMILES، بينما جزيء الأكسجين الثنائي هو O2. يرجى اختيار المركب المقصود."
            ),
            "did_you_mean": [m["name_ar"] for m in exact_matches],
            "candidates": exact_matches,
            "selected": None,
            "data_source": "القاموس المحلي الموثق",
            "connection_status": "تم الاسترجاع من القاموس المحلي دون الحاجة للاتصال الخارجي",
            "message": "تم العثور على عدة مركبات محتملة للمدخل؛ يرجى تأكيد المركب المقصود.",
        }

    if len(exact_matches) > 1:
        return {
            "found": True,
            "match_type": "exact_ambiguous",
            "is_ambiguous": True,
            "ambiguity_message": (
                f"المدخل ({raw_q}) يطابق {len(exact_matches)} مركبات مختلفة (مثل المصاوغات المشتركة في الصيغة الجزيئية). "
                "لا يختار النظام مركبا عشوائيا؛ يرجى تحديد المركب المطلوب لتأكيد هويته."
            ),
            "did_you_mean": [f"{m['name_ar']} ({m['name_en']})" for m in exact_matches],
            "candidates": exact_matches,
            "selected": None,
            "data_source": "القاموس المحلي الموثق",
            "connection_status": "تم الاسترجاع من القاموس المحلي دون الحاجة للاتصال الخارجي",
            "message": "يرجى اختيار المركب الصحيح من بين النتائج المطابقة.",
        }

    if len(exact_matches) == 1:
        primary = exact_matches[0]
        candidates = [primary]
        if primary["id"] == "molecular_oxygen":
            atom_o = next((dict(x) for x in LOCAL_COMPOUND_DICTIONARY if x["id"] == "atomic_oxygen"), None)
            if atom_o:
                candidates.append(atom_o)
        return {
            "found": True,
            "match_type": "exact",
            "is_ambiguous": False,
            "ambiguity_message": None,
            "did_you_mean": [],
            "candidates": candidates,
            "selected": primary,
            "data_source": primary.get("source", "القاموس المحلي الموثق"),
            "connection_status": "تم الاسترجاع من القاموس المحلي الموثق",
            "message": f"تم تحديد المركب بنجاح: {primary['name_ar']} ({primary['name_en']}).",
        }

    if search_type in ("auto", "smiles"):
        mol, canon_smi, _ = validate_smiles(raw_q)
        if mol is not None and any(c in raw_q for c in ("=", "#", "(", ")", "[", "]", "1", "2", "c", "n")):
            formula = get_molecular_formula(mol)
            custom_entry = {
                "id": "custom_smiles",
                "name_ar": f"مركب مدخل بصيغة SMILES ({formula})",
                "name_en": f"Custom SMILES Compound ({formula})",
                "synonyms_ar": [],
                "synonyms_en": [],
                "abbreviations": [],
                "brand_names": [],
                "iupac_name": "غير محدد محليا (مدخل مباشر بصيغة SMILES)",
                "molecular_formula": formula,
                "pubchem_cid": None,
                "chembl_id": None,
                "canonical_smiles": canon_smi,
                "info_ar": "تم التحقق من صحة بنية SMILES مباشرة باستخدام مكتبة RDKit.",
                "source": "إدخال SMILES مباشر تم التحقق منه عبر RDKit",
            }
            return {
                "found": True,
                "match_type": "direct_smiles",
                "is_ambiguous": False,
                "ambiguity_message": None,
                "did_you_mean": [],
                "candidates": [custom_entry],
                "selected": custom_entry,
                "data_source": custom_entry["source"],
                "connection_status": "تحقق محلي عبر RDKit",
                "message": "تم التحقق من صحة صيغة SMILES المدخلة.",
            }

    scored = []
    for item in LOCAL_COMPOUND_DICTIONARY:
        score = _compute_fuzzy_score(item, raw_q, norm_q)
        if score >= 0.56:
            c_item = dict(item)
            c_item["similarity_score"] = round(score, 3)
            scored.append((score, c_item))
    scored.sort(key=lambda x: x[0], reverse=True)
    fuzzy_candidates = [x[1] for x in scored[:5]]

    connection_status = "وضع عدم الاتصال (القاموس المحلي فقط)"
    if allow_network:
        try:
            online_hits = query_pubchem_online(raw_q)
            if online_hits:
                is_ambig = len(online_hits) > 1
                return {
                    "found": True,
                    "match_type": "pubchem_online",
                    "is_ambiguous": is_ambig,
                    "ambiguity_message": (
                        "تم العثور على أكثر من سجل مطابق في PubChem؛ يرجى اختيار المركب المقصود."
                        if is_ambig
                        else None
                    ),
                    "did_you_mean": [h["name_en"] for h in online_hits],
                    "candidates": online_hits,
                    "selected": None if is_ambig else online_hits[0],
                    "data_source": "قاعدة بيانات PubChem (متصل)",
                    "connection_status": "متصل بقاعدة بيانات PubChem",
                    "message": "تم استرجاع بيانات المركب من PubChem والتحقق من البنية عبر RDKit.",
                }
            connection_status = "متصل بـ PubChem (لم يعثر على تطابق مباشر، تم الرجوع للقاموس المحلي)"
        except Exception as exc:
            logger.warning("PubChem lookup failed for '%s': %s", raw_q, exc)
            connection_status = "تعذر الاتصال بـ PubChem - يعمل النظام تلقائيا عبر المسار البديل (القاموس المحلي)"

    if fuzzy_candidates:
        did_you_mean_list = [f"{c['name_ar']} ({c['name_en']} - {c['molecular_formula']})" for c in fuzzy_candidates]
        return {
            "found": True,
            "match_type": "approximate",
            "is_ambiguous": True,
            "ambiguity_message": (
                f"لم يتم العثور على تطابق حرفي تام للمدخل ({raw_q}). "
                "هل تقصد أحد المركبات القريبة التالية؟ يرجى اختيار المركب الصحيح لتأكيده."
            ),
            "did_you_mean": did_you_mean_list,
            "candidates": fuzzy_candidates,
            "selected": None,
            "data_source": "القاموس المحلي الموثق (بحث تقريبي وتصحيح أخطاء)",
            "connection_status": connection_status,
            "message": "تم العثور على اقتراحات قريبة عبر البحث التقريبي (Fuzzy Search).",
        }

    return {
        "found": False,
        "match_type": "none",
        "is_ambiguous": False,
        "ambiguity_message": None,
        "did_you_mean": ["الماء (Water - H2O)", "سيتاجليبتين (Sitagliptin)", "فيلداجليبتين (Vildagliptin)", "الأسبرين (Aspirin)"],
        "candidates": [],
        "selected": None,
        "data_source": "القاموس المحلي",
        "connection_status": connection_status,
        "message": f"لم يتم العثور على مركب مطابق للمدخل ({raw_q}) في القاموس المحلي أو المصادر المتاحة.",
    }

