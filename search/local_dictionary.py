import re

VALID_ELEMENTS = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
}

SUBSCRIPT_MAP = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")

SINGLE_LETTER_FORMULA_ELEMENTS = {"C", "H", "N", "O", "P", "S", "F", "K", "I", "B"}

KNOWN_NO_DIGIT_FORMULAS = {
    "O", "CO", "NO", "SO", "NaCl", "KCl", "HCl", "NaOH", "KOH", "HOOH", "LiCl", "MgO", "CaO", "FeO", "ZnO",
}


def replace_subscripts(text: str) -> str:
    if not text:
        return ""
    return text.translate(SUBSCRIPT_MAP)


def parse_chemical_formula(formula_str: str) -> dict[str, int] | None:
    if not formula_str or not isinstance(formula_str, str):
        return None
    cleaned = replace_subscripts(formula_str.strip()).replace(" ", "")
    if not cleaned or not re.fullmatch(r"[A-Za-z0-9]+", cleaned):
        return None

    candidate = cleaned
    if any(ch.isdigit() for ch in candidate) and candidate.islower():
        if all(ch.isdigit() or ch.upper() in SINGLE_LETTER_FORMULA_ELEMENTS for ch in candidate):
            candidate = candidate.upper()
    elif candidate.lower() in {"co2", "h2o", "o2", "h2o2", "c6h12o6", "nh3", "ch4", "so2", "no2", "h2s", "c2h6o", "c2h5oh", "ch3cooh"}:
        candidate = candidate.upper()

    tokens = re.findall(r"([A-Z][a-z]?)(\d*)", candidate)
    if not tokens:
        return None
    reconstructed = "".join(f"{el}{cnt}" for el, cnt in tokens)
    if reconstructed != candidate:
        return None

    counts: dict[str, int] = {}
    for elem, count_str in tokens:
        if elem not in VALID_ELEMENTS:
            return None
        count = int(count_str) if count_str else 1
        if count <= 0 or count > 500:
            return None
        counts[elem] = counts.get(elem, 0) + count

    return counts if counts else None


def format_hill_formula(counts: dict[str, int]) -> str:
    if not counts:
        return ""
    parts = []
    if "C" in counts:
        c_cnt = counts["C"]
        parts.append("C" if c_cnt == 1 else f"C{c_cnt}")
        if "H" in counts:
            h_cnt = counts["H"]
            parts.append("H" if h_cnt == 1 else f"H{h_cnt}")
        for elem in sorted(k for k in counts if k not in {"C", "H"}):
            cnt = counts[elem]
            parts.append(elem if cnt == 1 else f"{elem}{cnt}")
    else:
        for elem in sorted(counts.keys()):
            cnt = counts[elem]
            parts.append(elem if cnt == 1 else f"{elem}{cnt}")
    return "".join(parts)


def normalize_chemical_formula(formula_str: str) -> str | None:
    counts = parse_chemical_formula(formula_str)
    if not counts:
        return None
    return format_hill_formula(counts)


def formulas_match(formula_a: str, formula_b: str) -> bool:
    if not formula_a or not formula_b:
        return False
    clean_a = replace_subscripts(formula_a.strip()).replace(" ", "").upper()
    clean_b = replace_subscripts(formula_b.strip()).replace(" ", "").upper()
    if clean_a == clean_b:
        return True
    counts_a = parse_chemical_formula(formula_a)
    counts_b = parse_chemical_formula(formula_b)
    if counts_a is not None and counts_b is not None:
        return counts_a == counts_b
    return False


def is_chemical_formula_query(query: str) -> bool:
    if not query or not isinstance(query, str):
        return False
    cleaned = replace_subscripts(query.strip()).replace(" ", "")
    if not cleaned:
        return False
    if cleaned in KNOWN_NO_DIGIT_FORMULAS or cleaned.upper() in {"CO", "NO"}:
        return parse_chemical_formula(cleaned if cleaned in KNOWN_NO_DIGIT_FORMULAS else cleaned.upper()) is not None
    if not any(ch.isdigit() for ch in cleaned):
        return False
    return parse_chemical_formula(cleaned) is not None


LOCAL_COMPOUND_DICTIONARY = [
    {
        "key": "sitagliptin",
        "name_ar": "سيتاجليبتين",
        "name_en": "Sitagliptin",
        "iupac_name": "(3R)-3-amino-1-[3-(trifluoromethyl)-6,8-dihydro-5H-[1,2,4]triazolo[4,3-a]pyrazin-7-yl]-4-(2,4,5-trifluorophenyl)butan-1-one",
        "molecular_formula": "C16H15F6N5O",
        "canonical_smiles": "N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F",
        "pubchem_cid": 4369359,
        "chembl_id": "CHEMBL1422",
        "abbreviations": ["MK-0431"],
        "brand_names": ["Januvia", "جانوفيا"],
        "synonyms_ar": ["سيتاغليبتين", "سيتاجلبتين", "جانوفيا"],
        "synonyms_en": ["Sitagliptin", "MK-0431", "Januvia"],
        "info_ar": "مثبط انتقائي لإنزيم DPP4 يستخدم في ضبط سكر الدم لدى مرضى السكري من النوع الثاني.",
    },
    {
        "key": "vildagliptin",
        "name_ar": "فيلداجليبتين",
        "name_en": "Vildagliptin",
        "iupac_name": "(2S)-1-[2-[(3-hydroxy-1-adamantyl)amino]acetyl]pyrrolidine-2-carbonitrile",
        "molecular_formula": "C17H25N3O2",
        "canonical_smiles": "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2",
        "pubchem_cid": 6918537,
        "chembl_id": "CHEMBL142703",
        "abbreviations": ["LAF237"],
        "brand_names": ["Galvus", "جالفوس"],
        "synonyms_ar": ["فيلداغليبتين", "فيلداجلبتين", "جالفوس"],
        "synonyms_en": ["Vildagliptin", "LAF237", "Galvus"],
        "info_ar": "مثبط تساهمي عكسي لإنزيم DPP4 يحتوي على مجموعة سيانوبيروليدين وقاعدة أدامانتان.",
    },
    {
        "key": "saxagliptin",
        "name_ar": "ساكساجليبتين",
        "name_en": "Saxagliptin",
        "iupac_name": "(1S,3S,5S)-2-[(2S)-2-amino-2-(3-hydroxy-1-adamantyl)acetyl]-2-azabicyclo[3.1.0]hexane-3-carbonitrile",
        "molecular_formula": "C18H25N3O2",
        "canonical_smiles": "N#C[C@@H]1[C@@H]2C[C@H]2CN1C(=O)[C@H](N)C12CC3CC(CC(O)(C3)C1)C2",
        "pubchem_cid": 11243969,
        "chembl_id": "CHEMBL385517",
        "abbreviations": ["BMS-477118"],
        "brand_names": ["Onglyza", "أونجليزا"],
        "synonyms_ar": ["ساكساغليبتين", "ساكساجلبتين", "أونجليزا"],
        "synonyms_en": ["Saxagliptin", "BMS-477118", "Onglyza"],
        "info_ar": "مثبط قوي لإنزيم DPP4 يتميز بوجود حلقة ثلاثية مندمجة تزيد الثبات الأيضي.",
    },
    {
        "key": "linagliptin",
        "name_ar": "ليناجليبتين",
        "name_en": "Linagliptin",
        "iupac_name": "8-[(3R)-3-aminopiperidin-1-yl]-7-but-2-ynyl-3-methyl-1-[(4-methylquinazolin-2-yl)methyl]purine-2,6-dione",
        "molecular_formula": "C25H28N8O2",
        "canonical_smiles": "CC#CCn1c(N2CCC[C@@H](N)C2)nc2n(Cc3nc4ccccc4c(C)n3)c(=O)n(C)c(=O)c12",
        "pubchem_cid": 10096344,
        "chembl_id": "CHEMBL237500",
        "abbreviations": ["BI-1356"],
        "brand_names": ["Tradjenta", "Trajenta", "تراجينتا"],
        "synonyms_ar": ["ليناغليبتين", "ليناجلبتين", "تراجينتا"],
        "synonyms_en": ["Linagliptin", "BI-1356", "Tradjenta", "Trajenta"],
        "info_ar": "مثبط DPP4 مبني على نواة الزانثين، يتميز بفاعلية عالية جدا بتركيزات نانومولارية منخفضة.",
    },
    {
        "key": "alogliptin",
        "name_ar": "ألوجليبتين",
        "name_en": "Alogliptin",
        "iupac_name": "2-[[6-[(3R)-3-aminopiperidin-1-yl]-3-methyl-2,4-dioxopyrimidin-1-yl]methyl]benzonitrile",
        "molecular_formula": "C18H21N5O2",
        "canonical_smiles": "Cn1c(=O)cc(N2CCC[C@@H](N)C2)n(Cc2ccccc2C#N)c1=O",
        "pubchem_cid": 11450633,
        "chembl_id": "CHEMBL376359",
        "abbreviations": ["SYR-322"],
        "brand_names": ["Nesina", "نيسينا", "Vipidia"],
        "synonyms_ar": ["ألوغليبتين", "ألوجلبتين", "التراجليبتين", "الوجليبتين", "نيسينا"],
        "synonyms_en": ["Alogliptin", "SYR-322", "Nesina", "Vipidia"],
        "info_ar": "مثبط انتقائي لإنزيم DPP4 مشتق من نواة البيريميدين ديون.",
    },
    {
        "key": "carbon_dioxide",
        "name_ar": "ثاني أكسيد الكربون",
        "name_en": "Carbon Dioxide",
        "iupac_name": "carbon dioxide",
        "molecular_formula": "CO2",
        "canonical_smiles": "O=C=O",
        "pubchem_cid": 280,
        "chembl_id": "CHEMBL1231871",
        "abbreviations": ["CO2"],
        "brand_names": [],
        "synonyms_ar": ["غاز ثاني أكسيد الكربون", "ثنائي أكسيد الكربون", "ثاني اكسيد الكربون", "الجليد الجاف"],
        "synonyms_en": ["Carbon dioxide", "Carbonic anhydride", "Dry ice", "CO2"],
        "info_ar": "مركب كيميائي غازي يتكون من ذرة كربون مرتبطة برابطتين مزدوجتين مع ذرتي أكسجين (O=C=O)، وهو ناتج أساسي للتنفس الخلوي.",
    },
    {
        "key": "carbon_monoxide",
        "name_ar": "أول أكسيد الكربون",
        "name_en": "Carbon Monoxide",
        "iupac_name": "carbon monoxide",
        "molecular_formula": "CO",
        "canonical_smiles": "[C-]#[O+]",
        "pubchem_cid": 281,
        "chembl_id": "CHEMBL1231840",
        "abbreviations": ["CO"],
        "brand_names": [],
        "synonyms_ar": ["غاز أول أكسيد الكربون", "أحادي أكسيد الكربون", "اول اكسيد الكربون"],
        "synonyms_en": ["Carbon monoxide", "Carbonic oxide", "CO"],
        "info_ar": "غاز عديم اللون والرائحة يتكون من ذرة كربون وذرة أكسجين يرتبط بشدة بالهيموجلوبين.",
    },
    {
        "key": "water",
        "name_ar": "الماء",
        "name_en": "Water",
        "iupac_name": "oxidane",
        "molecular_formula": "H2O",
        "canonical_smiles": "O",
        "pubchem_cid": 962,
        "chembl_id": "CHEMBL1098659",
        "abbreviations": ["H2O"],
        "brand_names": [],
        "synonyms_ar": ["ماء", "الماء المقطر", "أكسيد الهيدروجين", "أوكسيدان"],
        "synonyms_en": ["Water", "Oxidane", "Dihydrogen monoxide", "Aqua", "H2O"],
        "info_ar": "المذيب الحيوي الأساسي في الأنظمة البيولوجية؛ يتكون من ذرتي هيدروجين وذرة أكسجين واحدة.",
    },
    {
        "key": "dioxygen",
        "name_ar": "الأكسجين الجزيئي",
        "name_en": "Molecular Oxygen",
        "iupac_name": "molecular oxygen",
        "molecular_formula": "O2",
        "canonical_smiles": "O=O",
        "pubchem_cid": 977,
        "chembl_id": "CHEMBL1234886",
        "abbreviations": ["O2"],
        "brand_names": [],
        "synonyms_ar": ["الأكسجين", "أكسجين", "ثنائي الأكسجين", "غاز الأكسجين"],
        "synonyms_en": ["Oxygen", "Dioxygen", "Molecular oxygen", "O2"],
        "info_ar": "جزيء الأكسجين الثنائي (O2) الضروري للتنفس الخلوي الهوائي وإنتاج الطاقة.",
    },
    {
        "key": "hydrogen_peroxide",
        "name_ar": "فوق أكسيد الهيدروجين",
        "name_en": "Hydrogen Peroxide",
        "iupac_name": "hydrogen peroxide",
        "molecular_formula": "H2O2",
        "canonical_smiles": "OO",
        "pubchem_cid": 784,
        "chembl_id": "CHEMBL718",
        "abbreviations": ["H2O2"],
        "brand_names": [],
        "synonyms_ar": ["ماء الأكسجين", "بيروكسيد الهيدروجين", "فوق اكسيد الهيدروجين"],
        "synonyms_en": ["Hydrogen peroxide", "Dioxidane", "Hydroperoxide", "H2O2"],
        "info_ar": "مركب مؤكسد شائع يعرف تجاريا باسم ماء الأكسجين (H2O2) ويحتوي على رابطة بيروكسيد أحادية O-O.",
    },
    {
        "key": "atp",
        "name_ar": "أدينوسين ثلاثي الفوسفات",
        "name_en": "Adenosine Triphosphate",
        "iupac_name": "[[(2R,3S,4R,5R)-5-(6-aminopurin-9-yl)-3,4-dihydroxyoxolan-2-yl]methoxy-hydroxyphosphoryl] phosphono hydrogen phosphate",
        "molecular_formula": "C10H16N5O13P3",
        "canonical_smiles": "Nc1ncnc2c1ncn2[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
        "pubchem_cid": 5957,
        "chembl_id": "CHEMBL14249",
        "abbreviations": ["ATP"],
        "brand_names": [],
        "synonyms_ar": ["ثلاثي فوسفات الأدينوسين", "عملة الطاقة الخلوية", "ادينوسين ثلاثي الفوسفات"],
        "synonyms_en": ["ATP", "Adenosine 5'-triphosphate", "Adenosine triphosphate"],
        "info_ar": "الناقل الكيميائي الرئيسي للطاقة الحيوية داخل الخلايا الحية.",
    },
    {
        "key": "nadh",
        "name_ar": "ثنائي نيوكليوتيد الأدينين والنيكوتيناميد المختزل",
        "name_en": "Nicotinamide Adenine Dinucleotide",
        "iupac_name": "[[(2R,3S,4R,5R)-5-(6-aminopurin-9-yl)-3,4-dihydroxyoxolan-2-yl]methoxy-hydroxyphosphoryl] [(2R,3S,4R,5R)-5-(3-carbamoyl-4H-pyridin-1-yl)-3,4-dihydroxyoxolan-2-yl]methyl phosphate",
        "molecular_formula": "C21H29N7O14P2",
        "canonical_smiles": "NC(=O)C1=CN([C@@H]2O[C@H](COP(=O)(O)OP(=O)(O)OC[C@H]3O[C@@H](n4cnc5c(N)ncnc54)[C@H](O)[C@@H]3O)[C@@H](O)[C@H]2O)C=CC1",
        "pubchem_cid": 439153,
        "chembl_id": "CHEMBL1234613",
        "abbreviations": ["NADH"],
        "brand_names": [],
        "synonyms_ar": ["مرافق الإنزيم NADH", "ناد إتش"],
        "synonyms_en": ["NADH", "Reduced nicotinamide adenine dinucleotide", "Beta-NADH"],
        "info_ar": "مرافق إنزيمي أساسي لتفاعلات الأكسدة والاختزال في سلسلة نقل الإلكترونات.",
    },
    {
        "key": "glucose",
        "name_ar": "الجلوكوز (سكر العنب)",
        "name_en": "D-Glucose",
        "iupac_name": "(3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol",
        "molecular_formula": "C6H12O6",
        "canonical_smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
        "pubchem_cid": 5793,
        "chembl_id": "CHEMBL1222250",
        "abbreviations": ["Glc"],
        "brand_names": ["Dextrose"],
        "synonyms_ar": ["جلوكوز", "دكستروز", "سكر الدم", "سكر العنب"],
        "synonyms_en": ["Glucose", "D-Glucose", "Dextrose", "Grape sugar"],
        "info_ar": "سكر أحادي سداسي الكربون (ألدوهكسوز) يشترك في الصيغة الجزيئية C6H12O6 مع الفركتوز والجلاكتوز.",
    },
    {
        "key": "fructose",
        "name_ar": "الفركتوز (سكر الفاكهة)",
        "name_en": "D-Fructose",
        "iupac_name": "(3S,4R,5R)-2-(hydroxymethyl)oxane-2,3,4,5-tetrol",
        "molecular_formula": "C6H12O6",
        "canonical_smiles": "OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@@H]1O",
        "pubchem_cid": 2723872,
        "chembl_id": "CHEMBL46740",
        "abbreviations": ["Fru"],
        "brand_names": ["Levulose"],
        "synonyms_ar": ["فركتوز", "سكر الفاكهة", "ليفولوز"],
        "synonyms_en": ["Fructose", "D-Fructose", "Fruit sugar", "Levulose"],
        "info_ar": "سكر أحادي كيتوني (كيتوهكسوز) يشترك في الصيغة الجزيئية C6H12O6 مع الجلوكوز والجلاكتوز.",
    },
    {
        "key": "galactose",
        "name_ar": "الجلاكتوز",
        "name_en": "D-Galactose",
        "iupac_name": "(3R,4S,5R,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol",
        "molecular_formula": "C6H12O6",
        "canonical_smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O",
        "pubchem_cid": 6036,
        "chembl_id": "CHEMBL15932",
        "abbreviations": ["Gal"],
        "brand_names": [],
        "synonyms_ar": ["جلاكتوز", "غالاكتوز"],
        "synonyms_en": ["Galactose", "D-Galactose"],
        "info_ar": "سكر أحادي سداسي الكربون يدخل في تركيب اللاكتوز ويشترك في الصيغة الجزيئية C6H12O6.",
    },
    {
        "key": "ethanol",
        "name_ar": "الإيثانول (الكحول الإيثيلي)",
        "name_en": "Ethanol",
        "iupac_name": "ethanol",
        "molecular_formula": "C2H6O",
        "canonical_smiles": "CCO",
        "pubchem_cid": 702,
        "chembl_id": "CHEMBL545",
        "abbreviations": ["EtOH", "C2H5OH"],
        "brand_names": [],
        "synonyms_ar": ["إيثانول", "الكحول الإيثيلي", "كحول طبي"],
        "synonyms_en": ["Ethanol", "Ethyl alcohol", "Grain alcohol", "C2H5OH"],
        "info_ar": "كحول أولي مذيب شائع يشترك في الصيغة الجزيئية C2H6O مع ثنائي ميثيل الإيثر.",
    },
    {
        "key": "dimethyl_ether",
        "name_ar": "ثنائي ميثيل الإيثر",
        "name_en": "Dimethyl Ether",
        "iupac_name": "methoxymethane",
        "molecular_formula": "C2H6O",
        "canonical_smiles": "COC",
        "pubchem_cid": 8254,
        "chembl_id": "CHEMBL1232461",
        "abbreviations": ["DME"],
        "brand_names": [],
        "synonyms_ar": ["إيثر ثنائي الميثيل", "ميثوكسي ميثان"],
        "synonyms_en": ["Dimethyl ether", "Methoxymethane", "Methyl ether"],
        "info_ar": "أبسط إيثر أليفاتي وهو مصاوغ بنيوي (Isomer) للإيثانول يشترك معه في الصيغة الجزيئية C2H6O.",
    },
    {
        "key": "ammonia",
        "name_ar": "الأمونيا (النشادر)",
        "name_en": "Ammonia",
        "iupac_name": "azane",
        "molecular_formula": "H3N",
        "canonical_smiles": "N",
        "pubchem_cid": 222,
        "chembl_id": "CHEMBL1160819",
        "abbreviations": ["NH3"],
        "brand_names": [],
        "synonyms_ar": ["الأمونيا", "النشادر", "غاز النشادر", "أزان"],
        "synonyms_en": ["Ammonia", "Azane", "NH3", "H3N"],
        "info_ar": "مركب لاعضوي أساسي يتكون من ذرة نيتروجين وثلاث ذرات هيدروجين (NH3).",
    },
    {
        "key": "methane",
        "name_ar": "الميثان",
        "name_en": "Methane",
        "iupac_name": "methane",
        "molecular_formula": "CH4",
        "canonical_smiles": "C",
        "pubchem_cid": 297,
        "chembl_id": "CHEMBL17564",
        "abbreviations": ["CH4"],
        "brand_names": [],
        "synonyms_ar": ["ميثان", "غاز الميثان", "غاز المستنقعات"],
        "synonyms_en": ["Methane", "Marsh gas", "CH4"],
        "info_ar": "أبسط هيدروكربون ألكاني يتكون من ذرة كربون واحدة وأربع ذرات هيدروجين (CH4).",
    },
    {
        "key": "metformin",
        "name_ar": "ميتفورمين",
        "name_en": "Metformin",
        "iupac_name": "3-(diaminomethylidene)-1,1-dimethylguanidine",
        "molecular_formula": "C4H11N5",
        "canonical_smiles": "CN(C)C(=N)NC(=N)N",
        "pubchem_cid": 4091,
        "chembl_id": "CHEMBL1431",
        "abbreviations": [],
        "brand_names": ["Glucophage", "جلوكوفاج"],
        "synonyms_ar": ["ميتفورمن", "جلوكوفاج", "منظم السكر"],
        "synonyms_en": ["Metformin", "Glucophage", "Dimethylbiguanide"],
        "info_ar": "الدواء الأول لعلاج السكري من النوع الثاني، يعمل عبر تقليل إنتاج الجلوكوز الكبدي وليس عبر تثبيط DPP4 مباشرة.",
    },
    {
        "key": "aspirin",
        "name_ar": "الأسبرين (حمض أسيتيل الساليسيليك)",
        "name_en": "Aspirin",
        "iupac_name": "2-acetyloxybenzoic acid",
        "molecular_formula": "C9H8O4",
        "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "pubchem_cid": 2244,
        "chembl_id": "CHEMBL25",
        "abbreviations": ["ASA"],
        "brand_names": ["Aspirin", "Jusprin", "جوسبرين"],
        "synonyms_ar": ["أسبرين", "اسبرين", "حمض أسيتيل الساليسيليك", "جوسبرين"],
        "synonyms_en": ["Aspirin", "Acetylsalicylic acid", "2-Acetoxybenzoic acid"],
        "info_ar": "مسكن ومضاد التهاب غير ستيرويدي ومضاد لتكدس الصفائح الدموية (مركب مرجعي للمقارنة خارج نطاق DPP4).",
    },
    {
        "key": "paracetamol",
        "name_ar": "باراسيتامول (أسيتامينوفين)",
        "name_en": "Paracetamol",
        "iupac_name": "N-(4-hydroxyphenyl)acetamide",
        "molecular_formula": "C8H9NO2",
        "canonical_smiles": "CC(=O)Nc1ccc(O)cc1",
        "pubchem_cid": 1983,
        "chembl_id": "CHEMBL112",
        "abbreviations": ["APAP"],
        "brand_names": ["Panadol", "Tylenol", "بنادول", "تايلينول", "فيفادول"],
        "synonyms_ar": ["باراسيتامول", "أسيتامينوفين", "بنادول", "فيفادول"],
        "synonyms_en": ["Paracetamol", "Acetaminophen", "Panadol", "Tylenol", "APAP"],
        "info_ar": "مسكن للألم وخافض للحرارة واسع الاستخدام في الممارسة الصيدلانية.",
    },
    {
        "key": "caffeine",
        "name_ar": "الكافيين",
        "name_en": "Caffeine",
        "iupac_name": "1,3,7-trimethylpurine-2,6-dione",
        "molecular_formula": "C8H10N4O2",
        "canonical_smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
        "pubchem_cid": 2519,
        "chembl_id": "CHEMBL113",
        "abbreviations": [],
        "brand_names": [],
        "synonyms_ar": ["كافيين", "قهوين"],
        "synonyms_en": ["Caffeine", "1,3,7-Trimethylxanthine", "Theine"],
        "info_ar": "شبه قلوي من عائلة الزانثين يعمل منبها للجهاز العصبي المركزي، يشترك جزئيا في نواة الزانثين مع دواء ليناجليبتين.",
    },
    {
        "key": "ibuprofen",
        "name_ar": "إيبوبروفين",
        "name_en": "Ibuprofen",
        "iupac_name": "2-[4-(2-methylpropyl)phenyl]propanoic acid",
        "molecular_formula": "C13H18O2",
        "canonical_smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1",
        "pubchem_cid": 3672,
        "chembl_id": "CHEMBL521",
        "abbreviations": [],
        "brand_names": ["Brufen", "Advil", "بروفين", "أدفيل"],
        "synonyms_ar": ["ايبوبروفين", "بروفين", "أدفيل"],
        "synonyms_en": ["Ibuprofen", "Brufen", "Advil", "Motrin"],
        "info_ar": "مضاد التهاب غير ستيرويدي شائع لتسكين الألم وخفض الحرارة.",
    },
]

AMBIGUOUS_SHORT_QUERIES = {
    "o": {
        "reason_ar": "الرمز الكيميائي 'O' قد يقصد به ذرة الأكسجين المنفردة أو جزيء الأكسجين الثنائي (O2)، كما يمثل في صيغة SMILES جزيء الماء (H2O). يرجى تحديد المركب المقصود بدقة.",
        "candidate_keys": ["dioxygen", "water"],
    },
    "c6h12o6": {
        "reason_ar": "الصيغة الكيميائية C6H12O6 تشترك فيها عدة مصاوغات سكرية سداسية (Isomers) مثل الجلوكوز والفركتوز والجلاكتوز. يرجى اختيار المركب المقصود من القائمة.",
        "candidate_keys": ["glucose", "fructose", "galactose"],
    },
    "c2h6o": {
        "reason_ar": "الصيغة الكيميائية C2H6O تشترك فيها مصاوغات بنيوية مختلفة مثل الإيثانول (CCO) وثنائي ميثيل الإيثر (COC). يرجى اختيار المركب المقصود.",
        "candidate_keys": ["ethanol", "dimethyl_ether"],
    },
}


def get_all_local_compounds() -> list[dict]:
    return [dict(c) for c in LOCAL_COMPOUND_DICTIONARY]

