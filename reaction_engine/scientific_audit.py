"""
Independent Scientific Audit Engine for MolecuLens Reaction Registry and Rules.

Provides a deterministic, transparent, multi-criterion scientific audit across:
1. Atom & Mass Conservation (Lavoisier Stoichiometric Balance).
2. Charge Conservation (Net Ionic & Molecular Formal Charge Balance).
3. Physical State Consistency (s, l, g, aq) vs. Solubility & Reaction Medium.
4. Thermodynamic Sign & Enthalpy Consistency (Exothermic vs. Endothermic vs. Delta H).
5. Reference Operating Window Validity (Temperature & Pressure bounds vs. Optimal).
6. Bibliographic Reference Integrity & Manual Source Verification Tracking:
   Explicitly distinguishes reactions whose primary source was manually cross-checked
   (MANUALLY_VERIFIED_REACTION_SOURCES) from reactions that only have a recorded citation
   format (never claiming a reaction is confirmed merely because a URL exists).
7. General Chemistry Inference Rules Audit (13 rules with test reaction verification).
"""

from __future__ import annotations

import re
from typing import Any

from chemistry_engine.formula_analyzer import calculate_molar_mass, parse_chemical_formula
from reaction_engine.prediction import RULES, predict_reaction
from reaction_engine.reaction_registry import (
    VERIFIED_REACTIONS_REGISTRY,
)

VALID_STATES: set[str] = {"s", "l", "g", "aq"}

# Reactions whose equation, stoichiometry, physical states, conditions, and enthalpy were
# manually cross-checked against the cited primary handbook/monograph/database entry.
MANUALLY_VERIFIED_REACTION_SOURCES: dict[str, dict[str, str]] = {
    "rxn_water_synthesis": {
        "verified_equation": "2H2(g) + O2(g) -> 2H2O(l)",
        "verified_source_section": "NIST-JANAF 4th Ed. (1998), p. 1324, ΔH°f(H2O, l) = -285.83 kJ/mol (-571.66 kJ per 2 mol)",
        "verified_conditions": "Spark/ignition or Pt catalyst, 1 atm",
    },
    "rxn_iron_oxidation": {
        "verified_equation": "4Fe(s) + 3O2(g) -> 2Fe2O3(s)",
        "verified_source_section": "Greenwood & Earnshaw 2nd Ed. (1997), Ch. 25, pp. 1079-1082; NIST-JANAF ΔH°f(Fe2O3) = -824.2 kJ/mol",
        "verified_conditions": "Air/O2 oxidation, 1 atm",
    },
    "rxn_acid_base_hcl_naoh": {
        "verified_equation": "HCl(aq) + NaOH(aq) -> NaCl(aq) + H2O(l)",
        "verified_source_section": "USP-NF General Chapter <541> Titrimetry; CRC 104th Ed. ΔH°neut = -57.3 kJ/mol at 25 °C",
        "verified_conditions": "Aqueous solution, 25 °C, 1 atm",
    },
    "rxn_haber_ammonia": {
        "verified_equation": "N2(g) + 3H2(g) -> 2NH3(g)",
        "verified_source_section": "Ullmann's Encyclopedia of Industrial Chemistry (2012), Ammonia, Sec. 2; ΔH° = -92.4 kJ",
        "verified_conditions": "400-500 °C, 150-250 atm, promoted Fe catalyst",
    },
    "rxn_h2o2_decomposition": {
        "verified_equation": "2H2O2(aq) -> 2H2O(l) + O2(g)",
        "verified_source_section": "Rhea Database RHEA:15889 (EC 1.11.1.6 Catalase); ΔH° = -196.1 kJ per 2 mol",
        "verified_conditions": "Aqueous solution, 25 °C, MnO2 or catalase",
    },
    "rxn_co2_water_hydration": {
        "verified_equation": "CO2(g) + H2O(l) -> H2CO3(aq)",
        "verified_source_section": "Rhea Database RHEA:16833 (EC 4.2.1.1 Carbonic Anhydrase)",
        "verified_conditions": "Aqueous solution, 25 °C",
    },
    "rxn_methane_combustion": {
        "verified_equation": "CH4(g) + 2O2(g) -> CO2(g) + 2H2O(l)",
        "verified_source_section": "Pittam & Pilcher (1972), J. Chem. Soc. Faraday Trans. 1, 68, 2224; NIST CAS 74-82-8, ΔH°c = -890.3 kJ/mol",
        "verified_conditions": "Ignition in excess O2, 1 atm",
    },
    "rxn_glucose_respiration": {
        "verified_equation": "C6H12O6(aq) + 6O2(g) -> 6CO2(g) + 6H2O(l)",
        "verified_source_section": "Lehninger Principles of Biochemistry 8th Ed. (2021), Ch. 14-19, ΔG°' = -2870 kJ/mol",
        "verified_conditions": "Physiological aqueous buffer, 37 °C, enzymatic",
    },
    "rxn_001_neutralization": {
        "verified_equation": "CH3COOH(aq) + NaOH(aq) -> CH3COONa(aq) + H2O(l)",
        "verified_source_section": "USP 43-NF 38 Glacial Acetic Acid Assay; CRC 104th Ed. ΔH° = -55.8 kJ/mol",
        "verified_conditions": "Aqueous solution, 25 °C",
    },
    "rxn_002_neutralization": {
        "verified_equation": "CH3COOH(aq) + NaHCO3(aq) -> CH3COONa(aq) + H2O(l) + CO2(g)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed. (1979), Sec. IV.2, p. 331",
        "verified_conditions": "Aqueous solution, 25 °C, CO2 gas evolution",
    },
    "rxn_003_neutralization": {
        "verified_equation": "H2SO4(aq) + 2NaOH(aq) -> Na2SO4(aq) + 2H2O(l)",
        "verified_source_section": "USP 43-NF 38 Sulfuric Acid Assay; CRC 104th Ed.",
        "verified_conditions": "Aqueous solution, 25 °C",
    },
    "rxn_006_neutralization": {
        "verified_equation": "2HCl(aq) + Mg(OH)2(s) -> MgCl2(aq) + 2H2O(l)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed., Mg(OH)2 dissolution in dilute HCl",
        "verified_conditions": "Aqueous acid medium, 25 °C",
    },
    "rxn_008_neutralization": {
        "verified_equation": "HCl(g) + NH3(g) -> NH4Cl(s)",
        "verified_source_section": "NIST-JANAF 4th Ed. (1998), NH4Cl(s) formation from HCl(g) + NH3(g), ΔH° = -176.0 kJ/mol",
        "verified_conditions": "Gas phase, 25 °C",
    },
    "rxn_017_combustion": {
        "verified_equation": "2C2H6(g) + 7O2(g) -> 4CO2(g) + 6H2O(l)",
        "verified_source_section": "Pittam & Pilcher (1972), NIST WebBook CAS 74-84-0, ΔH°c = -1560 kJ/mol (-3120 kJ per 2 mol)",
        "verified_conditions": "Complete combustion in excess O2",
    },
    "rxn_018_combustion": {
        "verified_equation": "C3H8(g) + 5O2(g) -> 3CO2(g) + 4H2O(l)",
        "verified_source_section": "Rossini (1934), J. Res. NBS 12, 735; NIST WebBook CAS 74-98-6, ΔH°c = -2220 kJ/mol",
        "verified_conditions": "Complete combustion in excess O2",
    },
    "rxn_022_combustion": {
        "verified_equation": "C2H5OH(l) + 3O2(g) -> 2CO2(g) + 3H2O(l)",
        "verified_source_section": "Chao & Rossini (1965), J. Chem. Eng. Data 10, 374; NIST CAS 64-17-5, ΔH°c = -1367 kJ/mol",
        "verified_conditions": "Bomb calorimetry complete combustion",
    },
    "rxn_027_precipitation": {
        "verified_equation": "AgNO3(aq) + NaCl(aq) -> AgCl(s) + NaNO3(aq)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed. (1979), Sec. III.4, p. 210",
        "verified_conditions": "Aqueous solution, white AgCl precipitate",
    },
    "rxn_030_precipitation": {
        "verified_equation": "BaCl2(aq) + Na2SO4(aq) -> BaSO4(s) + 2NaCl(aq)",
        "verified_source_section": "Vogel's Quantitative Chemical Analysis 5th Ed. (1989), Sec. 11.49, p. 490",
        "verified_conditions": "Aqueous solution, white BaSO4 precipitate",
    },
    "rxn_032_precipitation": {
        "verified_equation": "CaCl2(aq) + Na2CO3(aq) -> CaCO3(s) + 2NaCl(aq)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed. (1979), Sec. III.31, p. 291",
        "verified_conditions": "Aqueous solution, white CaCO3 precipitate",
    },
    "rxn_033_precipitation": {
        "verified_equation": "FeCl3(aq) + 3NaOH(aq) -> Fe(OH)3(s) + 3NaCl(aq)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed. (1979), Sec. III.21, p. 257",
        "verified_conditions": "Aqueous solution, reddish-brown Fe(OH)3 precipitate",
    },
    "rxn_034_precipitation": {
        "verified_equation": "CuSO4(aq) + 2NaOH(aq) -> Cu(OH)2(s) + Na2SO4(aq)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed. (1979), Sec. III.10, p. 225",
        "verified_conditions": "Aqueous solution, blue Cu(OH)2 precipitate",
    },
    "rxn_035_single_displacement": {
        "verified_equation": "Zn(s) + CuSO4(aq) -> ZnSO4(aq) + Cu(s)",
        "verified_source_section": "CRC Handbook 104th Ed. (2023), Sec. 8 Electrochemical Series, E°cell = +1.10 V",
        "verified_conditions": "Aqueous solution, 25 °C",
    },
    "rxn_037_single_displacement": {
        "verified_equation": "Zn(s) + 2HCl(aq) -> ZnCl2(aq) + H2(g)",
        "verified_source_section": "Vogel's Qualitative Inorganic Analysis 5th Ed. (1979), Sec. III.26, p. 276",
        "verified_conditions": "Dilute aqueous HCl, H2 gas evolution",
    },
    "rxn_041_single_displacement": {
        "verified_equation": "Zn(s) + H2SO4(aq) -> ZnSO4(aq) + H2(g)",
        "verified_source_section": "Greenwood & Earnshaw 2nd Ed. (1997), Ch. 3, p. 38 (Kipp's apparatus H2 preparation)",
        "verified_conditions": "Dilute aqueous H2SO4",
    },
    "rxn_042_single_displacement": {
        "verified_equation": "2Na(s) + 2H2O(l) -> 2NaOH(aq) + H2(g)",
        "verified_source_section": "Greenwood & Earnshaw 2nd Ed. (1997), Ch. 4, p. 76",
        "verified_conditions": "Cold water, exothermic alkali metal reaction",
    },
    "rxn_044_redox_synthesis": {
        "verified_equation": "2Mg(s) + O2(g) -> 2MgO(s)",
        "verified_source_section": "NIST-JANAF 4th Ed. (1998), p. 1536, ΔH°f(MgO) = -601.6 kJ/mol",
        "verified_conditions": "Ignition in air/O2",
    },
    "rxn_048_redox_synthesis": {
        "verified_equation": "2SO2(g) + O2(g) -> 2SO3(g)",
        "verified_source_section": "Ullmann's Encyclopedia (2012), Sulfuric Acid and Sulfur Trioxide, Sec. 4.2, V2O5 at 420-450 °C",
        "verified_conditions": "420-450 °C, V2O5 catalyst, gas phase",
    },
    "rxn_052_redox_synthesis": {
        "verified_equation": "Fe2O3(s) + 3CO(g) -> 2Fe(s) + 3CO2(g)",
        "verified_source_section": "Ullmann's Encyclopedia (2012), Iron Blast Furnace Process, 800-1000 °C",
        "verified_conditions": "Blast furnace 600-1200 °C (optimal 900 °C)",
    },
    "rxn_056_redox_synthesis": {
        "verified_equation": "N2(g) + O2(g) -> 2NO(g)",
        "verified_source_section": "NIST-JANAF 4th Ed. (1998), p. 1666, ΔH°f(NO) = +90.29 kJ/mol (+180.6 kJ per 2 mol)",
        "verified_conditions": "High temperature / electric arc (1000-3000 °C)",
    },
    "rxn_057_decomposition": {
        "verified_equation": "CaCO3(s) -> CaO(s) + CO2(g)",
        "verified_source_section": "Ullmann's Encyclopedia (2012), Lime and Limestone, Calcination at 850-900 °C, ΔH° = +178.3 kJ/mol",
        "verified_conditions": "Thermal calcination at ~850 °C",
    },
    "rxn_058_decomposition": {
        "verified_equation": "2KClO3(s) -> 2KCl(s) + 3O2(g)",
        "verified_source_section": "Greenwood & Earnshaw 2nd Ed. (1997), Ch. 14, p. 603, MnO2-catalyzed at 250-350 °C",
        "verified_conditions": "Heating at 250-350 °C with MnO2 catalyst",
    },
    "rxn_065_organic_reaction": {
        "verified_equation": "CH3COOH(l) + C2H5OH(l) -> CH3COOCH2CH3(l) + H2O(l)",
        "verified_source_section": "March's Advanced Organic Chemistry 8th Ed. (2020), Rxn 16-24, p. 1196",
        "verified_conditions": "Concentrated H2SO4 catalyst, heating",
    },
    "rxn_071_organic_reaction": {
        "verified_equation": "C6H6(l) + HNO3(l) -> C6H5NO2(l) + H2O(l)",
        "verified_source_section": "Vogel's Textbook of Practical Organic Chemistry 5th Ed. (1989), Exp. 6.11, p. 855 (50-60 °C)",
        "verified_conditions": "Mixed H2SO4/HNO3 acid, 40-65 °C (optimal 55 °C)",
    },
    "rxn_077_biochemical_reaction": {
        "verified_equation": "C10H16N5O13P3(aq) + H2O(l) -> C10H15N5O10P2(aq) + H3PO4(aq)",
        "verified_source_section": "Rhea Database RHEA:13065 (EC 3.6.1.3 ATPase), ΔG°' = -30.5 kJ/mol",
        "verified_conditions": "Physiological aqueous buffer, 37 °C, ATPase",
    },
    "rxn_078_biochemical_reaction": {
        "verified_equation": "CO(NH2)2(aq) + H2O(l) -> CO2(g) + 2NH3(aq)",
        "verified_source_section": "Rhea Database RHEA:16113 (EC 3.5.1.5 Urease)",
        "verified_conditions": "Aqueous buffer, 37 °C, urease",
    },
    "rxn_079_biochemical_reaction": {
        "verified_equation": "C12H22O11(aq) + H2O(l) -> 2C6H12O6(aq)",
        "verified_source_section": "Rhea Database RHEA:10068 (EC 3.2.1.26 Invertase / Sucrase)",
        "verified_conditions": "Aqueous medium, 37 °C, sucrase",
    },
}

_REGISTRY_IDS: set[str] = {r["reaction_id"] for r in VERIFIED_REACTIONS_REGISTRY}
REFERENCE_VERIFIED_IDS: set[str] = set(MANUALLY_VERIFIED_REACTION_SOURCES.keys()) & _REGISTRY_IDS

REFERENCE_CHECK_SCOPE_AR = (
    f"تم التحقق يدوياً من المصدر العلمي الأساسي والمعادلة والظروف لـ ({len(REFERENCE_VERIFIED_IDS)} / {len(_REGISTRY_IDS)}) تفاعلاً محورياً؛ "
    "أما بقية التفاعلات فمدققة آلياً بالكامل (حفظ الذرات، حفظ الشحنة، الكتلة المولية، اتساق الأطوار، وإشارة ΔH) "
    "مع تسجيل المرجع الببليوغرافي، ولا يُعرض أي تفاعل على أنه مراجَع الصفحة يدوياً لمجرد وجود رابط مرجعي."
)

RULE_REPRESENTATIVE_PROBES: list[tuple[str, list[str], list[str]]] = [
    ("combustion", ["C5H12", "O2"], ["CO2", "H2O"]),
    ("acid_base_neutralization", ["HBr", "LiOH"], ["LiBr", "H2O"]),
    ("acid_carbonate", ["CH3COOH", "K2CO3"], ["KCH3COO", "H2O", "CO2"]),
    ("acid_oxide", ["HBr", "MgO"], ["MgBr2", "H2O"]),
    ("acid_ammonia", ["HBr", "NH3"], ["NH4Br"]),
    ("acid_active_metal", ["H2SO4", "Mg"], ["MgSO4", "H2"]),
    ("alkali_water", ["K", "H2O"], ["KOH", "H2"]),
    ("metal_oxygen", ["Li", "O2"], ["Li2O"]),
    ("metal_halogen", ["Ca", "Cl2"], ["CaCl2"]),
    ("salt_salt_precipitation", ["Pb(NO3)2", "KI"], ["PbI2", "KNO3"]),
    ("salt_hydroxide_precipitation", ["AlCl3", "NaOH"], ["Al(OH)3", "NaCl"]),
    ("salt_metal_displacement", ["Zn", "Pb(NO3)2"], ["Zn(NO3)2", "Pb"]),
    ("esterification", ["HCOOH", "CH3OH"], ["C2H4O2", "H2O"]),
]


def _extract_formal_charge(formula: str) -> int:
    """Extracts net formal charge from a species formula (0 for neutral species)."""
    s = (formula or "").strip()
    if not s:
        return 0
    m = re.search(r"(\d*)([+⊕⁺]|[-⊖⁻])$", s)
    if m:
        mag = int(m.group(1)) if m.group(1) else 1
        sign = 1 if m.group(2) in "+⊕⁺" else -1
        return sign * mag
    return 0


def _side_net_charge(species_list: list[dict[str, Any]]) -> int:
    total_q = 0
    for sp in species_list:
        coeff = int(sp.get("coefficient", 1) or 1)
        q = int(sp.get("charge", _extract_formal_charge(sp.get("formula", ""))))
        total_q += coeff * q
    return total_q


def _count_side_atoms(species_list: list[dict[str, Any]]) -> dict[str, int]:
    totals: dict[str, int] = {}
    for sp in species_list:
        coeff = int(sp.get("coefficient", 1) or 1)
        atoms = parse_chemical_formula(sp.get("formula", ""))
        for el, count in atoms.items():
            totals[el] = totals.get(el, 0) + coeff * count
    return totals


def _side_molar_mass_sum(species_list: list[dict[str, Any]]) -> float:
    total = 0.0
    for sp in species_list:
        coeff = float(sp.get("coefficient", 1) or 1)
        mm = calculate_molar_mass(sp.get("formula", ""))
        total += coeff * mm
    return round(total, 3)


def get_reaction_verification_metadata(rxn_id: str) -> dict[str, Any]:
    """
    Returns transparent verification status for a reaction ID, strictly distinguishing
    manually source-verified reactions from format-checked bibliographic citations.
    """
    if rxn_id in MANUALLY_VERIFIED_REACTION_SOURCES:
        info = MANUALLY_VERIFIED_REACTION_SOURCES[rxn_id]
        return {
            "reference_manually_verified": True,
            "verification_status": "manually_verified_against_source",
            "verification_badge_ar": "موثق ومراجَع يدوياً من المصدر (Manually Verified Against Primary Source)",
            "verified_source_section": info["verified_source_section"],
            "verified_equation": info["verified_equation"],
            "verified_conditions": info["verified_conditions"],
        }
    return {
        "reference_manually_verified": False,
        "verification_status": "algorithmically_audited_citation_recorded",
        "verification_badge_ar": (
            "مرجع ببليوغرافي مسجل + تحقق برمجي صارم من حفظ الذرات والشحنة والأطوار (لم تُراجع الصفحة يدوياً بعد)"
        ),
        "verified_source_section": "",
        "verified_equation": "",
        "verified_conditions": "",
    }


def audit_single_reaction(rxn: dict[str, Any]) -> dict[str, Any]:
    """
    Runs an independent 6-point scientific audit on a single reaction record.
    """
    rxn_id = rxn.get("reaction_id", "unknown")
    reactants = rxn.get("reactants", [])
    products = rxn.get("products", [])
    conditions = rxn.get("conditions", {}) or {}
    spec = rxn.get("reference_conditions_spec", {}) or {}
    category = rxn.get("category", "")

    errors: list[str] = []
    warnings: list[str] = []
    checks: dict[str, bool] = {}

    # 1. Atom & Mass Conservation
    r_atoms = _count_side_atoms(reactants)
    p_atoms = _count_side_atoms(products)
    atoms_balanced = bool(r_atoms) and (r_atoms == p_atoms)
    r_mass = _side_molar_mass_sum(reactants)
    p_mass = _side_molar_mass_sum(products)
    mass_balanced = abs(r_mass - p_mass) <= 0.05
    checks["atom_and_mass_conservation"] = atoms_balanced and mass_balanced
    if not atoms_balanced:
        errors.append(f"[{rxn_id}] اختلال موازنة الذرات: المتفاعلات {r_atoms} مقابل النواتج {p_atoms}.")
    if not mass_balanced:
        errors.append(f"[{rxn_id}] اختلال حفظ الكتلة المولية: {r_mass} g مقابل {p_mass} g.")

    # 2. Charge Conservation
    r_charge = _side_net_charge(reactants)
    p_charge = _side_net_charge(products)
    charge_balanced = r_charge == p_charge
    checks["charge_conservation"] = charge_balanced
    if not charge_balanced:
        errors.append(f"[{rxn_id}] اختلال حفظ الشحنة الكهربائية: شحنة المتفاعلات ({r_charge}) مقابل شحنة النواتج ({p_charge}).")

    # 3. Physical State Validity & Phase Consistency
    states_ok = True
    cond_text = " ".join(str(v) for v in conditions.values())
    for side_name, sp_list in (("reactants", reactants), ("products", products)):
        for sp in sp_list:
            st = sp.get("state", "")
            formula = sp.get("formula", "")
            if st not in VALID_STATES:
                states_ok = False
                errors.append(f"[{rxn_id}] حالة فيزيائية غير صالحة ({st}) للمادة {formula} في {side_name}.")
            if category == "precipitation" and side_name == "products":
                if formula in {"AgCl", "BaSO4", "CaCO3", "Fe(OH)3", "Cu(OH)2"} and st != "s":
                    states_ok = False
                    errors.append(f"[{rxn_id}] الراسب {formula} يجب أن يكون في الحالة الصلبة (s) وليس ({st}).")
            if category in {"neutralization", "precipitation", "single_displacement"} and "طور غازي" not in cond_text:
                if formula in {"Mg(OH)2", "CuO", "CaCO3"} and st == "aq":
                    states_ok = False
                    errors.append(f"[{rxn_id}] المادة شحيحة الذوبان {formula} صُنفت خطأ كمحلول مائي (aq).")
    checks["physical_state_consistency"] = states_ok

    # 4. Thermodynamic Sign & Enthalpy Consistency
    thermo_ar = str(conditions.get("thermodynamics_ar", "") or "")
    thermo_ok = True
    has_quantitative_enthalpy = False
    delta_h_value: float | None = None
    if thermo_ar:
        m_num = re.search(r"([+-]?\d+(?:\.\d+)?)\s*kJ", thermo_ar, flags=re.IGNORECASE)
        if m_num:
            has_quantitative_enthalpy = True
            delta_h_value = float(m_num.group(1))
            if ("طارد" in thermo_ar or "باعث" in thermo_ar) and delta_h_value > 0:
                thermo_ok = False
                errors.append(f"[{rxn_id}] تناقض حراري: وُصف التفاعل بأنه طارد للحرارة بينما ΔH موجبة ({delta_h_value} kJ).")
            if "ماص" in thermo_ar and delta_h_value < 0:
                thermo_ok = False
                errors.append(f"[{rxn_id}] تناقض حراري: وُصف التفاعل بأنه ماص للحرارة بينما ΔH سالبة ({delta_h_value} kJ).")
    checks["thermodynamic_consistency"] = thermo_ok

    # 5. Reference Conditions Operating Window Validity
    window_ok = True
    if spec:
        t_min = float(spec.get("temp_c_min", 0.0))
        t_max = float(spec.get("temp_c_max", 0.0))
        t_opt = float(spec.get("temp_c_optimal", 0.0))
        p_min = float(spec.get("pressure_atm_min", 0.0))
        p_max = float(spec.get("pressure_atm_max", 0.0))
        p_opt = float(spec.get("pressure_atm_optimal", 0.0))
        if not (t_min <= t_opt <= t_max):
            window_ok = False
            errors.append(f"[{rxn_id}] نطاق الحرارة المرجعي غير متسق: [{t_min}, {t_opt}, {t_max}] °C.")
        if not (0.0 < p_min <= p_opt <= p_max):
            window_ok = False
            errors.append(f"[{rxn_id}] نطاق الضغط المرجعي غير متسق: [{p_min}, {p_opt}, {p_max}] atm.")
    else:
        window_ok = False
        errors.append(f"[{rxn_id}] يفتقر إلى مواصفات الظروف المرجعية (reference_conditions_spec).")
    checks["reference_window_validity"] = window_ok

    # 6. Bibliographic Reference Completeness & Manual Verification Status
    source_id = str(rxn.get("source_id", "") or "").strip()
    source_text = str(rxn.get("source", "") or "").strip()
    ref_url = str(rxn.get("reference_url", "") or "").strip()
    doc_type = str(rxn.get("documentation_type", "") or "").strip()
    ref_ok = (
        len(source_id) >= 5
        and len(source_text) >= 20
        and ref_url.startswith("https://")
        and doc_type in {"experimentally_documented", "literature_referenced", "theoretical_model"}
    )
    if not ref_ok:
        errors.append(f"[{rxn_id}] نقص في بيانات التوثيق المرجعي (source_id / source / reference_url).")
    checks["bibliographic_reference_integrity"] = ref_ok

    verif_meta = get_reaction_verification_metadata(rxn_id)

    return {
        "reaction_id": rxn_id,
        "reaction_name_ar": rxn.get("reaction_name_ar", ""),
        "category": category,
        "passed": len(errors) == 0,
        "checks": checks,
        "stoichiometric_mass_g": r_mass,
        "net_charge": r_charge,
        "atom_counts": r_atoms,
        "has_quantitative_enthalpy": has_quantitative_enthalpy,
        "delta_h_kj": delta_h_value,
        "documentation_type": doc_type,
        "is_theoretical_model": doc_type == "theoretical_model",
        "citation_exists": bool(rxn.get("source") and rxn.get("source_id")),
        "source_manually_verified": verif_meta["reference_manually_verified"],
        "reference_manually_verified": verif_meta["reference_manually_verified"],
        "verified_conditions_linked": bool(rxn.get("reference_conditions_spec")),
        "verification_status": verif_meta["verification_status"],
        "verification_badge_ar": verif_meta["verification_badge_ar"],
        "verified_source_section": verif_meta["verified_source_section"],
        "source_id": source_id,
        "source": source_text,
        "reference_url": ref_url,
        "errors": errors,
        "warnings": warnings,
    }


def audit_registry_and_rules() -> dict[str, Any]:
    """
    Performs a complete independent scientific audit across all registry reactions
    and all 13 general chemistry inference rules.
    """
    reaction_reports: list[dict[str, Any]] = []
    all_errors: list[str] = []
    source_ids: set[str] = set()
    sources: set[str] = set()
    doc_type_counts: dict[str, int] = {}
    category_counts: dict[str, int] = {}
    quantitative_enthalpy_count = 0

    for rxn in VERIFIED_REACTIONS_REGISTRY:
        rep = audit_single_reaction(rxn)
        reaction_reports.append(rep)
        all_errors.extend(rep["errors"])

        sid = rep["source_id"]
        src = rep["source"]
        if sid in source_ids:
            all_errors.append(f"[{rep['reaction_id']}] تكرار معرّف المرجع source_id: {sid}")
        source_ids.add(sid)
        if src in sources:
            all_errors.append(f"[{rep['reaction_id']}] تكرار نص المرجع العلمي source: {src}")
        sources.add(src)

        dt = rep["documentation_type"]
        doc_type_counts[dt] = doc_type_counts.get(dt, 0) + 1
        cat = rep["category"]
        category_counts[cat] = category_counts.get(cat, 0) + 1
        if rep["has_quantitative_enthalpy"]:
            quantitative_enthalpy_count += 1

    # Audit the 13 general-chemistry rules via non-registry probe reactions
    rule_reports: list[dict[str, Any]] = []
    for rule_id, reactant_formulas, expected_products in RULE_REPRESENTATIVE_PROBES:
        rule_meta = RULES.get(rule_id, {})
        inputs = [{"substance": {"symbol_or_formula": f}} for f in reactant_formulas]
        pred = predict_reaction(inputs)
        got_products = [p["formula"] for p in pred.get("products", [])]
        atom_cons = pred.get("atom_conservation", {}) or {}
        eq_unicode = pred.get("equation_unicode", "") or pred.get("equation_text", "")
        rule_passed = (
            bool(rule_meta)
            and pred.get("supported") is True
            and pred.get("mode") in {"rule_inferred", "documented"}
            and set(got_products) == set(expected_products)
            and bool(atom_cons.get("is_balanced"))
            and bool(rule_meta.get("source"))
            and str(rule_meta.get("reference_url", "")).startswith("https://")
        )
        if not rule_passed:
            all_errors.append(
                f"[Rule:{rule_id}] فشل تدقيق القاعدة على المتفاعلات {reactant_formulas}: النواتج {got_products} مقابل المتوقع {expected_products}."
            )
        rule_reports.append(
            {
                "rule_id": rule_id,
                "rule_name_ar": rule_meta.get("rule_name_ar", rule_id),
                "probe_reactants": reactant_formulas,
                "predicted_products": got_products,
                "balanced_equation": eq_unicode,
                "source_id": rule_meta.get("source_id", ""),
                "source": rule_meta.get("source", ""),
                "passed": rule_passed,
            }
        )

    passed_reactions = sum(1 for r in reaction_reports if r["passed"])
    passed_rules = sum(1 for r in rule_reports if r["passed"])
    manually_verified_count = sum(1 for r in reaction_reports if r["reference_manually_verified"])

    return {
        "overall_passed": len(all_errors) == 0,
        "total_reactions": len(reaction_reports),
        "passed_reactions": passed_reactions,
        "manually_verified_references": manually_verified_count,
        "algorithmically_verified_only_references": len(reaction_reports) - manually_verified_count,
        "reference_check_scope_ar": REFERENCE_CHECK_SCOPE_AR,
        "unique_source_ids": len(source_ids),
        "unique_bibliographic_sources": len(sources),
        "quantitative_enthalpy_reactions": quantitative_enthalpy_count,
        "documentation_type_breakdown": doc_type_counts,
        "category_breakdown": category_counts,
        "total_rules": len(RULES),
        "passed_rules": passed_rules,
        "errors": all_errors,
        "reaction_reports": reaction_reports,
        "rule_reports": rule_reports,
    }
