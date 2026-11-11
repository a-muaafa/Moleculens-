"""
Concentration and solution preparation calculator engine for MolecuLens.
Supports Molarity, Normality, ppm/ppb, Mass/Volume %, and Dilution calculations.
"""

from typing import Dict, Any, Optional
import math


class ConcentrationCalculator:
    """Calculates chemical solution concentrations and preparation recipes."""

    @staticmethod
    def calculate_molarity(
        mass_g: Optional[float] = None,
        molar_mass_g_mol: Optional[float] = None,
        volume_l: Optional[float] = None,
        molarity_m: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Calculates the missing parameter in M = m / (MW * V).
        Exactly one parameter of (mass_g, volume_l, molarity_m) should be None.
        molar_mass_g_mol is required and must be positive.
        """
        if molar_mass_g_mol is None or molar_mass_g_mol <= 0:
            return {"success": False, "error": "Molar mass must be greater than zero."}

        none_count = sum(1 for v in [mass_g, volume_l, molarity_m] if v is None)
        if none_count != 1:
            return {
                "success": False,
                "error": "Provide exactly two of mass, volume, and molarity to solve for the third.",
            }

        if mass_g is None:
            if volume_l <= 0 or molarity_m <= 0:
                return {"success": False, "error": "Volume and molarity must be positive."}
            mass_calc = molarity_m * molar_mass_g_mol * volume_l
            moles = molarity_m * volume_l
            return {
                "success": True,
                "solved_for": "mass_g",
                "mass_g": round(mass_calc, 4),
                "moles": round(moles, 6),
                "molarity_m": molarity_m,
                "volume_l": volume_l,
                "molar_mass_g_mol": molar_mass_g_mol,
                "formula": "m = M * MW * V",
                "recipe_en": f"Weigh {mass_calc:.4f} g of solute, dissolve in ~{volume_l * 0.7 * 1000:.1f} mL of distilled water, transfer to a {volume_l * 1000:.1f} mL volumetric flask, and dilute to the calibration mark.",
                "recipe_ar": f"قم بوزن {mass_calc:.4f} جم من المادة المذابة، وقم بإذابتها في حوالي {volume_l * 0.7 * 1000:.1f} مل من الماء المقطر، ثم انقل المحلول إلى دورق عياري سعة {volume_l * 1000:.1f} مل وأكمل الحجم حتى علامة المعايرة.",
            }

        elif molarity_m is None:
            if mass_g <= 0 or volume_l <= 0:
                return {"success": False, "error": "Mass and volume must be positive."}
            moles = mass_g / molar_mass_g_mol
            molarity_calc = moles / volume_l
            return {
                "success": True,
                "solved_for": "molarity_m",
                "molarity_m": round(molarity_calc, 6),
                "moles": round(moles, 6),
                "mass_g": mass_g,
                "volume_l": volume_l,
                "molar_mass_g_mol": molar_mass_g_mol,
                "formula": "M = (m / MW) / V",
                "recipe_en": f"A solution containing {mass_g} g solute in {volume_l * 1000:.1f} mL gives a concentration of {molarity_calc:.4f} M.",
                "recipe_ar": f"محلول يحتوي على {mass_g} جم من المادة المذابة في {volume_l * 1000:.1f} مل يعطي تركيزاً مقداره {molarity_calc:.4f} مولار.",
            }

        else:
            if mass_g <= 0 or molarity_m <= 0:
                return {"success": False, "error": "Mass and molarity must be positive."}
            moles = mass_g / molar_mass_g_mol
            volume_calc = moles / molarity_m
            return {
                "success": True,
                "solved_for": "volume_l",
                "volume_l": round(volume_calc, 4),
                "volume_ml": round(volume_calc * 1000.0, 2),
                "moles": round(moles, 6),
                "mass_g": mass_g,
                "molarity_m": molarity_m,
                "molar_mass_g_mol": molar_mass_g_mol,
                "formula": "V = (m / MW) / M",
                "recipe_en": f"Dissolve {mass_g} g solute to a final volume of {volume_calc * 1000:.1f} mL to yield a {molarity_m} M solution.",
                "recipe_ar": f"قم بإذابة {mass_g} جم من المادة المذابة لإكمال الحجم النهائي إلى {volume_calc * 1000:.1f} مل للحصول على محلول بتركيز {molarity_m} مولار.",
            }

    @staticmethod
    def calculate_normality(
        molarity_m: float,
        equivalents_per_mole: float,
    ) -> Dict[str, Any]:
        """Calculates Normality N = M * n_eq."""
        if molarity_m <= 0 or equivalents_per_mole <= 0:
            return {"success": False, "error": "Molarity and valence equivalents must be positive."}
        normality = molarity_m * equivalents_per_mole
        return {
            "success": True,
            "molarity_m": molarity_m,
            "equivalents_per_mole": equivalents_per_mole,
            "normality_n": round(normality, 6),
            "formula": "N = M * n",
            "explanation_en": f"Normality = {molarity_m} M * {equivalents_per_mole} eq/mol = {normality:.4f} N",
            "explanation_ar": f"العيارية (Normality) = {molarity_m} مولار * {equivalents_per_mole} مكافئ/مول = {normality:.4f} عياري (N)",
        }

    @staticmethod
    def calculate_ppm_ppb(
        solute_mg: float,
        solution_volume_l: float,
    ) -> Dict[str, Any]:
        """Calculates ppm (mg/L) and ppb (ug/L) in aqueous solutions."""
        if solute_mg <= 0 or solution_volume_l <= 0:
            return {"success": False, "error": "Solute mass and solution volume must be positive."}
        ppm = solute_mg / solution_volume_l
        ppb = ppm * 1000.0
        return {
            "success": True,
            "ppm": round(ppm, 4),
            "ppb": round(ppb, 2),
            "solute_mg": solute_mg,
            "solution_volume_l": solution_volume_l,
            "formula": "ppm = mg / L; ppb = ug / L",
            "explanation_en": f"{solute_mg} mg solute in {solution_volume_l} L gives {ppm:.2f} ppm ({ppb:.1f} ppb).",
            "explanation_ar": f"{solute_mg} مجم من المذاب في {solution_volume_l} لتر يعطي {ppm:.2f} جزء في المليون (ppm) و {ppb:.1f} جزء في المليار (ppb).",
        }

    @staticmethod
    def calculate_percentage(
        percent_type: str,
        solute_amount: float,
        solution_amount: float,
    ) -> Dict[str, Any]:
        """
        Calculates mass or volume percentage:
        - 'w/w': (mass solute g / mass solution g) * 100
        - 'w/v': (mass solute g / volume solution mL) * 100
        - 'v/v': (volume solute mL / volume solution mL) * 100
        """
        if solute_amount <= 0 or solution_amount <= 0:
            return {"success": False, "error": "Solute and solution quantities must be positive."}
        if solute_amount > solution_amount and percent_type in ["w/w", "v/v"]:
            return {"success": False, "error": "Solute amount cannot exceed total solution amount."}

        percentage = (solute_amount / solution_amount) * 100.0
        unit_map = {
            "w/w": ("g solute", "g solution", "% (w/w)"),
            "w/v": ("g solute", "mL solution", "% (w/v)"),
            "v/v": ("mL solute", "mL solution", "% (v/v)"),
        }
        units = unit_map.get(percent_type, ("units", "units", "%"))

        return {
            "success": True,
            "percent_type": percent_type,
            "percentage": round(percentage, 4),
            "solute_amount": solute_amount,
            "solution_amount": solution_amount,
            "formula": f"% ({percent_type}) = ({units[0]} / {units[1]}) * 100",
            "explanation_en": f"{solute_amount} {units[0]} in {solution_amount} {units[1]} = {percentage:.2f}{units[2]}",
            "explanation_ar": f"التركيز المئوي ({percent_type}) = {percentage:.2f}%",
        }

    @staticmethod
    def calculate_dilution(
        c1: Optional[float] = None,
        v1: Optional[float] = None,
        c2: Optional[float] = None,
        v2: Optional[float] = None,
    ) -> Dict[str, Any]:
        """
        Calculates dilution parameters using C1 * V1 = C2 * V2.
        Exactly one parameter of (c1, v1, c2, v2) must be None.
        Concentrations and volumes must be in consistent units.
        """
        params = [c1, v1, c2, v2]
        none_count = sum(1 for p in params if p is None)
        if none_count != 1:
            return {
                "success": False,
                "error": "Provide exactly three of C1, V1, C2, and V2 to calculate the fourth.",
            }

        if c1 is not None and c1 <= 0:
            return {"success": False, "error": "C1 must be positive."}
        if v1 is not None and v1 <= 0:
            return {"success": False, "error": "V1 must be positive."}
        if c2 is not None and c2 <= 0:
            return {"success": False, "error": "C2 must be positive."}
        if v2 is not None and v2 <= 0:
            return {"success": False, "error": "V2 must be positive."}

        if v1 is None:
            if c1 <= c2:
                return {
                    "success": False,
                    "error": "Stock concentration C1 must be strictly greater than target concentration C2.",
                }
            v1_calc = (c2 * v2) / c1
            solvent_added = v2 - v1_calc
            return {
                "success": True,
                "solved_for": "v1",
                "c1": c1,
                "v1": round(v1_calc, 4),
                "c2": c2,
                "v2": v2,
                "solvent_to_add": round(solvent_added, 4),
                "formula": "V1 = (C2 * V2) / C1",
                "recipe_en": f"Pipette {v1_calc:.2f} mL of stock solution (C1={c1}) and add {solvent_added:.2f} mL of solvent to reach final volume {v2:.2f} mL (C2={c2}).",
                "recipe_ar": f"اسحب {v1_calc:.2f} مل من المحلول المركز (C1={c1}) وأضف {solvent_added:.2f} مل من المذيب للوصول إلى الحجم النهائي {v2:.2f} مل (C2={c2}).",
            }

        elif c2 is None:
            if v1 >= v2:
                return {
                    "success": False,
                    "error": "Final volume V2 must be strictly greater than initial volume V1 for dilution.",
                }
            c2_calc = (c1 * v1) / v2
            solvent_added = v2 - v1
            return {
                "success": True,
                "solved_for": "c2",
                "c1": c1,
                "v1": v1,
                "c2": round(c2_calc, 4),
                "v2": v2,
                "solvent_to_add": round(solvent_added, 4),
                "formula": "C2 = (C1 * V1) / V2",
                "recipe_en": f"Diluting {v1:.2f} mL of stock ({c1}) to {v2:.2f} mL yields a final concentration of {c2_calc:.4f}.",
                "recipe_ar": f"تخفيف {v1:.2f} مل من المحلول المركز ({c1}) إلى {v2:.2f} مل ينتج تركيزاً نهائياً قدره {c2_calc:.4f}.",
            }

        elif v2 is None:
            if c1 <= c2:
                return {
                    "success": False,
                    "error": "Stock concentration C1 must be strictly greater than target concentration C2.",
                }
            v2_calc = (c1 * v1) / c2
            solvent_added = v2_calc - v1
            return {
                "success": True,
                "solved_for": "v2",
                "c1": c1,
                "v1": v1,
                "c2": c2,
                "v2": round(v2_calc, 4),
                "solvent_to_add": round(solvent_added, 4),
                "formula": "V2 = (C1 * V1) / C2",
                "recipe_en": f"Dilute {v1:.2f} mL of stock ({c1}) by adding {solvent_added:.2f} mL solvent to achieve a final volume of {v2_calc:.2f} mL at {c2}.",
                "recipe_ar": f"خفف {v1:.2f} مل من المحلول المركز ({c1}) بإضافة {solvent_added:.2f} مل من المذيب للوصول إلى حجم نهائي {v2_calc:.2f} مل بتركيز {c2}.",
            }

        else:
            v1_calc = v1
            c1_calc = (c2 * v2) / v1
            return {
                "success": True,
                "solved_for": "c1",
                "c1": round(c1_calc, 4),
                "v1": v1,
                "c2": c2,
                "v2": v2,
                "formula": "C1 = (C2 * V2) / V1",
                "recipe_en": f"Required stock concentration is {c1_calc:.4f}.",
                "recipe_ar": f"التركيز المطلوب للمحلول المركز هو {c1_calc:.4f}.",
            }

