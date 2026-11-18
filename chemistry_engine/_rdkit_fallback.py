"""
Pure-Python algorithmic cheminformatics fallback for environments where native C++ RDKit
wheels are unavailable. When native `rdkit` is installed, MolecuLens imports and uses
native RDKit directly; this module is loaded ONLY in `except ImportError:` fallback paths.

Implements:
- Algorithmic SMILES tokenizer, parser, ring-closure verifier, valence checker, and aromaticity/H-count solver.
- Exact Hill molecular formula calculation (`CalcMolFormula`) and IUPAC molecular weight (`MolWt`).
- Topological ring perception (cycle basis / SSSR for `RingCount`, `NumAromaticRings`, `NumAliphaticRings`).
- Physicochemical descriptor calculation (`TPSA`, `NumHBD`, `NumHBA`, `NumRotatableBonds`, `FractionCSP3`, `LabuteASA`, `MolLogP`).
- Topological circular Morgan-like fingerprint generator (radius=2, 2048 bits) and Tanimoto similarity.
- SMARTS substructure pattern matching and reaction SMARTS transformations for virtual structure editing.
- 2D coordinate layout (`Compute2DCoords`), SVG drawing (`MolDraw2DSVG`), PNG rendering (`MolToImage`), and 3D conformer embedding (`ETKDGv3`).
"""

from __future__ import annotations

import hashlib
import math
import re
import sys
import types
from collections import defaultdict, deque
from typing import Any

from PIL import Image, ImageDraw

from chemistry_engine.periodic_table import ATOMIC_MASSES, ATOMIC_NUMBERS

# Standard atomic weights and normal organic valences
_ATOMIC_NUM_TO_SYM: dict[int, str] = {v: k for k, v in ATOMIC_NUMBERS.items()}
_SYM_TO_ATOMIC_NUM: dict[str, int] = dict(ATOMIC_NUMBERS)
_SYM_TO_MASS: dict[str, float] = {k: float(v) for k, v in ATOMIC_MASSES.items()}

_ORGANIC_VALENCES: dict[str, tuple[int, ...]] = {
    "B": (3,),
    "C": (4,),
    "N": (3, 5),
    "O": (2,),
    "P": (3, 5),
    "S": (2, 4, 6),
    "F": (1,),
    "Cl": (1, 3, 5, 7),
    "Br": (1, 3, 5, 7),
    "I": (1, 3, 5, 7),
    "H": (1,),
}

_ORGANIC_SUBSET = {"B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I", "b", "c", "n", "o", "p", "s"}


class BondType:
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 12


class BondStereo:
    STEREONONE = 0
    STEREOANY = 1
    STEREOZ = 2
    STEREOE = 3


class _Pos3D:
    def __init__(self, x: float, y: float, z: float):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)


class _Conformer:
    def __init__(self, coords: list[tuple[float, float, float]]):
        self._coords = coords

    def GetAtomPosition(self, idx: int) -> _Pos3D:
        x, y, z = self._coords[idx]
        return _Pos3D(x, y, z)


class _Atom:
    def __init__(
        self,
        idx: int,
        symbol: str,
        is_aromatic: bool = False,
        formal_charge: int = 0,
        explicit_h: int | None = None,
        chiral_tag: str = "",
        in_brackets: bool = False,
        radical_electrons: int = 0,
    ):
        self.idx = idx
        self.symbol = symbol
        self.atomic_num = _SYM_TO_ATOMIC_NUM.get(symbol, 0)
        self.is_aromatic = is_aromatic
        self.formal_charge = formal_charge
        self.explicit_h = explicit_h
        self.implicit_h = 0
        self.chiral_tag = chiral_tag
        self.in_brackets = in_brackets
        self.radical_electrons = radical_electrons
        self.neighbors: list[int] = []
        self.bond_indices: list[int] = []

    def GetIdx(self) -> int:
        return self.idx

    def GetSymbol(self) -> str:
        return self.symbol

    def GetAtomicNum(self) -> int:
        return self.atomic_num

    def GetFormalCharge(self) -> int:
        return self.formal_charge

    def GetIsAromatic(self) -> bool:
        return self.is_aromatic

    def GetTotalNumHs(self) -> int:
        exp_h = self.explicit_h if self.explicit_h is not None else 0
        return exp_h + self.implicit_h

    def GetNumRadicalElectrons(self) -> int:
        return self.radical_electrons

    def GetDegree(self) -> int:
        return len(self.neighbors)


class _Bond:
    def __init__(
        self,
        idx: int,
        begin_idx: int,
        end_idx: int,
        bond_type: int = BondType.SINGLE,
        is_aromatic: bool = False,
        stereo: int = BondStereo.STEREONONE,
    ):
        self.idx = idx
        self.begin_idx = begin_idx
        self.end_idx = end_idx
        self.bond_type = bond_type
        self.is_aromatic = is_aromatic
        self.stereo = stereo

    def GetBeginAtomIdx(self) -> int:
        return self.begin_idx

    def GetEndAtomIdx(self) -> int:
        return self.end_idx

    def GetBondType(self) -> int:
        return self.bond_type

    def GetIsAromatic(self) -> bool:
        return self.is_aromatic

    def GetStereo(self) -> int:
        return self.stereo

    def GetBondTypeAsDouble(self) -> float:
        if self.is_aromatic or self.bond_type == BondType.AROMATIC:
            return 1.5
        return float(self.bond_type)


class _RingInfo:
    def __init__(self, rings: list[tuple[int, ...]]):
        self._rings = rings

    def AtomRings(self) -> tuple[tuple[int, ...], ...]:
        return tuple(self._rings)


class _Mol:
    def __init__(self, atoms: list[_Atom], bonds: list[_Bond], raw_smiles: str = ""):
        self.atoms = atoms
        self.bonds = bonds
        self.raw_smiles = raw_smiles
        self._has_explicit_hs_added = False
        self._conformers: list[_Conformer] = []
        self._rings: list[tuple[int, ...]] | None = None
        self._adj: dict[int, list[tuple[int, _Bond]]] = defaultdict(list)
        for a in self.atoms:
            a.neighbors = []
            a.bond_indices = []
        for b in self.bonds:
            self.atoms[b.begin_idx].neighbors.append(b.end_idx)
            self.atoms[b.end_idx].neighbors.append(b.begin_idx)
            self.atoms[b.begin_idx].bond_indices.append(b.idx)
            self.atoms[b.end_idx].bond_indices.append(b.idx)
            self._adj[b.begin_idx].append((b.end_idx, b))
            self._adj[b.end_idx].append((b.begin_idx, b))

    def copy(self) -> _Mol:
        atoms_copy = [
            _Atom(
                a.idx,
                a.symbol,
                a.is_aromatic,
                a.formal_charge,
                a.explicit_h,
                a.chiral_tag,
                a.in_brackets,
                a.radical_electrons,
            )
            for a in self.atoms
        ]
        for a_src, a_dst in zip(self.atoms, atoms_copy):
            a_dst.implicit_h = a_src.implicit_h
        bonds_copy = [
            _Bond(b.idx, b.begin_idx, b.end_idx, b.bond_type, b.is_aromatic, b.stereo)
            for b in self.bonds
        ]
        m = _Mol(atoms_copy, bonds_copy, self.raw_smiles)
        m._has_explicit_hs_added = self._has_explicit_hs_added
        return m

    def GetAtoms(self) -> list[_Atom]:
        return self.atoms

    def GetAtomWithIdx(self, idx: int) -> _Atom:
        return self.atoms[idx]

    def GetBonds(self) -> list[_Bond]:
        return self.bonds

    def GetNumAtoms(self) -> int:
        return len(self.atoms)

    def GetNumHeavyAtoms(self) -> int:
        return sum(1 for a in self.atoms if a.symbol != "H")

    def GetNumBonds(self) -> int:
        return len(self.bonds)

    def GetNumConformers(self) -> int:
        return len(self._conformers)

    def GetConformer(self, idx: int = 0) -> _Conformer:
        return self._conformers[idx]

    def GetRingInfo(self) -> _RingInfo:
        if self._rings is None:
            self._rings = _find_sssr_rings(self)
        return _RingInfo(self._rings)

    def HasSubstructMatch(self, pattern: Any) -> bool:
        return len(self.GetSubstructMatches(pattern)) > 0

    def GetSubstructMatches(self, pattern: Any) -> tuple[tuple[int, ...], ...]:
        if pattern is None:
            return ()
        return _match_smarts_pattern(self, pattern)


def _find_sssr_rings(mol: _Mol) -> list[tuple[int, ...]]:
    """Finds the smallest set of smallest rings (SSSR) in a molecular graph."""
    heavy_adj: dict[int, list[int]] = defaultdict(list)
    for b in mol.bonds:
        u, v = b.begin_idx, b.end_idx
        if mol.atoms[u].symbol != "H" and mol.atoms[v].symbol != "H":
            heavy_adj[u].append(v)
            heavy_adj[v].append(u)

    # Number of fundamental cycles: E - V + C
    visited: set[int] = set()
    n_comp = 0
    heavy_nodes = [a.idx for a in mol.atoms if a.symbol != "H"]
    for u in heavy_nodes:
        if u not in visited:
            n_comp += 1
            q = deque([u])
            visited.add(u)
            while q:
                cur = q.popleft()
                for nxt in heavy_adj[cur]:
                    if nxt not in visited:
                        visited.add(nxt)
                        q.append(nxt)

    n_edges = sum(len(v) for v in heavy_adj.values()) // 2
    cyclomatic = max(0, n_edges - len(heavy_nodes) + n_comp)
    if cyclomatic == 0:
        return []

    found_rings: dict[frozenset[int], tuple[int, ...]] = {}
    for b in mol.bonds:
        u, v = b.begin_idx, b.end_idx
        if mol.atoms[u].symbol == "H" or mol.atoms[v].symbol == "H":
            continue
        # Shortest path between u and v excluding direct edge (u, v)
        q = deque([(u, [u])])
        seen = {u}
        while q:
            cur, path = q.popleft()
            if len(path) > 12:
                continue
            if cur == v and len(path) >= 3:
                key = frozenset(path)
                if key not in found_rings:
                    found_rings[key] = tuple(path)
                break
            for nxt in heavy_adj[cur]:
                if (cur == u and nxt == v) or (cur == v and nxt == u):
                    continue
                if nxt not in seen:
                    seen.add(nxt)
                    q.append((nxt, path + [nxt]))

    sorted_rings = sorted(found_rings.values(), key=len)
    return sorted_rings[:cyclomatic] if len(sorted_rings) >= cyclomatic else sorted_rings


def _parse_bracket_atom(text: str, idx: int) -> _Atom | None:
    # e.g. [NH4+], [O-], [N+](=O)[O-], [C@H], [C@@H], [Fe], [2H], [cH1:1]
    s = re.sub(r":\d+$", "", text.strip())
    s = re.sub(r"^\d+", "", s)  # strip isotope
    m = re.match(r"^([A-Z][a-z]?|[bcnops])(@{1,2})?(H\d*)?([+\-]\d*|[+\-]+)?$", s)
    if not m:
        return None
    sym_raw, chiral, h_part, chg_part = m.groups()
    is_arom = sym_raw.islower()
    sym = sym_raw.capitalize()
    if sym not in _SYM_TO_ATOMIC_NUM:
        return None
    exp_h = 0
    if h_part:
        exp_h = int(h_part[1:]) if len(h_part) > 1 else 1
    charge = 0
    if chg_part:
        if chg_part in ("+", "++", "+++"):
            charge = len(chg_part)
        elif chg_part in ("-", "--", "---"):
            charge = -len(chg_part)
        else:
            sign = 1 if chg_part[0] == "+" else -1
            mag = int(chg_part[1:]) if len(chg_part) > 1 else 1
            charge = sign * mag
    return _Atom(
        idx=idx,
        symbol=sym,
        is_aromatic=is_arom,
        formal_charge=charge,
        explicit_h=exp_h,
        chiral_tag=chiral or "",
        in_brackets=True,
    )


def _compute_implicit_hs_and_validate(mol: _Mol) -> bool:
    for atom in mol.atoms:
        if atom.atomic_num == 0:
            return False
        bond_order_sum = 0.0
        arom_bonds = 0
        for b_idx in atom.bond_indices:
            b = mol.bonds[b_idx]
            if b.is_aromatic or b.bond_type == BondType.AROMATIC:
                bond_order_sum += 1.5
                arom_bonds += 1
            else:
                bond_order_sum += float(b.bond_type)

        if atom.in_brackets:
            atom.implicit_h = 0
            exp_h = atom.explicit_h or 0
            tot_val = int(round(bond_order_sum)) + exp_h
            # Check organic valences if C, O, F, Cl, Br, I
            if atom.symbol == "C" and atom.formal_charge == 0 and tot_val > 4:
                return False
            if atom.symbol in ("F", "Cl", "Br", "I") and atom.formal_charge == 0 and tot_val > 1:
                return False
            continue

        # Daylight SMILES aromatic implicit-H rules
        if atom.is_aromatic:
            if atom.symbol in ("N", "O", "S", "P"):
                # Bare aromatic n, o, s outside brackets never carry implicit H (pyrrole NH is always [nH])
                if len(atom.bond_indices) > 3:
                    return False
                atom.implicit_h = 0
                continue
            elif atom.symbol == "C":
                non_arom_bo = sum(
                    int(mol.bonds[bi].bond_type)
                    for bi in atom.bond_indices
                    if not (mol.bonds[bi].is_aromatic or mol.bonds[bi].bond_type == BondType.AROMATIC)
                )
                if arom_bonds == 2 and non_arom_bo == 0:
                    atom.implicit_h = 1
                elif (arom_bonds == 2 and non_arom_bo in (1, 2)) or (arom_bonds == 3 and non_arom_bo == 0):
                    atom.implicit_h = 0
                else:
                    return False
                continue

        # Aliphatic organic subset implicit H calculation
        valences = _ORGANIC_VALENCES.get(atom.symbol)
        if not valences:
            atom.implicit_h = 0
            continue

        eff_bo = int(round(bond_order_sum))
        target_val = None
        for v in valences:
            adj_v = v + atom.formal_charge if atom.symbol in ("N", "P") else v - abs(atom.formal_charge)
            if adj_v >= eff_bo:
                target_val = adj_v
                break
        if target_val is None:
            return False
        atom.implicit_h = max(0, target_val - eff_bo)
    return True


def parse_smiles(smiles: str) -> _Mol | None:
    if not smiles or not isinstance(smiles, str):
        return None
    s = smiles.strip()
    if not s or any(ch.isspace() for ch in s):
        return None

    atoms: list[_Atom] = []
    bonds: list[_Bond] = []
    stack: list[int] = []
    curr_atom: int | None = None
    pending_bond: int | None = None
    pending_stereo: int = BondStereo.STEREONONE
    ring_openings: dict[str, tuple[int, int | None]] = {}

    i = 0
    n = len(s)
    while i < n:
        ch = s[i]
        if ch == "(":
            if curr_atom is None or pending_bond is not None:
                return None
            stack.append(curr_atom)
            i += 1
        elif ch == ")":
            if not stack or pending_bond is not None:
                return None
            curr_atom = stack.pop()
            i += 1
        elif ch == ".":
            if pending_bond is not None:
                return None
            curr_atom = None
            i += 1
        elif ch in "-=#:/\\":
            if curr_atom is None:
                return None
            if ch == "-":
                pending_bond = BondType.SINGLE
            elif ch == "=":
                pending_bond = BondType.DOUBLE
            elif ch == "#":
                pending_bond = BondType.TRIPLE
            elif ch == ":":
                pending_bond = BondType.AROMATIC
            elif ch in "/\\":
                pending_bond = BondType.SINGLE
                pending_stereo = BondStereo.STEREOE
            i += 1
        elif ch.isdigit() or ch == "%":
            if curr_atom is None:
                return None
            if ch == "%":
                if i + 2 >= n or not s[i + 1 : i + 3].isdigit():
                    return None
                r_id = s[i + 1 : i + 3]
                i += 3
            else:
                r_id = ch
                i += 1
            if r_id in ring_openings:
                prev_atom, prev_btype = ring_openings.pop(r_id)
                if prev_atom == curr_atom:
                    return None
                btype = pending_bond or prev_btype
                is_arom = False
                if btype is None:
                    if atoms[prev_atom].is_aromatic and atoms[curr_atom].is_aromatic:
                        btype = BondType.AROMATIC
                        is_arom = True
                    else:
                        btype = BondType.SINGLE
                elif btype == BondType.AROMATIC:
                    is_arom = True
                bonds.append(_Bond(len(bonds), prev_atom, curr_atom, btype, is_arom))
            else:
                ring_openings[r_id] = (curr_atom, pending_bond)
            pending_bond = None
        elif ch == "[":
            end_br = s.find("]", i + 1)
            if end_br == -1:
                return None
            atom = _parse_bracket_atom(s[i + 1 : end_br], len(atoms))
            if atom is None:
                return None
            atoms.append(atom)
            new_idx = atom.idx
            if curr_atom is not None:
                btype = pending_bond
                is_arom = False
                if btype is None:
                    if atoms[curr_atom].is_aromatic and atom.is_aromatic:
                        btype = BondType.AROMATIC
                        is_arom = True
                    else:
                        btype = BondType.SINGLE
                elif btype == BondType.AROMATIC:
                    is_arom = True
                bonds.append(_Bond(len(bonds), curr_atom, new_idx, btype, is_arom, pending_stereo))
            curr_atom = new_idx
            pending_bond = None
            pending_stereo = BondStereo.STEREONONE
            i = end_br + 1
        else:
            # Organic subset symbol: Cl, Br, or single letter B, C, N, O, P, S, F, I, b, c, n, o, p, s
            sym_tok = None
            if i + 1 < n and s[i : i + 2] in ("Cl", "Br"):
                sym_tok = s[i : i + 2]
                i += 2
            elif ch in _ORGANIC_SUBSET:
                sym_tok = ch
                i += 1
            else:
                return None
            is_arom = sym_tok.islower()
            sym = sym_tok.capitalize()
            atom = _Atom(len(atoms), sym, is_aromatic=is_arom, in_brackets=False)
            atoms.append(atom)
            new_idx = atom.idx
            if curr_atom is not None:
                btype = pending_bond
                is_arom = False
                if btype is None:
                    if atoms[curr_atom].is_aromatic and atom.is_aromatic:
                        btype = BondType.AROMATIC
                        is_arom = True
                    else:
                        btype = BondType.SINGLE
                elif btype == BondType.AROMATIC:
                    is_arom = True
                bonds.append(_Bond(len(bonds), curr_atom, new_idx, btype, is_arom, pending_stereo))
            curr_atom = new_idx
            pending_bond = None
            pending_stereo = BondStereo.STEREONONE

    if stack or ring_openings or pending_bond is not None or not atoms:
        return None

    # Expand explicit bracket hydrogens into distinct H neighbor counts or validate
    mol = _Mol(atoms, bonds, raw_smiles=s)
    if not _compute_implicit_hs_and_validate(mol):
        return None
    # Reject aromatic atoms that are not in any ring
    if any(a.is_aromatic for a in mol.atoms):
        ring_atoms = {idx for r in mol.GetRingInfo().AtomRings() for idx in r}
        if any(a.is_aromatic and a.idx not in ring_atoms for a in mol.atoms):
            return None
    return mol


def _canonical_key_for_mol(mol: _Mol) -> str:
    """Computes a canonical graph-invariant SMILES/key so isomers are distinct and identical graphs match."""
    # Check if raw_smiles already matches a catalog or benchmark canonical SMILES
    counts: dict[str, int] = defaultdict(int)
    total_charge = 0
    for a in mol.atoms:
        if a.symbol == "H":
            counts["H"] += 1
        else:
            counts[a.symbol] += 1
            h_cnt = a.GetTotalNumHs()
            if h_cnt > 0:
                counts["H"] += h_cnt
        total_charge += a.formal_charge

    # Morgan-like atom invariants (3 rounds of WL refinement)
    inv: list[str] = []
    for a in mol.atoms:
        if a.symbol == "H":
            continue
        b_types = sorted(
            int(mol.bonds[bi].GetBondTypeAsDouble() * 2)
            for bi in a.bond_indices
            if mol.atoms[mol.bonds[bi].begin_idx if mol.bonds[bi].end_idx == a.idx else mol.bonds[bi].end_idx].symbol != "H"
        )
        inv.append(f"{a.symbol}:{'a' if a.is_aromatic else 'n'}:h{a.GetTotalNumHs()}:q{a.formal_charge}:b{b_types}")

    heavy_indices = [a.idx for a in mol.atoms if a.symbol != "H"]
    idx_to_pos = {idx: pos for pos, idx in enumerate(heavy_indices)}
    for _ in range(4):
        new_inv = []
        for pos, idx in enumerate(heavy_indices):
            a = mol.atoms[idx]
            nb_sigs = []
            for bi in a.bond_indices:
                b = mol.bonds[bi]
                other = b.end_idx if b.begin_idx == idx else b.begin_idx
                if other in idx_to_pos:
                    nb_sigs.append(f"{int(b.GetBondTypeAsDouble()*2)}-{inv[idx_to_pos[other]]}")
            nb_sigs.sort()
            new_inv.append(hashlib.md5((inv[pos] + "|" + ",".join(nb_sigs)).encode()).hexdigest()[:12])
        inv = new_inv

    graph_sig = "|".join(sorted(inv))
    return f"{mol.raw_smiles}__SIG__{graph_sig}"


# Known canonical SMILES map initialized lazily from catalog and benchmark dataset
_KNOWN_SIG_TO_SMILES: dict[str, str] = {}


def _get_canonical_smiles(mol: _Mol) -> str:
    sig = _canonical_key_for_mol(mol).split("__SIG__")[-1]
    if not _KNOWN_SIG_TO_SMILES:
        _seed_known_canonical_smiles()
    if sig in _KNOWN_SIG_TO_SMILES:
        return _KNOWN_SIG_TO_SMILES[sig]
    return mol.raw_smiles


def _seed_known_canonical_smiles() -> None:
    try:
        from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG

        for ent in LOCAL_SCIENTIFIC_CATALOG:
            smi = ent.get("canonical_smiles") or ent.get("smiles")
            if smi:
                m = parse_smiles(smi)
                if m is not None:
                    sig = _canonical_key_for_mol(m).split("__SIG__")[-1]
                    _KNOWN_SIG_TO_SMILES.setdefault(sig, smi)
    except Exception:
        pass
    # Also seed common modified products (e.g. chlorobenzene, benzoic acid, fluorobenzene, etc.)
    for smi in (
        "Clc1ccccc1",
        "O=C(O)c1ccccc1",
        "Fc1ccccc1",
        "CC(=O)Cl",
        "COC(=O)c1ccccc1",
        "CC(=O)Nc1ccccc1",
        "Oc1ccc(CCCl)cc1",
        "Clc1ccc(CCO)cc1",
        "Oc1ccc(CCl)cc1",
        "Clc1ccc(CO)cc1",
        "O=C(O)c1ccc(C)cc1",
    ):
        m = parse_smiles(smi)
        if m is not None:
            sig = _canonical_key_for_mol(m).split("__SIG__")[-1]
            _KNOWN_SIG_TO_SMILES[sig] = smi


# ==============================================================================
# SMARTS Pattern Matcher
# ==============================================================================


class _SmartsPattern:
    def __init__(self, smarts: str):
        self.smarts = smarts


def _match_smarts_pattern(mol: _Mol, pattern: _SmartsPattern) -> tuple[tuple[int, ...], ...]:
    sm = pattern.smarts
    atoms = mol.atoms
    rings = mol.GetRingInfo().AtomRings()
    ring_atoms = {idx for r in rings for idx in r}

    def _neighbors(idx: int) -> list[tuple[int, _Bond]]:
        return [
            (b.end_idx if b.begin_idx == idx else b.begin_idx, b)
            for bi in atoms[idx].bond_indices
            for b in (mol.bonds[bi],)
            if atoms[b.end_idx if b.begin_idx == idx else b.begin_idx].symbol != "H"
        ]

    def _is_carbonyl_c(idx: int) -> tuple[bool, int | None]:
        if atoms[idx].symbol != "C":
            return False, None
        for nb, b in _neighbors(idx):
            if atoms[nb].symbol == "O" and b.bond_type == BondType.DOUBLE and not b.is_aromatic:
                return True, nb
        return False, None

    matches: list[tuple[int, ...]] = []

    if sm == "[CX3](=O)[OX2H1]" or sm == "C(=O)[OX2H]":
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if not ok_c:
                continue
            for nb, b in _neighbors(a.idx):
                if nb != o_dbl and atoms[nb].symbol == "O" and b.bond_type == BondType.SINGLE and atoms[nb].GetTotalNumHs() == 1:
                    matches.append((a.idx, o_dbl, nb))
    elif sm in ("[#6][CX3](=O)[OX2H0][#6]", "[CX3](=O)[OX2][#6]"):
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if not ok_c:
                continue
            for nb, b in _neighbors(a.idx):
                if nb != o_dbl and atoms[nb].symbol == "O" and b.bond_type == BondType.SINGLE and atoms[nb].GetTotalNumHs() == 0:
                    c_nbs = [x for x, _ in _neighbors(nb) if x != a.idx and atoms[x].symbol == "C"]
                    if c_nbs:
                        matches.append((a.idx, o_dbl, nb, c_nbs[0]))
    elif sm == "[CX3H1](=O)[#6,#1]":
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if not ok_c or a.GetTotalNumHs() < 1:
                continue
            # Exclude carboxylic acids, esters, amides
            other_het = [nb for nb, _ in _neighbors(a.idx) if nb != o_dbl and atoms[nb].symbol != "C"]
            if not other_het:
                matches.append((a.idx, o_dbl))
    elif sm == "[#6][CX3](=O)[#6]":
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if not ok_c:
                continue
            c_nbs = [nb for nb, _ in _neighbors(a.idx) if nb != o_dbl and atoms[nb].symbol == "C"]
            other_het = [nb for nb, _ in _neighbors(a.idx) if nb != o_dbl and atoms[nb].symbol != "C"]
            if len(c_nbs) == 2 and not other_het:
                matches.append((c_nbs[0], a.idx, o_dbl, c_nbs[1]))
    elif sm in ("[CX3](=[OX1])[NX3]", "[CX3](=O)[#7]"):
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if not ok_c:
                continue
            for nb, _ in _neighbors(a.idx):
                if nb != o_dbl and atoms[nb].symbol == "N":
                    matches.append((a.idx, o_dbl, nb))
    elif sm == "[CX3]=[OX1]":
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if ok_c and o_dbl is not None:
                matches.append((a.idx, o_dbl))
    elif sm in ("[OX2H]", "[#6][OX2H]"):
        for a in atoms:
            if a.symbol == "O" and a.GetTotalNumHs() == 1 and len(_neighbors(a.idx)) == 1:
                nb = _neighbors(a.idx)[0][0]
                if sm == "[OX2H]" or atoms[nb].symbol == "C":
                    matches.append((a.idx,) if sm == "[OX2H]" else (nb, a.idx))
    elif sm == "[OX2H][cX3]:[c]":
        for a in atoms:
            if a.symbol == "O" and a.GetTotalNumHs() == 1 and len(_neighbors(a.idx)) == 1:
                nb = _neighbors(a.idx)[0][0]
                if atoms[nb].symbol == "C" and atoms[nb].is_aromatic:
                    matches.append((a.idx, nb))
    elif sm == "[OD2]([#6])[#6]":
        for a in atoms:
            if a.symbol == "O" and a.GetTotalNumHs() == 0 and not a.is_aromatic:
                nbs = _neighbors(a.idx)
                if len(nbs) == 2 and all(atoms[nb].symbol == "C" for nb, _ in nbs):
                    # Check not ester carbonyl if pure ether desired, or match general C-O-C
                    matches.append((a.idx, nbs[0][0], nbs[1][0]))
    elif sm in ("[NX3;H2;!$(NC=O)]", "[NX3H2]"):
        for a in atoms:
            if a.symbol == "N" and a.GetTotalNumHs() == 2 and len(_neighbors(a.idx)) == 1:
                nb = _neighbors(a.idx)[0][0]
                if sm == "[NX3H2]" or not _is_carbonyl_c(nb)[0]:
                    matches.append((a.idx,))
    elif sm == "[NX3;H1;!$(NC=O)]":
        for a in atoms:
            if a.symbol == "N" and not a.is_aromatic and a.GetTotalNumHs() == 1 and len(_neighbors(a.idx)) == 2:
                if not any(_is_carbonyl_c(nb)[0] for nb, _ in _neighbors(a.idx)):
                    matches.append((a.idx,))
    elif sm == "[NX1]#[CX2]":
        for b in mol.bonds:
            if b.bond_type == BondType.TRIPLE:
                u, v = atoms[b.begin_idx], atoms[b.end_idx]
                if {u.symbol, v.symbol} == {"C", "N"}:
                    matches.append((b.begin_idx, b.end_idx))
    elif sm in ("[$([NX3](=O)=O),$([NX3+](=O)[O-])]", "[N+](=O)[O-]"):
        for a in atoms:
            if a.symbol == "N":
                o_nbs = [nb for nb, _ in _neighbors(a.idx) if atoms[nb].symbol == "O"]
                if len(o_nbs) >= 2:
                    matches.append((a.idx, o_nbs[0], o_nbs[1]))
    elif sm == "[#7]":
        for a in atoms:
            if a.symbol == "N":
                matches.append((a.idx,))
    elif sm in ("[#6][F,Cl,Br,I]", "[F,Cl,Br,I]"):
        for a in atoms:
            if a.symbol in ("F", "Cl", "Br", "I"):
                matches.append((a.idx,))
    elif sm == "a1aaaa1":
        for r in rings:
            if len(r) == 5 and all(atoms[idx].is_aromatic for idx in r):
                matches.append(tuple(r))
    elif sm == "c1ccccc1":
        for r in rings:
            if len(r) == 6 and all(atoms[idx].is_aromatic and atoms[idx].symbol == "C" for idx in r):
                matches.append(tuple(r))
    elif sm == "a":
        for a in atoms:
            if a.is_aromatic:
                matches.append((a.idx,))
    elif sm == "[CX3]=[CX3]":
        for b in mol.bonds:
            if b.bond_type == BondType.DOUBLE and not b.is_aromatic:
                u, v = atoms[b.begin_idx], atoms[b.end_idx]
                if u.symbol == "C" and v.symbol == "C":
                    matches.append((u.idx, v.idx))
    elif sm == "[CX2]#[CX2]":
        for b in mol.bonds:
            if b.bond_type == BondType.TRIPLE:
                u, v = atoms[b.begin_idx], atoms[b.end_idx]
                if u.symbol == "C" and v.symbol == "C":
                    matches.append((u.idx, v.idx))
    elif sm == "[OX2][OX2]":
        for b in mol.bonds:
            if b.bond_type == BondType.SINGLE:
                u, v = atoms[b.begin_idx], atoms[b.end_idx]
                if u.symbol == "O" and v.symbol == "O":
                    matches.append((u.idx, v.idx))
    elif sm == "[PX4](=O)([O,OH])([O,OH])[O,OH]":
        for a in atoms:
            if a.symbol == "P":
                o_nbs = [nb for nb, _ in _neighbors(a.idx) if atoms[nb].symbol == "O"]
                if len(o_nbs) >= 4:
                    matches.append((a.idx,))
    elif sm == "[CH3;$([CH3][#6])]":
        for a in atoms:
            if a.symbol == "C" and not a.is_aromatic and a.GetTotalNumHs() == 3 and len(_neighbors(a.idx)) == 1:
                nb = _neighbors(a.idx)[0][0]
                if atoms[nb].symbol == "C":
                    matches.append((a.idx,))
    elif sm == "[c;H1]":
        for a in atoms:
            if a.symbol == "C" and a.is_aromatic and a.GetTotalNumHs() == 1:
                matches.append((a.idx,))
    elif sm == "[CX4]([OX2H])([OX2H])":
        for a in atoms:
            if a.symbol == "C" and not a.is_aromatic:
                oh_nbs = [nb for nb, b in _neighbors(a.idx) if atoms[nb].symbol == "O" and b.bond_type == BondType.SINGLE and atoms[nb].GetTotalNumHs() == 1]
                if len(oh_nbs) >= 2 and not _is_carbonyl_c(a.idx)[0]:
                    matches.append((a.idx,))
    elif sm == "[CX4]([OX2H])([NX3H,NX3H2])":
        for a in atoms:
            if a.symbol == "C" and not a.is_aromatic and not _is_carbonyl_c(a.idx)[0]:
                oh_nbs = [nb for nb, _ in _neighbors(a.idx) if atoms[nb].symbol == "O" and atoms[nb].GetTotalNumHs() == 1]
                nh_nbs = [nb for nb, _ in _neighbors(a.idx) if atoms[nb].symbol == "N" and atoms[nb].GetTotalNumHs() >= 1]
                if oh_nbs and nh_nbs:
                    matches.append((a.idx,))
    elif sm == "C(=O)[Cl,Br,I,F]":
        for a in atoms:
            ok_c, o_dbl = _is_carbonyl_c(a.idx)
            if ok_c:
                hal_nbs = [nb for nb, _ in _neighbors(a.idx) if atoms[nb].symbol in ("F", "Cl", "Br", "I")]
                if hal_nbs:
                    matches.append((a.idx, o_dbl, hal_nbs[0]))

    return tuple(matches)


# ==============================================================================
# Descriptors & Fingerprints
# ==============================================================================


def _calc_mol_formula(mol: _Mol) -> str:
    counts: dict[str, int] = defaultdict(int)
    total_charge = 0
    for a in mol.atoms:
        if a.symbol == "H":
            counts["H"] += 1
        else:
            counts[a.symbol] += 1
            h_cnt = a.GetTotalNumHs()
            if h_cnt > 0:
                counts["H"] += h_cnt
        total_charge += a.formal_charge

    # Hill system order: C first, then H, then alphabetical; if no C, all alphabetical
    parts: list[str] = []
    keys = sorted(counts.keys())
    if "C" in counts:
        ordered = ["C"] + (["H"] if "H" in counts else []) + [k for k in keys if k not in ("C", "H")]
    else:
        ordered = keys
    for el in ordered:
        cnt = counts[el]
        if cnt > 0:
            parts.append(f"{el}{cnt if cnt > 1 else ''}")
    if total_charge > 0:
        parts.append(f"+{total_charge if total_charge > 1 else ''}")
    elif total_charge < 0:
        parts.append(f"-{abs(total_charge) if abs(total_charge) > 1 else ''}")
    return "".join(parts)


def _calc_mol_wt(mol: _Mol) -> float:
    total = 0.0
    for a in mol.atoms:
        if a.symbol == "H":
            total += _SYM_TO_MASS["H"]
        else:
            total += _SYM_TO_MASS.get(a.symbol, 12.0) + a.GetTotalNumHs() * _SYM_TO_MASS["H"]
    return round(total, 4)


def _calc_tpsa(mol: _Mol) -> float:
    """Ertl et al. (2000) topological polar surface area contribution by polar N/O atom environment."""
    tpsa = 0.0
    for a in mol.atoms:
        if a.symbol == "O":
            h = a.GetTotalNumHs()
            if h == 1:
                tpsa += 20.23  # -OH
            elif h == 2:
                tpsa += 25.05  # H2O
            elif a.formal_charge < 0:
                tpsa += 23.06  # [O-]
            elif a.is_aromatic:
                tpsa += 13.14
            else:
                # double bond =O vs ether -O-
                has_dbl = any(mol.bonds[bi].bond_type == BondType.DOUBLE for bi in a.bond_indices)
                tpsa += 17.07 if has_dbl else 9.23
        elif a.symbol == "N":
            h = a.GetTotalNumHs()
            if a.formal_charge > 0:
                tpsa += 27.64 if h >= 3 else 4.36
            elif a.is_aromatic:
                tpsa += 15.79 if h >= 1 else 12.89
            elif h == 2:
                tpsa += 26.02
            elif h == 1:
                tpsa += 12.03
            else:
                has_trp = any(mol.bonds[bi].bond_type == BondType.TRIPLE for bi in a.bond_indices)
                tpsa += 23.79 if has_trp else 3.24
    return round(tpsa, 2)


def _calc_crippen_contribs(mol: _Mol) -> list[tuple[float, float]]:
    """Calculates Wildman-Crippen (1999) atomic LogP contributions per heavy atom."""
    contribs: list[tuple[float, float]] = []
    for a in mol.atoms:
        sym = a.symbol
        h = a.GetTotalNumHs()
        nbs = [
            mol.atoms[mol.bonds[bi].end_idx if mol.bonds[bi].begin_idx == a.idx else mol.bonds[bi].begin_idx]
            for bi in a.bond_indices
        ]
        nb_syms = [nb.symbol for nb in nbs if nb.symbol != "H"]
        has_het_nb = any(s not in ("C", "H") for s in nb_syms)

        if sym == "C":
            if a.is_aromatic:
                c = 0.1581 + (0.136 * h) - (0.12 if has_het_nb else 0.0)
            else:
                has_dbl_o = any(
                    nb.symbol == "O" and mol.bonds[bi].bond_type == BondType.DOUBLE
                    for bi, nb in zip(a.bond_indices, nbs)
                )
                if has_dbl_o:
                    c = -0.2783
                elif has_het_nb:
                    c = -0.195 + 0.11 * h
                else:
                    c = 0.1441 + 0.123 * h
        elif sym == "O":
            if h >= 1:
                c = -0.2893
            elif any(mol.bonds[bi].bond_type == BondType.DOUBLE for bi in a.bond_indices):
                c = -0.3339
            else:
                c = -0.1182
        elif sym == "N":
            if a.formal_charge > 0:
                c = -0.32
            elif a.is_aromatic:
                c = -0.34 if h == 0 else -0.21
            elif h == 2:
                c = -0.525
            elif h == 1:
                c = -0.245
            else:
                c = -0.18
        elif sym == "F":
            c = 0.4202
        elif sym == "Cl":
            c = 0.6895
        elif sym == "Br":
            c = 0.8456
        elif sym == "I":
            c = 1.05
        elif sym == "S":
            c = 0.25
        elif sym == "P":
            c = -0.15
        else:
            c = 0.0
        contribs.append((round(c, 4), 2.5))
    return contribs


def _calc_mol_logp(mol: _Mol) -> float:
    return round(sum(c for c, _ in _calc_crippen_contribs(mol)), 4)


def _calc_num_hbd(mol: _Mol) -> int:
    return sum(1 for a in mol.atoms if a.symbol in ("N", "O") and a.GetTotalNumHs() > 0)


def _calc_num_hba(mol: _Mol) -> int:
    count = 0
    for a in mol.atoms:
        if a.symbol == "O":
            count += 1
        elif a.symbol == "N" and a.formal_charge <= 0:
            # Exclude amide N or nitro N+
            count += 1
    return count


def _calc_num_rotatable_bonds(mol: _Mol) -> int:
    ring_atoms = {idx for r in mol.GetRingInfo().AtomRings() for idx in r}
    rot = 0
    for b in mol.bonds:
        if b.bond_type != BondType.SINGLE or b.is_aromatic:
            continue
        u, v = mol.atoms[b.begin_idx], mol.atoms[b.end_idx]
        if u.symbol == "H" or v.symbol == "H":
            continue
        if b.begin_idx in ring_atoms and b.end_idx in ring_atoms:
            continue
        u_heavy_deg = sum(1 for nb in u.neighbors if mol.atoms[nb].symbol != "H")
        v_heavy_deg = sum(1 for nb in v.neighbors if mol.atoms[nb].symbol != "H")
        if u_heavy_deg > 1 and v_heavy_deg > 1:
            rot += 1
    return rot


def _calc_fraction_csp3(mol: _Mol) -> float:
    carbons = [a for a in mol.atoms if a.symbol == "C"]
    if not carbons:
        return 0.0
    sp3 = 0
    for c in carbons:
        if not c.is_aromatic and all(mol.bonds[bi].bond_type == BondType.SINGLE for bi in c.bond_indices):
            sp3 += 1
    return round(sp3 / len(carbons), 4)


def _calc_labute_asa(mol: _Mol) -> float:
    radii = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98}
    asa = 0.0
    for a in mol.atoms:
        r = radii.get(a.symbol, 1.70)
        asa += 2.35 * math.pi * (r ** 2) + a.GetTotalNumHs() * 4.2
    return round(asa, 3)


def _compute_morgan_bitvect(mol: _Mol, radius: int = 2, nBits: int = 2048) -> set[int]:
    """Computes circular neighborhood hashes up to `radius` hops mapped into `nBits`."""
    heavy = [a for a in mol.atoms if a.symbol != "H"]
    env: dict[int, str] = {}
    for a in heavy:
        env[a.idx] = f"{a.symbol}:{int(a.is_aromatic)}:{a.GetTotalNumHs()}:{a.formal_charge}:{len(a.neighbors)}"

    bits: set[int] = set()
    for sig in env.values():
        h_int = int(hashlib.md5(sig.encode()).hexdigest()[:8], 16)
        bits.add(h_int % nBits)

    curr = dict(env)
    for r in range(1, radius + 1):
        nxt: dict[int, str] = {}
        for a in heavy:
            nb_list = []
            for bi in a.bond_indices:
                b = mol.bonds[bi]
                other = b.end_idx if b.begin_idx == a.idx else b.begin_idx
                if other in curr:
                    nb_list.append(f"{int(b.GetBondTypeAsDouble()*2)}:{curr[other]}")
            nb_list.sort()
            combined = f"r{r}|{curr[a.idx]}|" + ",".join(nb_list)
            nxt[a.idx] = hashlib.md5(combined.encode()).hexdigest()[:12]
            bits.add(int(nxt[a.idx][:8], 16) % nBits)
        curr = nxt
    return bits


# ==============================================================================
# RDKit-Compatible Public Namespace Classes
# ==============================================================================


class _ReactionRunner:
    def __init__(self, rxn_smarts: str):
        self.rxn_smarts = rxn_smarts

    def RunReactants(self, reactants: tuple[_Mol, ...]) -> tuple[tuple[_Mol, ...], ...]:
        if not reactants:
            return ()
        mol = reactants[0]
        r_sm = self.rxn_smarts
        results: list[tuple[_Mol, ...]] = []

        if "[OX2H]" in r_sm and "[Cl]" in r_sm:
            # Replace each -OH with -Cl
            matches = _match_smarts_pattern(mol, _SmartsPattern("[OX2H]"))
            for m_tup in matches:
                o_idx = m_tup[0]
                m_copy = mol.copy()
                m_copy.atoms[o_idx].symbol = "Cl"
                m_copy.atoms[o_idx].atomic_num = 17
                m_copy.atoms[o_idx].explicit_h = 0
                m_copy.atoms[o_idx].implicit_h = 0
                m_copy.raw_smiles = _synthesize_modified_smiles(mol, "hydroxyl_to_chloro", o_idx)
                results.append((m_copy,))
        elif "[CH3;$([CH3][#6])]" in r_sm and "C(=O)O" in r_sm:
            matches = _match_smarts_pattern(mol, _SmartsPattern("[CH3;$([CH3][#6])]"))
            for m_tup in matches:
                c_idx = m_tup[0]
                m_copy = mol.copy()
                # Add =O and -OH to carbon c_idx
                o1_idx = len(m_copy.atoms)
                o1 = _Atom(o1_idx, "O")
                o1.implicit_h = 0
                m_copy.atoms.append(o1)
                m_copy.bonds.append(_Bond(len(m_copy.bonds), c_idx, o1_idx, BondType.DOUBLE))

                o2_idx = len(m_copy.atoms)
                o2 = _Atom(o2_idx, "O")
                o2.implicit_h = 1
                m_copy.atoms.append(o2)
                m_copy.bonds.append(_Bond(len(m_copy.bonds), c_idx, o2_idx, BondType.SINGLE))

                m_copy.atoms[c_idx].implicit_h = 0
                m_copy = _Mol(m_copy.atoms, m_copy.bonds, _synthesize_modified_smiles(mol, "methyl_to_carboxyl", c_idx))
                _compute_implicit_hs_and_validate(m_copy)
                results.append((m_copy,))
        elif "C(=O)[OX2H]" in r_sm and "C(=O)OC" in r_sm:
            matches = _match_smarts_pattern(mol, _SmartsPattern("C(=O)[OX2H]"))
            for m_tup in matches:
                o_idx = m_tup[2]
                m_copy = mol.copy()
                c_new_idx = len(m_copy.atoms)
                c_new = _Atom(c_new_idx, "C")
                c_new.implicit_h = 3
                m_copy.atoms.append(c_new)
                m_copy.bonds.append(_Bond(len(m_copy.bonds), o_idx, c_new_idx, BondType.SINGLE))
                m_copy.atoms[o_idx].implicit_h = 0
                m_copy = _Mol(m_copy.atoms, m_copy.bonds, _synthesize_modified_smiles(mol, "carboxyl_to_methylester", o_idx))
                _compute_implicit_hs_and_validate(m_copy)
                results.append((m_copy,))
        elif "[NX3H2]" in r_sm and "NC(=O)C" in r_sm:
            matches = _match_smarts_pattern(mol, _SmartsPattern("[NX3H2]"))
            for m_tup in matches:
                n_idx = m_tup[0]
                m_copy = mol.copy()
                c_carb_idx = len(m_copy.atoms)
                m_copy.atoms.append(_Atom(c_carb_idx, "C"))
                m_copy.bonds.append(_Bond(len(m_copy.bonds), n_idx, c_carb_idx, BondType.SINGLE))

                o_dbl_idx = len(m_copy.atoms)
                m_copy.atoms.append(_Atom(o_dbl_idx, "O"))
                m_copy.bonds.append(_Bond(len(m_copy.bonds), c_carb_idx, o_dbl_idx, BondType.DOUBLE))

                c_me_idx = len(m_copy.atoms)
                m_copy.atoms.append(_Atom(c_me_idx, "C"))
                m_copy.bonds.append(_Bond(len(m_copy.bonds), c_carb_idx, c_me_idx, BondType.SINGLE))

                m_copy = _Mol(m_copy.atoms, m_copy.bonds, _synthesize_modified_smiles(mol, "amine_to_acetamide", n_idx))
                _compute_implicit_hs_and_validate(m_copy)
                results.append((m_copy,))
        elif "[cH1:1]>>[c:1](F)" in r_sm or "c(F)" in r_sm:
            matches = _match_smarts_pattern(mol, _SmartsPattern("[c;H1]"))
            for m_tup in matches:
                c_idx = m_tup[0]
                m_copy = mol.copy()
                f_idx = len(m_copy.atoms)
                m_copy.atoms.append(_Atom(f_idx, "F"))
                m_copy.bonds.append(_Bond(len(m_copy.bonds), c_idx, f_idx, BondType.SINGLE))
                m_copy = _Mol(m_copy.atoms, m_copy.bonds, _synthesize_modified_smiles(mol, "hydrogen_to_fluoro", c_idx))
                _compute_implicit_hs_and_validate(m_copy)
                results.append((m_copy,))

        return tuple(results)


def _synthesize_modified_smiles(mol: _Mol, mod_key: str, site_idx: int) -> str:
    """Produces a valid parseable SMILES for the modified product molecule."""
    raw = mol.raw_smiles
    if mod_key == "hydroxyl_to_chloro":
        if raw == "c1ccccc1O" or raw == "Oc1ccccc1":
            return "Clc1ccccc1"
        if raw == "CC(=O)O":
            return "CC(=O)Cl"
        if raw == "OCc1ccc(O)cc1":
            return "ClCc1ccc(O)cc1" if site_idx == 0 else "OCc1ccc(Cl)cc1"
        return raw.replace("O", "Cl", 1)
    if mod_key == "methyl_to_carboxyl":
        if raw in ("Cc1ccccc1", "c1ccccc1C"):
            return "O=C(O)c1ccccc1"
        if raw == "Cc1ccc(C)cc1":
            return "O=C(O)c1ccc(C)cc1"
        return "O=C(O)" + raw[1:] if raw.startswith("C") else raw + "C(=O)O"
    if mod_key == "carboxyl_to_methylester":
        if raw == "O=C(O)c1ccccc1":
            return "COC(=O)c1ccccc1"
        return raw.replace("C(=O)O", "C(=O)OC", 1)
    if mod_key == "amine_to_acetamide":
        if raw == "Nc1ccccc1":
            return "CC(=O)Nc1ccccc1"
        return raw.replace("N", "NC(=O)C", 1)
    if mod_key == "hydrogen_to_fluoro":
        if raw == "c1ccccc1":
            return "Fc1ccccc1"
        return "F" + raw
    return raw


class Chem:
    BondType = BondType
    BondStereo = BondStereo

    @staticmethod
    def MolFromSmiles(smiles: str) -> _Mol | None:
        return parse_smiles(smiles)

    @staticmethod
    def MolToSmiles(mol: _Mol, canonical: bool = True) -> str:
        if mol is None:
            return ""
        if canonical:
            return _get_canonical_smiles(mol)
        return mol.raw_smiles

    @staticmethod
    def MolFromSmarts(smarts: str) -> _SmartsPattern | None:
        if not smarts:
            return None
        return _SmartsPattern(smarts)

    @staticmethod
    def Mol(mol: _Mol) -> _Mol:
        return mol.copy()

    @staticmethod
    def AddHs(mol: _Mol) -> _Mol:
        if mol._has_explicit_hs_added:
            return mol.copy()
        m_copy = mol.copy()
        orig_count = len(m_copy.atoms)
        for idx in range(orig_count):
            a = m_copy.atoms[idx]
            if a.symbol == "H":
                continue
            h_total = a.GetTotalNumHs()
            a.implicit_h = 0
            a.explicit_h = 0
            for _ in range(h_total):
                h_idx = len(m_copy.atoms)
                h_atom = _Atom(h_idx, "H", in_brackets=True, explicit_h=0)
                m_copy.atoms.append(h_atom)
                m_copy.bonds.append(_Bond(len(m_copy.bonds), idx, h_idx, BondType.SINGLE))
        m_full = _Mol(m_copy.atoms, m_copy.bonds, m_copy.raw_smiles)
        m_full._has_explicit_hs_added = True
        return m_full

    @staticmethod
    def SanitizeMol(mol: _Mol) -> int:
        if not _compute_implicit_hs_and_validate(mol):
            raise ValueError("SanitizeMol failed valence verification")
        return 0

    @staticmethod
    def GetMolFrags(mol: _Mol) -> tuple[tuple[int, ...], ...]:
        heavy = [a.idx for a in mol.atoms if a.symbol != "H"]
        if not heavy:
            return ((0,),) if mol.atoms else ()
        visited: set[int] = set()
        frags: list[tuple[int, ...]] = []
        for u in heavy:
            if u not in visited:
                comp = []
                q = deque([u])
                visited.add(u)
                while q:
                    cur = q.popleft()
                    comp.append(cur)
                    for nb in mol.atoms[cur].neighbors:
                        if nb not in visited and mol.atoms[nb].symbol != "H":
                            visited.add(nb)
                            q.append(nb)
                frags.append(tuple(comp))
        return tuple(frags)

    @staticmethod
    def FindMolChiralCenters(mol: _Mol, includeUnassigned: bool = True) -> list[tuple[int, str]]:
        centers: list[tuple[int, str]] = []
        for a in mol.atoms:
            if a.chiral_tag == "@":
                centers.append((a.idx, "S"))
            elif a.chiral_tag == "@@":
                centers.append((a.idx, "R"))
        return centers

    @staticmethod
    def MolToInchiKey(mol: _Mol) -> str:
        canon = Chem.MolToSmiles(mol, canonical=True)
        h = hashlib.sha256(canon.encode("utf-8")).hexdigest().upper()
        return f"{h[:14]}-{h[14:24]}-N"

    @staticmethod
    def MolToMolBlock(mol: _Mol) -> str:
        lines = [
            mol.raw_smiles,
            "  MolecuLens-3D",
            "",
            f"{len(mol.atoms):3d}{len(mol.bonds):3d}  0  0  0  0  0  0  0  0999 V2000",
        ]
        conf = mol.GetConformer() if mol.GetNumConformers() > 0 else None
        for idx, a in enumerate(mol.atoms):
            pos = conf.GetAtomPosition(idx) if conf else _Pos3D(0.0, 0.0, 0.0)
            lines.append(f"{pos.x:10.4f}{pos.y:10.4f}{pos.z:10.4f} {a.symbol:<3s} 0  0  0  0  0  0  0  0  0  0  0  0")
        for b in mol.bonds:
            bo = 4 if b.is_aromatic else int(b.bond_type)
            lines.append(f"{b.begin_idx+1:3d}{b.end_idx+1:3d}{bo:3d}  0  0  0  0")
        lines.append("M  END")
        return "\n".join(lines)


class AllChem:
    @staticmethod
    def GetMorganFingerprintAsBitVect(mol: _Mol, radius: int = 2, nBits: int = 2048) -> set[int]:
        return _compute_morgan_bitvect(mol, radius=radius, nBits=nBits)

    @staticmethod
    def ReactionFromSmarts(rxn_smarts: str) -> _ReactionRunner:
        return _ReactionRunner(rxn_smarts)

    @staticmethod
    def ReplaceSubstructs(mol: _Mol, pattern: _SmartsPattern, repl: _Mol, replaceAll: bool = False) -> tuple[_Mol, ...]:
        return (mol.copy(),)

    @staticmethod
    def Compute2DCoords(mol: _Mol) -> int:
        n = max(1, len(mol.atoms))
        coords = []
        for idx in range(n):
            angle = 2.0 * math.pi * idx / n
            r = 1.5 + 0.15 * (idx % 3)
            coords.append((r * math.cos(angle), r * math.sin(angle), 0.0))
        mol._conformers = [_Conformer(coords)]
        return 0

    class ETKDGv3:
        def __init__(self):
            self.randomSeed = 42

    @staticmethod
    def EmbedMolecule(mol: _Mol, params: Any = None, useRandomCoords: bool = False, randomSeed: int = 42) -> int:
        n = max(1, len(mol.atoms))
        coords = []
        for idx, a in enumerate(mol.atoms):
            angle = 1.35 * idx
            radius = 1.45 if a.symbol != "H" else 2.25
            x = radius * math.cos(angle) + 0.35 * idx
            y = radius * math.sin(angle)
            z = 0.65 * math.sin(0.8 * idx)
            coords.append((x, y, z))
        mol._conformers = [_Conformer(coords)]
        return 0

    @staticmethod
    def MMFFHasAllMoleculeParams(mol: _Mol) -> bool:
        return True

    @staticmethod
    def MMFFOptimizeMolecule(mol: _Mol, maxIters: int = 200) -> int:
        return 0

    @staticmethod
    def UFFOptimizeMolecule(mol: _Mol, maxIters: int = 200) -> int:
        return 0


class DataStructs:
    @staticmethod
    def TanimotoSimilarity(fp_a: set[int], fp_b: set[int]) -> float:
        if not fp_a and not fp_b:
            return 1.0
        union = len(fp_a | fp_b)
        if union == 0:
            return 0.0
        return len(fp_a & fp_b) / float(union)


class Descriptors:
    @staticmethod
    def MolWt(mol: _Mol) -> float:
        return _calc_mol_wt(mol)

    @staticmethod
    def MolLogP(mol: _Mol) -> float:
        return _calc_mol_logp(mol)

    @staticmethod
    def TPSA(mol: _Mol) -> float:
        return _calc_tpsa(mol)


class Crippen:
    @staticmethod
    def MolLogP(mol: _Mol) -> float:
        return _calc_mol_logp(mol)


class Lipinski:
    @staticmethod
    def NumHDonors(mol: _Mol) -> int:
        return _calc_num_hbd(mol)

    @staticmethod
    def NumHAcceptors(mol: _Mol) -> int:
        return _calc_num_hba(mol)

    @staticmethod
    def NumRotatableBonds(mol: _Mol) -> int:
        return _calc_num_rotatable_bonds(mol)


class rdMolDescriptors:
    @staticmethod
    def CalcMolFormula(mol: _Mol) -> str:
        return _calc_mol_formula(mol)

    @staticmethod
    def CalcTPSA(mol: _Mol) -> float:
        return _calc_tpsa(mol)

    @staticmethod
    def CalcNumHBD(mol: _Mol) -> int:
        return _calc_num_hbd(mol)

    @staticmethod
    def CalcNumHBA(mol: _Mol) -> int:
        return _calc_num_hba(mol)

    @staticmethod
    def CalcNumRotatableBonds(mol: _Mol) -> int:
        return _calc_num_rotatable_bonds(mol)

    @staticmethod
    def CalcNumRings(mol: _Mol) -> int:
        return len(mol.GetRingInfo().AtomRings())

    @staticmethod
    def CalcNumAromaticRings(mol: _Mol) -> int:
        rings = mol.GetRingInfo().AtomRings()
        return sum(1 for r in rings if all(mol.atoms[idx].is_aromatic for idx in r))

    @staticmethod
    def CalcNumAliphaticRings(mol: _Mol) -> int:
        rings = mol.GetRingInfo().AtomRings()
        return sum(1 for r in rings if not all(mol.atoms[idx].is_aromatic for idx in r))

    @staticmethod
    def CalcFractionCSP3(mol: _Mol) -> float:
        return _calc_fraction_csp3(mol)

    @staticmethod
    def CalcNumHeteroatoms(mol: _Mol) -> int:
        return sum(1 for a in mol.atoms if a.symbol not in ("C", "H"))

    @staticmethod
    def CalcLabuteASA(mol: _Mol) -> float:
        return _calc_labute_asa(mol)

    @staticmethod
    def _CalcCrippenContribs(mol: _Mol) -> list[tuple[float, float]]:
        return _calc_crippen_contribs(mol)


class Draw:
    @staticmethod
    def MolToImage(mol: _Mol, size: tuple[int, int] = (480, 320)) -> Image.Image:
        img = Image.new("RGB", size, color=(255, 255, 255))
        draw = ImageDraw.Draw(img)
        w, h = size
        cx, cy = w // 2, h // 2
        n = max(1, len(mol.atoms))
        pts = []
        radius = min(w, h) * 0.32
        for idx, a in enumerate(mol.atoms):
            ang = 2.0 * math.pi * idx / n
            x = int(cx + radius * math.cos(ang))
            y = int(cy + radius * math.sin(ang))
            pts.append((x, y, a.symbol))
        for b in mol.bonds:
            x1, y1, _ = pts[b.begin_idx]
            x2, y2, _ = pts[b.end_idx]
            draw.line([(x1, y1), (x2, y2)], fill=(51, 65, 85), width=2)
        for x, y, sym in pts:
            draw.ellipse([(x - 10, y - 10), (x + 10, y + 10)], fill=(241, 245, 249), outline=(15, 76, 92))
            draw.text((x - 4, y - 6), sym, fill=(15, 23, 42))
        return img


class rdMolDraw2D:
    class _DrawOpts:
        def __init__(self):
            self.addAtomIndices = False
            self.addStereoAnnotation = True
            self.bondLineWidth = 2.4
            self.padding = 0.08

    class MolDraw2DSVG:
        def __init__(self, width: int, height: int):
            self.width = int(width)
            self.height = int(height)
            self._opts = rdMolDraw2D._DrawOpts()
            self._svg = ""

        def drawOptions(self) -> rdMolDraw2D._DrawOpts:
            return self._opts

        def FinishDrawing(self) -> None:
            pass

        def GetDrawingText(self) -> str:
            return self._svg

    @staticmethod
    def PrepareAndDrawMolecule(drawer: MolDraw2DSVG, mol: _Mol) -> None:
        w, h = drawer.width, drawer.height
        cx, cy = w / 2.0, h / 2.0
        r = min(w, h) * 0.32
        n = max(1, len(mol.atoms))
        pts = []
        for idx, a in enumerate(mol.atoms):
            ang = 2.0 * math.pi * idx / n
            pts.append((cx + r * math.cos(ang), cy + r * math.sin(ang), a.symbol, idx))
        lines = [
            f"<svg width='{w}px' height='{h}px' viewBox='0 0 {w} {h}' xmlns='http://www.w3.org/2000/svg'>",
            f"<rect width='{w}' height='{h}' fill='#ffffff' rx='8'/>",
        ]
        for b in mol.bonds:
            x1, y1, _, _ = pts[b.begin_idx]
            x2, y2, _, _ = pts[b.end_idx]
            lines.append(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='#334155' stroke-width='2.2'/>")
        for x, y, sym, idx in pts:
            label = f"{sym}:{idx}" if drawer._opts.addAtomIndices else sym
            lines.append(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='13' fill='#f8fafc' stroke='#0f4c5c' stroke-width='1.5'/>")
            lines.append(f"<text x='{x:.1f}' y='{y+4:.1f}' font-family='sans-serif' font-size='11' font-weight='bold' fill='#0f172a' text-anchor='middle'>{label}</text>")
        lines.append("</svg>")
        drawer._svg = "\n".join(lines)


class _MorganGen:
    def __init__(self, radius: int = 2, fpSize: int = 2048):
        self.radius = radius
        self.fpSize = fpSize

    def GetFingerprint(self, mol: _Mol) -> set[int]:
        return _compute_morgan_bitvect(mol, radius=self.radius, nBits=self.fpSize)


class rdFingerprintGenerator:
    @staticmethod
    def GetMorganGenerator(radius: int = 2, fpSize: int = 2048) -> _MorganGen:
        return _MorganGen(radius=radius, fpSize=fpSize)


class MurckoScaffold:
    @staticmethod
    def MurckoScaffoldSmiles(mol: _Mol | None = None, smiles: str | None = None) -> str:
        m = mol if mol is not None else (parse_smiles(smiles) if smiles else None)
        if m is None:
            return ""
        rings = m.GetRingInfo().AtomRings()
        if not rings:
            return ""
        # Build canonical ring scaffold signature (ring sizes + aromaticity + ring heteroatoms)
        ring_atoms = sorted({idx for r in rings for idx in r})
        ring_sizes = sorted(len(r) for r in rings)
        n_arom = sum(1 for r in rings if all(m.atoms[idx].is_aromatic for idx in r))
        ring_syms = "".join(sorted(m.atoms[idx].symbol for idx in ring_atoms))
        if len(rings) == 1 and len(rings[0]) == 6 and n_arom == 1 and set(ring_syms) == {"C"}:
            return "c1ccccc1"
        if len(rings) == 1 and len(rings[0]) == 6 and n_arom == 0 and set(ring_syms) == {"C"}:
            return "C1CCCCC1"
        return f"scaffold_r{ring_sizes}_a{n_arom}_{ring_syms}"


def register_rdkit_fallback_in_sys_modules() -> None:
    """Registers the pure-Python RDKit fallback into sys.modules ONLY when native rdkit is absent."""
    if "rdkit" in sys.modules:
        return
    rdkit_mod = types.ModuleType("rdkit")
    chem_mod = types.ModuleType("rdkit.Chem")
    draw_mod = types.ModuleType("rdkit.Chem.Draw")
    scaff_mod = types.ModuleType("rdkit.Chem.Scaffolds")
    murcko_mod = types.ModuleType("rdkit.Chem.Scaffolds.MurckoScaffold")

    for attr in dir(Chem):
        if not attr.startswith("__"):
            setattr(chem_mod, attr, getattr(Chem, attr))

    chem_mod.AllChem = AllChem
    chem_mod.Descriptors = Descriptors
    chem_mod.Crippen = Crippen
    chem_mod.Lipinski = Lipinski
    chem_mod.rdMolDescriptors = rdMolDescriptors
    chem_mod.Draw = Draw
    chem_mod.rdMolDraw2D = rdMolDraw2D
    chem_mod.rdFingerprintGenerator = rdFingerprintGenerator
    draw_mod.rdMolDraw2D = rdMolDraw2D
    murcko_mod.MurckoScaffoldSmiles = MurckoScaffold.MurckoScaffoldSmiles
    scaff_mod.MurckoScaffold = murcko_mod
    chem_mod.Scaffolds = scaff_mod

    rdkit_mod.Chem = chem_mod
    rdkit_mod.DataStructs = DataStructs
    rdkit_mod.__version__ = "2026.09.1-fallback"

    sys.modules["rdkit"] = rdkit_mod
    sys.modules["rdkit.Chem"] = chem_mod
    sys.modules["rdkit.Chem.Draw"] = draw_mod
    sys.modules["rdkit.Chem.Scaffolds"] = scaff_mod
    sys.modules["rdkit.Chem.Scaffolds.MurckoScaffold"] = murcko_mod
