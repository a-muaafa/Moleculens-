from typing import Any
from config.settings import OP_SCIENTIFIC_RETRIEVAL
from core.schemas import build_provenance_tag


def build_biology_section(entity: dict[str, Any]) -> dict[str, Any]:
    bio = dict(entity.get("biology_details", {}))
    ids = entity.get("identifiers", {})
    uniprot_id = bio.get("uniprot_id") or ids.get("uniprot_id", "غير متوفر")
    pdb_ids = ids.get("pdb_ids", [])
    rhea_ids = ids.get("rhea_ids", [])

    pdb_links = [
        {"pdb_id": pid, "url": f"https://www.rcsb.org/structure/{pid}", "viewer_3d_url": f"https://www.rcsb.org/3d-view/{pid}"}
        for pid in pdb_ids
    ]

    return {
        "uniprot_id": uniprot_id,
        "uniprot_url": f"https://www.uniprot.org/uniprotkb/{uniprot_id}/entry" if uniprot_id and uniprot_id != "غير متوفر" else None,
        "gene_name": bio.get("gene_name", "غير محدد"),
        "organism": bio.get("organism", "Homo sapiens"),
        "sequence_length": bio.get("sequence_length", "غير محدد"),
        "subcellular_location": bio.get("subcellular_location", "خلوي / خارج خلوي"),
        "biological_function_ar": bio.get("biological_function_ar", "يؤدي وظيفة حيوية متخصصة في الكائن الحي."),
        "pdb_structures": bio.get("pdb_structures", pdb_ids),
        "pdb_links": pdb_links,
        "rhea_ids": rhea_ids,
        "provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=["UniProt", "PDB", "Rhea"],
            notes_ar="بيانات التسلسل البروتيني والوظيفة والتراكيب ثلاثية الأبعاد مسترجعة من UniProt و RCSB PDB و Rhea.",
        ),
    }

