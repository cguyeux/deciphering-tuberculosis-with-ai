"""P6.2 -- ESM Atlas : similarite fonctionnelle par embedding (SAE features).

Modalite orthogonale a HHpred (P3.1, profil-profil) et Foldseek (P3.2, structure 3D) :
similarite dans l'espace des representations apprises par ESMC. Controle positif
obligatoire (ML2204, orthologue M. leprae du meme domaine DUF3073) avant d'interpreter
toute absence de signal sur Rv0810c comme negative.
"""
import json
import math
import sys
from pathlib import Path

SRC = "/home/christophe/docs/codes/claude_plugins/bio_population_genetics/skills/esm-atlas-cli/src"
sys.path.insert(0, SRC)
from esm_atlas_cli import EsmAtlasClient  # noqa: E402

RV0810C = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSGTGTDRLDGDGPSDDDSWNDEDDWRR"
ML2204 = "MGRGRAKAKQTKVARELKYSSPQTDFQRLQRELSSTGAADPGQLDGDDRVSEDSWDEDAWRR"
# Temoin negatif : RpmG2/Rv0634B (50S ribosomal L33, 55 aa), deja calibre en P3.4 comme
# temoin de taille comparable sans rapport biologique avec Rv0810c.
RPMG2 = "MASSTDVRPKITLACEVCKHRNYITKKNRRNDPDRLELKKFCPNCGKHQAHRETR"

OUT = Path(__file__).resolve().parent.parent / "résultats" / "p6_2_esm_atlas.json"


def sparse_cosine(pa, pb):
    """Cosine similarity entre deux vecteurs SAE creux au format {indices:[[...]], values:[...]}."""
    ia = pa["indices"][0]
    va = pa["values"]
    ib = pb["indices"][0]
    vb = pb["values"]
    da = dict(zip(ia, va))
    db = dict(zip(ib, vb))
    common = set(da) & set(db)
    dot = sum(da[i] * db[i] for i in common)
    norm_a = math.sqrt(sum(v * v for v in va))
    norm_b = math.sqrt(sum(v * v for v in vb))
    if norm_a == 0 or norm_b == 0:
        return 0.0
    return dot / (norm_a * norm_b)


def top_features(data, k=10):
    feats = data.get("sae_features") or []
    out = []
    for f in feats[:k]:
        out.append({
            "feature_index": f.get("feature_index"),
            "label": f.get("label"),
            "description": (f.get("description") or "")[:300],
        })
    return out


def main():
    client = EsmAtlasClient()
    result = {}

    # --- Controle positif : ML2204 (M. leprae, meme domaine DUF3073) ---
    ml2204 = client.lookup_sequence(ML2204)
    result["control_ML2204"] = {
        "protein_hash": ml2204.get("protein_hash"),
        "source": ml2204.get("source"),
        "accession": ml2204.get("accession"),
        "sequence_length": ml2204.get("sequence_length"),
        "top_features": top_features(ml2204),
    }

    # --- Similarity search a partir de la sequence Rv0810c : doit retrouver ML2204 ---
    try:
        sim_from_rv0810c = client.similarity_search(RV0810C, topk_results=25)
    except Exception as e:
        sim_from_rv0810c = {"error": f"{type(e).__name__}: {e}"}
    result["similarity_search_Rv0810c"] = sim_from_rv0810c

    try:
        sim_from_ml2204 = client.similarity_search(ML2204, topk_results=25)
    except Exception as e:
        sim_from_ml2204 = {"error": f"{type(e).__name__}: {e}"}
    result["similarity_search_ML2204_control"] = sim_from_ml2204

    # --- Temoin negatif : RpmG2 (55 aa, taille comparable, sans rapport biologique) ---
    rpmg2 = client.lookup_sequence(RPMG2)
    result["control_RpmG2_negative"] = {
        "protein_hash": rpmg2.get("protein_hash"),
        "source": rpmg2.get("source"),
        "sequence_length": rpmg2.get("sequence_length"),
        "top_features": top_features(rpmg2),
    }

    # --- Lookup Rv0810c lui-meme ---
    rv0810c = client.lookup_sequence(RV0810C)
    result["Rv0810c"] = {
        "protein_hash": rv0810c.get("protein_hash"),
        "source": rv0810c.get("source"),
        "accession": rv0810c.get("accession"),
        "sequence_length": rv0810c.get("sequence_length"),
        "cluster_rep_protein_hash": rv0810c.get("cluster_rep_protein_hash"),
        "top_features": top_features(rv0810c, k=20),
    }

    # --- Cluster info si disponible ---
    rep = rv0810c.get("cluster_rep_protein_hash") or rv0810c.get("protein_hash")
    try:
        cl = client.cluster(rep)
        result["cluster_Rv0810c"] = cl
    except Exception as e:
        result["cluster_Rv0810c"] = {"error": f"{type(e).__name__}: {e}"}

    # --- Meme chose pour le hash propre de Rv0810c si different du representant ---
    try:
        cl_self = client.cluster(rv0810c.get("protein_hash"))
        result["cluster_Rv0810c_self_hash"] = cl_self
    except Exception as e:
        result["cluster_Rv0810c_self_hash"] = {"error": f"{type(e).__name__}: {e}"}

    # --- Similarite cosinus directe sur le vecteur SAE pooled (protein_activations) ---
    cos_self = sparse_cosine(rv0810c["protein_activations"], rv0810c["protein_activations"])
    cos_ml2204 = sparse_cosine(rv0810c["protein_activations"], ml2204["protein_activations"])
    cos_rpmg2 = sparse_cosine(rv0810c["protein_activations"], rpmg2["protein_activations"])
    result["cosine_similarity"] = {
        "Rv0810c_vs_self": cos_self,
        "Rv0810c_vs_ML2204_positive_control": cos_ml2204,
        "Rv0810c_vs_RpmG2_negative_control": cos_rpmg2,
    }

    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"Ecrit : {OUT}")

    print("\n=== SIMILARITE COSINUS (vecteur SAE pooled, protein_activations) ===")
    print(f"Rv0810c vs lui-meme (sanite)      : {cos_self:.4f}")
    print(f"Rv0810c vs ML2204 (temoin positif) : {cos_ml2204:.4f}")
    print(f"Rv0810c vs RpmG2  (temoin negatif) : {cos_rpmg2:.4f}")

    # Resume a l'ecran
    print("\n=== CONTROLE POSITIF (ML2204) ===")
    print("hash:", result["control_ML2204"]["protein_hash"], "source:", result["control_ML2204"]["source"])
    for f in result["control_ML2204"]["top_features"][:5]:
        print(" -", f["feature_index"], f["label"])

    print("\n=== Rv0810c ===")
    print("hash:", result["Rv0810c"]["protein_hash"], "source:", result["Rv0810c"]["source"])
    for f in result["Rv0810c"]["top_features"][:10]:
        print(" -", f["feature_index"], f["label"])

    print("\n=== SIMILARITY SEARCH depuis Rv0810c (topk) ===")
    print(json.dumps(sim_from_rv0810c, indent=2, ensure_ascii=False)[:3000])

    print("\n=== SIMILARITY SEARCH depuis ML2204 (controle) ===")
    print(json.dumps(sim_from_ml2204, indent=2, ensure_ascii=False)[:3000])


if __name__ == "__main__":
    main()
