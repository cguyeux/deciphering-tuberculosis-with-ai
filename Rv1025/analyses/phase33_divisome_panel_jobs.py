#!/usr/bin/env python3
"""phase33_divisome_panel_jobs.py — P4.3 : panel divisome élargi, jobs Boltz-2.

P4.1 a exclu un complexe binaire direct Rv1025-DivIC (AF3, ipTM max 0,25, PAE min 14,7 A).
Cette piste teste cinq autres partenaires potentiels du divisome mycobactérien : FtsZ
(Rv2150c, l'anneau contractile lui-même), FtsW (Rv2154c, polymérase de peptidoglycane),
FtsI/PbpB (Rv2163c, transpeptidase), FtsK (Rv2748c, translocase d'ADN) et SepF (Rv2147c,
ancrage membranaire de FtsZ chez les Actinobactéries, où FtsA est absent).

CALIBRATION OBLIGATOIRE (garde-fou fixé dans pistes.md avant tout run) : Boltz-2 n'est
PAS AF3, ses scores ne sont pas comparables en absolu. Le contrôle positif DivIC-FtsQ
(interface divisome connue, déjà validée en AF3 : ipTM 0,38 reproductible, PAE min 5,1 A)
est donc refait ICI, dans Boltz, pour établir le plafond d'ipTM que CE système atteint sur
une VRAIE interface. Toute paire du panel est lue relativement à ce plafond, jamais en
absolu, et la métrique décisive reste le PAE inter-chaînes minimum et sa reproductibilité
sur plusieurs graines — pas l'ipTM seul (cf. `~/.claude/knowledge/tuberculosis.md`,
"tester une interaction par AF-Multimer/AF3 : toujours un contrôle positif du même système").

Rv1025 réutilise le MSA déjà calculé par AF3 (résultats/af3_out/.../unpaired_msa_chains_a.a3m).
Le contrôle positif DivIC-FtsQ réutilise SES DEUX MSA réels (déjà calculés par le job AF3
correspondant), ce qui lui donne un repli individuel confiant des deux côtés. Les 5 partenaires
du panel n'ont PAS de MSA propre (aucun serveur MSA local) : Boltz les traite en single-sequence,
ce qui dégrade leur repli individuel. ASYMÉTRIE CONSCIENTE : le plafond mesuré par le contrôle
est donc optimiste par rapport aux conditions du panel — l'interprétation devra vérifier le
repli individuel des partenaires (pLDDT hors interface) avant de lire un ipTM bas comme une
absence d'interaction plutôt que comme un simple repli raté faute de MSA.

Lit   : UniProt (séquences), résultats/af3_out/.../unpaired_msa_chains_a.a3m (MSA Rv1025)
Écrit : résultats/phase33_divisome_panel/<paire>.yaml (x6), run_divisome_panel.sh
"""
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "résultats/phase33_divisome_panel"
OUT.mkdir(parents=True, exist_ok=True)
RV1025_MSA = (ROOT / "résultats/af3_out/fold_rv1025_divic_test/msas/"
             "fold_rv1025_divic_test_unpaired_msa_chains_a.a3m")
assert RV1025_MSA.exists(), f"MSA Rv1025 introuvable : {RV1025_MSA}"

# MSA du contrôle positif : le job AF3 divIC-ftsQ a déjà calculé les DEUX côtés — les
# réutiliser rend la calibration fidèle (repli individuel confiant des deux chaînes),
# ce que le panel lui-même NE POURRA PAS avoir (aucun MSA local pour FtsZ/FtsW/PbpB/
# FtsK/SepF, pas de serveur MSA disponible dans cet environnement). Asymétrie
# consciente : le plafond mesuré par le contrôle est donc un plafond OPTIMISTE,
# atteint dans de meilleures conditions que le panel — à rappeler à l'interprétation,
# pas seulement en calibrant l'ipTM mais en regardant si les CHAÎNES SEULES (repli
# individuel) sont déjà confiantes dans le panel avant de lire l'interface.
POSCTRL_MSA_A = (ROOT / "résultats/af3_out/fold_divic_ftsq_posctrl/msas/"
                 "fold_divic_ftsq_posctrl_unpaired_msa_chains_a.a3m")
POSCTRL_MSA_B = (ROOT / "résultats/af3_out/fold_divic_ftsq_posctrl/msas/"
                 "fold_divic_ftsq_posctrl_unpaired_msa_chains_b.a3m")
assert POSCTRL_MSA_A.exists() and POSCTRL_MSA_B.exists(), "MSA du contrôle positif manquants"

RV1025 = ("VVTRQLGRAPRGVLAIAYRCPNGEPGVVKTAPRLPDGTPFPTLYYLTHPVLTAAASRLETTGLMREMNRRLGQDAELAAAYRR"
          "AHESYLSERDALEPLGTTVSAGGMPDRVKCLHVLIAHSLAKGPGLNPFGDEALALLAAEPRTAATLVAGQWR")

# Contrôle positif : séquences déjà utilisées pour le job AF3 divIC-ftsQ (extraites de
# résultats/af3_out/fold_divic_ftsq_posctrl/..._job_request.json), reprises telles quelles.
DIVIC = ("MPEAKRPESKRRSPASRPGKAGDSVRGGRATKPSAKPSTPAPHASRKTTRTPHEHIVEPIKRAITESVEKRSEQRLGFTARR"
         "AAILAAVVCVLTLTIARPVRTYFAQRAEMEQLAATEAMLRRQIADLEEQQVKLADPAYIAAQARERLGFVMPGDIPFQVQLP"
         "STPLAPPQPGSDAATATNNEPWYTALWHTIADDPHLPPAAPPAPEPGRPGPLPPASPNPEQPGG")
FTSQ = ("MTEHNEDPQIERVADDAADEEAVTEPLATESKDEPAEHPEFEGPRRRARRERAERRAAQARATAIEQARRAAKRRARGQIVS"
        "EQNPAKPAARGVVRGLKALLATVVLAVVGIGLGLALYFTPAMSAREIVIIGIGAVSREEVLDAARVRPATPLLQIDTQQVAD"
        "RVATIRRVASARVQRQYPSALRITIVERVPVVVKDFSDGPHLFDRDGVDFATDPPPPALPYFDVDNPGPSDPTTKAALQVLT"
        "ALHPEVASQVGRIAAPSVASITLTLADGRVVIWGTTDRCEEKAEKLAALLTQPGRTYDVSSPDLPTVK")

PANEL = {"ftsZ_Rv2150c": "P9WN95", "ftsW_Rv2154c": "P9WN97", "pbpB_ftsI_Rv2163c": "L0T911",
         "ftsK_Rv2748c": "P9WNA3", "sepF_Rv2147c": "P9WGJ5"}


def fetch(acc):
    r = subprocess.run(["curl", "-sS", "--max-time", "60",
                        f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"],
                       capture_output=True, text=True)
    lines = r.stdout.splitlines()
    seq = "".join(l.strip() for l in lines if l and not l.startswith(">"))
    assert seq, f"séquence vide pour {acc} (échec réseau ?)"
    return seq


def yaml_pair(seq_a, seq_b, msa_a=None, msa_b=None):
    a_block = f"      msa: {msa_a}\n" if msa_a else "      msa: empty\n"
    b_block = f"      msa: {msa_b}\n" if msa_b else "      msa: empty\n"
    return (f"version: 1\nsequences:\n"
            f"  - protein:\n      id: A\n      sequence: {seq_a}\n{a_block}"
            f"  - protein:\n      id: B\n      sequence: {seq_b}\n{b_block}")


written = []

# Contrôle positif AVEC MSA réel des deux côtés (cf. note ci-dessus) : plafond optimiste,
# à ne pas comparer terme à terme aux ipTM du panel sans regarder aussi le repli individuel.
p = OUT / "posctrl_divIC_ftsQ.yaml"
p.write_text(yaml_pair(DIVIC, FTSQ, msa_a=str(POSCTRL_MSA_A), msa_b=str(POSCTRL_MSA_B)))
written.append(p.name)

for label, acc in PANEL.items():
    seq = fetch(acc)
    p = OUT / f"rv1025_{label}.yaml"
    p.write_text(yaml_pair(RV1025, seq, msa_a=str(RV1025_MSA)))
    written.append(p.name)
    print(f"{label:<20} {acc}  {len(seq):>4} aa  -> {p.name}")

runner = OUT / "run_divisome_panel.sh"
runner.write_text(
    "#!/bin/bash\n"
    "# Panel divisome (P4.3) — 6 jobs (1 contrôle positif + 5 partenaires), séquentiel.\n"
    "# Suit le garde-fou du skill boltz : CPU exclusif, plusieurs heures par job.\n"
    "set -e\n"
    f"cd {OUT}\n"
    + "\n".join(
        f'echo "=== {n} ==="; ~/venvs/boltz/bin/boltz predict {n} --out_dir . '
        f'--accelerator cpu --devices 1 --output_format mmcif --diffusion_samples 1 '
        f'--num_workers 4' for n in written)
    + "\necho DONE\n"
)
runner.chmod(0o755)
print(f"\n{len(written)} jobs écrits dans {OUT.relative_to(ROOT)}/")
print(f"Lanceur : {runner.relative_to(ROOT)}")
