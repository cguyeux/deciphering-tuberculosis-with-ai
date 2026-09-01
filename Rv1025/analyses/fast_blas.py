"""fast_blas.py — garantir une BLAS rapide avant tout calcul matriciel lourd.

PROBLÈME (mesuré le 2026-07-31 en lançant phase26_dca.py) : le python système de
cette machine a un numpy lié à la BLAS de RÉFÉRENCE netlib (`cblas`), mono-thread.
Benchmark matmul 1200^3 : **0,4 GFLOPS** contre **15,4 GFLOPS** pour un numpy
installé par pip (qui embarque `openblas64` dans sa wheel), soit un facteur ~38.
Concrètement, un DCA de 155 positions sur 8 700 séquences passe de ~30 s à ~15 min.
Le piège est silencieux : le calcul est JUSTE, il est seulement 38x trop lent, donc
rien ne signale l'anomalie — on croit que l'analyse est intrinsèquement coûteuse et
on la sous-dimensionne.

SOLUTION : aucun paquet système à installer, aucun sudo. Il suffit d'exécuter le
script avec un interpréteur dont le numpy vient de pip. `ensure()` le fait tout
seul en se relançant dans le premier interpréteur rapide trouvé.

Usage, en tête de script AVANT les gros calculs :

    from fast_blas import ensure
    ensure()                 # se relance dans un interpréteur OpenBLAS si besoin
    import numpy as np
"""
import os
import subprocess
import sys
from pathlib import Path

# Interpréteurs candidats, du plus spécifique au plus général. Un venv créé par pip
# suffit : les wheels numpy de PyPI embarquent OpenBLAS.
CANDIDATES = [
    Path.home() / "venvs/boltz/bin/python",
    Path.home() / "venvs/sci/bin/python",
    Path("/usr/bin/python3"),
]
_GUARD = "FAST_BLAS_REEXEC"
_PROBE = ("import numpy;"
          "print(numpy.show_config('dicts')['Build Dependencies']['blas']['name'])")


def blas_name(python=None):
    """Nom de la BLAS derrière numpy, pour l'interpréteur courant ou un autre."""
    if python is None:
        import numpy as np
        try:
            return str(np.show_config("dicts")["Build Dependencies"]["blas"]["name"]).lower()
        except Exception:
            return "inconnue"
    try:
        out = subprocess.run([str(python), "-c", _PROBE], capture_output=True,
                             text=True, timeout=60)
        return out.stdout.strip().lower() if out.returncode == 0 else ""
    except Exception:
        return ""


def is_fast(name):
    return any(k in name for k in ("openblas", "mkl", "accelerate", "blis"))


def ensure(verbose=True):
    """Si la BLAS courante est lente, se relancer dans un interpréteur rapide.

    Ne fait rien si la BLAS est déjà rapide, si aucun candidat n'est disponible
    (le calcul tourne alors quand même, en avertissant), ou si on est déjà dans
    la relance (garde anti-boucle par variable d'environnement).
    """
    if os.environ.get(_GUARD):
        return
    current = blas_name()
    if is_fast(current):
        if verbose:
            print(f"[fast_blas] BLAS = {current} : rapide, on continue.")
        return
    # NE PAS comparer les chemins RÉSOLUS : `python -m venv` crée un SYMLINK vers
    # l'interpréteur système, donc un venv se résout sur le même binaire que nous.
    # Ce qui distingue un venv n'est pas son binaire mais son site-packages (et donc
    # son numpy). Comparer les chemins tels quels suffit ; la boucle est déjà bloquée
    # par la variable de garde. (Bug vécu le 2026-07-31 : le venv était écarté à tort.)
    for cand in CANDIDATES:
        if not cand.exists() or str(cand) == sys.executable:
            continue
        name = blas_name(cand)
        if is_fast(name):
            if verbose:
                print(f"[fast_blas] BLAS courante = {current} (lente, mono-thread) ; "
                      f"relance avec {cand} ({name}).")
            os.environ[_GUARD] = "1"
            os.execv(str(cand), [str(cand), *sys.argv])
    if verbose:
        print(f"[fast_blas] ATTENTION : BLAS = {current} (netlib, mono-thread, ~0,4 GFLOPS) "
              f"et aucun interpréteur OpenBLAS trouvé. Le calcul sera JUSTE mais ~38x plus "
              f"lent. Remède : python -m venv ~/venvs/sci && ~/venvs/sci/bin/pip install numpy")


def benchmark(n=1200):
    """GFLOPS effectifs sur un produit matriciel n^3 (diagnostic)."""
    import time

    import numpy as np
    A = np.random.rand(n, n)
    B = np.random.rand(n, n)
    _ = A @ B                               # chauffe
    t = time.time()
    _ = A @ B
    return 2 * n ** 3 / (time.time() - t) / 1e9


if __name__ == "__main__":
    print(f"interpréteur : {sys.executable}")
    print(f"BLAS         : {blas_name()}")
    print(f"performance  : {benchmark():.1f} GFLOPS")
