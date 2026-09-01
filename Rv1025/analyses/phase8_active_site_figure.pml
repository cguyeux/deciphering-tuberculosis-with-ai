# P6 - Figure du site actif métallique de Rv1025 (DUF501), modèle holo Fe (AF3).
# Rendu headless : pymol -cq phase8_active_site_figure.pml

load /home/christophe/docs/codes/mtbc/Rv1025/résultats/af3_metal_out/fold_rv1025_holo_fe/fold_rv1025_holo_fe_model_2.cif, rv

bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set antialias, 2
set cartoon_fancy_helices, 1
set valence, 0
hide everything

# --- protéine ---
show cartoon, rv and chain A and polymer
color grey70, rv and chain A
set cartoon_transparency, 0.55

# --- ion métallique ---
select metal, resn FE
show spheres, metal
color orange, metal
set sphere_scale, 0.45, metal

# --- triade coordinante ---
select triad, chain A and resi 59+113+115
show sticks, triad
set stick_radius, 0.20, triad
color slate, triad and elem C
util.cnc triad

# --- 2e sphère (poche) ---
select shell2, chain A and resi 92+112
show sticks, shell2
set stick_radius, 0.12, shell2
color palegreen, shell2 and elem C
util.cnc shell2

# --- liaisons de coordination (pointillés) ---
distance c1, metal, chain A and resi 113 and name SG
distance c2, metal, chain A and resi 115 and name ND1
distance c3, metal, chain A and resi 59 and name OE2
hide labels, c1 c2 c3
color grey30, c1 c2 c3
set dash_radius, 0.06
set dash_gap, 0.25

# --- labels résidus ---
set label_size, 22
set label_color, black
set label_outline_color, white
set label_font_id, 7
label chain A and resi 113 and name CA, "Cys113"
label chain A and resi 115 and name CA, "His115"
label chain A and resi 59  and name CA, "Glu59"
label metal, "Fe"

# ======== Vue 1 : gros plan du site ========
orient triad or metal
zoom triad or metal, 3.5
turn y, 15
ray 2200, 1800
png /home/christophe/docs/codes/mtbc/Rv1025/article/figures/active_site_closeup.png, dpi=600

# ======== Vue 2 : vue d'ensemble du repli ========
set cartoon_transparency, 0.15
spectrum count, rainbow, rv and chain A and name CA
ray 2000, 1800
png /home/christophe/docs/codes/mtbc/Rv1025/article/figures/fold_overview.png, dpi=600
