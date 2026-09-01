# phase15_pocket_figure.pml — P5.1 druggability figure (PyMOL headless)
# Renders two panels from the apo AF model + fpocket alpha spheres:
#   panelA = whole protein, both pockets located (druggable Pocket 1 + metal Pocket 3)
#   panelB = close-up of the conserved metal cavity (Pocket 3) with the Cys113/His115/Glu59 triad
# Run: pymol -cq phase15_pocket_figure.pml
# fpocket STP alpha spheres: chain C, resi = pocket index (1..7).

reinitialize
bg_color white
set ray_opaque_background, 0
set ray_shadows, 0
set specular, 0.15
set surface_quality, 1

load résultats/druggability/AF-P96375-F1_out/AF-P96375-F1_out.pdb, m
hide everything

select prot,  m and chain A and polymer
select poc1,  m and chain C and resi 1
select poc3,  m and chain C and resi 3
select triad, m and chain A and resi 59+113+115

# --- protein ---
show cartoon, prot
color grey70, prot
set cartoon_transparency, 0.1, prot

# --- pockets as alpha-sphere clouds ---
show spheres, poc1
set sphere_scale, 0.45, poc1
color marine, poc1
show spheres, poc3
set sphere_scale, 0.45, poc3
color orange, poc3

# --- triad sticks ---
show sticks, triad
set stick_radius, 0.22, triad
color yellow, triad and elem C
util.cnc triad

# ============ Panel A: overall, semi-transparent surface ============
set transparency, 0.55
show surface, prot
color grey85, prot
orient prot
turn y, 25
turn x, 10
ray 1700, 1500
png /tmp/claude-1000/-home-christophe-docs-codes-mtbc-Rv1025/d9a98c91-4927-4b5f-9291-2769bad55c2b/scratchpad/panelA.png, dpi=300

# ============ Panel B: close-up of the metal cavity ============
hide surface, prot
set cartoon_transparency, 0.65, prot
set label_size, 22
set label_color, black
set label_outline_color, white
set label_font_id, 7
label m and chain A and resi 113 and name CB, "Cys113"
label m and chain A and resi 115 and name CB, "His115"
label m and chain A and resi 59 and name CB, "Glu59"
orient triad
zoom triad, 5
ray 1700, 1500
png /tmp/claude-1000/-home-christophe-docs-codes-mtbc-Rv1025/d9a98c91-4927-4b5f-9291-2769bad55c2b/scratchpad/panelB.png, dpi=300
