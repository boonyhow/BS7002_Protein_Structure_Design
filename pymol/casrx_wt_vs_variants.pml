# =============================================================================
# CasRx wild-type (Boltz control) vs deletion variants D7 and D9
#
#   pymol -cq casrx_wt_vs_variants.pml
#
# Files needed in the working directory:
#   cas13d_boltz2_structure_eval002_bin_ctrl_CasRx_crRNA.pdb   (WT, 967 aa)
#   cas13d_boltz2_structure_eval002_bin_D7_r0_crRNA.pdb        (D7, 799 aa)
#   cas13d_boltz2_structure_eval002_bin_D9_r0_crRNA.pdb        (D9, 656 aa)
#   9M31.cif                                                    (optional, experimental)
#
# Writes: *_WT_D7.png and *_WT_D9.png
#
# COLOUR SCHEME
#   NTD        purple     HEPN-1  red      HEPN-2  orange
#   Helical-1  yellow     Helical-2  green  crRNA  grey / slate
#   Wild-type gets the saturated shade, the variant a pale shade of the SAME
#   hue, so an overlaid pair reads as one domain in two tones.
#   Catalytic residues follow their own domain: HEPN-1 sticks are red,
#   HEPN-2 sticks are orange. Nothing is black.
#
# NUMBERING
#   Each model is numbered from 1, and the variants are shorter, so the same
#   residue has a different number in each. Ranges below were derived by
#   Needleman-Wunsch mapping of each model back onto CasRx.
#
#   domain       WT          D7          D9
#   NTD           1-180       1-154       1-108
#   HEPN-1      181-310     155-274     109-211
#   Helical-1   311-470     275-409     212-316
#   Helical-2   471-790     410-677     317-548
#   HEPN-2      791-967     678-799     549-656
#   catalytic   239,244     203,208     153,158     (HEPN-1 R-X4-H)
#               858,863     706,711     576,581     (HEPN-2 R-X4-H)
#
#   Domain boundaries are APPROXIMATE, anchored on the catalytic motifs.
# =============================================================================

reinitialize
set ray_shadows, 0
set antialias, 1
set depth_cue, 0
set orthoscopic, 1
set ray_opaque_background, 0
set cartoon_transparency, 0
set cartoon_ring_mode, 3
set cartoon_ring_finder, 1
set dash_width, 3
set dash_gap, 0.4
set label_size, 20
set sphere_scale, 0.4
set stick_radius, 0.25
bg_color white

# ---- palette ----------------------------------------------------------------
# saturated = wild-type, pale = variant, same hue
set_color wt_ntd,    [0.49, 0.22, 0.64]
set_color wt_hepn1,  [0.78, 0.13, 0.13]
set_color wt_hel1,   [0.90, 0.72, 0.05]
set_color wt_hel2,   [0.13, 0.55, 0.28]
set_color wt_hepn2,  [0.95, 0.52, 0.10]
set_color wt_rna,    [0.35, 0.35, 0.40]

set_color va_ntd,    [0.78, 0.63, 0.88]
set_color va_hepn1,  [0.97, 0.60, 0.58]
set_color va_hel1,   [0.99, 0.91, 0.55]
set_color va_hel2,   [0.60, 0.85, 0.65]
set_color va_hepn2,  [0.99, 0.80, 0.55]
set_color va_rna,    [0.70, 0.72, 0.80]

# =============================================================================
python
from pymol import cmd

MODELS = {
  "wt": dict(file="cas13d_boltz2_structure_eval002_bin_ctrl_CasRx_crRNA.pdb",
             dom=dict(ntd=(1,180), hepn1=(181,310), hel1=(311,470),
                      hel2=(471,790), hepn2=(791,967)),
             cat1=(239,244), cat2=(858,863), tone="wt"),
  "d7": dict(file="cas13d_boltz2_structure_eval002_bin_D7_r0_crRNA.pdb",
             dom=dict(ntd=(1,154), hepn1=(155,274), hel1=(275,409),
                      hel2=(410,677), hepn2=(678,799)),
             cat1=(203,208), cat2=(706,711), tone="va"),
  "d9": dict(file="cas13d_boltz2_structure_eval002_bin_D9_r0_crRNA.pdb",
             dom=dict(ntd=(1,108), hepn1=(109,211), hel1=(212,316),
                      hel2=(317,548), hepn2=(549,656)),
             cat1=(153,158), cat2=(576,581), tone="va"),
}

def setup(tag):
    m = MODELS[tag]
    cmd.load(m["file"], tag)
    cmd.remove(f"{tag} and not polymer")
    cmd.remove(f"{tag} and hydrogens")
    t = m["tone"]
    for dom, (a, b) in m["dom"].items():
        cmd.select(f"{tag}_{dom}", f"{tag} and chain A and resi {a}-{b}")
        cmd.color(f"{t}_{dom}", f"{tag}_{dom}")
    cmd.select(f"{tag}_prot", f"{tag} and chain A")
    cmd.select(f"{tag}_rna",  f"{tag} and chain B")
    cmd.color(f"{t}_rna", f"{tag}_rna")
    # catalytic residues take their own domain's colour
    cmd.select(f"{tag}_cat1", f"{tag} and chain A and resi {m['cat1'][0]}+{m['cat1'][1]}")
    cmd.select(f"{tag}_cat2", f"{tag} and chain A and resi {m['cat2'][0]}+{m['cat2'][1]}")
    cmd.color(f"{t}_hepn1", f"{tag}_cat1")
    cmd.color(f"{t}_hepn2", f"{tag}_cat2")
    cmd.select(f"{tag}_cat", f"{tag}_cat1 or {tag}_cat2")

for tag in MODELS:
    setup(tag)

def pair(var, suffix):
    """Four figures for one wild-type / variant comparison."""
    m, w = MODELS[var], MODELS["wt"]
    cmd.disable("all"); cmd.enable("wt"); cmd.enable(var)
    cmd.cealign("wt_prot", f"{var}_prot")

    # --- fig 1: overlay, cartoon ------------------------------------------
    cmd.hide("everything")
    cmd.show("cartoon", f"wt_prot or {var}_prot or wt_rna or {var}_rna")
    cmd.orient("wt_prot")
    cmd.png(f"fig1_overlay_{suffix}.png", width=1600, height=1200, dpi=200, ray=1)

    # --- fig 2: variant alone, domains labelled ---------------------------
    cmd.disable("wt")
    cmd.orient(f"{var}_prot")
    cmd.png(f"fig2_domains_{suffix}.png", width=1600, height=1200, dpi=200, ray=1)
    cmd.enable("wt")

    # --- fig 3: crRNA channel ---------------------------------------------
    cmd.select("wt_pocket",  f"byres (wt_prot within 4.5 of wt_rna)")
    cmd.select(f"{var}_pocket", f"byres ({var}_prot within 4.5 of {var}_rna)")
    cmd.hide("everything")
    cmd.show("cartoon", f"wt_prot or {var}_prot")
    cmd.set("cartoon_transparency", 0.85, "wt_prot")
    cmd.set("cartoon_transparency", 0.85, f"{var}_prot")
    cmd.show("sticks", f"wt_pocket or {var}_pocket")
    cmd.show("cartoon", f"wt_rna or {var}_rna")
    cmd.orient(f"wt_rna or {var}_rna")
    cmd.png(f"fig3_crRNA_pocket_{suffix}.png", width=1600, height=1200, dpi=200, ray=1)
    cmd.set("cartoon_transparency", 0, "wt_prot")
    cmd.set("cartoon_transparency", 0, f"{var}_prot")

    # --- fig 4: active site, HEPN-coloured --------------------------------
    cmd.hide("everything")
    cmd.show("cartoon", f"wt_prot or {var}_prot")
    cmd.set("cartoon_transparency", 0.9, "wt_prot")
    cmd.set("cartoon_transparency", 0.9, f"{var}_prot")
    cmd.show("sticks",  f"wt_cat or {var}_cat")
    cmd.show("spheres", f"(wt_cat or {var}_cat) and name CA")
    cmd.distance(f"cat_wt_{suffix}",
                 f"wt_prot and resi {w['cat1'][0]} and name CA",
                 f"wt_prot and resi {w['cat2'][0]} and name CA")
    cmd.distance(f"cat_{var}_{suffix}",
                 f"{var}_prot and resi {m['cat1'][0]} and name CA",
                 f"{var}_prot and resi {m['cat2'][0]} and name CA")
    cmd.color("wt_hepn1",  f"cat_wt_{suffix}")
    cmd.color("va_hepn1",  f"cat_{var}_{suffix}")
    cmd.orient(f"wt_cat or {var}_cat")
    cmd.zoom(f"wt_cat or {var}_cat", 10)
    cmd.png(f"fig4_active_site_{suffix}.png", width=1600, height=1200, dpi=200, ray=1)
    cmd.set("cartoon_transparency", 0, "wt_prot")
    cmd.set("cartoon_transparency", 0, f"{var}_prot")

    d_wt = cmd.get_distance(f"wt_prot and resi {w['cat1'][0]} and name CA",
                            f"wt_prot and resi {w['cat2'][0]} and name CA")
    d_va = cmd.get_distance(f"{var}_prot and resi {m['cat1'][0]} and name CA",
                            f"{var}_prot and resi {m['cat2'][0]} and name CA")
    print(f"=== {suffix} ===")
    print(f"  catalytic Ca-Ca   WT {d_wt:.2f} A   {var.upper()} {d_va:.2f} A   "
          f"diff {abs(d_wt-d_va):.2f} A")
    print(f"  crRNA contacts    WT {cmd.count_atoms('wt_pocket and name CA')}   "
          f"{var.upper()} {cmd.count_atoms(f'{var}_pocket and name CA')}")
    print(f"  all-atom RMSD     {cmd.align(f'{var}_prot', 'wt_prot', cycles=0)[0]:.2f} A")

pair("d7", "WT_D7")
pair("d9", "WT_D9")

# --- fig 5: where the deletions fall, painted on wild-type -----------------
DEL = {
 "WT_D7": "41-48+62-68+88-94+177-190+372-389+446-452+577-583+633-644+739-748+768-829+952-967",
 "WT_D9": ("39-82+87-96+156-162+170-194+260-265+301-307+345-351+371-398+433-452+"
           "564-584+632-645+687-695+740-757+765-830+900-911+951-967"),
}
for suffix, sel in DEL.items():
    cmd.disable("all"); cmd.enable("wt")
    cmd.hide("everything")
    cmd.show("cartoon", "wt_prot or wt_rna")
    cmd.color("grey80", "wt_prot")
    cmd.color("wt_rna", "wt_rna")
    cmd.select("deleted", f"wt_prot and resi {sel}")
    cmd.color("wt_hepn1", "deleted")
    cmd.show("cartoon", "wt_cat")
    cmd.show("spheres", "wt_cat and name CA")
    cmd.color("wt_hepn2", "wt_cat")
    cmd.orient("wt_prot")
    cmd.png(f"fig5_deleted_regions_{suffix}.png", width=1600, height=1200, dpi=200, ray=1)
    print(f"  {suffix}: {cmd.count_atoms('deleted and name CA')} residues deleted")

print("done - wrote fig1..fig5 for WT_D7 and WT_D9")
python end
