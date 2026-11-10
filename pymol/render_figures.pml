# =============================================================================
#  Structure figures: wild-type CasRx vs deletion variants D7 and D9,
#  with the experimental complex PDB 9M31.
#
#    pymol -cq render_figures.pml
#
#  Everything runs through cmd.* inside a python block. Bare PyMOL commands
#  silently rename short object names (load f, x -> object "x_"), after which
#  every selection matches nothing and the renders come out blank.
#
#  Files needed alongside this script:
#    9M31.cif
#    cas13d_boltz2_structure_eval002_bin_ctrl_CasRx_crRNA.pdb
#    cas13d_boltz2_structure_eval002_bin_D7_r0_crRNA.pdb
#    cas13d_boltz2_structure_eval002_bin_D9_r0_crRNA.pdb
#
#  Domain colours: NTD purple, HEPN-1 red, Helical-1 yellow, Helical-2 green,
#  HEPN-2 orange. Wild type saturated, variant pale, same hue.
#  Catalytic residues take their own domain's colour.
# =============================================================================

python
from pymol import cmd

cmd.reinitialize()
for k, v in dict(max_threads=1, ray_shadows=0, antialias=1, depth_cue=0,
                 orthoscopic=1, ray_opaque_background=0, cartoon_ring_mode=3,
                 cartoon_ring_finder=1, sphere_scale=0.4, stick_radius=0.25,
                 dash_width=3, label_size=18).items():
    cmd.set(k, v)
cmd.bg_color("white")

# ---- palette --------------------------------------------------------------
for n, rgb in [("wt_ntd",[.49,.22,.64]), ("wt_hepn1",[.78,.13,.13]),
               ("wt_hel1",[.90,.72,.05]), ("wt_hel2",[.13,.55,.28]),
               ("wt_hepn2",[.95,.52,.10]), ("wt_rna",[.35,.35,.40]),
               ("va_ntd",[.78,.63,.88]), ("va_hepn1",[.97,.60,.58]),
               ("va_hel1",[.99,.91,.55]), ("va_hel2",[.60,.85,.65]),
               ("va_hepn2",[.99,.80,.55]), ("va_rna",[.70,.72,.80])]:
    cmd.set_color(n, rgb)

# ---- models: domain ranges derived by NW mapping of each model onto CasRx --
MODELS = {
 "wt": dict(f="cas13d_boltz2_structure_eval002_bin_ctrl_CasRx_crRNA.pdb", tone="wt",
            dom=dict(ntd=(1,180), hepn1=(181,310), hel1=(311,470),
                     hel2=(471,790), hepn2=(791,967)),
            cat1=(239,244), cat2=(858,863)),
 "d7": dict(f="cas13d_boltz2_structure_eval002_bin_D7_r0_crRNA.pdb", tone="va",
            dom=dict(ntd=(1,154), hepn1=(155,274), hel1=(275,409),
                     hel2=(410,677), hepn2=(678,799)),
            cat1=(203,208), cat2=(706,711)),
 "d9": dict(f="cas13d_boltz2_structure_eval002_bin_D9_r0_crRNA.pdb", tone="va",
            dom=dict(ntd=(1,108), hepn1=(109,211), hel1=(212,316),
                     hel2=(317,548), hepn2=(549,656)),
            cat1=(153,158), cat2=(576,581)),
}
# segments removed from CasRx, in CasRx/9M31 numbering
DELETED = {
 "d7": "41-48+62-68+88-94+177-190+372-389+446-452+577-583+633-644+739-748+768-829+952-967",
 "d9": ("39-82+87-96+156-162+170-194+260-265+301-307+345-351+371-398+433-452+"
        "564-584+632-645+687-695+740-757+765-830+900-911+951-967"),
}

def setup(tag):
    m = MODELS[tag]; t = m["tone"]
    cmd.load(m["f"], tag)
    cmd.remove(f"{tag} and not polymer")
    cmd.remove(f"{tag} and hydrogens")
    cmd.dss(tag)                      # cartoons need secondary structure assigned
    cmd.select(f"{tag}_prot", f"{tag} and chain A")
    cmd.select(f"{tag}_rna",  f"{tag} and chain B")
    for d, (a, b) in m["dom"].items():
        cmd.select(f"{tag}_{d}", f"{tag}_prot and resi {a}-{b}")
        cmd.color(f"{t}_{d}", f"{tag}_{d}")
    cmd.color(f"{t}_rna", f"{tag}_rna")
    cmd.select(f"{tag}_cat1", f"{tag}_prot and resi {m['cat1'][0]}+{m['cat1'][1]}")
    cmd.select(f"{tag}_cat2", f"{tag}_prot and resi {m['cat2'][0]}+{m['cat2'][1]}")
    cmd.color(f"{t}_hepn1", f"{tag}_cat1")
    cmd.color(f"{t}_hepn2", f"{tag}_cat2")
    cmd.select(f"{tag}_cat", f"{tag}_cat1 or {tag}_cat2")
    print(f"{tag}: {cmd.count_atoms(f'{tag}_prot and name CA')} residues, "
          f"{cmd.count_atoms(f'{tag}_rna and name P')} nt")

cmd.load("9M31.cif", "xtal")
cmd.remove("xtal and not polymer")
cmd.dss("xtal")
for tag in MODELS:
    setup(tag)

W, H = 1600, 1200
def shot(name):
    cmd.png(name, width=W, height=H, dpi=200, ray=1)
    print("  wrote", name)

def rebuild_base():
    cmd.hide("everything")
    cmd.rebuild()

# ===========================================================================
#  Figure 1 — wild type alone, domains labelled
# ===========================================================================
cmd.disable("all"); cmd.enable("wt")
rebuild_base()
cmd.show("cartoon", "wt_prot or wt_rna")
cmd.orient("wt_prot")
shot("figS1_wildtype_domains.png")

# ===========================================================================
#  Figures 2a/2b — each variant overlaid on wild type
# ===========================================================================
for var in ("d7", "d9"):
    cmd.disable("all"); cmd.enable("wt"); cmd.enable(var)
    cmd.cealign("wt_prot", f"{var}_prot")
    rebuild_base()
    cmd.show("cartoon", f"wt_prot or {var}_prot or wt_rna or {var}_rna")
    cmd.orient("wt_prot")
    shot(f"figS2_overlay_{var}.png")

    # --- active site close-up, catalytic residues in their domain colours ---
    rebuild_base()
    cmd.show("cartoon", f"wt_prot or {var}_prot")
    cmd.set("cartoon_transparency", 0.88, "wt_prot")
    cmd.set("cartoon_transparency", 0.88, f"{var}_prot")
    cmd.show("sticks",  f"wt_cat or {var}_cat")
    cmd.show("spheres", f"(wt_cat or {var}_cat) and name CA")
    m, w = MODELS[var], MODELS["wt"]
    cmd.distance("dwt", f"wt_prot and resi {w['cat1'][0]} and name CA",
                        f"wt_prot and resi {w['cat2'][0]} and name CA")
    cmd.distance("dva", f"{var}_prot and resi {m['cat1'][0]} and name CA",
                        f"{var}_prot and resi {m['cat2'][0]} and name CA")
    cmd.color("wt_hepn1", "dwt"); cmd.color("va_hepn1", "dva")
    d_wt = cmd.get_distance(f"wt_prot and resi {w['cat1'][0]} and name CA",
                            f"wt_prot and resi {w['cat2'][0]} and name CA")
    d_va = cmd.get_distance(f"{var}_prot and resi {m['cat1'][0]} and name CA",
                            f"{var}_prot and resi {m['cat2'][0]} and name CA")
    print(f"  {var}: catalytic Ca-Ca  WT {d_wt:.2f} A   {var.upper()} {d_va:.2f} A")
    cmd.orient(f"wt_cat or {var}_cat")
    cmd.zoom(f"wt_cat or {var}_cat", 12)
    shot(f"figS3_active_site_{var}.png")
    cmd.delete("dwt"); cmd.delete("dva")
    cmd.set("cartoon_transparency", 0, "wt_prot")
    cmd.set("cartoon_transparency", 0, f"{var}_prot")

    # --- crRNA channel --------------------------------------------------
    cmd.select("wt_pocket", "byres (wt_prot within 4.5 of wt_rna)")
    cmd.select(f"{var}_pocket", f"byres ({var}_prot within 4.5 of {var}_rna)")
    rebuild_base()
    cmd.show("cartoon", f"wt_prot or {var}_prot")
    cmd.set("cartoon_transparency", 0.82, "wt_prot")
    cmd.set("cartoon_transparency", 0.82, f"{var}_prot")
    cmd.show("sticks", f"wt_pocket or {var}_pocket")
    cmd.show("cartoon", f"wt_rna or {var}_rna")
    cmd.orient(f"wt_rna or {var}_rna")
    print(f"  {var}: crRNA contacts  WT {cmd.count_atoms('wt_pocket and name CA')}  "
          f"{var.upper()} {cmd.count_atoms(f'{var}_pocket and name CA')}")
    shot(f"figS4_crRNA_pocket_{var}.png")
    cmd.set("cartoon_transparency", 0, "wt_prot")
    cmd.set("cartoon_transparency", 0, f"{var}_prot")

# ===========================================================================
#  Figure 3 — deleted segments painted on wild type
# ===========================================================================
for var, sel in DELETED.items():
    cmd.disable("all"); cmd.enable("wt")
    rebuild_base()
    cmd.show("cartoon", "wt_prot or wt_rna")
    cmd.color("grey80", "wt_prot")
    cmd.color("wt_rna", "wt_rna")
    cmd.select("deleted", f"wt_prot and resi {sel}")
    cmd.color("wt_hepn1", "deleted")
    cmd.show("spheres", "wt_cat and name CA")
    cmd.color("wt_hepn2", "wt_cat")
    cmd.orient("wt_prot")
    print(f"  {var}: {cmd.count_atoms('deleted and name CA')} residues deleted")
    shot(f"figS5_deleted_{var}.png")

# ===========================================================================
#  Figure 4 — prediction against the experimental complex
# ===========================================================================
for var in ("d7", "d9"):
    cmd.disable("all"); cmd.enable("xtal"); cmd.enable(var)
    cmd.cealign("xtal and chain A", f"{var}_prot")
    rebuild_base()
    cmd.show("cartoon", f"xtal or {var}_prot or {var}_rna")
    cmd.color("grey70", "xtal and chain A")
    cmd.color("wt_rna", "xtal and chain B")
    cmd.orient("xtal and chain A")
    shot(f"figS6_vs_9M31_{var}.png")

print("all figures written")
python end
