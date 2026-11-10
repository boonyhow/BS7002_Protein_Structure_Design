python
from pymol import cmd
cmd.reinitialize()
for k,v in dict(max_threads=1, ray_shadows=0, antialias=1, depth_cue=0, orthoscopic=1,
                ray_opaque_background=0, cartoon_ring_mode=3, cartoon_ring_finder=1,
                sphere_scale=0.42, stick_radius=0.28, dash_width=4, dash_gap=0.35,
                dash_radius=0.08, label_size=22, label_color="black").items():
    cmd.set(k,v)
cmd.bg_color("white")
cmd.set_color("c_hepn1",[.78,.13,.13])      # HEPN-1 red
cmd.set_color("c_hepn2",[.95,.52,.10])      # HEPN-2 orange
cmd.set_color("pale",   [.88,.88,.90])      # everything else
cmd.set_color("rna_wt", [.20,.42,.65])      # wild-type guide, blue
cmd.set_color("rna_va", [.85,.35,.12])      # variant guide, orange-red

M={"wt":("cas13d_boltz2_structure_eval002_bin_ctrl_CasRx_crRNA.pdb",(239,244),(858,863)),
   "d7":("cas13d_boltz2_structure_eval002_bin_D7_r0_crRNA.pdb",(203,208),(706,711)),
   "d9":("cas13d_boltz2_structure_eval002_bin_D9_r0_crRNA.pdb",(153,158),(576,581))}
for t,(f,c1,c2) in M.items():
    cmd.load(f,t); cmd.remove(f"{t} and not polymer"); cmd.remove(f"{t} and hydrogens"); cmd.dss(t)
    cmd.select(f"{t}_prot",f"{t} and chain A"); cmd.select(f"{t}_rna",f"{t} and chain B")
    cmd.select(f"{t}_c1",f"{t}_prot and resi {c1[0]}+{c1[1]}")
    cmd.select(f"{t}_c2",f"{t}_prot and resi {c2[0]}+{c2[1]}")
    cmd.select(f"{t}_cat",f"{t}_c1 or {t}_c2")

# ---------- ACTIVE SITE: grey ribbon, only catalytic residues coloured -------
for var in ("d7","d9"):
    cmd.disable("all"); cmd.enable("wt"); cmd.enable(var)
    cmd.cealign("wt_prot",f"{var}_prot")
    cmd.hide("everything"); cmd.rebuild()
    cmd.show("cartoon",f"wt_prot or {var}_prot")
    cmd.color("pale",f"wt_prot or {var}_prot")
    cmd.set("cartoon_transparency",0.80,"wt_prot")
    cmd.set("cartoon_transparency",0.80,f"{var}_prot")
    for tag,alpha in (("wt",1.0),(var,1.0)):
        cmd.show("sticks",f"{tag}_cat"); cmd.show("spheres",f"{tag}_cat and name CA")
        cmd.color("c_hepn1",f"{tag}_c1"); cmd.color("c_hepn2",f"{tag}_c2")
    w=M["wt"]; m=M[var]
    cmd.distance("dw",f"wt_prot and resi {w[1][0]} and name CA",f"wt_prot and resi {w[2][0]} and name CA")
    cmd.distance("dv",f"{var}_prot and resi {m[1][0]} and name CA",f"{var}_prot and resi {m[2][0]} and name CA")
    cmd.hide("labels")                       # numbers go in the caption, not on the image
    cmd.color("grey40","dw"); cmd.color("black","dv")
    cmd.orient(f"wt_cat or {var}_cat"); cmd.zoom(f"wt_cat or {var}_cat",14)
    cmd.png(f"figF_active_{var}.png",width=1500,height=1200,dpi=200,ray=1)
    print("wrote active",var)
    cmd.delete("dw"); cmd.delete("dv")
    cmd.set("cartoon_transparency",0,"wt_prot"); cmd.set("cartoon_transparency",0,f"{var}_prot")

# ---------- crRNA: protein decoloured, only the two guides in colour --------
for var in ("d7","d9"):
    cmd.disable("all"); cmd.enable("wt"); cmd.enable(var)
    cmd.cealign("wt_prot",f"{var}_prot")
    cmd.hide("everything"); cmd.rebuild()
    cmd.show("cartoon",f"wt_prot or {var}_prot")
    cmd.color("pale",f"wt_prot or {var}_prot")
    cmd.set("cartoon_transparency",0.90,"wt_prot")
    cmd.set("cartoon_transparency",0.90,f"{var}_prot")
    cmd.show("cartoon","wt_rna"); cmd.show("cartoon",f"{var}_rna")
    cmd.color("rna_wt","wt_rna"); cmd.color("rna_va",f"{var}_rna")
    cmd.set("cartoon_nucleic_acid_mode",4)
    cmd.orient("wt_rna")
    cmd.png(f"figF_crRNA_{var}.png",width=1500,height=1200,dpi=200,ray=1)
    print("wrote crRNA",var)
    cmd.set("cartoon_transparency",0,"wt_prot"); cmd.set("cartoon_transparency",0,f"{var}_prot")
print("done")
python end
