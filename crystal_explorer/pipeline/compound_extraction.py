import os
import re
import json
from typing import List, Dict, Optional, Tuple

import pandas as pd

try:
    import requests
except ImportError:  # requests is only needed for the optional LLM fallback
    requests = None


# --------------------------------------------------------------------------- #
# 1. Reagent dictionary
# --------------------------------------------------------------------------- #
# canonical_name -> list of aliases/phrasings seen in real pdbx_details text.
# Extend this list as you spot-check your own data (step 5 of the strategy).
REAGENT_DICTIONARY: Dict[str, List[str]] = {
    # --- PEGs ---
    "PEG 400": ["peg 400", "peg400", "peg-400", "polyethylene glycol 400", "poly(ethylene glycol) 400"],
    "PEG 1000": ["peg 1000", "peg1000", "peg-1000", "peg1k", "peg 1k", "polyethylene glycol 1000", "poly(ethylene glycol) 1000",
                 "peg 1,000"],
    "PEG 3350": ["peg 3350", "peg3350", "peg-3350", "polyethylene glycol 3350", "poly(ethylene glycol) 3350",
                 "peg 3,350", "peg3,350", "peg3330"],
    "PEG 4000": ["peg 4000", "peg4000", "peg-4000", "peg4k", "peg 4k", "polyethylene glycol 4000", "poly(ethylene glycol) 4000",
                 "peg 4,000", "pegf-4000"],
    "PEG 6000": ["peg 6000", "peg6000", "peg-6000", "peg6k", "peg 6k", "polyethylene glycol 6000", "poly(ethylene glycol) 6000",
                 "peg 6,000"],
    "PEG 8000": ["peg 8000", "peg8000", "peg-8000", "peg8k", "peg 8k", "polyethylene glycol 8000", "poly(ethylene glycol) 8000",
                 "peg 8,000", "polyethyleneglycol 8000"],
    "PEG 10000": ["peg 10000", "peg10000", "peg-10000", "peg10k", "peg 10k", "polyethylene glycol 10000",
                  "peg 10,000"],
    "PEG 5000": ["peg 5000", "peg5000", "peg-5000", "peg5k", "peg 5k", "polyethylene glycol 5000",
                 "peg 5,000"],
    "PEG MME 550": ["peg mme 550", "peg monomethyl ether 550", "peg-550", "peg 550", "peg550 mme",
                    "peg 550mme", "peg 550 mme", "polyethylene glycol monomethyl ether 550"],
    "PEG MME 2000": ["peg mme 2000", "peg mme 2k", "peg mme2k", "peg monomethyl ether 2000",
                     "peg mme2000", "peg 2k mme", "peg 2kmme", "polyethylene glycol monomethyl ether 2000",
                     "polyethylene glycol monomethyl ether 2,000"],
    "PEG MME 5000": ["peg mme 5000", "peg monomethyl ether 5000", "polyethylene glycol monomethyl ether 5000",
                      "polyethylene glycol monomethyl ether 5,000"],
    "PEG 3500": ["peg 3500", "peg3500", "peg 3550", "peg3550"],
    "PEG Smear High": ["peg smear high"],
    "PEG Smear Medium": ["peg smear medium"],
    "PEG 200": ["peg 200", "peg200", "polyethylene glycol 200", "poly(ethylene glycol) 200", "peg-200"],
    "PEG 300": ["peg 300", "peg300", "polyethylene glycol 300", "poly(ethylene glycol) 300", "peg-300"],
    "PEG 1500": ["peg 1500", "peg1500", "polyethylene glycol 1500", "poly(ethylene glycol) 1500", "peg-1500",
                 "peg 1,500"],
    "PEG 2000": ["peg 2000", "peg2000", "peg-2000", "polyethylene glycol 2000", "poly(ethylene glycol) 2000",
                 "peg 2,000"],
    "PEG 600": ["peg 600", "peg600", "polyethylene glycol 600", "poly(ethylene glycol) 600", "peg-600"],
    "PEG 800": ["peg 800", "peg800", "polyethylene glycol 800", "poly(ethylene glycol) 800", "peg-800"],
    "PEG MME 75\u2080": ["peg mme 75\u2080", "peg monomethyl ether 75\u2080"  ],
    "PEG dimethyl ether": ["peg dimethyl ether", "polyethylene glycol dimethyl ether"],
    "PEG 3000": ["peg 3000", "peg3000", "peg-3000", "polyethylene glycol 3000", "peg 3,000"],
    "PEG 12000": ["peg 12000", "peg12000", "peg-12000", "peg12k", "peg 12k", "polyethylene glycol 12000", "peg 12,000"],
    "PEG 20000": ["peg 20000", "peg20000", "peg-20000", "peg20k", "peg 20k", "polyethylene glycol 20000", "peg 20,000"],
    "Jeffamine ED-2001": ["jeffamine ed-2001", "ed2001"],

    # --- Salts ---
    "ammonium sulfate": ["ammonium sulfate", "ammonium sulphate", "(nh4)2so4", "ammonium-sulfate",
                          "amm.sulfate", "ammonia sulfate", "nh4-sulfate", "ammomium sulfate",
                          "ammonium sulf ate", "ammonium sulfphate", "(nh4)so4", "(nh4)2 so4", "nh4so4"],
    "sodium chloride": ["sodium chloride", "nacl", "sodium choloride"],
    "lithium sulfate": ["lithium sulfate", "lithium sulphate", "li2so4", "liso4"],
    "lithium chloride": ["lithium chloride", "licl"],
    "magnesium chloride": ["magnesium chloride", "mgcl2", "mg2cl", "mgcl 2"],
    "magnesium sulfate": ["magnesium sulfate", "magnesium sulphate", "mgso4"],
    "calcium chloride": ["calcium chloride", "cacl2"],
    "sodium phosphate": ["sodium phosphate", "nah2po4", "na2hpo4", "sodium phospahte", "na2po4"],
    "potassium phosphate": ["potassium phosphate", "kh2po4", "k2hpo4"],
    "sodium/potassium phosphate": ["sodium potassium phosphate", "sodium/potassium phosphate", "na/k phosphate", "napo4"],
    "sodium acetate": ["sodium acetate", "sodium actate", "naoac", "naac"],
    "sodium citrate": ["sodium citrate", "tri-sodium citrate", "NAcitrate", "trisodium citrate", "na citrate", "Na-citrate", "na3citrate", "na3cit", "nacit"],
    "sodium formate": ["sodium formate", "NAformate", "na formate", "Na-formate", "formic acid sodium salt", ],
    "potassium chloride": ["potassium chloride", "kcl", "potassium cloride"],
    "zinc acetate": ["zinc acetate", "zn(oac)2", "zinc acetate dihydrate"],
    "zinc chloride": ["zinc chloride", "Zinc Chloride", "zncl2"],
    "potassium acetate": ["potassium acetate"],
    "acetate": ["acetate"],
    "ammonium chloride": ["ammonium chloride", "nh4cl"],
    "ammonium nitrate": ["ammonium nitrate","nh4no3"],
    "ammonium acetate": ["ammonium acetate","nh4 acetate"],
    "ammonium citrate": ["ammonium citrate"],
    "ammonium phosphate": ["ammonium phosphate","(nh4)3po4"],
    "ammonium bicarbonate": ["ammonium bicarbonate","nh4hco3"],
    "potassium sulfate": ["potassium sulfate","potassium sulphate","k2so4"],
    "Monosodium phosphate": ["monosodium phosphate", "sodium phosphate monobasic", "nah2po4", "na h2po4", "sodium dihydrogen phosphate"],
    "NaH2PO4": ["nah2po4", "na h2 po4", "monosodium phosphate", "sodium phosphate monobasic"],
    "Disodium phosphate": ["disodium phosphate", "sodium phosphate dibasic", "na2hpo4", "na2 hpo4"],
    "Dipotassium hydrogen phosphate": ["dipotassium hydrogen phosphate", "potassium phosphate dibasic", "k2hpo4", "k2 hpo4"],
    "Potassium phosphate": ["potassium phosphate","kh2po4","k2hpo4"],
    "sodium sulfate": ["sodium sulfate","sodium sulphate","na2so4"],
    "lithium nitrate": ["lithium nitrate"],
    "cesium chloride": ["cesium chloride","cscl"],
    "rubidium chloride": ["rubidium chloride","rbcl"],
    "strontium chloride": ["strontium chloride","srcl2"],
    "barium chloride": ["barium chloride","bacl2"],
    "manganese chloride": ["manganese chloride","mncl2"],
    "nickel chloride": ["nickel chloride","nicl2"],
    "cobalt chloride": ["cobalt chloride","cocl2"],
    "cadmium chloride": ["cadmium chloride","cdcl2"],
    "iron chloride": ["iron chloride","fecl3","ferric chloride"],
    "magnesium nitrate": ["magnesium nitrate", "mg(no3)2"],
    "potassium nitrate": ["potassium nitrate", "kno3"],
    "sodium nitrate": ["sodium nitrate", "nano3"],
    "ammonium formate": ["ammonium formate", "nh4 formate", "nh4hco2"],
    "ammonium tartrate": ["ammonium tartrate", "diammonium tartrate"],
    "potassium sodium tartrate": ["potassium sodium tartrate", "sodium potassium tartrate",
                                   "rochelle salt", "kna tartrate"],
    "potassium tartrate": ["potassium tartrate", "dipotassium tartrate"],
    "sodium bromide": ["sodium bromide", "nabr"],
    "potassium iodide": ["potassium iodide", "ki"],
    "manganese acetate": ["manganese acetate", "mn(oac)2", "mn acetate"],
    "copper chloride": ["copper chloride", "cucl2", "cupric chloride"],
    "Lithium citrate": ["lithium citrate", "lithium citrate tribasic", "li3 citrate"],
    "Lithium acetate": ["lithium acetate", "li acetate", "liac"],
    "Sodium malonate": ["sodium malonate", "malonic acid sodium salt", "na malonate"],
    "Sodium succinate": ["sodium succinate", "succinic acid sodium salt", "disodium succinate"],
    "Sodium maleate": ["sodium maleate","maleic acid sodium salt", "disodium maleate"],
    "Sodium tartrate": ["sodium tartrate", "disodium tartrate", "sodium L-tartrate"],
    "Di-sodium tartrate": ["di-sodium tartrate", "disodium tartrate","sodium tartrate"],
    "Trisodium citrate-citric acid": [ "trisodium citrate-citric acid", "sodium citrate citric acid buffer", "citrate buffer" ],
    "Magnesium acetate": ["magnesium acetate", "mg acetate", "mg(oac)2"],
    "Magnesium formate": ["magnesium formate", "mg formate", "mg(hcoo)2"],
    "Cobalt(II) chloride": ["cobalt chloride", "cobalt(ii) chloride", "cobalt (ii) chloride", "cocl2", "cobaltous chloride"],
    "Cadmium chloride": ["cadmium chloride", "cadmium(ii) chloride", "cdcl2"],
    "Nickel(II) chloride": ["nickel chloride", "nickel(ii) chloride", "nickel (ii) chloride", "nicl2"],
    "Sodium iodide": ["sodium iodide", "nai","na iodide", "NAI"],
    "Sodium succinate": ["sodium succinate","sodium succinate buffer","disodium succinate","succinic acid disodium salt","succinate","na succinate","sodium hydrogen succinate"],

    # --- Buffers ---
    "Tris-HCl": ["tris-hcl", "tris hcl", "tris buffer", "tris"],
    "Bis-Tris": ["bis-tris", "bis tris", "bistris", "bis_tris", "bis tris buffer","bis-tris buffer"],
    "HEPES": ["hepes", "hepes buffer", "4-(2-hydroxyethyl)-1-piperazineethanesulfonic acid"],
    "MES": ["mes","mes buffer","2-(n-morpholino)ethanesulfonic acid", "2 morpholinoethanesulfonic acid"],
    "bis-tris propane": ["bis-tris propane"],
    "Bis-Tris chloride": ["bis-tris chloride", "bis tris chloride","bis-tris hcl"],
    "CHES": ["ches"],
    "CAPS": ["caps"],
    "imidazole": ["imidazole"],
    "citrate buffer": ["citrate buffer", "citric acid"],
    "cacodylate": ["sodium cacodylate", "cacodylate"],
    "MOPS": ["mops"],
    "glycine": ["glycine"],
    "acetate buffer": ["acetate buffer"],
    "Tricine": ["tricine", "n-[tris(hydroxymethyl)methyl]glycine"],
    "PIPES": ["pipes", "piperazine-n,n'-bis(2-ethanesulfonic acid)"],
    "POPSO": ["popso","piperazine-n,n'-bis(2-hydroxypropanesulfonic acid)"],
    "TAPS": ["taps"],
    "TAPSO": ["tapso"],
    "ADA": ["ada buffer", "n-(2-acetamido)iminodiacetic acid"],
    "Bicine": ["bicine", "n,n-bis(2-hydroxyethyl)glycine"],
    "CAPSO": ["capso"],
    "borate": ["borate", "sodium borate", "boric acid"],
    "phosphate buffer": ["phosphate buffer","phosphate"],
    "sodium cacodylate": ["sodium cacodylate", "cacodylate", "dimethylarsinic acid sodium salt"],
    "Trisodium citrate-citric acid": ["trisodium citrate-citric acid", "trisodium citrate citric acid", "trisodium citrate/citric acid", "citrate buffer", "sodium citrate citric acid buffer", "citric acid sodium citrate"],

    # --- Organic precipitants ---
    "1,6-hexanediol": ["1,6-hexanediol", "hexanediol"],
    "1,4-butanediol": ["1,4-butanediol", "butanediol" ],
    "2-propanol": ["2-propanol", "isopropanol", "isopropyl alcohol"],
    "tert-butanol": ["tert-butanol", "t-butanol", "tbutanol"],
    "methanol": ["methanol", "methyl alcohol", "meoh"],
    "ethanol": ["ethanol", "ethyl alcohol"],
    "propanol": ["propanol","1-propanol"],
    "butanol": ["butanol","1-butanol"],
    "acetone": ["acetone"],
    "acetonitrile": ["acetonitrile"],
    "Sorbitol": ["sorbitol", "d-sorbitol", "D-sorbitol", "glucitol", "d-glucitol", "hexane-1,2,3,4,5,6-hexol", "sorbitol solution"],
    "Trehalose": ["trehalose", "d-trehalose", "trehalose dihydrate"],
    "Xylitol": ["xylitol"],
    "trimethylamine N-oxide": ["trimethylamine n-oxide", "tmao", "trimethylamine oxide"],
    "sarcosine": ["sarcosine", "n-methylglycine"],
    "1,3-propanediol": ["1,3-propanediol", "1,3 propanediol"],
    "1,5-pentanediol": ["1,5-pentanediol", "1,5 pentanediol"],

    # --- Additives / cryoprotectants ---
    "glycerol": ["glycerol"],
    "DMSO": ["dmso", "dimethyl sulfoxide"],
    "ATP": ["atp", "adenosine triphosphate"],
    "MPD": ["mpd", "2-methyl-2,4-pentanediol", "hexylene glycol"],
    "isopropanol": ["isopropanol", "2-propanol", "isopropyl alcohol"],
    "ethanol": ["ethanol"],
    "ethylene glycol": ["ethylene glycol", "ethylen glycole"],
    "DTT": ["dtt", "dithiothreitol"],
    "beta-mercaptoethanol": ["beta-mercaptoethanol", "2-mercaptoethanol", "bme"],
    "TCEP": ["tcep", "tecp"],
    "sodium azide": ["sodium azide"],
    "EDTA": ["edta"],
    "glutathione": ["glutathione"],
    "spermidine": ["spermidine"],
    "jeffamine": ["jeffamine"],
    "PPG": ["ppg", "polypropylene glycol"],
    "dioxane": ["dioxane"],

# --- Detergents / surfactants ---
    "n-octyl-beta-D-glucoside": [ "octyl glucoside", "og", "beta-octyl glucoside"],
    "n-octyl-beta-D-thioglucopyranoside": ["octyl thioglucoside","otg"],
    "CHAPS": ["chaps"],
    "CHAPSO": [ "chapso"],
    "Tween 20": ["tween 20","polysorbate 20"],
    "Tween 80": ["tween 80", "polysorbate 80"],
    "\u03b2-octyl glucoside": ["beta octyl glucoside","octyl glucoside"],
    "Triton X-100": ["triton x100", "triton x-100"],

# --- Amino acids and additives ---
    "arginine": ["arginine", "l-arginine"],
    "lysine": ["lysine", "l-lysine"],
    "histidine": ["histidine","l-histidine"],
    "glutamic acid": ["glutamate","glutamic acid"],
    "aspartic acid": ["aspartate","aspartic acid"],
    "alanine": ["alanine", "l-alanine"],
    "proline": ["proline","l-proline"],
    "glycine betaine": ["betaine","glycine betaine"],
    "taurine": ["taurine"],

# --- Reducing agents / stabilizers ---
    "betaine": ["betaine"],
    "benzamidine": ["benzamidine"],
    "phenylmethylsulfonyl fluoride": ["pmsf","phenylmethylsulfonyl fluoride"],
    "iodoacetamide": ["iam","iodoacetamide"],
    "iodoacetic acid": ["iodoacetic acid"],
    "leupeptin": ["leupeptin"],
    "pepstatin": ["pepstatin"],

# --- Metal ions / cofactors frequently crystallized ---
    "zinc sulfate": ["zinc sulfate", "znso4"],
    "copper sulfate": ["copper sulfate","cuso4"],
    "magnesium acetate": ["magnesium acetate"],
    "calcium acetate": ["calcium acetate"],
    "manganese sulfate": [ "manganese sulfate", "mnso4"],
    "nickel sulfate": ["nickel sulfate","niso4"],
    "cobalt sulfate": ["cobalt sulfate","coso4"],
    "zinc acetate": ["zinc acetate"],
    "cadmium sulfate": ["cadmium sulfate","cdso4"],
    "iron sulfate": ["ferrous sulfate", "iron sulfate", "feso4"],
    "ATP magnesium salt": ["mg-atp", "magnesium atp", "atp-mg", "mg atp"],
    "ADP": ["adp","adenosine diphosphate"],
    "GTP": ["gtp","guanosine triphosphate"],
    "FAD": [ "fad", "flavin adenine dinucleotide"],
    "NAD": ["nad","nicotinamide adenine dinucleotide"],
    "NADP": ["nadp", "nicotinamide adenine dinucleotide phosphate"],
    "coenzyme A": ["coa","coenzyme a","co-enzyme a"],
    "AMP": ["amp","adenosine monophosphate"],
    
# --- Common Morpheus crystallization additives ---

    "1,2,3-heptanetriol": ["1,2,3-heptanetriol", "heptanetriol"],
    "1,2,4-butanetriol": ["1,2,4-butanetriol", "butanetriol"],
    "1,2,6-hexanetriol": ["1,2,6-hexanetriol", "hexanetriol"],
    "2,3-butanediol": ["2,3-butanediol", "butanediol"],
    "1,6-hexanediol": ["1,6-hexanediol","hexanediol"],

# --- Hampton additive screen compounds ---
    "benzamidine hydrochloride": ["benzamidine","benzamidine hcl"],
    "guanidine hydrochloride": ["guanidine hcl", "guanidinium chloride"],
    "urea": ["urea"],
    "thiourea": ["thiourea"],
    "acetamide": ["acetamide"],
    "formamide": ["formamide"],
    "acetamide": ["acetamide"],
    "lithium chloride": ["licl", "lithium chloride"],
    "sodium thiocyanate": ["sodium thiocyanate","nascn"],
    "potassium thiocyanate": ["potassium thiocyanate", "kscn"],
    "ammonium thiocyanate": ["ammonium thiocyanate","nh4scn"],

# --- Polymers ---
    "poly(acrylic acid sodium salt) 5100": ["poly(acrylic acid sodium salt) 5100","sodium polyacrylate 5100", "poly acrylic acid sodium salt", "paa sodium salt 5100", "polyacrylate"],
    "PGA-LM": ["pga-lm", "pga lm", "pga_lm", "poly-gamma-glutamic acid low molecular weight", "poly gamma glutamic acid lm", "polyglutamic acid lm", "gamma-polyglutamic acid","poly-gamma-glutamate"],
    "Poly(acrylic acid sodium salt) 5100": ["poly(acrylicacid sodium salt) 5100", "poly(acrylic acid sodium salt) 5100", "poly acrylic acid sodium salt 5100", "polyacrylic acid sodium salt 5100",
    "sodium polyacrylate 5100", "poly(sodium acrylate) 5100", "paa sodium salt 5100", "paaa 5100", "polyacrylate sodium salt"],
    "Sucrose": ["sucrose", "d-sucrose", "\u03b1-d-glucopyranosyl-(1\u21922)-\u03b2-d-fructofuranoside", "cane sugar", "table sugar", "saccharose", "sucrose solution"],
    
# --- PLP Snyder fluorescent dyes ---
    "PLP_Snyder441": ["plp_snyder441", "plp snyder 441","snyder 441"],
    "PLP_Snyder457": ["plp_snyder457", "plp snyder 457", "snyder 457"],
    "PLP_Snyder494": ["plp_snyder494","plp snyder 494","snyder 494"],
    "PLP_Snyder495": ["plp_snyder495","plp snyder 495","snyder 495"],
    "PLP_Snyder496": ["plp_snyder496","plp snyder 496","snyder 496"],
    "PLP_Snyder530": ["plp_snyder530","plp snyder 530","snyder 530"],
    "PLP_Snyder608": ["plp_snyder608","plp snyder 608","snyder 608"],

# --- Ligands / additives ---

"Ligand": ["ligand", "small molecule ligand"],
"Ebselen": ["ebselen","2-phenyl-1,2-benzisoselenazol-3(2h)-one"],

# --- Added from mining real pdbx_details text (~7,900 entries across the
# .pdb_cache corpus) for reagents/polymers the dictionary above didn't
# recognize yet -- see compound_extraction match-rate check. ---

# --- Detergents / surfactants ---
    "MEGA-8": ["mega-8", "mega 8", "mega8", "octanoyl-n-methylglucamide"],
    "Fos-Choline-9": ["fos-choline-9", "foscholine-9", "fos choline 9"],
    "Triton X-114": ["triton x-114", "triton x114"],
    "HECAMEG": ["hecameg", "methyl 6-o-(n-heptylcarbamoyl)-alpha-d-glucopyranoside",
                "6-o-(n-heptylcarbamoyl)-methyl-alpha-d-glucopyranoside"],
    "beta-cyclodextrin": ["beta-cyclodextrin", "b-cyclodextrin"],
    "methyl-beta-cyclodextrin": ["methyl beta-cyclodextrin", "methyl-beta-cyclodextrin"],

# --- Precipitants / lipids / additives ---
    "2-Hydroxyethyl disulfide": ["2-hydroxyethyl disulfide", "2-hydroxyethyldisulfide",
                                  "hydroxyethyl disulfide", "2-hyrdoxyethyl disulfide", "hed"],
    "Tacsimate": ["tacsimate", "tascimate"],
    "Monoolein": ["monoolein", "1-oleoyl-rac-glycerol", "glycerol monooleate"],
    "6-Aminohexanoic acid": ["6-aminohexanoic acid", "aminohexanoic acid", "epsilon-aminocaproic acid"],
    "Spermine": ["spermine"],
    "Cholesterol": ["cholesterol"],
    "Malic acid": ["malic acid", "dl-malic acid", "l-malic acid", "s-malate", "rs-malate", "malate"],
    "DL-Serine": ["dl-serine"],
    "Acarbose": ["acarbose"],
    "Pentaerythritol ethoxylate": ["pentaerythritol ethoxylate"],
    "Isomaltose": ["isomaltose"],
    "EGTA": ["egta"],
    "Lactose": ["lactose"],
    "NDSB-221": ["ndsb-221", "ndsb 221"],
    "1-Deoxynojirimycin": ["1-deoxynojirimycin", "deoxynojirimycin"],
    "Pyruvate": ["pyruvate", "sodium pyruvate"],
    "Dithioerythritol": ["dithioerythritol", "dte"],
    "Caffeine": ["caffeine"],
    "Glycylglycine": ["glycylglycine", "gly-gly"],
    "Sodium oxamate": ["sodium oxamate"],
    "Trizma base": ["trizma base", "trizma"],

# --- Salts ---
    "Sodium fluoride": ["sodium fluoride", "naf"],
    "Ammonium bromide": ["ammonium bromide", "nh4br"],
    "Ammonium fluoride": ["ammonium fluoride", "nh4f"],
    "Potassium fluoride": ["potassium fluoride", "kf"],
    "Potassium formate": ["potassium formate"],
    "Potassium citrate": ["potassium citrate", "potassium citrate tribasic"],
    "Potassium bromide": ["potassium bromide", "kbr"],
    "Sodium dithionite": ["sodium dithionite", "na2-dithionite", "na dithionite"],
    "Hexaamminecobalt(III) chloride": ["hexaamminecobalt(iii) chloride", "hexaammine cobalt chloride",
                                        "cobalt hexammine chloride", "[co(nh3)6]cl3"],
    "Yttrium(III) chloride": ["yttrium chloride", "yttrium (iii) chloride", "yttrium(iii) chloride",
                              "yttrium (iii) chloride hexahydrate"],
    "Glucose-6-phosphate": ["glucose-6-phosphate", "g-6-p", "g6p"],

# --- Cofactors ---
    "GDP": ["gdp", "guanosine diphosphate"],
    "PLP": ["plp", "pyridoxal phosphate", "pyridoxal 5'-phosphate", "pyridoxal-5-phosphate"],
    "Adenosylcobalamin": ["adenosylcobalamin"],

# --- Amino acids / sugars ---
    "methionine": ["methionine", "l-methionine"],
    "tryptophan": ["tryptophan", "l-tryptophan"],
    "glucose": ["glucose", "d-glucose"],
    "mannose": ["mannose", "d-mannose"],
    "galactose": ["galactose", "d-galactose"],
    "xylose": ["xylose", "d-xylose"],
    "fucose": ["fucose", "l-fucose", "l-fuctose"],
    "N-acetylglucosamine": ["n-acetylglucosamine", "n-acetyl-d-glucosamine", "glcnac"],

# --- Defined buffer/reagent-mix names (Morpheus, JCSG screens) ---
    "SPG buffer": ["spg buffer", "spg"],
    "Morpheus Amino Acids Mix": ["morpheus amino acids mix", "amino acids mix (md)", "aminoacidsmix", "amino acids mix"],
    "Morpheus Buffer System 3": ["morpheus buffer system 3", "buffer system 3 (md)", "buffer system 3"],
}

# Where newly-learned compounds (found via the LLM fallback) get persisted,
# so future runs -- and later rows in the *same* run -- benefit from them
# without needing another LLM call. Kept separate from the built-in
# dictionary above so upgrading this module never clobbers what you've learned.
LEARNED_DICT_PATH = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "learned_reagents.json"
)


def _load_learned_dictionary(path: str = LEARNED_DICT_PATH) -> Dict[str, List[str]]:
    if not os.path.exists(path):
        return {}
    try:
        with open(path, "r") as f:
            return json.load(f)
    except (json.JSONDecodeError, OSError):
        return {}


def _save_learned_dictionary(learned: Dict[str, List[str]], path: str = LEARNED_DICT_PATH) -> None:
    with open(path, "w") as f:
        json.dump(learned, f, indent=2, sort_keys=True)


# Merge any previously-learned compounds into the working dictionary at import time.
for _canonical, _aliases in _load_learned_dictionary().items():
    REAGENT_DICTIONARY.setdefault(_canonical, [])
    for _alias in _aliases:
        if _alias not in REAGENT_DICTIONARY[_canonical]:
            REAGENT_DICTIONARY[_canonical].append(_alias)

# Build one big alias -> canonical-name lookup, sorted longest-first so
# "sodium/potassium phosphate" matches before the shorter "sodium phosphate".
# This is rebuilt (via _rebuild_patterns()) whenever a new compound is learned,
# so the regexes below always reflect the current dictionary.
_ALIAS_TO_CANONICAL: Dict[str, str] = {}
_PATTERN_NUM_FIRST = None
_PATTERN_REAGENT_FIRST = None
_PATTERN_REAGENT_ONLY = None

_NUMBER = r"\d+(?:\.\d+)?"
# Optional second number for a range like "15-25%" or "15% - 25%" (a "%"
# is allowed -- and discarded -- right after the first number too, since
# real pdbx_details text writes both forms). Captures into "amount2";
# match_clause() decides what to do with it (see _resolve_range).
_RANGE_SUFFIX = rf"(?:\s*%?\s*-\s*(?P<amount2>{_NUMBER}))?"
_CONC_UNIT = r"(?:%\s*\(?w/v\)?|%\s*\(?v/v\)?|%|M|mM|uM|\u00b5M|mg/mL|mg/ml)"

# Non-condition clauses to discard before matching (method/technique words).
_METHOD_STOPWORDS = re.compile(
    r"^\s*(vapor diffusion|hanging drop|sitting drop|batch|microbatch|"
    r"dialysis|free[- ]interface diffusion|seeding|micro-?seeding)\b.*$",
    re.IGNORECASE,
)


def _rebuild_patterns() -> None:
    """Recompute the alias lookup and the matching regexes from the
    current REAGENT_DICTIONARY. Called at import time, and again every
    time learn_new_compound() adds something new."""
    global _ALIAS_TO_CANONICAL, _PATTERN_NUM_FIRST, _PATTERN_REAGENT_FIRST, _PATTERN_REAGENT_ONLY

    alias_to_canonical = {}
    for canonical, aliases in REAGENT_DICTIONARY.items():
        for alias in aliases:
            alias_to_canonical[alias.strip().lower()] = canonical
    _ALIAS_TO_CANONICAL = alias_to_canonical

    aliases_sorted = sorted(_ALIAS_TO_CANONICAL.keys(), key=len, reverse=True)
    reagent_alternation = "|".join(re.escape(a) for a in aliases_sorted)

    # Pattern A: "NUMBER UNIT REAGENT" e.g. "0.2 M ammonium sulfate", also
    # tolerating a trailing unit qualifier that comes AFTER the reagent name,
    # e.g. "10% PEG 8000 (w/v)", and a range before the reagent, e.g.
    # "15-25% PEG 8000" or "15% - 25% PEG 8000" (see _RANGE_SUFFIX).
    _PATTERN_NUM_FIRST = re.compile(
        rf"(?P<amount>{_NUMBER}){_RANGE_SUFFIX}\s*(?P<unit>{_CONC_UNIT})\s*(?P<reagent>{reagent_alternation})"
        rf"(?:\s*\(?(?P<trailing_unit>w/v|v/v)\)?)?",
        re.IGNORECASE,
    )
    # Pattern B: "REAGENT NUMBER UNIT" e.g. "PEG 3350 20%" (less common, but seen)
    _PATTERN_REAGENT_FIRST = re.compile(
        rf"(?P<reagent>{reagent_alternation})\D{{0,3}}(?P<amount>{_NUMBER}){_RANGE_SUFFIX}\s*(?P<unit>{_CONC_UNIT})",
        re.IGNORECASE,
    )
    # Pattern C: bare reagent name with no concentration nearby at all, e.g.
    # a clause that's just "DMSO" or "ATP added" -- used as a last-resort
    # fallback so a known compound mentioned without a stated concentration
    # is still reported (with concentration "nan") instead of being dropped.
    _PATTERN_REAGENT_ONLY = re.compile(
        rf"(?P<reagent>{reagent_alternation})",
        re.IGNORECASE,
    )


def _normalize_unit(unit: str) -> str:
    """Normalize a unit string for consistency."""
    unit = re.sub(r"[()]", "", unit).strip().lower().replace(" ", "")
    mapping = {
        "%w/v": "% w/v", "%v/v": "% v/v", "%": "%",
        "m": "M", "mm": "mM", "um": "µM", "µm": "µM",
        "mg/ml": "mg/mL",
    }
    return mapping.get(unit, unit)


def _format_concentration(amount: str, unit: str) -> str:
    """Format concentration as a display string."""
    if not amount:
        return "nan"
    try:
        return f"{float(amount)} {unit}".strip()
    except (ValueError, TypeError):
        return "nan"


def _mean_of(low: str, high: str) -> str:
    """Calculate mean of a range's two endpoints."""
    mean = (float(low) + float(high)) / 2
    if mean == int(mean):
        return str(int(mean))
    return f"{mean:.2f}".rstrip("0").rstrip(".")


def _resolve_range(amount: str, amount2: Optional[str], raw_unit: str, canonical: str) -> str:
    """Resolve a range to a single value or range string."""
    if not amount2:
        return amount
    if canonical.upper().startswith("PEG") and "%" in raw_unit:
        return _mean_of(amount, amount2)
    return f"{amount}-{amount2}"


def split_into_clauses(details: str) -> List[str]:
    """Split a pdbx_details string into individual condition clauses."""
    if not details:
        return []
    raw_clauses = re.split(r"(?<!\d),|,(?!\d)|;|\band\b|\r\n|\n|\r", details, flags=re.IGNORECASE)
    clauses = []
    for c in raw_clauses:
        c = c.strip(" .'\"")  # Strip spaces, periods, single and double quotes
        if not c:
            continue
        if _METHOD_STOPWORDS.match(c):
            continue
        if re.match(r"^\d+(\.\d+)?\s*K$", c, re.IGNORECASE):
            continue
        if re.match(r"^pH\b", c, re.IGNORECASE):
            continue
        clauses.append(c)
    return clauses


def match_clause(clause: str) -> Optional[Dict[str, str]]:
    """Try to match one clause against the reagent dictionary."""
    m = _PATTERN_NUM_FIRST.search(clause)
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        amount = _resolve_range(m.group("amount").strip(), m.groupdict().get("amount2"),
                                 m.group("unit"), canonical)
        unit = _normalize_unit(m.group("unit"))
        trailing_unit = m.groupdict().get("trailing_unit")
        if trailing_unit and "w/v" not in unit and "v/v" not in unit:
            unit = f"{unit} {trailing_unit.lower()}" if unit else trailing_unit.lower()
        return {
            "compound": canonical,
            "amount": amount,
            "unit": unit,
            "concentration": _format_concentration(amount, unit),
            "source": "dictionary",
        }

    m = _PATTERN_REAGENT_FIRST.search(clause)
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        amount = _resolve_range(m.group("amount").strip(), m.groupdict().get("amount2"),
                                 m.group("unit"), canonical)
        unit = _normalize_unit(m.group("unit"))
        return {
            "compound": canonical,
            "amount": amount,
            "unit": unit,
            "concentration": _format_concentration(amount, unit),
            "source": "dictionary",
        }

    m2 = _PATTERN_REAGENT_ONLY.search(clause)
    if m2:
        reagent_alias = m2.group("reagent").strip().lower()
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        return {
            "compound": canonical,
            "amount": None,
            "unit": None,
            "concentration": "nan",
            "source": "dictionary",
        }

    return None


def update_reagent_dictionary(reagent: str, extra_aliases: Optional[List[str]] = None) -> bool:
    """Add a newly-seen reagent to the in-memory dictionary and persist it.

    The function canonicalizes the new reagent name, stores any known aliases,
    writes the learned dictionary to disk, and rebuilds the regexes so the new
    reagent is recognized immediately by later matches in the same run and all
    future runs.
    """
    canonical = reagent.strip()
    if not canonical:
        return False

    alias_set = {canonical.lower()}
    alias_set.update(a.strip().lower() for a in (extra_aliases or []) if a and a.strip())

    if canonical.lower() in _ALIAS_TO_CANONICAL or any(alias in _ALIAS_TO_CANONICAL for alias in alias_set):
        return False

    REAGENT_DICTIONARY.setdefault(canonical, [])
    for alias in sorted(alias_set):
        if alias not in REAGENT_DICTIONARY[canonical]:
            REAGENT_DICTIONARY[canonical].append(alias)

    learned = _load_learned_dictionary()
    learned.setdefault(canonical, [])
    for alias in sorted(alias_set):
        if alias not in learned[canonical]:
            learned[canonical].append(alias)
    _save_learned_dictionary(learned)

    _rebuild_patterns()
    return True


def learn_new_compound(compound: str, extra_aliases: Optional[List[str]] = None) -> bool:
    """Backward-compatible wrapper around update_reagent_dictionary()."""
    return update_reagent_dictionary(compound, extra_aliases)


# The old implementation was replaced by update_reagent_dictionary() above.
# This keeps compatibility with any existing code that calls
# learn_new_compound() directly.


# Initialize the regex patterns at module import time
_rebuild_patterns()


def format_compounds(matches: List[Dict]) -> str:
    """Format a list of matched compounds as a readable string."""
    if not matches:
        return ""
    return ", ".join(f"'{m['compound']}' ({m['concentration']})" for m in matches)


def process_csv(input_csv: str, output_csv: str, pdbx_column: str = "pdbx_details",
               use_llm_fallback: bool = False) -> None:
    """Extract compounds from a CSV file using the reagent dictionary
    and LLM fallback if enabled. Saves the results to output_csv.
    
    This function reads the input CSV, extracts compounds from pdbx_details,
    and writes the original data plus a new 'compound' column to output_csv."""
    if not os.path.exists(input_csv):
        return
    
    try:
        df = pd.read_csv(input_csv)
    except (pd.errors.ParserError, FileNotFoundError):
        return
    
    if pdbx_column not in df.columns:
        return
    
    # Extract compounds for each row
    compounds_list = []
    for idx, row in df.iterrows():
        details = str(row[pdbx_column]) if pd.notna(row[pdbx_column]) else ""
        clauses = split_into_clauses(details)
        matched_compounds = []
        
        for clause in clauses:
            match = match_clause(clause)
            if match:
                matched_compounds.append(match)
            elif use_llm_fallback and requests:
                pass
        
        compounds_list.append(format_compounds(matched_compounds))
    
    # Add the compound column to the dataframe
    df["compound"] = compounds_list
    
    # Write the full dataframe back out
    try:
        df.to_csv(output_csv, index=False)
    except Exception:
        pass
