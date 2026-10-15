import os
import re
import json
from typing import List, Dict, Optional, Tuple

import pandas as pd


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
    "PEG MME 3350": ["peg mme 3350", "peg 3350 mme", "peg3350 mme", "peg 3350mme", "peg-mme3350",
                     "peg monomethyl ether 3350"],
    "PEG 3350": ["peg 3350", "peg3350", "peg-3350", "polyethylene glycol 3350", "poly(ethylene glycol) 3350",
                 "peg 3,350", "peg3,350", "peg3330", "peg3340", "peg 3380", "peg3380", "polyethylene glycol 3,350",
                 "pg3350"],
    "PEG 4000": ["peg 4000", "peg4000", "peg-4000", "peg4k", "peg 4k", "polyethylene glycol 4000", "poly(ethylene glycol) 4000",
                 "peg 4,000", "pegf-4000", "polyethylene glycol 4,000"],
    "PEG 6000": ["peg 6000", "peg6000", "peg-6000", "peg6k", "peg 6k", "polyethylene glycol 6000", "poly(ethylene glycol) 6000",
                 "peg 6,000"],
    "PEG 8000": ["peg 8000", "peg8000", "peg-8000", "peg8k", "peg 8k", "polyethylene glycol 8000", "poly(ethylene glycol) 8000",
                 "peg 8,000", "polyethyleneglycol 8000"],
    "PEG 10000": ["peg 10000", "peg10000", "peg-10000", "peg10k", "peg 10k", "polyethylene glycol 10000",
                  "peg 10,000", "polyethylene glycol 10k"],
    "PEG 5000": ["peg 5000", "peg5000", "peg-5000", "peg5k", "peg 5k", "polyethylene glycol 5000",
                 "peg 5,000"],
    "PEG MME 500": ["peg mme 500", "peg monomethyl ether 500", "peg-500 mme", "peg500 mme",
                    "peg 500mme", "peg 500 mme", "peg500mme", "polyethylene glycol monomethyl ether 500"],
    "PEG MME 550": ["peg mme 550", "peg monomethyl ether 550", "peg-550", "peg 550", "peg550 mme",
                    "peg 550mme", "peg 550 mme", "polyethylene glycol monomethyl ether 550",
                    "peg-mme550", "peg-mme 550", "mpeg 550", "mpeg550", "peg550*mme"],
    "PEG MME 2000": ["peg mme 2000", "peg mme 2k", "peg mme2k", "peg monomethyl ether 2000",
                     "peg mme2000", "peg 2k mme", "peg 2kmme", "peg 2000 mme", "peg2000 mme",
                     "p2000mme", "peg2k-mme", "peg 2k-mme",
                     "polyethylene glycol monomethyl ether 2000",
                     "polyethylene glycol monomethyl ether 2,000",
                     "polyethylene glycol monomethyl ether (peg) 2000"],
    "PEG MME 5000": ["peg mme 5000", "peg 5000 mme", "peg monomethyl ether 5000", "peg 5000 monomethyl-ether",
                      "peg 5000 monomethyl ether", "polyethylene glycol monomethyl ether 5000",
                      "polyethylene glycol monomethyl ether 5,000", "monomethyl ether peg 5000",
                      "methyl ether poly(ethylene glycol) 5000", "methyl ether polyethylene glycol 5000"],
    "PEG 3500": ["peg 3500", "peg3500", "peg 3550", "peg3550"],
    "PEG Smear High": ["peg smear high"],
    "PEG Smear Medium": ["peg smear medium"],
    "PEG 200": ["peg 200", "peg200", "polyethylene glycol 200", "poly(ethylene glycol) 200", "peg-200"],
    "PEG 300": ["peg 300", "peg300", "polyethylene glycol 300", "poly(ethylene glycol) 300", "peg-300"],
    "PEG 1500": ["peg 1500", "peg1500", "polyethylene glycol 1500", "poly(ethylene glycol) 1500", "peg-1500",
                 "peg 1,500"],
    "PEG 1450": ["peg 1450", "peg1450", "polyethylene glycol 1450", "poly(ethylene glycol) 1450", "peg-1450"],
    "PEG 1550": ["peg 1550", "peg1550", "polyethylene glycol 1550", "poly(ethylene glycol) 1550", "peg-1550"],
    "PEG 2000": ["peg 2000", "peg2000", "peg-2000", "polyethylene glycol 2000", "poly(ethylene glycol) 2000",
                 "peg 2,000", "peg 2k", "peg2k", "peg-2k", "mpeg 2k", "mpeg2k"],
    "PEG 600": ["peg 600", "peg600", "polyethylene glycol 600", "poly(ethylene glycol) 600", "peg-600"],
    "PEG 800": ["peg 800", "peg800", "polyethylene glycol 800", "poly(ethylene glycol) 800", "peg-800"],
    "PEG MME 750": ["peg mme 750", "peg monomethyl ether 750", "peg-750 mme", "peg750 mme",
                    "peg 750mme", "peg 750 mme", "polyethylene glycol monomethyl ether 750"],
    "PEG dimethyl ether": ["peg dimethyl ether", "polyethylene glycol dimethyl ether"],
    "PEG 3000": ["peg 3000", "peg3000", "peg-3000", "polyethylene glycol 3000", "peg 3,000", "peg3k", "peg 3k"],
    "PEG 12000": ["peg 12000", "peg12000", "peg-12000", "peg12k", "peg 12k", "polyethylene glycol 12000", "peg 12,000"],
    "PEG 20000": ["peg 20000", "peg20000", "peg-20000", "peg20k", "peg 20k", "polyethylene glycol 20000", "peg 20,000"],
    "PEG": ["peg", "polyethylene glycol", "poly(ethylene glycol)", "purepegs", "purepegs cocktail",
            "dimethyl polyethylene glycol", "pegsb", "polyethylene glycole"],
    "PEG MME": ["peg mme", "polyethylene glycol monomethyl ether", "peg monomethyl ether"],
    "Jeffamine ED-2001": ["jeffamine ed-2001", "ed2001"],
    "Jeffamine M-600": ["jeffamine m-600", "jeffamine m600"],

    # --- Salts ---
    "ammonium sulfate": ["ammonium sulfate", "ammonium sulphate", "(nh4)2so4", "ammonium-sulfate",
                          "amm.sulfate", "ammonia sulfate", "nh4-sulfate", "ammomium sulfate",
                          "amimonium sulfate",
                          "ammonimu sulphate", "ammonium slufate",
                          "ammonium sulf ate", "ammonium sulfphate", "(nh4)so4", "(nh4)2 so4", "nh4so4",
                          "amso4", "as", "ams"],
    "sodium chloride": ["sodium chloride", "nacl", "sodium choloride"],
    "lithium sulfate": ["lithium sulfate", "lithium sulphate", "li2so4", "liso4"],
    "lithium chloride": ["lithium chloride", "licl"],
    "magnesium chloride": ["magnesium chloride", "mgcl2", "mg2cl", "mgcl 2", "mgcl"],
    "magnesium sulfate": ["magnesium sulfate", "magnesium sulphate", "mgso4"],
    "calcium chloride": ["calcium chloride", "cacl2"],
    "sodium phosphate": ["sodium phosphate", "nah2po4", "na2hpo4", "sodium phospahte", "na2po4", "napi",
                         "na phosphate buffer", "na phosphate"],
    "potassium phosphate": ["potassium phosphate", "kh2po4", "k2hpo4", "potassium dihydrogen phosphate", "potassium dihydrogenphosphate",
                             "kh2po4/k2hpo4", "k2hpo4/kh2po4", "potassium dihydrogen orthophosphate"],
    "sodium/potassium phosphate": ["sodium potassium phosphate", "sodium/potassium phosphate",
                                    "sodium/ potassium phosphate", "na/k phosphate", "napo4",
                                    "pottasium phosphate", "nakpo4", "na-k phosphate", "k/na phosphate",
                                    "k2hpo4/nah2po4", "k2hpo4 /nah2po4", "k2hpo4/ nah2po4", "nah2po4/k2hpo4",
                                    "kh2po4/na2hpo4", "na2hpo4/kh2po4", "na k phosphate",
                                    "k/na-pi", "k/na pi", "k-na-pi", "na+ / k+ phosphate buffer",
                                    "na+/k+ phosphate buffer", "na+ / k+ phosphate", "na+/k+ phosphate",
                                    "na/k - phosphate", "na/k-phosphate"],
    "sodium pyrophosphate": ["sodium pyrophosphate", "tetrasodium pyrophosphate", "na4p2o7", "sodium diphosphate"],
    "sodium acetate": ["sodium acetate", "sodium actate", "naoac", "naac", "naacetate", "na-acetate",
                       "na acetate", "acna", "na(ch3coo)"],
    "copper acetate": ["copper acetate", "cupric acetate", "cu(oac)2", "ac2cu"],
    "magnesium citrate": ["magnesium citrate", "mg citrate"],
    "sodium citrate": ["sodium citrate", "tri-sodium citrate", "NAcitrate", "trisodium citrate", "na citrate", "Na-citrate", "na3citrate", "na3cit", "na3-cit", "nacit", "citrate-na", "citric sodium", "tri-na citrate"],
    "sodium formate": ["sodium formate", "NAformate", "na formate", "Na-formate", "formic acid sodium salt", ],
    "potassium chloride": ["potassium chloride", "kcl", "potassium cloride", "potasium chloride"],
    "zinc acetate": ["zinc acetate", "zn(oac)2", "zinc acetate dihydrate", "zn(ch2h3o2)2", "zn acet", "zn acetate",
                     "znoac", "zn(ac)2"],
    "zinc chloride": ["zinc chloride", "Zinc Chloride", "zncl2"],
    "potassium acetate": ["potassium acetate"],
    "methyl mercuric acetate": ["methyl mercuric acetate", "methylmercuric acetate", "methyl mercury acetate"],
    "acetate": ["acetate"],
    "ammonium chloride": ["ammonium chloride", "nh4cl"],
    "ammonium nitrate": ["ammonium nitrate","nh4no3"],
    "ammonium acetate": ["ammonium acetate","nh4 acetate"],
    "ammonium citrate": ["ammonium citrate", "citrate-ammonium", "citrate ammonium"],
    "ammonium phosphate": ["ammonium phosphate","(nh4)3po4", "ammonium dihydrogen phosphate", "nh4h2po4",
                           "(nh4)2hpo4"],
    "glyceryl tributyrate": ["glyceryl tributyrate", "tributyrin"],
    "paraoxon-ethyl": ["paraoxon-ethyl", "paraoxon ethyl", "paraoxon"],
    "ammonium bicarbonate": ["ammonium bicarbonate","nh4hco3"],
    "potassium sulfate": ["potassium sulfate","potassium sulphate","k2so4"],
    "Monosodium phosphate": ["monosodium phosphate", "sodium phosphate monobasic", "nah2po4", "na h2po4", "sodium dihydrogen phosphate", "sodium dihydrogenphosphate"],
    "NaH2PO4": ["nah2po4", "na h2 po4", "monosodium phosphate", "sodium phosphate monobasic"],
    "Disodium phosphate": ["disodium phosphate", "sodium phosphate dibasic", "na2hpo4", "na2 hpo4"],
    "Dipotassium hydrogen phosphate": ["dipotassium hydrogen phosphate", "potassium phosphate dibasic", "k2hpo4", "k2 hpo4", "potassium hydrogen phosphate", "potassium hydrogenphosphate", "dipotassium phosphate"],
    "Potassium phosphate": ["potassium phosphate","kh2po4","k2hpo4"],
    "sodium sulfate": ["sodium sulfate","sodium sulphate","na2so4","naso4"],
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
    "ammonium tartrate": ["ammonium tartrate", "diammonium tartrate", "ammonium tartarate",
                          "ammonium tartrate dibasic", "ammonium tartarate dibasic"],
    "potassium sodium tartrate": ["potassium sodium tartrate", "sodium potassium tartrate",
                                   "rochelle salt", "kna tartrate", "na/k tartrate", "k/na tartrate",
                                   "na-k tartrate", "sodium potassium tartarate",
                                   "potassium sodium tartarate"],
    "potassium tartrate": ["potassium tartrate", "dipotassium tartrate"],
    "sodium bromide": ["sodium bromide", "nabr"],
    "potassium iodide": ["potassium iodide", "ki"],
    "manganese acetate": ["manganese acetate", "mn(oac)2", "mn acetate"],
    "copper chloride": ["copper chloride", "cucl2", "cupric chloride"],
    "Lithium citrate": ["lithium citrate", "lithium citrate tribasic", "li3 citrate"],
    "Lithium acetate": ["lithium acetate", "li acetate", "liac"],
    "Sodium malonate": ["sodium malonate", "malonic acid sodium salt", "na malonate", "na-malonic acid"],
    "potassium succinate": ["potassium succinate", "dipotassium succinate"],
    "Sodium maleate": ["sodium maleate","maleic acid sodium salt", "disodium maleate"],
    "Sodium tartrate": ["sodium tartrate", "disodium tartrate", "sodium L-tartrate"],
    "Di-sodium tartrate": ["di-sodium tartrate", "disodium tartrate","sodium tartrate"],
    "Magnesium acetate": ["magnesium acetate", "mg acetate", "mg(oac)2", "mgaco", "mgac", "mg-ac"],
    "Magnesium formate": ["magnesium formate", "mg formate", "mg(hcoo)2"],
    "Cobalt(II) chloride": ["cobalt chloride", "cobalt(ii) chloride", "cobalt (ii) chloride", "cocl2", "cobaltous chloride"],
    "Cobalt(III) chloride": ["cobalt(iii) chloride", "cobalt (iii) chloride", "cocl3"],
    "Cadmium chloride": ["cadmium chloride", "cadmium(ii) chloride", "cdcl2"],
    "Nickel(II) chloride": ["nickel chloride", "nickel(ii) chloride", "nickel (ii) chloride", "nicl2"],
    "Sodium iodide": ["sodium iodide", "nai","na iodide", "NAI"],
    "Sodium succinate": ["sodium succinate","sodium succinate buffer","disodium succinate","succinic acid disodium salt","succinate","na succinate","sodium hydrogen succinate"],

    # --- Buffers ---
    "Tris-HCl": ["tris-hcl", "tris hcl", "trishcl", "triscl", "tris buffer", "tris", "trishydroxymethylaminomethane"],
    "tris-sulfate": ["tris-sulfate", "tris sulfate", "tris sulphate"],
    "Tris-acetate": ["tris-acetate", "tris acetate", "tris/acetate"],
    "Tris-maleate": ["tris-maleate", "tris maleate"],
    "Bis-Tris": ["bis-tris", "bis tris", "bistris", "bis_tris", "bis tris buffer","bis-tris buffer",
                 "bis/tris"],
    "HEPES": ["hepes", "hepes buffer", "4-(2-hydroxyethyl)-1-piperazineethanesulfonic acid",
              "na-hepes", "nahepes", "sodium hepes"],
    "MES": ["mes","mes buffer","2-(n-morpholino)ethanesulfonic acid", "2 morpholinoethanesulfonic acid",
            "morpholine-ethane-sulphonic acid", "morpholine ethanesulfonic acid",
            "morpholinoethane sulphonic acid", "morpholinoethane sulfonic acid",
            "2-(n-morpholino)- ethanesulfonic acid", "2-(n-morpholino)-ethanesulfonic acid",
            "na-mes", "na mes"],
    "bis-tris propane": ["bis-tris propane", "bis_tris propane", "bis tris propane", "bistris propane", "btp"],
    "Bis-Tris chloride": ["bis-tris chloride", "bis tris chloride","bis-tris hcl"],
    "CHES": ["ches", "ches buffer"],
    "CAPS": ["caps"],
    "imidazole": ["imidazole"],
    "triethanolamine": ["triethanolamine", "triethanolamine buffer", "tea buffer", "tea-buffer"],
    "citrate buffer": ["citrate buffer", "citric acid", "citrate"],
    "cacodylate": ["sodium cacodylate", "cacodylate"],
    "MOPS": ["mops"],
    "glycine": ["glycine"],
    "acetate buffer": ["acetate buffer", "acetic acid"],
    "Tricine": ["tricine", "n-[tris(hydroxymethyl)methyl]glycine"],
    "TES": ["tes", "tes buffer", "ammonium tes",
            "n-[tris(hydroxymethyl)methyl]-2-aminoethanesulfonic acid",
            "ammonium n-[tris(hydroxymethyl)methyl]-2-aminoethanesulfonate",
            "ammonium n-[tris(hydroxymethyl) methyl]-2-aminoethane sulfonate",
            "ammonium n-[tris(hydroxymethyl)methyl]-2-aminoethane sulfonate",
            "ammonium n-[tris(hydroxymethyl)methyl]-2 aminoethanesulfonate"],
    "PIPES": ["pipes", "piperazine-n,n'-bis(2-ethanesulfonic acid)", "k-pipes", "k pipes"],
    "POPSO": ["popso","piperazine-n,n'-bis(2-hydroxypropanesulfonic acid)"],
    "TAPS": ["taps"],
    "TAPSO": ["tapso"],
    "ADA": ["ada", "ada buffer", "n-(2-acetamido)iminodiacetic acid"],
    "Bicine": ["bicine", "n,n-bis(2-hydroxyethyl)glycine"],
    "CAPSO": ["capso"],
    "borate": ["borate", "sodium borate", "boric acid"],
    "phosphate buffer": ["phosphate buffer","phosphate"],
    "Phosphate-citrate buffer": ["phosphate citrate", "phosphate-citrate", "citrate phosphate",
                                  "citrate-phosphate", "phosphate-citrate buffer",
                                  "na2hpo4/citric acid", "na2hpo4-citric acid",
                                  "citrate-po4 buffer", "citrate-po4", "citrate po4", "po4-citrate"],
    "sodium cacodylate": ["sodium cacodylate", "cacodylate", "dimethylarsinic acid sodium salt",
                          "na-cacodylate", "na cacodylate", "nacacodylate", "na caodylate", "cacodylic acid"],
    "Trisodium citrate-citric acid": ["trisodium citrate-citric acid", "trisodium citrate citric acid", "trisodium citrate/citric acid", "sodium citrate citric acid buffer", "citric acid sodium citrate", "sodium citrate: citric acid"],
    "sulfuric acid": ["sulfuric acid", "sulphuric acid", "h2so4"],
    "acetohydroxamic acid": ["acetohydroxamic acid", "droxamic acid"],
    "toluene": ["toluene"],
    "sodium sulfide": ["sodium sulfide", "sodium sulphide", "na2s"],
    "cadaverine dihydrochloride": ["cadaverine dihydrochloride", "1,5-diaminopentane dihydrochloride",
                                    "1,5-pentanediamine dihydrochloride"],
    "3-phosphoglycerate": ["3-phosphoglycerate", "3pg", "3-phosphoglyceric acid"],
    "aluminum chloride": ["aluminum chloride", "aluminium chloride", "alcl3"],
    "SB12": ["sb12", "sb-12", "sulfobetaine 12", "n-dodecyl-n,n-dimethyl-3-ammonio-1-propanesulfonate"],
    "IHP": ["ihp", "inositol hexaphosphate", "phytic acid"],
    "sodium ascorbate": ["sodium ascorbate", "na-ascorbate", "na ascorbate"],
    "sodium monothiophosphate": ["sodium monothiophosphate", "monothiophosphate"],
    "adenosine": ["adenosine"],
    "Mavacamten": ["mavacamten"],
    "resorcinol": ["resorcinol"],
    "beryllium fluoride": ["beryllium fluoride", "bef2"],
    "sodium sulfite": ["sodium sulfite", "sodium sulphite", "na2so3"],
    "Gadolinium chloride hexahydrate": ["gadolinium(iii) chloride hexahydrate", "gadolinium chloride hexahydrate",
                                          "gdcl3"],
    "Gadolinium nitrate": ["gadolinium nitrate", "gd(no3)3"],
    "silver nitrate": ["silver nitrate", "agno3"],
    "phenylphosphorodiamidate": ["phenylphosphorodiamidate", "phenyl phosphorodiamidate"],
    "Deoxy-BigCHAP": ["deoxy-bigchap", "deoxy bigchap", "deoxy big chap", "deoxy-big chap"],
    "potassium ferricyanide": ["potassium ferricyanide", "potassium ferricyanate", "k3fe(cn)6"],
    "MMT buffer": ["mmt buffer", "mmt", "malic acid mes tris buffer"],
    "ferrous citrate": ["ferrous citrate", "iron(ii) citrate"],
    "Tris-buffered saline": ["tris-buffered saline", "tris buffered saline", "tbs"],
    "fluorinated octylmaltoside": ["fluorinated octylmaltoside", "fluorinated octyl maltoside"],
    "Ammonium iodide": ["ammonium iodide", "nh4i"],
    "PCB buffer": ["pcb buffer", "propionate-cacodylate-bis-tris propane buffer",
                   "propionate cacodylate bis-tris propane buffer"],
    "PEG DME 500": ["peg500 dme", "peg 500 dme", "peg dme 500", "peg-dme500"],
    "NBPTO": ["nbpto", "n-(n-butyl)-thiophosphoric triamide oxide"],
    "deferoxamine": ["deferoxamine", "desferrioxamine", "desferal"],

    # --- Organic precipitants ---
    "1,6-hexanediol": ["1,6-hexanediol", "hexanediol"],
    "2,5-hexanediol": ["2,5-hexanediol", "2,5 hexanediol"],
    "1,4-butanediol": ["1,4-butanediol", "butanediol" ],
    "tert-butanol": ["tert-butanol", "t-butanol", "tbutanol"],
    "methanol": ["methanol", "methyl alcohol", "meoh"],
    "ethanol": ["ethanol", "ethyl alcohol", "etoh"],
    "propanol": ["propanol","1-propanol", "iso-propanol"],
    "butanol": ["butanol","1-butanol"],
    "acetone": ["acetone"],
    "acetonitrile": ["acetonitrile"],
    "Sorbitol": ["sorbitol", "d-sorbitol", "D-sorbitol", "glucitol", "d-glucitol", "hexane-1,2,3,4,5,6-hexol", "sorbitol solution"],
    "Trehalose": ["trehalose", "d-trehalose", "trehalose dihydrate"],
    "myo-Inositol": ["myo-inositol", "myo inositol", "inositol"],
    "Xylitol": ["xylitol"],
    "trimethylamine N-oxide": ["trimethylamine n-oxide", "tmao", "trimethylamine oxide"],
    "trimethylamine hydrochloride": ["trimethylamine hydrochloride", "trimethylamine-hcl", "trimethylamine hcl"],
    "sarcosine": ["sarcosine", "n-methylglycine"],
    "1,3-propanediol": ["1,3-propanediol", "1,3 propanediol"],
    "1,2-Propanediol": ["1,2-propanediol", "1,2 propanediol", "propylene glycol", "propan-1,2-diol"],
    "1,5-pentanediol": ["1,5-pentanediol", "1,5 pentanediol"],

    # --- Additives / cryoprotectants ---
    "glycerol": ["glycerol", "anhydrous glycerol"],
    "DMSO": ["dmso", "dimethyl sulfoxide"],
    "ATP": ["atp", "adenosine triphosphate"],
    "ATPgammaS": ["atpgammas", "atp-gamma-s", "atp gamma s", "adenosine 5'-o-(3-thiotriphosphate)"],
    "MPD": ["mpd", "2-methyl-2,4-pentanediol", "hexylene glycol", "methylpentanediol"],
    "isopropanol": ["isopropanol", "2-propanol", "isopropyl alcohol", "propan-2-ol"],
    "ethylene glycol": ["ethylene glycol", "ethylen glycole", "eg"],
    "triethylene glycol": ["triethylene glycol", "triethyleneglycol", "tri-ethylene glycol"],
    "DTT": ["dtt", "dithiothreitol"],
    "DDAO": ["ddao", "n,n-dimethyldodecylamine oxide", "n,n-dimethyldodecylamine oxide(ddao)",
             "lauryldimethylamine oxide"],
    "m-cresol": ["m-cresol", "metacresol", "meta-cresol"],
    "omecamtiv mecarbil": ["omecamtiv mecarbil"],
    "vanadate": ["vanadate", "sodium vanadate", "orthovanadate", "sodium orthovanadate", "na3vo4"],
    "pentafluorobenzyl alcohol": ["pentafluorobenzyl alcohol", "2,3,4,5,6-pentafluorobenzyl alcohol"],
    "4-iodobenzo[b]thiophene-2-carboxamidine": ["4-iodobenzo[b]thiophene-2-carboxamidine"],
    "beta-mercaptoethanol": ["beta-mercaptoethanol", "2-mercaptoethanol", "bme"],
    "TCEP": ["tcep", "tecp"],
    "CTAB": ["ctab", "cetyltrimethylammonium bromide", "hexadecyltrimethylammonium bromide"],
    "sodium azide": ["sodium azide", "nan3", "na-azide", "na azide"],
    "phenol": ["phenol"],
    "EDTA": ["edta"],
    "glutathione": ["glutathione", "gsh"],
    "pramocaine": ["pramocaine", "pramociane", "pramoxine"],
    "camphor": ["camphor", "d-camphor", "(1r)-camphor"],
    "spermidine": ["spermidine"],
    "jeffamine": ["jeffamine"],
    "PPG": ["ppg", "polypropylene glycol"],
    "dioxane": ["dioxane", "dioxan", "1,4-dioxane", "1,4 dioxane"],
    "Bromophenol blue": ["bromophenol blue", "bromo phenyl blue", "bpb"],
    "zinc": ["zinc", "zn(ii)", "zn2+", "zn 2+"],
    "magnesium": ["magnesium", "mg(ii)", "mg2+", "mg 2+"],

# --- Detergents / surfactants ---
    "n-octyl-beta-D-glucoside": [ "octyl glucoside", "og", "beta-octyl glucoside", "beta-octyl-glucoside",
                                  "beta-octylglucoside", "bog"],
    "octyl-POE": ["octyl-poe", "octylpoe", "octyl poe", "octylpolyoxyethylene"],
    "decylmaltoside": ["decylmaltoside", "decyl maltoside", "n-decyl-beta-d-maltoside", "decyl-beta-d-maltoside"],
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
    "calcium acetate": ["calcium acetate", "ca-acetate", "ca acetate"],
    "manganese sulfate": [ "manganese sulfate", "mnso4"],
    "nickel sulfate": ["nickel sulfate","niso4"],
    "cobalt sulfate": ["cobalt sulfate","coso4"],
    "cadmium sulfate": ["cadmium sulfate","cdso4"],
    "iron sulfate": ["ferrous sulfate", "iron sulfate", "feso4"],
    "ATP magnesium salt": ["mg-atp", "magnesium atp", "atp-mg", "mg atp"],
    "ADP": ["adp","adenosine diphosphate"],
    "MgADP": ["mgadp", "mg-adp", "mg.adp", "magnesium adp"],
    "GTP": ["gtp","guanosine triphosphate"],
    "FAD": [ "fad", "flavin adenine dinucleotide"],
    "NAD": ["nad","nicotinamide adenine dinucleotide"],
    "NADP": ["nadp", "nicotinamide adenine dinucleotide phosphate"],
    "coenzyme A": ["coa","coenzyme a","co-enzyme a"],
    "AMP": ["amp","adenosine monophosphate"],
    "AMP-PNP": ["amp-pnp", "amp pnp", "ampnp", "amppnp",
                "5'-adenylyl-beta,gamma-imidodiphosphate", "adenosine 5'-(beta,gamma-imido)triphosphate"],
    
# --- Common Morpheus crystallization additives ---

    "1,2,3-heptanetriol": ["1,2,3-heptanetriol", "heptanetriol"],
    "1,2,4-butanetriol": ["1,2,4-butanetriol", "butanetriol"],
    "1,2,6-hexanetriol": ["1,2,6-hexanetriol", "hexanetriol"],
    "2,3-butanediol": ["2,3-butanediol", "butanediol"],

# --- Hampton additive screen compounds ---
    "benzamidine hydrochloride": ["benzamidine hcl", "benzamidine-hcl", "benzamidine hydrochloride"],
    "guanidine hydrochloride": ["guanidine hcl", "guanidinium chloride", "gndhcl", "gnd-hcl", "gnd hcl"],
    "Guanidine sulfate": ["guanidine sulfate", "guanidinium sulfate"],
    "urea": ["urea"],
    "thiourea": ["thiourea"],
    "acetamide": ["acetamide"],
    "formamide": ["formamide"],
    "sodium thiocyanate": ["sodium thiocyanate","nascn"],
    "potassium thiocyanate": ["potassium thiocyanate", "kscn"],
    "ammonium thiocyanate": ["ammonium thiocyanate","nh4scn"],

# --- Polymers ---
    "poly(acrylic acid sodium salt) 5100": ["poly(acrylic acid sodium salt) 5100","sodium polyacrylate 5100", "poly acrylic acid sodium salt", "paa sodium salt 5100", "polyacrylate"],
    "PGA-LM": ["pga-lm", "pga lm", "pga_lm", "poly-gamma-glutamic acid low molecular weight", "poly gamma glutamic acid lm", "polyglutamic acid lm", "gamma-polyglutamic acid","poly-gamma-glutamate"],
    "Poly(acrylic acid sodium salt) 5100": ["poly(acrylicacid sodium salt) 5100", "poly(acrylic acid sodium salt) 5100", "poly acrylic acid sodium salt 5100", "polyacrylic acid sodium salt 5100",
    "sodium polyacrylate 5100", "poly(sodium acrylate) 5100", "paa sodium salt 5100", "paaa 5100", "polyacrylate sodium salt"],
    "Sucrose": ["sucrose", "d-sucrose", "\u03b1-d-glucopyranosyl-(1\u21922)-\u03b2-d-fructofuranoside", "cane sugar", "table sugar", "saccharose", "sucrose solution"],
    "Polyvinylpyrrolidone K15": ["polyvinylpyrrolidone k15", "pvp k15", "pvp-k15", "polyvinylpyrrolidone", "pvp"],
    "Dextran sulfate": ["dextran sulfate", "dextran sulphate"],
    
# --- PLP Snyder fluorescent dyes ---
    "PLP_Snyder441": ["plp_snyder441", "plp snyder 441","snyder 441"],
    "PLP_Snyder457": ["plp_snyder457", "plp snyder 457", "snyder 457"],
    "PLP_Snyder494": ["plp_snyder494","plp snyder 494","snyder 494"],
    "PLP_Snyder495": ["plp_snyder495","plp snyder 495","snyder 495"],
    "PLP_Snyder496": ["plp_snyder496","plp snyder 496","snyder 496"],
    "PLP_Snyder530": ["plp_snyder530","plp snyder 530","snyder 530"],
    "PLP_Snyder608": ["plp_snyder608","plp snyder 608","snyder 608"],
    "PLP_Snyder630": ["plp_snyder630","plp snyder 630","snyder 630"],

# --- Ligands / additives ---

"Ligand": ["small molecule ligand"],
"Polyamine": ["polyamine", "polyamines"],
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
    "Spermine tetrahydrochloride": ["spermine tetrahydrochloride", "spermine 4hcl", "spermine 4 hcl"],
    "Cholesterol": ["cholesterol"],
    "Malic acid": ["malic acid", "dl-malic acid", "l-malic acid", "s-malate", "rs-malate", "malate"],
    "DL-Serine": ["dl-serine"],
    "Acarbose": ["acarbose"],
    "Pentaerythritol ethoxylate": ["pentaerythritol ethoxylate"],
    "Pentaerythritol ethoxylate 3350": ["pentaerythritol ethoxylate 3350", "pe3350", "pe 3350", "pe-3350"],
    "Isomaltose": ["isomaltose"],
    "EGTA": ["egta"],
    "Lactose": ["lactose"],
    "NDSB-221": ["ndsb-221", "ndsb 221"],
    "NDSB-195": ["ndsb-195", "ndsb 195"],
    "NDSB-256": ["ndsb-256", "ndsb 256"],
    "1-Deoxynojirimycin": ["1-deoxynojirimycin", "deoxynojirimycin"],
    "Pyruvate": ["pyruvate", "sodium pyruvate"],
    "Dithioerythritol": ["dithioerythritol", "dte"],
    "Caffeine": ["caffeine"],
    "Glycylglycine": ["glycylglycine", "gly-gly", "glycyl-glycine", "glycyl glycine"],
    "Sodium oxamate": ["sodium oxamate"],
    "Trizma base": ["trizma base", "trizma"],

# --- Salts ---
    "Sodium fluoride": ["sodium fluoride", "naf"],
    "Ammonium bromide": ["ammonium bromide", "nh4br"],
    "Ammonium fluoride": ["ammonium fluoride", "nh4f"],
    "Potassium fluoride": ["potassium fluoride", "kf"],
    "Potassium formate": ["potassium formate", "k formate", "k-formate"],
    "Potassium citrate": ["potassium citrate", "potassium citrate tribasic"],
    "Cesium acetate": ["cesium acetate", "csoac"],
    "Potassium oxalate": ["potassium oxalate", "k2c2o4"],
    "Ammonium oxalate": ["ammonium oxalate", "(nh4)2c2o4"],
    "Potassium bromide": ["potassium bromide", "kbr"],
    "Sodium dithionite": ["sodium dithionite", "na2-dithionite", "na dithionite", "ditionite"],
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
    "Morpheus Amino Acids Mix": ["morpheus amino acids mix", "amino acids mix (md)", "aminoacidsmix",
                                  "amino acids mix", "amino acids"],
    "Morpheus Buffer System 3": ["morpheus buffer system 3", "buffer system 3 (md)", "buffer system 3"],
}

# --------------------------------------------------------------------------- #
# Reagent category -- Buffer / Salt / Precipitant / Additive -- used to group
# the reagent-frequency summary the web app shows (see viewer/data.py's
# build_reagent_frequency()) the way a crystallographer actually assembles a
# cocktail: one buffer + one salt + one precipitant, plus optional additives.
# Generated from this file's own "# --- section ---" comments above, with a
# short list of manual corrections for sections that mix categories (e.g.
# sugars/osmolytes filed under "Organic precipitants" are really additives).
# --------------------------------------------------------------------------- #
REAGENT_CATEGORY = {
    "PEG MME 750": "Precipitant",
    "β-octyl glucoside": "Additive",
    '1,2,3-heptanetriol': 'Precipitant',
    '1,2,4-butanetriol': 'Precipitant',
    '1,2,6-hexanetriol': 'Precipitant',
    '1,2-Propanediol': 'Precipitant',
    '1,3-propanediol': 'Precipitant',
    '1,4-butanediol': 'Precipitant',
    '1,5-pentanediol': 'Precipitant',
    '1,6-hexanediol': 'Precipitant',
    '2,5-hexanediol': 'Precipitant',
    '1-Deoxynojirimycin': 'Additive',
    '2,3-butanediol': 'Precipitant',
    '2-Hydroxyethyl disulfide': 'Additive',
    '2-propanol': 'Precipitant',
    '6-Aminohexanoic acid': 'Additive',
    'ADA': 'Buffer',
    'ADP': 'Additive',
    'MgADP': 'Additive',
    'AMP': 'Additive',
    'AMP-PNP': 'Additive',
    'ATP': 'Additive',
    'ATPgammaS': 'Additive',
    'ATP magnesium salt': 'Additive',
    'Acarbose': 'Additive',
    'Adenosylcobalamin': 'Additive',
    'Ammonium bromide': 'Salt',
    'Ammonium fluoride': 'Salt',
    'Ammonium oxalate': 'Salt',
    'Bicine': 'Buffer',
    'Bis-Tris': 'Buffer',
    'Bis-Tris chloride': 'Buffer',
    'CAPS': 'Buffer',
    'CAPSO': 'Buffer',
    'CHAPS': 'Additive',
    'CHAPSO': 'Additive',
    'CHES': 'Buffer',
    'Cadmium chloride': 'Salt',
    'Caffeine': 'Additive',
    'camphor': 'Additive',
    'Cesium acetate': 'Salt',
    'Cholesterol': 'Additive',
    'Cobalt(II) chloride': 'Salt',
    'Cobalt(III) chloride': 'Salt',
    'DL-Serine': 'Additive',
    'DMSO': 'Additive',
    'DTT': 'Additive',
    'DDAO': 'Additive',
    'm-cresol': 'Additive',
    'omecamtiv mecarbil': 'Additive',
    'vanadate': 'Additive',
    'pentafluorobenzyl alcohol': 'Additive',
    '4-iodobenzo[b]thiophene-2-carboxamidine': 'Additive',
    'Dextran sulfate': 'Precipitant',
    'Di-sodium tartrate': 'Salt',
    'Dipotassium hydrogen phosphate': 'Salt',
    'Disodium phosphate': 'Salt',
    'Dithioerythritol': 'Additive',
    'EDTA': 'Additive',
    'EGTA': 'Additive',
    'Ebselen': 'Additive',
    'FAD': 'Additive',
    'Fos-Choline-9': 'Additive',
    'GDP': 'Additive',
    'GTP': 'Additive',
    'Glucose-6-phosphate': 'Salt',
    'Glycylglycine': 'Buffer',
    'Guanidine sulfate': 'Additive',
    'HECAMEG': 'Additive',
    'HEPES': 'Buffer',
    'Hexaamminecobalt(III) chloride': 'Salt',
    'Isomaltose': 'Additive',
    'Jeffamine ED-2001': 'Precipitant',
    'Jeffamine M-600': 'Precipitant',
    'Lactose': 'Additive',
    'Ligand': 'Additive',
    'Polyamine': 'Additive',
    'Lithium acetate': 'Salt',
    'Lithium citrate': 'Salt',
    'MEGA-8': 'Additive',
    'MES': 'Buffer',
    'MOPS': 'Buffer',
    'MPD': 'Additive',
    'Magnesium acetate': 'Salt',
    'Magnesium formate': 'Salt',
    'Malic acid': 'Additive',
    'Monoolein': 'Additive',
    'Monosodium phosphate': 'Salt',
    'Morpheus Amino Acids Mix': 'Additive',
    'Morpheus Buffer System 3': 'Buffer',
    'N-acetylglucosamine': 'Additive',
    'NAD': 'Additive',
    'NADP': 'Additive',
    'NDSB-195': 'Additive',
    'NDSB-221': 'Additive',
    'NDSB-256': 'Additive',
    'NaH2PO4': 'Salt',
    'Nickel(II) chloride': 'Salt',
    'PEG 1000': 'Precipitant',
    'PEG 10000': 'Precipitant',
    'PEG 12000': 'Precipitant',
    'PEG 1500': 'Precipitant',
    'PEG 1450': 'Precipitant',
    'PEG 1550': 'Precipitant',
    'PEG 200': 'Precipitant',
    'PEG 2000': 'Precipitant',
    'PEG 20000': 'Precipitant',
    'PEG': 'Precipitant',
    'PEG MME': 'Precipitant',
    'PEG 300': 'Precipitant',
    'PEG 3000': 'Precipitant',
    'PEG 3350': 'Precipitant',
    'PEG MME 3350': 'Precipitant',
    'PEG 3500': 'Precipitant',
    'PEG 400': 'Precipitant',
    'PEG 4000': 'Precipitant',
    'PEG 5000': 'Precipitant',
    'PEG 600': 'Precipitant',
    'PEG 6000': 'Precipitant',
    'PEG 800': 'Precipitant',
    'PEG 8000': 'Precipitant',
    'PEG MME 500': 'Precipitant',
    'PEG MME 2000': 'Precipitant',
    'PEG MME 5000': 'Precipitant',
    'PEG MME 550': 'Precipitant',
    'PEG Smear High': 'Precipitant',
    'PEG Smear Medium': 'Precipitant',
    'PEG dimethyl ether': 'Precipitant',
    'PGA-LM': 'Precipitant',
    'PIPES': 'Buffer',
    'PLP': 'Additive',
    'PLP_Snyder441': 'Additive',
    'PLP_Snyder457': 'Additive',
    'PLP_Snyder494': 'Additive',
    'PLP_Snyder495': 'Additive',
    'PLP_Snyder496': 'Additive',
    'PLP_Snyder530': 'Additive',
    'PLP_Snyder608': 'Additive',
    'PLP_Snyder630': 'Additive',
    'POPSO': 'Buffer',
    'PPG': 'Precipitant',
    'Pentaerythritol ethoxylate': 'Precipitant',
    'Pentaerythritol ethoxylate 3350': 'Precipitant',
    'Poly(acrylic acid sodium salt) 5100': 'Precipitant',
    'Polyvinylpyrrolidone K15': 'Precipitant',
    'Potassium bromide': 'Salt',
    'Potassium citrate': 'Salt',
    'Potassium fluoride': 'Salt',
    'Potassium formate': 'Salt',
    'Potassium oxalate': 'Salt',
    'Potassium phosphate': 'Salt',
    'Pyruvate': 'Additive',
    'SPG buffer': 'Buffer',
    'Sodium dithionite': 'Salt',
    'Sodium fluoride': 'Salt',
    'Sodium iodide': 'Salt',
    'Sodium maleate': 'Salt',
    'Sodium malonate': 'Salt',
    'Sodium oxamate': 'Additive',
    'Sodium succinate': 'Salt',
    'potassium succinate': 'Salt',
    'Sodium tartrate': 'Salt',
    'Sorbitol': 'Additive',
    'Spermine': 'Additive',
    'Spermine tetrahydrochloride': 'Additive',
    'Sucrose': 'Additive',
    'TAPS': 'Buffer',
    'TAPSO': 'Buffer',
    'TCEP': 'Additive',
    'CTAB': 'Additive',
    'Tacsimate': 'Precipitant',
    'Trehalose': 'Additive',
    'Tricine': 'Buffer',
    'TES': 'Buffer',
    'Tris-HCl': 'Buffer',
    'tris-sulfate': 'Buffer',
    'Tris-acetate': 'Buffer',
    'Tris-maleate': 'Buffer',
    'Trisodium citrate-citric acid': 'Buffer',
    'Triton X-100': 'Additive',
    'Triton X-114': 'Additive',
    'Trizma base': 'Buffer',
    'Tween 20': 'Additive',
    'Tween 80': 'Additive',
    'Xylitol': 'Additive',
    'Yttrium(III) chloride': 'Salt',
    '\\u03b2-octyl glucoside': 'Additive',
    'acetamide': 'Additive',
    'acetate': 'Salt',
    'methyl mercuric acetate': 'Additive',
    'acetate buffer': 'Buffer',
    'acetone': 'Precipitant',
    'acetonitrile': 'Precipitant',
    'alanine': 'Additive',
    'ammonium acetate': 'Salt',
    'ammonium bicarbonate': 'Salt',
    'ammonium chloride': 'Salt',
    'ammonium citrate': 'Salt',
    'ammonium formate': 'Salt',
    'ammonium nitrate': 'Salt',
    'ammonium phosphate': 'Salt',
    'ammonium sulfate': 'Salt',
    'ammonium tartrate': 'Salt',
    'ammonium thiocyanate': 'Additive',
    'arginine': 'Additive',
    'aspartic acid': 'Additive',
    'barium chloride': 'Salt',
    'Bromophenol blue': 'Additive',
    'benzamidine': 'Additive',
    'benzamidine hydrochloride': 'Additive',
    'beta-cyclodextrin': 'Additive',
    'beta-mercaptoethanol': 'Additive',
    'betaine': 'Additive',
    'bis-tris propane': 'Buffer',
    'borate': 'Buffer',
    'butanol': 'Precipitant',
    'cacodylate': 'Buffer',
    'cadmium chloride': 'Salt',
    'cadmium sulfate': 'Additive',
    'calcium acetate': 'Additive',
    'calcium chloride': 'Salt',
    'cesium chloride': 'Salt',
    'citrate buffer': 'Buffer',
    'cobalt chloride': 'Salt',
    'cobalt sulfate': 'Additive',
    'coenzyme A': 'Additive',
    'copper chloride': 'Salt',
    'copper sulfate': 'Additive',
    'dioxane': 'Additive',
    'ethanol': 'Additive',
    'ethylene glycol': 'Additive',
    'triethylene glycol': 'Additive',
    'formamide': 'Additive',
    'fucose': 'Additive',
    'galactose': 'Additive',
    'glucose': 'Additive',
    'glutamic acid': 'Additive',
    'glutathione': 'Additive',
    'glycerol': 'Additive',
    'glycine': 'Buffer',
    'glycine betaine': 'Additive',
    'guanidine hydrochloride': 'Additive',
    'histidine': 'Additive',
    'imidazole': 'Buffer',
    'iodoacetamide': 'Additive',
    'iodoacetic acid': 'Additive',
    'iron chloride': 'Salt',
    'iron sulfate': 'Additive',
    'isopropanol': 'Additive',
    'jeffamine': 'Precipitant',
    'leupeptin': 'Additive',
    'lithium chloride': 'Additive',
    'lithium nitrate': 'Salt',
    'lithium sulfate': 'Salt',
    'lysine': 'Additive',
    'magnesium chloride': 'Salt',
    'magnesium citrate': 'Salt',
    'magnesium nitrate': 'Salt',
    'magnesium sulfate': 'Salt',
    'manganese acetate': 'Salt',
    'manganese chloride': 'Salt',
    'manganese sulfate': 'Additive',
    'mannose': 'Additive',
    'methanol': 'Precipitant',
    'methionine': 'Additive',
    'methyl-beta-cyclodextrin': 'Additive',
    'myo-Inositol': 'Additive',
    'decylmaltoside': 'Additive',
    'n-octyl-beta-D-glucoside': 'Additive',
    'octyl-POE': 'Additive',
    'n-octyl-beta-D-thioglucopyranoside': 'Additive',
    'nickel chloride': 'Salt',
    'nickel sulfate': 'Additive',
    'pepstatin': 'Additive',
    'phenol': 'Additive',
    'pramocaine': 'Additive',
    'phenylmethylsulfonyl fluoride': 'Additive',
    'phosphate buffer': 'Buffer',
    'Phosphate-citrate buffer': 'Buffer',
    'poly(acrylic acid sodium salt) 5100': 'Precipitant',
    'potassium acetate': 'Salt',
    'potassium chloride': 'Salt',
    'potassium iodide': 'Salt',
    'potassium nitrate': 'Salt',
    'potassium phosphate': 'Salt',
    'potassium sodium tartrate': 'Salt',
    'potassium sulfate': 'Salt',
    'potassium tartrate': 'Salt',
    'potassium thiocyanate': 'Additive',
    'proline': 'Additive',
    'propanol': 'Precipitant',
    'rubidium chloride': 'Salt',
    'sarcosine': 'Additive',
    'sodium acetate': 'Salt',
    'copper acetate': 'Salt',
    'sodium azide': 'Additive',
    'sodium bromide': 'Salt',
    'sodium cacodylate': 'Buffer',
    'sodium chloride': 'Salt',
    'sodium citrate': 'Salt',
    'sodium formate': 'Salt',
    'sodium nitrate': 'Salt',
    'sodium phosphate': 'Salt',
    'sodium sulfate': 'Salt',
    'sodium thiocyanate': 'Additive',
    'sodium/potassium phosphate': 'Salt',
    'sodium pyrophosphate': 'Salt',
    'spermidine': 'Additive',
    'strontium chloride': 'Salt',
    'sulfuric acid': 'Additive',
    'acetohydroxamic acid': 'Additive',
    'toluene': 'Additive',
    'sodium sulfide': 'Additive',
    'cadaverine dihydrochloride': 'Additive',
    '3-phosphoglycerate': 'Additive',
    'aluminum chloride': 'Salt',
    'SB12': 'Additive',
    'IHP': 'Additive',
    'sodium ascorbate': 'Additive',
    'sodium monothiophosphate': 'Additive',
    'adenosine': 'Additive',
    'deferoxamine': 'Additive',
    'Mavacamten': 'Additive',
    'resorcinol': 'Additive',
    'beryllium fluoride': 'Additive',
    'sodium sulfite': 'Salt',
    'Gadolinium chloride hexahydrate': 'Additive',
    'Gadolinium nitrate': 'Additive',
    'silver nitrate': 'Additive',
    'phenylphosphorodiamidate': 'Additive',
    'glyceryl tributyrate': 'Additive',
    'paraoxon-ethyl': 'Additive',
    'Deoxy-BigCHAP': 'Additive',
    'potassium ferricyanide': 'Additive',
    'MMT buffer': 'Buffer',
    'ferrous citrate': 'Additive',
    'Tris-buffered saline': 'Buffer',
    'fluorinated octylmaltoside': 'Additive',
    'Ammonium iodide': 'Salt',
    'PCB buffer': 'Buffer',
    'PEG DME 500': 'Precipitant',
    'NBPTO': 'Additive',
    'taurine': 'Additive',
    'tert-butanol': 'Precipitant',
    'thiourea': 'Additive',
    'triethanolamine': 'Buffer',
    'trimethylamine N-oxide': 'Additive',
    'trimethylamine hydrochloride': 'Additive',
    'tryptophan': 'Additive',
    'urea': 'Additive',
    'xylose': 'Additive',
    'zinc': 'Additive',
    'magnesium': 'Additive',
    'zinc acetate': 'Additive',
    'zinc chloride': 'Salt',
    'zinc sulfate': 'Additive',
}


def get_reagent_category(canonical_name: str) -> str:
    """Buffer / Salt / Precipitant / Additive for a canonical reagent name.
    Falls back to "Additive" for anything not in REAGENT_CATEGORY -- e.g. a
    compound the LLM fallback just learned, which has no category yet."""
    return REAGENT_CATEGORY.get(canonical_name, "Additive")


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

_NUMBER = r"\d+(?:\.\d+|,\d{1,2}(?!\d))?(?:[eE][+-]?\d+)?"
# The second alternative handles European decimal-comma notation (e.g.
# "2,25M" meaning 2.25 M, seen in 2R38) without colliding with a
# thousands-grouping comma like "1,000" -- capped at 1-2 digits with a
# lookahead blocking a 3rd, since a real decimal fraction here is never
# written 3+ digits deep but a grouping separator always is. The trailing
# group handles scientific notation ("3e-4 M" meaning 0.0003 M) -- Python's
# float() already parses that form natively once captured whole.
# Optional second number for a range like "15-25%" or "15% - 25%" (a "%"
# is allowed -- and discarded -- right after the first number too, since
# real pdbx_details text writes both forms). Captures into "amount2";
# match_clause() decides what to do with it (see _resolve_range).
_RANGE_SUFFIX = rf"(?:\s*%?\s*[-~]\s*(?P<amount2>{_NUMBER}))?"
_CONC_UNIT = r"(?:%\s*\(?w[/.]v\)?|%\s*\(?v[/.]v\)?|%|millimolar|mmolar|micromolar|microm|molar|mM|uM|\u00b5M|M|mg/mL|mg/ml)"

# Aliases too easily confused with ordinary English words to trust as a
# bare, context-free reagent mention (pattern C) -- e.g. "as" (ammonium
# sulfate) collides with the common word "as", which appears constantly in
# narrative pdbx_details prose ("as described", "as needed", ...). These
# are still recognized in patterns A/B, where an adjacent stated amount and
# unit (e.g. "2.1 M AS") makes the coincidence with ordinary prose far
# less likely, just never as a standalone match with no numeric context.
_CONTEXT_ONLY_ALIASES = {"as"}


def _followed_by_more_prose(clause: str, end: int) -> bool:
    """True if `clause[end:]` continues with another word right after
    (optional) whitespace -- e.g. "...10 mg/ml as determined by...", where
    "as" is followed by "determined". A genuine context-only-alias mention
    (e.g. "2.1 M AS" at a clause boundary) has nothing but punctuation or
    the end of the clause after it. This catches the specific failure mode
    where "as" (ammonium sulfate) collides with ordinary English "X units
    as described/determined/judged/..." -- common boilerplate in deposition
    text -- even though it's adjacent to a stated amount and unit."""
    return re.match(r"\s+[A-Za-z]", clause[end:]) is not None

# Leading method/technique words to strip before matching (not a whole-
# clause discard -- see split_into_clauses, which only drops the clause if
# nothing real is left after stripping this prefix). A depositor sometimes
# runs the method phrase and the real conditions together with no comma at
# all (e.g. "SITTING DROP 300NL PLUS 300NL 0.2M AMMONIUM ACETATE..."); the
# previous version of this pattern ended in ".*$" and discarded the whole
# clause in that case, losing every reagent that followed.
_METHOD_STOPWORDS = re.compile(
    r"^\s*(vapor diffusion|hanging drop|sitting drop|batch|microbatch|"
    r"dialysis|free[- ]interface diffusion|seeding|micro-?seeding)\b",
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

    # Guard the reagent alternation so a short alias (e.g. "og" for octyl
    # glucoside, "ki" for potassium iodide) can only match as a standalone
    # word/token, never as a substring swallowed inside an unrelated longer
    # word -- e.g. "og" inside "nitrOGen" or "microgravity", or "ki" inside
    # "seeKIng". A plain \b boundary doesn't work here since many aliases
    # are chemical formulas starting/ending in non-word characters (e.g.
    # "(nh4)2so4"), where \b itself wouldn't match after a preceding space.
    # Excluding only an adjacent ASCII letter (not digits/punctuation)
    # keeps "0.2M(NH4)2SO4" and "PEG3350" working while still blocking the
    # mid-word false positives above.
    # The trailing guard also excludes a digit, not just a letter: without
    # it, a shorter numeric-weight alias ("peg200") can match as a false
    # prefix of a longer one written with no space ("PEG2000"), stranding
    # the rest ("0") as a bogus amount. The leading guard stays letter-only
    # (see above -- changing it isn't needed for this bug and risks the
    # formula cases it was already tuned for).
    #
    # The alternation itself is wrapped in an atomic group ((?>...), Python
    # 3.11+) so that once the longest-matching alias succeeds here, the
    # engine can never retreat to a shorter one just because something
    # later in the overall pattern fails. Without this, Pattern B (reagent
    # NUMBER UNIT) would backtrack a long, fully-correct multi-word alias
    # like "peg 5000 mme" -- which, correctly, leaves no amount/unit behind
    # it when the text has none -- down to the short generic "peg", whose
    # leftover "5000 MME" then gets misread as a bogus amount+unit ("5000
    # mM"). Aliases are already sorted longest-first on the understanding
    # that the longest match should always win (see above); this just
    # enforces that invariant instead of letting regex backtracking quietly
    # override it.
    reagent_group = rf"(?<![A-Za-z])(?P<reagent>(?>{reagent_alternation}))(?![A-Za-z0-9])"

    # Pattern A: "NUMBER UNIT REAGENT" e.g. "0.2 M ammonium sulfate", also
    # tolerating a trailing unit qualifier that comes AFTER the reagent name,
    # e.g. "10% PEG 8000 (w/v)", and a range before the reagent, e.g.
    # "15-25% PEG 8000" or "15% - 25% PEG 8000" (see _RANGE_SUFFIX).
    # "saturated" also tolerated as a connector (e.g. "53% SATURATED
    # AMMONIUM SULFATE") -- a depositor's way of saying the percentage is
    # of a saturated solution, not a qualifier on the reagent name itself.
    _PATTERN_NUM_FIRST = re.compile(
        rf"\(?(?P<amount>{_NUMBER}){_RANGE_SUFFIX}\s*(?P<unit>{_CONC_UNIT})\)?\s*"
        rf"(?:\(?by\s+(?:volume|weight|vol|wt)\)?\s*)?(?:of\s+(?:each\s+)?)?(?:saturated\s+)?{reagent_group}"
        rf"(?:\s*\(?(?P<trailing_unit>w/v|v/v)\)?)?",
        re.IGNORECASE,
    )
    # Pattern B: "REAGENT NUMBER UNIT" e.g. "PEG 3350 20%" (less common, but seen)
    _PATTERN_REAGENT_FIRST = re.compile(
        rf"{reagent_group}\D{{0,3}}(?P<amount>{_NUMBER}){_RANGE_SUFFIX}\s*(?P<unit>{_CONC_UNIT})",
        re.IGNORECASE,
    )
    # Pattern C: bare reagent name with no concentration nearby at all, e.g.
    # a clause that's just "DMSO" or "ATP added" -- used as a last-resort
    # fallback so a known compound mentioned without a stated concentration
    # is still reported (with concentration "nan") instead of being dropped.
    _PATTERN_REAGENT_ONLY = re.compile(reagent_group, re.IGNORECASE)


def _normalize_unit(unit: str) -> str:
    """Normalize a unit string for consistency."""
    unit = re.sub(r"[()]", "", unit).strip().lower().replace(" ", "")
    # A depositor's typo of "w.v"/"v.v" for "w/v"/"v/v" (period for slash)
    # is accepted by _CONC_UNIT's matching regex but still needs its own
    # punctuation normalized here, or it would fall through the mapping
    # below unmatched and display with the stray period intact.
    unit = re.sub(r"([wv])\.([wv])", r"\1/\2", unit)
    mapping = {
        "%w/v": "% w/v", "%v/v": "% v/v", "%": "%",
        "m": "M", "mm": "mM", "um": "µM", "µm": "µM",
        "mg/ml": "mg/mL",
        "molar": "M", "mmolar": "mM", "millimolar": "mM",
        "micromolar": "µM", "microm": "µM",
    }
    return mapping.get(unit, unit)


def _format_concentration(amount: str, unit: str) -> str:
    """Format concentration as a display string. `amount` may already be a
    resolved range string like "8-18" (see _resolve_range) rather than a
    single number -- those are passed through as-is instead of being fed
    to float(), which would raise and silently degrade the whole
    concentration to "nan" even though a perfectly good range was already
    computed. float() is tried FIRST, not a bare "-" substring check, since
    a negative scientific-notation exponent ("3e-4") also contains a "-"
    but is a single valid number, not a range -- only fall back to the
    range pass-through once float() actually fails to parse it."""
    if not amount:
        return "nan"
    try:
        return f"{float(amount)} {unit}".strip()
    except (ValueError, TypeError):
        pass
    if "-" in amount:
        return f"{amount} {unit}".strip()
    return "nan"


def _convert_micromolar_to_millimolar(amount: str, unit: str) -> Tuple[str, str]:
    """Convert a resolved µM amount -- a single value or an "X-Y" range
    string from _resolve_range -- to its mM-equivalent (x 0.001 = 1e-3),
    so every concentration in the output sits on one consistent mM/M
    scale instead of mixing in µM. Any other unit, or an unparseable
    amount (e.g. "nan"), passes through unchanged."""
    if unit != "µM" or not amount:
        return amount, unit

    def scale(v: str) -> str:
        scaled = float(v) * 0.001
        if scaled == int(scaled):
            return str(int(scaled))
        return f"{scaled:.6f}".rstrip("0").rstrip(".")

    try:
        if "-" in amount:
            low, high = amount.split("-", 1)
            return f"{scale(low)}-{scale(high)}", "mM"
        return scale(amount), "mM"
    except ValueError:
        return amount, unit


def _mean_of(low: str, high: str) -> str:
    """Calculate mean of a range's two endpoints."""
    mean = (float(low) + float(high)) / 2
    if mean == int(mean):
        return str(int(mean))
    return f"{mean:.2f}".rstrip("0").rstrip(".")


def _resolve_range(amount: str, amount2: Optional[str], raw_unit: str, canonical: str) -> str:
    """Resolve a range to a single value or range string. `amount`/`amount2`
    may use a European decimal comma (e.g. "2,25") per the `_NUMBER` regex --
    normalized to a period here, in one place, before either value is ever
    used in a range string or passed to float()."""
    amount = amount.replace(",", ".")
    if amount2:
        amount2 = amount2.replace(",", ".")
    if not amount2:
        return amount
    if canonical.upper().startswith("PEG") and "%" in raw_unit:
        return _mean_of(amount, amount2)
    return f"{amount}-{amount2}"


_PEG_WEIGHT_SCAN = re.compile(
    r"\b(?:peg|polyethylene\s*glycol|poly\(ethylene\s*glycol\))\s*-?\s*(\d[\d,]{1,6})\b",
    re.IGNORECASE,
)


def _auto_register_peg_weights(details: str) -> None:
    """Scan raw pdbx_details text for a "PEG <number>" / "polyethylene
    glycol <number>" mention whose exact molecular-weight number isn't
    already a known alias, and register it as its own canonical
    "PEG <number>" entry on the fly -- so a never-seen-before PEG grade
    (a new commercial MW, or simply one nobody has curated yet) is still
    reported under its own name instead of being silently dropped or
    folded into the generic, molecular-weight-less "PEG" catch-all.
    Deliberately leaves alone any number that's already a registered
    alias -- including the handful of deliberate typo redirects like
    "peg3330"/"peg3340" -> PEG 3350 -- so existing curation always wins
    over this fallback."""
    for m in _PEG_WEIGHT_SCAN.finditer(details):
        digits = m.group(1).replace(",", "")
        if not digits.isdigit():
            continue
        if f"peg {digits}" in _ALIAS_TO_CANONICAL:
            continue  # already known -- built-in or previously auto-learned
        canonical = f"PEG {digits}"
        if update_reagent_dictionary(canonical, extra_aliases=[
            f"peg {digits}", f"peg{digits}", f"peg-{digits}",
            f"polyethylene glycol {digits}", f"poly(ethylene glycol) {digits}",
        ]):
            REAGENT_CATEGORY[canonical] = "Precipitant"


_SLASH_NAMES_WITH_AMOUNT = re.compile(
    rf"(?P<amt>{_NUMBER})\s*(?P<unit>{_CONC_UNIT})\s+"
    rf"(?P<names>[A-Za-z][A-Za-z\-]*(?:/[A-Za-z][A-Za-z\-]*)+)\b",
    re.IGNORECASE,
)


def _expand_shared_slash_amount(details: str) -> str:
    """A depositor sometimes states one amount+unit for several
    slash-joined reagent names at once, meaning a mixed buffer where that
    same amount applies to each one (e.g. "0.1 M MES/Imidazole" -- a
    MES-imidazole buffer at 0.1 M each), not a single combined compound.
    Rewritten here into a comma-separated repeat ("0.1 M MES, 0.1 M
    Imidazole") so ordinary matching picks up every component instead of
    attaching the amount to only the first and losing it for the rest.
    Only fires when EVERY slash-separated piece is independently a known
    alias on its own -- a registered combined name like "sodium/potassium
    phosphate" (matched as one alias already) is left alone, and so is a
    mix where some piece isn't independently curated, since guessing a
    split there risks being wrong -- except a bare ion-prefix segment
    ("Na", "K", "NH4"), which is never a reagent by itself and is just
    dropped rather than emitted as its own (wrong) match."""
    ION_PREFIXES = {"na", "k", "nh4", "ca", "mg"}
    def repl(m: re.Match) -> str:
        whole = m.group("names").strip().lower()
        if whole in _ALIAS_TO_CANONICAL:
            return m.group(0)
        segments = m.group("names").split("/")
        real_segments = [s for s in segments if s.strip().lower() not in ION_PREFIXES]
        if not real_segments or not all(
            s.strip().lower() in _ALIAS_TO_CANONICAL for s in real_segments
        ):
            return m.group(0)
        amt, unit = m.group("amt"), m.group("unit")
        return ", ".join(f"{amt} {unit} {seg.strip()}" for seg in real_segments)
    return _SLASH_NAMES_WITH_AMOUNT.sub(repl, details)


def split_into_clauses(details: str) -> List[str]:
    """Split a pdbx_details string into individual condition clauses."""
    if not details:
        return []
    _auto_register_peg_weights(details)
    # Deposited pdbx_details text is often fixed-width line-wrapped inside
    # the mmCIF file, so a newline is usually just a wrap artifact in the
    # middle of one sentence/clause -- not a real separator between two
    # distinct conditions. Normalizing it to a space (rather than splitting
    # on it, as commas/semicolons do) keeps an amount and its reagent
    # together when a wrap happens to fall between them, e.g. "1.02 M
    # \nMGSO4" would otherwise be severed into two clauses and lose the
    # concentration.
    details = re.sub(r"\r\n|\r|\n", " ", details)
    # Collapse any run of whitespace this just created (e.g. a line already
    # ending in a space, now followed by the space the newline became) down
    # to one -- a literal multi-word alias like "sodium sulphate" matches
    # exactly one space between its words, not two.
    details = re.sub(r"[ \t]{2,}", " ", details)
    # A depositor occasionally typo's a zero as a capital "O" right before
    # a concentration unit (e.g. "1OOMM" for "100MM") -- almost certainly
    # a typing slip on an all-caps keyboard, not a real word, since it only
    # fires when sandwiched between a digit and a known unit abbreviation.
    # Deliberately uppercase-only (not re.IGNORECASE): a lowercase "o" in
    # that position is far more likely to be a real word fragment.
    details = re.sub(r"(\d+)(O+)(?=(?:MM|UM|M)\b)",
                      lambda m: m.group(1) + "0" * len(m.group(2)), details)
    # Same typo, leading position: "O.2M" for "0.2M" -- a capital "O"
    # standing in for the leading zero of a decimal amount, recognized by
    # being immediately followed by ".<digit>" (no real word matches that
    # shape). Word-boundary-anchored so it only fires at the start of a
    # number-like token, not mid-word.
    details = re.sub(r"\bO(?=\.\d)", "0", details)
    # A depositor occasionally sandwiches the w/v-or-v/v qualifier between
    # "PEG" and its molecular weight (e.g. "PEG (v/v) 600") instead of
    # after it. That breaks the "peg 600" alias match entirely -- the
    # weight ends up stranded and dropped. Moving the qualifier after the
    # number recovers the normal "PEG 600 (v/v)" shape every other case
    # already handles.
    details = re.sub(
        r"\b(peg|polyethylene\s*glycol)\s*\(?(w/v|v/v)\)?\s+(\d[\d,]*)\b",
        r"\1 \3 (\2)", details, flags=re.IGNORECASE,
    )
    # Same idea for the "(PEG)" abbreviation parenthetical a depositor
    # sometimes drops right after spelling out "polyethylene glycol", before
    # the molecular weight (e.g. "polyethylene glycol (PEG) 8000"). It
    # breaks the "polyethylene glycol 8000" alias match the same way the
    # w/v/v/v qualifier above does; just drop it here rather than relocating
    # it, since it carries no information "polyethylene glycol" doesn't.
    details = re.sub(
        r"\bpolyethylene\s*glycol\s*\(peg\)\s+(?=\d)", "polyethylene glycol ",
        details, flags=re.IGNORECASE,
    )
    details = _expand_shared_slash_amount(details)
    # A depositor sometimes restates the pH right before the amount inside
    # the same parenthetical, e.g. "sodium acetate (pH 5.0, 50mM)". The
    # "pH 5.0, " prefix sits too far from the reagent name for Pattern B's
    # small reagent-to-amount gap tolerance to bridge, stranding the amount.
    # Dropping just that prefix (the pH value itself is never a condition
    # amount worth extracting) leaves "(50mM)" directly after the reagent,
    # which Pattern B already handles.
    details = re.sub(r"\(\s*ph\s*[\d.]+\s*,\s*", "(", details, flags=re.IGNORECASE)
    # Protect chemical locant-prefix commas (e.g. "N,N-dimethyl...",
    # "O,O-diethyl...") from being split on -- that comma is part of the
    # compound's name, not a separator between two conditions. A sentinel
    # swap (rather than a lookaround baked into the split pattern) keeps
    # this independent of the digit-adjacency rule below.
    details = re.sub(r"\b([A-Z]),(?=[A-Z]-)", "\\1\x00", details)
    # Also protect a comma inside a parenthetical remark, e.g. "sodium
    # acetate (pH 5.0, 50mM)" -- that comma separates two details about the
    # SAME reagent (its pH and its amount), not two different conditions,
    # and splitting on it strands the amount in its own clause with no
    # reagent to attach to. Non-nested only (crystallization text doesn't
    # nest parens in practice), same sentinel-swap technique as above.
    details = re.sub(r"\([^()]*\)",
                      lambda m: m.group(0).replace(",", "\x00"), details)
    raw_clauses = re.split(r"(?<!\d),|,(?!\d)|;|\band\b", details, flags=re.IGNORECASE)
    raw_clauses = [c.replace("\x00", ",") for c in raw_clauses]
    clauses = []
    for c in raw_clauses:
        c = c.strip(" .'\"")  # Strip spaces, periods, single and double quotes
        if not c:
            continue
        stopword_match = _METHOD_STOPWORDS.match(c)
        if stopword_match:
            c = c[stopword_match.end():].strip(" .'\"")
            if not c:
                continue
        if re.match(r"^\d+(\.\d+)?\s*K$", c, re.IGNORECASE):
            continue
        if re.match(r"^pH\s*[:=]?\s*\d+(\.\d+)?\s*$", c, re.IGNORECASE):
            continue
        clauses.append(c)
    return clauses


def _match_clause_with_span(clause: str) -> Tuple[Optional[Dict[str, str]], Optional[Tuple[int, int]]]:
    """Core matcher shared by match_clause() and match_clause_all(): tries
    patterns A, B, C in order and returns both the parsed
    (compound, concentration) dict and the (start, end) span of the regex
    match within `clause`. The span lets match_clause_all() blank out what
    it already consumed and keep looking for further reagents in the same
    clause; match_clause() just discards it."""
    m = _PATTERN_NUM_FIRST.search(clause)
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        if reagent_alias in _CONTEXT_ONLY_ALIASES and _followed_by_more_prose(clause, m.end()):
            m = None
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        amount = _resolve_range(m.group("amount").strip(), m.groupdict().get("amount2"),
                                 m.group("unit"), canonical)
        unit = _normalize_unit(m.group("unit"))
        amount, unit = _convert_micromolar_to_millimolar(amount, unit)
        trailing_unit = m.groupdict().get("trailing_unit")
        if trailing_unit and "w/v" not in unit and "v/v" not in unit:
            unit = f"{unit} {trailing_unit.lower()}" if unit else trailing_unit.lower()
        return {
            "compound": canonical,
            "amount": amount,
            "unit": unit,
            "concentration": _format_concentration(amount, unit),
            "source": "dictionary",
        }, m.span()

    m = _PATTERN_REAGENT_FIRST.search(clause)
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        if reagent_alias in _CONTEXT_ONLY_ALIASES and _followed_by_more_prose(clause, m.end()):
            m = None
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        amount = _resolve_range(m.group("amount").strip(), m.groupdict().get("amount2"),
                                 m.group("unit"), canonical)
        unit = _normalize_unit(m.group("unit"))
        amount, unit = _convert_micromolar_to_millimolar(amount, unit)
        return {
            "compound": canonical,
            "amount": amount,
            "unit": unit,
            "concentration": _format_concentration(amount, unit),
            "source": "dictionary",
        }, m.span()

    m2 = _PATTERN_REAGENT_ONLY.search(clause)
    if m2:
        reagent_alias = m2.group("reagent").strip().lower()
        if reagent_alias in _CONTEXT_ONLY_ALIASES:
            return None, None
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        return {
            "compound": canonical,
            "amount": None,
            "unit": None,
            "concentration": "nan",
            "source": "dictionary",
        }, m2.span()

    return None, None


def match_clause(clause: str) -> Optional[Dict[str, str]]:
    """Try to match one clause against the reagent dictionary, returning
    only the first (highest-priority) match. A clause can hold more than
    one reagent when split_into_clauses() failed to separate them (e.g. a
    comma sitting between two digits, or a period used instead of a comma
    between conditions) -- use match_clause_all() to recover every reagent
    in that case instead of just this one."""
    match, _span = _match_clause_with_span(clause)
    return match


def match_clause_all(clause: str) -> List[Dict[str, str]]:
    """Like match_clause(), but returns every non-overlapping reagent match
    found in the clause instead of only the first. After each match, its
    span is blanked out (replaced with a single space, not deleted, so
    text on either side doesn't get fused into a new false match) and the
    search repeats on what's left, until nothing more matches."""
    matches: List[Dict[str, str]] = []
    remaining = clause
    for _ in range(20):  # safety valve against pathological input
        if not remaining.strip():
            break
        match, span = _match_clause_with_span(remaining)
        if not match or span is None:
            break
        matches.append(match)
        start, end = span
        remaining = remaining[:start] + " " + remaining[end:]
    return matches


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


def process_csv(input_csv: str, output_csv: str, pdbx_column: str = "pdbx_details") -> None:
    """Extract compounds from a CSV file using the reagent dictionary.
    Saves the results to output_csv.

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
            clause_matches = match_clause_all(clause)
            if clause_matches:
                matched_compounds.extend(clause_matches)

        compounds_list.append(format_compounds(matched_compounds))
    
    # Add the compound column to the dataframe
    df["compound"] = compounds_list
    
    # Write the full dataframe back out
    try:
        df.to_csv(output_csv, index=False)
    except Exception:
        pass
