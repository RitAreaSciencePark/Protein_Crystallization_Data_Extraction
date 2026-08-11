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
    "PEG 1000": ["peg 1000", "peg1000", "peg-1000", "peg1k", "peg 1k", "polyethylene glycol 1000", "poly(ethylene glycol) 1000", "peg 1,000"],
    "PEG 3350": ["peg 3350", "peg3350", "peg-3350", "polyethylene glycol 3350", "poly(ethylene glycol) 3350"],  
    "PEG 4000": ["peg 4000", "peg4000", "peg-4000", "peg4k", "peg 4k", "polyethylene glycol 4000", "poly(ethylene glycol) 4000"],
    "PEG 6000": ["peg 6000", "peg6000", "peg-6000", "peg6k", "peg 6k", "polyethylene glycol 6000", "poly(ethylene glycol) 6000"],
    "PEG 8000": ["peg 8000", "peg8000", "peg-8000", "peg8k", "peg 8k", "polyethylene glycol 8000", "poly(ethylene glycol) 8000"],
    "PEG 10000": ["peg 10000", "peg10000", "peg-10000", "peg10k", "peg 10k", "polyethylene glycol 10000"],
    "PEG 5000": ["peg 5000", "peg5000", "peg-5000", "peg5k", "peg 5k", "polyethylene glycol 5000"],
    "PEG MME 550": ["peg mme 550", "peg monomethyl ether 550"],
    "PEG MME 2000": ["peg mme 2000", "peg mme 2k", "peg mme2k", "peg monomethyl ether 2000"],
    "PEG MME 5000": ["peg mme 5000", "peg monomethyl ether 5000", "polyethylene glycol monomethyl ether 5,000", "Polyethylene glycol monomethyl ether 5000"],
    "PEG 200": ["peg 200", "peg200", "polyethylene glycol 200", "poly(ethylene glycol) 200", "peg-200"],
    "PEG 300": ["peg 300", "peg300", "polyethylene glycol 300", "poly(ethylene glycol) 300", "peg-300"],
    "PEG 1500": ["peg 1500", "peg1500", "polyethylene glycol 1500", "poly(ethylene glycol) 1500", "peg-1500"],
    "PEG 2000": ["peg 2000", "peg2000", "peg-2000", "polyethylene glycol 2000", "poly(ethylene glycol) 2000"],
    "PEG 600": ["peg 600", "peg600", "polyethylene glycol 600", "poly(ethylene glycol) 600", "peg-600"],
    "PEG 800": ["peg 800", "peg800", "polyethylene glycol 800", "poly(ethylene glycol) 800", "peg-800"],
    "PEG MME 75\u2080": ["peg mme 75\u2080", "peg monomethyl ether 75\u2080"  ],
    "PEG dimethyl ether": ["peg dimethyl ether", "polyethylene glycol dimethyl ether"],
    "PEG 3000": ["peg 3000", "peg3000", "peg-3000", "polyethylene glycol 3000"],
    "PEG 12000": ["peg 12000", "peg12000", "peg-12000", "peg12k", "peg 12k", "polyethylene glycol 12000"],
    "PEG 20000": ["peg 20000", "peg20000", "peg-20000", "peg20k", "peg 20k", "polyethylene glycol 20000"],
    "Jeffamine ED-2001": ["jeffamine ed-2001", "ed2001"],

    # --- Salts ---
    "ammonium sulfate": ["ammonium sulfate", "ammonium sulphate", "(nh4)2so4"],
    "sodium chloride": ["sodium chloride", "nacl"],
    "lithium sulfate": ["lithium sulfate", "lithium sulphate", "li2so4"],
    "lithium chloride": ["lithium chloride", "licl"],
    "magnesium chloride": ["magnesium chloride", "mgcl2"],
    "magnesium sulfate": ["magnesium sulfate", "magnesium sulphate", "mgso4"],
    "calcium chloride": ["calcium chloride", "cacl2"],
    "sodium phosphate": ["sodium phosphate", "nah2po4", "na2hpo4"],
    "potassium phosphate": ["potassium phosphate", "kh2po4", "k2hpo4"],
    "sodium/potassium phosphate": ["sodium potassium phosphate", "sodium/potassium phosphate", "na/k phosphate", "napo4"],
    "sodium acetate": ["sodium acetate"],
    "sodium citrate": ["sodium citrate", "tri-sodium citrate", "NAcitrate", "trisodium citrate", "na citrate, Na-citrate"],
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
    "Monosodium phosphate": ["monosodium phosphate", "sodium phosphate monobasic"],
    "NaH2PO4": ["nah2po4", "na h2 po4", "sodium dihydrogenphosphate", "monosodium phosphate", "sodium phosphate monobasic"], 
    "Sodium dihydrogen phosphate": ["nah2po4", "na h2po4", "sodium dihydrogen phosphate", "monosodium phosphate", "sodium phosphate monobasic"],
    "Potassium dihydrogen phosphate": ["kh2po4", "k h2po4", "potassium hydrogenphosphate'", "potassium dihydrogen phosphate", "monopotassium phosphate", "potassium phosphate monobasic"],
    "Disodium phosphate": ["disodium phosphate", "sodium phosphate dibasic", "na2hpo4", "na2 hpo4"],
    "Dipotassium hydrogen phosphate": ["dipotassium hydrogen phosphate", "potassium phosphate dibasic", "k2hpo4", "k2 hpo4"],
    "Sodium/Potassium phosphate": ["nah2po4/kh2po4", "nah2po4/kh2po4", "nah2po4 / kh2po4", "nah2po4/kh2po4", "nah2po4-kh2po4", "sodium/potassium phosphate", "sodium potassium phosphate",
    "mixed phosphate buffer", "na/k phosphate buffer", "na/k phosphate"],
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
                                   "rochelle salt", "kna tartrate", "K/Na tartrate tetrahydrate"],
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
    "Cobalt(II) chloride": ["cobalt chloride", "cobalt(ii) chloride", "cocl2", "cobaltous chloride"],
    "Cadmium chloride": ["cadmium chloride", "cadmium(ii) chloride", "cdcl2"],
    "Nickel(II) chloride": ["nickel chloride", "nickel(ii) chloride", "nicl2"],
    "Sodium iodide": ["sodium iodide", "nai","na iodide", "NAI"],
    "Sodium succinate": ["sodium succinate","sodium succinate buffer","disodium succinate","succinic acid disodium salt","succinate","na succinate","sodium hydrogen succinate"],
    "Cobalt(II) chloride": ["cobalt(ii) chloride", "cobalt chloride", "cocl2", "cocl₂", "cobaltous chloride"],
    "Cobalt(III) chloride": ["cobalt(iii) chloride", "cocl3", "cocl₃"],

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
    "MIB buffer": ["mib buffer", "mib", "MIB buffer"],

    "TCEP": ["tcep", "tris(2-carboxyethyl)phosphine", "tris(2-carboxyethyl)phosphine hydrochloride", "tris carboxyethyl phosphine", "tcep-hcl", "tcep hydrochloride"],

    # --- Organic precipitants ---
    "1,6-hexanediol": ["1,6-hexanediol", "hexanediol"],
    "1,4-butanediol": ["1,4-butanediol", "butanediol" ],
    "2-propanol": ["2-propanol", "isopropanol", "isopropyl alcohol"],
    "tert-butanol": ["tert-butanol", "t-butanol", "tbutanol"],
    "methanol": ["methanol", "methyl alcohol"],
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
    "ethylene glycol": ["ethylene glycol"],
    "DTT": ["dtt", "dithiothreitol"],
    "beta-mercaptoethanol": ["beta-mercaptoethanol", "2-mercaptoethanol", "bme"],
    "TCEP": ["tcep"],
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
    "coenzyme A": ["coa","coenzyme a"],
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
    "Tryptone": ["tryptone", "casein tryptic digest", "pancreatic digest of casein", "enzymatic digest of casein", "tryptone peptone", "bacto tryptone"],
    
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

_NUMBER = r"\d+(?:\.\d+)?(?:\s*-\s*\d+(?:\.\d+)?)?"
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
    # e.g. "10% PEG 8000 (w/v)".
    _PATTERN_NUM_FIRST = re.compile(
        rf"(?P<amount>{_NUMBER})\s*(?P<unit>{_CONC_UNIT})\s*(?P<reagent>{reagent_alternation})"
        rf"(?:\s*\(?(?P<trailing_unit>w/v|v/v)\)?)?",
        re.IGNORECASE,
    )
    # Pattern B: "REAGENT NUMBER UNIT" e.g. "PEG 3350 20%" (less common, but seen)
    _PATTERN_REAGENT_FIRST = re.compile(
        rf"(?P<reagent>{reagent_alternation})\D{{0,3}}(?P<amount>{_NUMBER})\s*(?P<unit>{_CONC_UNIT})",
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


def learn_new_compound(compound: str, extra_aliases: Optional[List[str]] = None) -> bool:
    """
    Add a newly-discovered compound to the dictionary (in memory) and
    persist it to LEARNED_DICT_PATH, then rebuild the matching regexes so
    it's immediately usable by the dictionary/regex pass -- for the rest
    of the current run, and for every run after this one.

    Returns True if the compound was actually new (and got added), False
    if it was already known (nothing to do).
    """
    canonical = compound.strip()
    alias = canonical.lower()
    if not canonical or alias in _ALIAS_TO_CANONICAL:
        return False  # already known under this name or some other alias

    aliases = [alias] + [a.strip().lower() for a in (extra_aliases or [])]
    REAGENT_DICTIONARY.setdefault(canonical, [])
    for a in aliases:
        if a not in REAGENT_DICTIONARY[canonical]:
            REAGENT_DICTIONARY[canonical].append(a)

    learned = _load_learned_dictionary()
    learned.setdefault(canonical, [])
    for a in aliases:
        if a not in learned[canonical]:
            learned[canonical].append(a)
    _save_learned_dictionary(learned)

    _rebuild_patterns()
    return True


_rebuild_patterns()  # build the initial regexes from built-in + learned dictionary


def _normalize_unit(unit: str) -> str:
    unit = re.sub(r"[()]", "", unit).strip().lower().replace(" ", "")
    mapping = {
        "%w/v": "% w/v", "%v/v": "% v/v", "%": "%",
        "m": "M", "mm": "mM", "um": "\u00b5M", "\u00b5m": "\u00b5M",
        "mg/ml": "mg/mL",
    }
    return mapping.get(unit, unit)


def split_into_clauses(details: str) -> List[str]:
    """Split a pdbx_details string into individual condition clauses,
    dropping clauses that are just method/technique descriptions."""
    if not details:
        return []
    raw_clauses = re.split(r"[,;]|\band\b|\r\n|\n|\r", details, flags=re.IGNORECASE)
    clauses = []
    for c in raw_clauses:
        c = c.strip(" .")
        if not c:
            continue
        if _METHOD_STOPWORDS.match(c):
            continue
        if re.match(r"^\d+(\.\d+)?\s*K$", c, re.IGNORECASE):  # bare temperature clause
            continue
        if re.match(r"^pH\b", c, re.IGNORECASE):  # bare "pH 7.5" clause, not a compound
            continue
        clauses.append(c)
    return clauses


def match_clause(clause: str) -> Optional[Dict[str, str]]:
    """Try to match one clause against the reagent dictionary. Returns a
    dict with compound/amount/unit/concentration, or None if no dictionary
    reagent was found in this clause at all.

    Tries, in order:
      A. "NUMBER UNIT REAGENT" (e.g. "0.2 M ammonium sulfate"), also
         catching a trailing unit qualifier after the reagent name (e.g.
         "10% PEG 8000 (w/v)").
      B. "REAGENT NUMBER UNIT" (e.g. "PEG 3350 20%")
      C. a bare reagent name with no concentration nearby (e.g. just
         "DMSO") -- reported with concentration "nan" rather than dropped.
    """
    m = _PATTERN_NUM_FIRST.search(clause)
    if m:
        reagent_alias = m.group("reagent").strip().lower()
        canonical = _ALIAS_TO_CANONICAL.get(reagent_alias, reagent_alias)
        amount = m.group("amount").strip()
        unit = _normalize_unit(m.group("unit"))

        trailing_unit = m.groupdict().get("trailing_unit")
        if trailing_unit and "w/v" not in unit and "v/v" not in unit:
            # e.g. base unit "%" + trailing "w/v" -> "% w/v"
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
        amount = m.group("amount").strip()
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

def _format_concentration(amount: str, unit: str) -> str:
    """
    Format concentration by converting molarity units to M.

    Examples
    --------
    200 mM  -> 0.2 M
    500 µM  -> 0.0005 M
    10 nM   -> 1e-08 M
    2 M     -> 2 M
    10% w/v -> 10% w/v
    5 mg/mL -> 5 mg/mL
    """

    if not amount or not unit:
        return f"{amount} {unit}".strip()

    unit = unit.strip()

    # Keep percentage units unchanged
    if unit.startswith("%"):
        return f"{amount}{unit}"

    try:
        value = float(amount)
    except ValueError:
        return f"{amount} {unit}"

    unit_lower = unit.lower().replace("μ", "u").replace("µ", "u")

    conversion = { "m": 1, "mm": 1e-3, "um": 1e-6, "nm": 1e-9, "pm": 1e-12,}
    if unit_lower in conversion:
        value *= conversion[unit_lower]

        # Remove trailing zeros
        value_str = f"{value:.12g}"

        return f"{value_str} M"

    # Leave all other units unchanged
    return f"{amount} {unit}"



def extract_compounds(details: Optional[str]) -> Tuple[List[Dict[str, str]], List[str]]:
    """
    Extract all recognizable (compound, concentration) pairs from a
    pdbx_details string.

    Returns
    -------
    (matches, unmatched_clauses)
        matches: list of dicts (compound/amount/unit/concentration/source)
        unmatched_clauses: clauses that had a number+unit pattern (or looked
        like a real condition) but didn't hit the dictionary -- candidates
        for the LLM fallback.
    """
    if not details or not isinstance(details, str):
        return [], []

    matches = []
    unmatched = []
    for clause in split_into_clauses(details):
        result = match_clause(clause)
        if result:
            matches.append(result)
        else:
            # only worth sending to the LLM fallback if it actually looks
            # like it contains a concentration (avoids wasting calls on
            # stray fragments like "crystals appeared after 3 days")
            if re.search(_NUMBER + r"\s*" + _CONC_UNIT, clause, re.IGNORECASE):
                unmatched.append(clause)

    return matches, unmatched


# --------------------------------------------------------------------------- #
# 2. Optional LLM fallback for clauses the dictionary can't resolve
# --------------------------------------------------------------------------- #
_LLM_SYSTEM_PROMPT = (
    "You extract chemical compound names and concentrations from short "
    "crystallization-condition text fragments. Given one or more text "
    "fragments, return ONLY a JSON array (no prose, no markdown fences) "
    "where each element is {\"compound\": str, \"amount\": str, \"unit\": str}. "
    "If a fragment has no extractable compound+concentration, omit it. "
    "Normalize compound names to their common chemical name."
)


def llm_extract_compounds(
    unmatched_clauses: List[str],
    api_key: Optional[str] = None,
    model: str = "claude-sonnet-5",
) -> List[Dict[str, str]]:
    """
    Send clauses the dictionary/regex pass couldn't resolve to Claude for
    structured extraction. Only call this on the residual unmatched text,
    not on every row -- that's what keeps this affordable at scale.

    Requires `requests` and an Anthropic API key (passed in, or read from
    the ANTHROPIC_API_KEY environment variable).
    """
    if not unmatched_clauses:
        return []
    if requests is None:
        raise RuntimeError("The 'requests' package is required for the LLM fallback.")

    api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise RuntimeError(
            "No Anthropic API key provided. Pass api_key=... or set "
            "the ANTHROPIC_API_KEY environment variable."
        )

    user_message = "Fragments:\n" + "\n".join(f"- {c}" for c in unmatched_clauses)

    resp = requests.post(
        "https://api.anthropic.com/v1/messages",
        headers={
            "x-api-key": api_key,
            "anthropic-version": "2023-06-01",
            "content-type": "application/json",
        },
        json={
            "model": model,
            "max_tokens": 1024,
            "system": _LLM_SYSTEM_PROMPT,
            "messages": [{"role": "user", "content": user_message}],
        },
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()

    text_blocks = [b["text"] for b in data.get("content", []) if b.get("type") == "text"]
    raw_text = "\n".join(text_blocks).strip()
    raw_text = re.sub(r"^```(json)?|```$", "", raw_text.strip(), flags=re.MULTILINE).strip()

    try:
        parsed = json.loads(raw_text)
    except json.JSONDecodeError:
        return []

    results = []
    for item in parsed:
        compound = item.get("compound")
        amount = item.get("amount")
        unit = _normalize_unit(item.get("unit", "")) if item.get("unit") else ""
        if compound:
            results.append({
                "compound": compound,
                "amount": amount or "",
                "unit": unit,
                "concentration": _format_concentration(amount or "", unit) if unit else str(amount or ""),
                "source": "llm",
            })
    return results


# --------------------------------------------------------------------------- #
# 3. CSV pipeline
# --------------------------------------------------------------------------- #
def format_compound_string(matches: List[Dict[str, str]]) -> Optional[str]:
    """Format a list of matches into a single string:
    'compound1' (concentration1), 'compound2' (concentration2), ...
    Returns None if there are no matches."""
    if not matches:
        return None
    return ", ".join(f"'{m['compound']}' ({m['concentration']})" for m in matches)


def process_csv(
    input_csv: str,
    output_csv: str,
    pdbx_column: str = "pdbx_details",
    use_llm_fallback: bool = False,
    api_key: Optional[str] = None,
) -> str:
    """
    Read `input_csv` (e.g. your pipeline's Output.csv), extract compounds
    from `pdbx_column` for every row, and write the SAME rows back out with
    one added `compound` column formatted as:

        'compound1' (concentration1), 'compound2' (concentration2), ...

    (or blank if nothing was found for that row). One row in -> one row
    out; no amount/unit/concentration/source columns are added.
    """
    df = pd.read_csv(input_csv)
    if pdbx_column not in df.columns:
        raise KeyError(f"Column '{pdbx_column}' not found. Available: {list(df.columns)}")

    compound_strings = []
    for _, row in df.iterrows():
        details = row.get(pdbx_column)
        matches, unmatched = extract_compounds(details)

        if use_llm_fallback and unmatched:
            llm_matches = llm_extract_compounds(unmatched, api_key=api_key)
            for m in llm_matches:
                if learn_new_compound(m["compound"]):
                    print(f"Learned new compound: '{m['compound']}' "
                          f"(saved to {LEARNED_DICT_PATH})")
            matches.extend(llm_matches)

        compound_strings.append(format_compound_string(matches))

    df["compound"] = compound_strings
    df.to_csv(output_csv, index=False)
    print(f"Wrote {len(df)} rows to {output_csv}")
    return output_csv


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Extract compounds/concentrations from pdbx_details.")
    p.add_argument("input_csv", help="Input CSV (e.g. Output.csv) containing a pdbx_details column.")
    p.add_argument("output_csv", help="Path to write the expanded compound-level CSV to.")
    p.add_argument("--llm-fallback", action="store_true",
                   help="Use Claude to resolve clauses the dictionary/regex pass misses.")
    p.add_argument("--api-key", default=None, help="Anthropic API key (or set ANTHROPIC_API_KEY).")
    args = p.parse_args()

    process_csv(args.input_csv, args.output_csv,
                use_llm_fallback=args.llm_fallback, api_key=args.api_key)