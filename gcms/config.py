"""
Static knowledge base used by the re-identification engine.

Editing these lists is the standard way to extend the tool to new
matrices or analytical protocols.
"""
from __future__ import annotations

import re

# ---------------------------------------------------------------------------
# Implausible substrings — any candidate matching one of these is rejected
# while walking the NIST hit list.
# ---------------------------------------------------------------------------
IMPLAUSIBLE_SUBSTRINGS: list[str] = [
    # Drugs / pharmaceuticals
    "Floxuridine", "Cidofovir", "Niflumic", "Debrisoquine",
    "Ergotaman", "Ergotamine", "Ergocristine", "Ergocornine",
    "Ergosine", "Bromocriptine", "Pergolide", "Cabergoline",
    # Industrial / plasticisers / phosphite stabilisers
    "Tris(2,4-di-tert-butylphenyl) phosphate",
    "tert-butylphenyl) phosphate",
    "Triphenyl phosphate",
    # Synthetic dyes / scintillators / oxadiazoles
    "oxadiazole", "1,3,4-oxadiazol", "biphenylyl",
    "fluoren", "quinolizine", "phenanthridine",
    # Halogenated synthetics
    "Bis(bromomethyl)", "bromomethyl)-1,1'-biphenyl",
    "1-iodo", "2-iodo", "3-iodo", ", iodo-",
    "Iodoethane",
    "1,54-dibromo", ", dibromo-", ", dichloro-",
    # Silyl artefacts (only the obviously synthetic ones — TMS suffixes
    # on real compounds are stripped post-hoc, see strip_tms())
    "Disilapentane", "Disilane",
    "Silacyclohex", "Silacyclopent",
    "Trioxa-5-aza-1-silabicyclo",
    # Triglycerides whose RI cannot match this oven program
    "Triarachine", "Tristearin", "Tripalmitin",
    "Trimyristin", "Trilinolein", "Trilaurin", "Tricaprin",
    # Synthetic indane
    "1H-Indene, 1-hexadecyl",
]

# ---------------------------------------------------------------------------
# Contaminant substrings — kept only as a last-resort fallback identity and
# tagged as `is_contaminant_final = True` so they can be excluded from
# biological diversity statistics down-stream.
# ---------------------------------------------------------------------------
CONTAMINANT_SUBSTRINGS: list[str] = [
    "phthalate", "Phthalic",
    "Benzenedicarboxylic acid, bis",
    "Benzenedicarboxylic acid, di",
    "Benzenedicarboxylic acid, decyl",
    # Cyclic siloxanes (column bleed)
    "Cyclotrisiloxane", "Cyclotetrasiloxane",
    "Cyclopentasiloxane", "Cyclohexasiloxane",
    "Cycloheptasiloxane", "Cyclooctasiloxane",
    "Cyclononasiloxane", "Cyclodecasiloxane",
    "Cycloundecasiloxane", "Cyclododecasiloxane",
    "Octamethylcyclo", "Decamethylcyclo", "Dodecamethylcyclo",
    "Tetradecamethyl", "Hexadecamethyl", "Octadecamethyl",
    "Heptasiloxane", "Trisiloxane", "Disiloxane", "polysiloxane",
    "1,4-Bis(trimethylsilyl)-1,3-butadiyne",
    # Squalene (column conditioning)
    "Squalene",
    # Acrylates / methacrylates
    "2-Propenoic acid, 2-methyl-, octyl ester",
    "methacrylate",
    # Perfluorinated derivatising-agent residues
    "pentafluoropropionate",
    "heptafluorobutyrate",
]


# ---------------------------------------------------------------------------
# Compound-class regular-expression rules. Order matters — more specific
# patterns evaluated first.
# ---------------------------------------------------------------------------
CLASS_RULES: list[tuple[str, re.Pattern]] = [
    ("Diketopiperazine",     re.compile(r"Pyrrolo\[1,2-a\]pyrazine|Piperazine-1,4-dione|"
                                         r"2,5-piperazinedione|Diketopiperazine", re.I)),
    ("Pyrazine",             re.compile(r"\bpyrazine\b(?!.*Pyrrolo)", re.I)),
    ("Pyranone/Furanone",    re.compile(r"\bpyranone|\bpyrone|furanone|\bfuran-2", re.I)),
    ("Indole/Pyrrole",       re.compile(r"\bindol|\bpyrrole(?!.*pyrazine)", re.I)),
    ("Phenolic",             re.compile(r"\bphenol\b|2,4-Di-tert-butylphenol|catechol|"
                                         r"hydroxybenzene|hydroxyphenyl|trihydroxybenzaldehyde|"
                                         r"4,4'-butylidenebis", re.I)),
    ("Sterol",               re.compile(r"sterol|cholesta|ergosta|stigmast|sitosterol|squalen", re.I)),
    ("Mono/Di-glyceride",    re.compile(r"Monopalmitin|Monostearin|Monoolein|Dipalmitin|Distearin|"
                                         r"hydroxy-1-\(hydroxymethyl\)ethyl ester|"
                                         r"2,3-dihydroxypropyl ester|"
                                         r"glyceryl|palmitoyl-glycerol|stearoyl-glycerol|"
                                         r"\bMonoglyceride\b|\bDiglyceride\b", re.I)),
    ("Fatty acid amide",     re.compile(r"Docosenamide|Hexadecanamide|Octadecanamide|"
                                         r"Tetradecanamide|Decanamide|Dodecanamide|"
                                         r"Oleamide|Erucamide|Palmitamide|Stearamide|"
                                         r"\b\w*anamide\b|\b\w*enamide\b", re.I)),
    ("Fatty acid methyl ester", re.compile(r"acid,\s*methyl ester|"
                                            r"\bMethyl\s+(?:stearate|palmitate|oleate|laurate|"
                                            r"myristate|linoleate|hexadec|octadec|tetradec|dodec)|"
                                            r"hexadecanoate|octadecanoate|tetradecanoate|"
                                            r"dodecanoate", re.I)),
    ("Fatty acid ethyl ester",  re.compile(r"acid,\s*ethyl ester|\bEthyl\s+\w+ate|"
                                            r"ethyl\s+(?:hexadec|octadec|tetradec|dodec)", re.I)),
    ("Hydroxy fatty acid",   re.compile(r"hydroxy(?:dodec|tetradec|hexadec|octadec)anoic|"
                                         r"hydroxy.*anoic acid|"
                                         r"\b\d+-hydroxy.*acid", re.I)),
    ("Fatty acid (free)",    re.compile(r"^\s*Hexadecanoic acid|^\s*n-Hexadecanoic|"
                                         r"^\s*Octadecanoic acid|^\s*Tetradecanoic acid|"
                                         r"^\s*Dodecanoic acid|^\s*Decanoic acid|"
                                         r"^\s*Octanoic acid|^\s*Nonanoic acid|"
                                         r"^\s*Tetracosanoic acid|^\s*Docosanoic acid|"
                                         r"^\s*Eicosanoic acid|^\s*Heptadecanoic acid|"
                                         r"\bPalmitic Acid|\bStearic acid|\bOleic Acid|"
                                         r"\bLinoleic acid|\bMyristic acid|\bLauric acid|"
                                         r"\bArachidic acid|\bBehenic acid|\bLignoceric acid|"
                                         r"\bCapric acid|\bCaprylic acid|\bCaproic acid|"
                                         r"\bButanoic acid,\s*\d+-methyl|"
                                         r"\bPropanoic acid,\s*\d+-methyl|"
                                         r"\bPentanoic acid|\bHexanoic acid", re.I)),
    ("Wax ester",            re.compile(r"docosyl ester|hexadecyl ester|octadecyl ester|"
                                         r"tetracosyl ester|wax", re.I)),
    ("Fatty alcohol",        re.compile(r"^1-(?:Octanol|Nonanol|Decanol|Undecanol|Dodecanol|"
                                         r"Tridecanol|Tetradecanol|Pentadecanol|Hexadecanol|"
                                         r"Heptadecanol|Octadecanol|Nonadecanol|Eicosanol|"
                                         r"Docosanol|Tetracosanol|Hexacosanol|Octacosanol)\b|"
                                         r"^(?:Octanol|Nonanol|Decanol|Undecanol|Dodecanol|"
                                         r"Tridecanol|Tetradecanol|Pentadecanol|Hexadecanol|"
                                         r"Heptadecanol|Octadecanol|Nonadecanol|Eicosanol|"
                                         r"Docosanol|Tetracosanol|Hexacosanol|Octacosanol)\b|"
                                         r"Behenic alcohol|cetyl alcohol|stearyl alcohol|"
                                         r"lauryl alcohol|myristyl alcohol|"
                                         r"\b\w*adecan-1-ol\b|\b\w*ocosan-1-ol\b|"
                                         r"\b\w*acosan-1-ol\b|"
                                         r"\bDecanol,\s*\d+-hexyl|\bDecanol,\s*\d+-methyl|"
                                         r"^\s*1-Decanol,\s*2-hexyl-", re.I)),
    ("Branched alkane",      re.compile(r"\b(?:Heptane|Octane|Nonane|Decane|Undecane|Dodecane|"
                                         r"Tridecane|Tetradecane|Pentadecane|Hexadecane|"
                                         r"Heptadecane|Octadecane|Nonadecane|Eicosane|"
                                         r"Heneicosane|Docosane|Tricosane|Tetracosane|"
                                         r"Pentacosane|Hexacosane|Heptacosane|Octacosane|"
                                         r"Nonacosane|Triacontane|Hentriacontane|Dotriacontane),"
                                         r"\s*[\d,]+-(?:methyl|ethyl|propyl|trimethyl|"
                                         r"dimethyl|tetramethyl|pentamethyl)|"
                                         r"\bisoprenoid|\bsqualane|\bpristane|\bphytane|"
                                         r"\bfarnesane", re.I)),
    ("Alkene",               re.compile(r"^\s*1-(?:Decene|Undecene|Dodecene|Tetradecene|"
                                         r"Hexadecene|Octadecene|Eicosene|Docosene)|"
                                         r"\b\w+ene,\s*\(?[EZ]\)|"
                                         r"\bTetradecene\b|\bTridecene\b", re.I)),
    ("n-Alkane",             re.compile(r"^\s*(?:Heptane|Octane|Nonane|Decane|Undecane|"
                                         r"Dodecane|Tridecane|Tetradecane|Pentadecane|"
                                         r"Hexadecane|Heptadecane|Octadecane|Nonadecane|"
                                         r"Eicosane|Heneicosane|Docosane|Tricosane|"
                                         r"Tetracosane|Pentacosane|Hexacosane|Heptacosane|"
                                         r"Octacosane|Nonacosane|Triacontane|Hentriacontane|"
                                         r"Dotriacontane|Tritriacontane|Tetratriacontane|"
                                         r"Pentatriacontane|Hexatriacontane|Tetratetracontane)"
                                         r"\s*$|^\s*n-\w+ane\s*$", re.I)),
    ("Cycloalkane",          re.compile(r"^\s*Cyclo(?:hexane|heptane|octane|nonane|decane|"
                                         r"undecane|dodecane|tridecane|tetradecane)|"
                                         r"\bCyclohexane,\s*\d", re.I)),
    ("Aromatic",             re.compile(r"\bbenzene\b|\btoluene\b|\bxylene\b|naphthal|"
                                         r"\bbiphenyl(?!yl)|\bIndan\b|\bIndan,", re.I)),
    ("Aldehyde",             re.compile(r"\b\w+anal\b|\baldehyde\b|carbaldehyde|carboxaldehyde", re.I)),
    ("Ketone",               re.compile(r"\b\w+anone\b|\b\w+enone\b|methylketone|"
                                         r"trimethyl-cyclohex-2-enone|2-Hydroxy-3,5,5-trimethyl|"
                                         r"oxaspiro.*dione", re.I)),
    ("Ether",                re.compile(r"\b\w+ether\b|^\s*Ether,|\boxide\b", re.I)),
    ("Thiophene/Sulfur",     re.compile(r"\bthiophene|\bsulfide|\bsulfonyl|disulfide|"
                                         r"methylthio|methanethiol", re.I)),
    ("Amino-acid derivative",re.compile(r"\bvalin|\bleucin|\bisoleucin|\bglycin|\balanin|"
                                         r"prolin|tyrosin|tryptophan|phenylalanin|"
                                         r"L-Valine|L-Leucine|L-Proline|L-Alanine", re.I)),
    ("Saccharide/Polyol",    re.compile(r"glucose|fructose|maltose|trehalose|xylose|ribose|"
                                         r"glycerol|mannitol|sorbitol|inositol", re.I)),
    ("Pyrimidine/Purine",    re.compile(r"pyrimidine|adenine|cytosine|thymine|uracil|guanine|"
                                         r"\buric acid\b|xanthine|hypoxanthine|barbital", re.I)),
    ("Quinoxaline",          re.compile(r"\bquinoxalin|quinazoline|quinoline|"
                                         r"\bQuinolinone\b", re.I)),
    ("Diol/Polyol",          re.compile(r"\bbutanediol\b|\bpropanediol\b|\bpentanediol\b|"
                                         r"\bhexanediol\b|\bethanediol\b", re.I)),
    ("Amine",                re.compile(r"\bdiethanolamine\b|\bethanolamine\b|"
                                         r"\bN-Ethyl|\bbutanamine\b|piperazine,?\s*N-|"
                                         r"\bpiperidine\b|\bpyrrolidine\b", re.I)),
    ("Ester (other)",        re.compile(r"\b\w+yl\s+\w+anoate\b|\b\w+yl\s+\w+enoate\b|"
                                         r"\b\w+yl\s+formate\b|\b\w+yl\s+acetate\b|"
                                         r"\bhexanoate\b|\bnonanoate\b|\bisopentyl|"
                                         r"\bbutyl\s+ester|\bpropyl\s+ester|\boctyl\s+ester", re.I)),
    ("Phenylpropanoid/Aromatic acid", re.compile(r"\bcinnamic|\bhydrocinnamic|\bbenzoic acid|"
                                                  r"\bbenzeneacetic|\bphenylacetic|"
                                                  r"\bphenylpropanoic|\bferulic|\bcaffeic", re.I)),
    ("Terpene",              re.compile(r"\bcamphor\b|\bcaryophyll|\bbicyclo\[7\.2\.0\]undec|"
                                         r"\bisoprene\b|\blimonen|\bpinene\b|\bbisabolen", re.I)),
]


# ---------------------------------------------------------------------------
# Light synonym normalisation applied to final names so synonymous library
# entries are counted as the same metabolite for diversity statistics.
# ---------------------------------------------------------------------------
SYNONYM_MAP: dict[str, str] = {
    "n-Hexadecanoic acid": "Hexadecanoic acid",
    "Palmitic Acid": "Palmitic acid",
    "Hexadecanoic acid, 2-hydroxy-1-(hydroxymethyl)ethyl ester": "1-Monopalmitin",
    "Octadecanoic acid, 2,3-dihydroxypropyl ester": "Glycerol monostearate",
}
