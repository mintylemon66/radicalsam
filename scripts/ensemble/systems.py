"""Radical SAM systems used to train the ensemble, with the partner family that counts as a positive.

Each system is one radical SAM anchor enzyme and a protein family known to work with it. A candidate
family scored against that anchor is a positive when it is that partner and a negative otherwise.

The futalosine system is deliberately absent from TRAIN. The model is applied to MqnE and DUF4229 after
training, so putting futalosine anywhere in the training set would make that application circular.

Every accession below was checked at InterPro and against primary literature on 17 October 2026.
Anchors are NCBIfam rather than Pfam because every radical SAM enzyme carries PF04055, so a Pfam anchor
cannot name a specific enzyme.
"""

TRAIN = {
    "biotin": {
        "anchor": "ncbifam-TIGR00433", "anchor_name": "BioB",
        "partner": "pfam-PF13500", "partner_name": "BioD",
    },
    "lipoyl": {
        "anchor": "ncbifam-TIGR00510", "anchor_name": "LipA",
        "partner": "pfam-PF21948", "partner_name": "LipB",
    },
    "pqq": {
        "anchor": "ncbifam-TIGR02109", "anchor_name": "PqqE",
        "partner": "pfam-PF05402", "partner_name": "PqqD",
    },
    "mycofactocin": {
        "anchor": "ncbifam-TIGR03962", "anchor_name": "MftC",
        "partner": "pfam-PF26520", "partner_name": "MftB",
    },
    "molybdopterin": {
        "anchor": "ncbifam-TIGR02666", "anchor_name": "MoaA",
        "partner": "pfam-PF01967", "partner_name": "MoaC",
    },
    "thiamine": {
        "anchor": "ncbifam-TIGR02351", "anchor_name": "ThiH",
        "partner": "pfam-PF05690", "partner_name": "ThiG",
    },
    "pfl": {
        "anchor": "ncbifam-TIGR02493", "anchor_name": "PflA",
        "partner": "ncbifam-TIGR01255", "partner_name": "PflB",
    },
    "queuosine": {
        "anchor": "ncbifam-TIGR03365", "anchor_name": "QueE",
        "partner": "pfam-PF06508", "partner_name": "QueC",
    },
}

# scored only after the model is fit, never trained on
HELDOUT = {
    "futalosine": {
        "anchor": "ncbifam-TIGR03700", "anchor_name": "MqnE",
        "partner": "pfam-PF02621", "partner_name": "MqnA/MqnD",
        "candidate_of_interest": "pfam-PF14012", "candidate_name": "DUF4229",
    },
}

# Checked and rejected, with the reason, so nobody adds them back later:
#   anSME (TIGR03942 / PF00884)   one maturase serves several sulfatases and they are usually not
#                                 genomically linked, so the co-location signal is absent by construction
#   RimO (TIGR01125 / PF00164)    the partner is the substrate, and S12 is a ribosomal protein present in
#                                 nearly every genome, so neighborhood co-occurrence carries no information
#   FO synthase (CofG / CofH)     fused into one polypeptide in 72 to 81 percent of genomes, which leaves
#                                 no partner gene and makes gene distance undefined
#   hydrogenase (HydE/HydF)       dropped after the pathway mapping: hydE and hydF have no reviewed
#                                 Swiss-Prot entries, and the only gene that resolved, hydG, maps to
#                                 PF25601 "AAA+ ATPase lid domain", a generic structural domain
REJECTED = ("anSME", "RimO", "FO synthase", "hydrogenase")
