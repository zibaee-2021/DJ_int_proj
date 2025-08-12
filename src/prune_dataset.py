"""
DEFINITIONS:

    "IDENTICAL SEQUENCES":
        - '100_IDENT_100_COV' = 100 % pident AND (100 % qcov OR 100 % tcov).

    "HOMOLOGOUS SEQUENCES":
        - '30_20_90' = evalue < 1e-3; pident >= 30.0; alnlen >= 20; qcov >= 0.9; tcov >= 0.9.

    'IDENTICAL STRUCTURES":
        - LT1RMSD_1TMS = RMSD < 1.0 OR TM-SCORE == 1.0.

======================================================================================================================
DATASET 1.0:

    - ALL HETEROMERIC PDB-CHAINS, BUT ONLY ONE (RANDOMLY-SELECTED) CHAIN FROM THE HOMOMERIC PDBS.
    - INCLUDES ONLY MULTIMODEL PDB-CHAINS.
    - INCLUDES ONLY PDB-CHAINS WITH 3 OR MORE RESIDUES.
    - ** DOES NOT REMOVE IDENTICAL PDB-CHAINS ** (I.E. 100 % pident AND 100% SAME STRUCTURE).

======================================================================================================================
DATASET 1.1:

     - AS DATASET 1.0, PLUS:

            - ANY ADDITIONAL MODIFICATIONS AS REQUESTED BY DAVID.
            - SUBSET OF SINGLE-MODEL PDB-CHAINS WHICH HAVE:
                - 30_20_90 HOMOLOGOUS TO OTHER PDB-CHAINS, BUT WITH NON-IDENTICAL STRUCTURES.
                    - WHEN COMPARING STRUCTURES OF SINGLE-MODEL TO MULTI-MODEL, MEAN COORDINATES OF MULTI-MODEL USED.

======================================================================================================================
DATASET 1.1.1:

    - AS DATASET 1.1, PLUS:

            - ANY ADDITIONAL MODIFICATIONS AS REQUESTED BY DAVID.
            - PRUNED TO REMOVE PDB-CHAINS THAT HAVE "IDENTICAL SEQUENCES" AND "IDENTICAL STRUCTURES":
                    - "IDENTICAL STRUCTURES" FOR MULTI-MODEL PDB-CHAINS, TO BE DONE VIA GENERATING
                      GAUSSIAN DISTRUBUTIONS OF MODEL COORDINATES AND CALCULATING SIMILARITY BY KL DIVERGENCE.

======================================================================================================================
DATASET 1.0.1:

    - AS DATASET 1.0, PLUS .. (SAME AS 1.1.1)

======================================================================================================================


-----------------------------------------------------------------------------------------------------------------
---------------- REMOVE PDBID-CHAINS SMALLER THAN 3 ALPHA-CARBONS: ----------------
-----------------------------------------------------------------------------------------------------------------
3 PIDC with only 1 CA:

    2AIZ_B has 1 CA
    1GAC_D has 1 CA
    1GAC_C has 1 CA

8 PIDC with 2 CAs:

    1GAC_A has 2 CA
    1GAC_B has 2 CA
    1WCO_A has 2 CA
    2K1Q_B has 2 CA
    2M9P_B has 2 CA
    2M9Q_B has 2 CA
    2MX6_B has 2 CA
    3CYS_B has 2 CA

3 PIDC with 3 CAs:

    2KID_B has 3 CA
    2RUI_B has 3 CA
    1CFA_B has 3 CA



-----------------------------------------------------------------------------------------------------------------
---------------- KEEP ONLY ONE PDB-CHAIN IF MORE THAN 1 HAVE SAME SEQUENCE *AND* SAME STRUCTURE: ----------------
-----------------------------------------------------------------------------------------------------------------

"SAME SEQUENCE": DEFINED AS 100 % PIDENT WITH 100 % COVERAGE:
"SAME STRUCTURE": DEFINED AS RMSD < 1.0 OR TM-SCORE = 1.0  (IF EITHER IS TRUE, THE STRUCTURES ARE ASSUMMED TO BE SAME.)

HOWEVER, FOR 'DATASET 1.0', I WILL NOT PRUNE ACCORDING TO THIS.

"""