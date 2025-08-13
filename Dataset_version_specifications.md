============
DEFINITIONS:
============
                                                    
    - "IDENTICAL PDB-CHAINS":
        - Both of following must be TRUE:
            - "IDENTICAL SEQUENCES"
            - "IDENTICAL STRUCTURES"

    - "IDENTICAL SEQUENCES"
        - Both of following must be TRUE:
            - 100 % pident
            - 100 % qcov - OR - 100 % tcov

    - "IDENTICAL STRUCTURES":
        - Either of following must be TRUE:
            - mean RMSD < 1.0
            - mean TM-score == 1.0

    - "HOMOLOGOUS":
        - All of following must be TRUE:
            - evalue < 1e-3 
            - pident >= 30.0 
            - alnlen >= 20 
            - qcov >= 0.9
            - tcov >= 0.9    

     - EXTREME RMSD:
        - Both of following must be TRUE:
            - mean RMSD > 10 ANGSTROMS
            - stdev RMSD < 0.01
    
    - EXTREME TM-SCORE:
        - Both of following must be TRUE:
            - mean TM-score < 0.2
            - stdev TM-score < 0.1 

============
DATASET 1.0:
============
                                                                                    # PDB-CHAINS DROPPED
    - INCLUDE ALL HETEROMERIC PDB-CHAINS:                                           0
    - INCLUDE ONLY ONE (RANDOMLY-SELECTED) CHAIN FROM THE HOMOMERIC PDBS:           812
    - INCLUDES ONLY MULTIMODEL PDB-CHAINS:                                          468
    - INCLUDES ONLY PDB-CHAINS WITH 3 OR MORE RESIDUES:                             11
    - EXCLUDES ANY THAT HAVE "EXTREME RMSD" AND "EXTREME TM-SCORE":                 0

    NOTE: I HAVE NOT LOOKED FOR OR REMOVED "IDENTICAL PDB-CHAINS" YET, 
          SO, RECORDS WITH IDENTICAL SEQUENCE-STRUCTURES BUT UNDER DIFFERENT PDB IDENTIFIERS CANNOT BE RULED OUT.

======================================================================================================================
DATASET 1.1:

    - SAME AS DATASET 1.0, PLUS:

        - ANY MODIFICATIONS/ADDITIONS REQUESTED BY DAVID.
        - REMOVE PDB-CHAINS THAT HAVE "IDENTICAL PDB-CHAINS":
                - METHOD/AGORITHM FOR QUANTIFYING "IDENTICAL STRUCTURES" WHEN DEALING WITH MULTI-MODEL PDB-CHAINS 
                  HAS NOT BEEN CHOSEN YET, BUT SOMETHING MORE SOPHISTICATED THAN PURELY COMPARING MEAN COORDINATES.
                  (E.G. GAUSSIAN DISTRUBUTIONS OF MODEL COORDINATES AND CALCULATING SIMILARITY BY KL DIVERGENCE?)

======================================================================================================================
DATASET 1.2:

    - SAME AS DATASET 1.1, PLUS:

        - ANY ADDITIONAL MODIFICATIONS AS REQUESTED BY DAVID.
        - SUBSET OF SINGLE-MODEL PDB-CHAINS WHICH ARE "HOMOLOGOUS" TO OTHER PDB-CHAINS, 
          BUT HAVE NON-IDENTICAL STRUCTURES.
        - WHEN COMPARING STRUCTURES OF SINGLE-MODEL TO MULTI-MODEL, MEAN COORDINATES OF MULTI-MODELS TO BE USED.

======================================================================================================================
