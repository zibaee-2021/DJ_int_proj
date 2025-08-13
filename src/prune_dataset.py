"""-----------------------------------------------------------------------------------------------------------------
---------------- REMOVE PDBID-CHAINS SMALLER THAN 3 ALPHA-CARBONS: ----------------
-----------------------------------------------------------------------------------------------------------------
3 PIDC with only 1 CA:

    1GAC_C has 1 CA
    1GAC_D has 1 CA
    2AIZ_B has 1 CA

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

    1CFA_B has 3 CA
    2KID_B has 3 CA
    2RUI_B has 3 CA
-----------------------------------------------------------------------------------------------------------------
---------------- KEEP ONLY ONE PDB-CHAIN IF MORE THAN 1 HAVE SAME SEQUENCE *AND* SAME STRUCTURE: ----------------
-----------------------------------------------------------------------------------------------------------------

"SAME SEQUENCE": DEFINED AS 100 % PIDENT WITH 100 % COVERAGE:
"SAME STRUCTURE": DEFINED AS RMSD < 1.0 OR TM-SCORE = 1.0  (IF EITHER IS TRUE, THE STRUCTURES ARE ASSUMMED TO BE SAME.)

HOWEVER, FOR 'DATASET 1.0', I WILL NOT PRUNE ACCORDING TO THIS.

"""
import os

if __name__ == "__main__":
    with open(os.path.join('..','data','NMR','multimodel_lists', 'multimod_2713_hetallchains_hom1chain.lst'), 'r') as f:
        pidc_2713 = f.read().splitlines()

    pidc_2713.sort()

    rp_datasets_dir = os.path.join('..','data','NMR', 'datasets')
    os.makedirs(rp_datasets_dir, exist_ok=True)
    pidc_2702 = []
    for pidc in pidc_2713:
        if pidc in ['2AIZ_B', '1GAC_D', '1GAC_C', '1GAC_A', '1GAC_B', '1WCO_A',
                   '2K1Q_B', '2M9P_B', '2M9Q_B', '2MX6_B', '3CYS_B']:
            continue
        else:
            pidc_2702.append(pidc)

