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
import pandas as pd

if __name__ == "__main__":

    # - INCLUDES ONLY PDB-CHAINS WITH 3 OR MORE RESIDUES:

    _1_AA_pidc = ['1GAC_C', '1GAC_D', '2AIZ_B']
    _2_AA_pidc = ['1GAC_A', '1GAC_B', '1WCO_A', '2K1Q_B', '2M9P_B', '2M9Q_B', '2MX6_B', '3CYS_B']
    _3_AA_pidc = ['1CFA_B', '2KID_B', '2RUI_B']

    with open(os.path.join('..','data','NMR','multimodel_lists', 'multimod_2713_hetallchains_hom1chain.lst'), 'r') as f:
        pidc_2713 = sorted(f.read().splitlines())
    rp_datasets_lists_dir = os.path.join('..','data','NMR', 'datasets', 'lists_specs')
    os.makedirs(rp_datasets_lists_dir, exist_ok=True)
    rp_dataset_1p0_lst = os.path.join(rp_datasets_lists_dir, '2702_pidc.lst')

    pidc_2702 = [pidc for pidc in pidc_2713 if pidc not in _1_AA_pidc + _2_AA_pidc]
    with open(rp_dataset_1p0_lst, 'w') as f:
        f.write('\n'.join(pidc_2702))

 # - EXCLUDES ANY THAT HAVE "EXTREME RMSD" AND "EXTREME TM-SCORE":

    rp_mm_2713_pidc = os.path.join('..', 'data', 'NMR', 'stats', 'multimod_2713_hetallchains_hom1chain', 'mm_2713_pidc.csv')
    mm_2713_pdf = pd.read_csv(rp_mm_2713_pidc)

    # 1. Apply size filter:
    mask_pidc = mm_2713_pdf['CA_count'] < 3
    # print(mask_pidc.sum())
    # print(mm_2713_pdf[mask_pidc])

    # 2. Apply structural filters:
    mask_structural = (
        (mm_2713_pdf['meanRMSD'] > 10) &
        (mm_2713_pdf['stdevRMSD'] < 0.01) &
        (mm_2713_pdf['meanTMS'] < 0.2) &
        (mm_2713_pdf['stdevTMS'] < 0.1)
    )
    # print(mask_structural.sum())
    # print(mm_2713_pdf[mask_structural])
    # print(mm_2713_pdf['meanTMS'].describe())
    # print(mm_2713_pdf['stdevTMS'].describe())
    # print(mm_2713_pdf['meanRMSD'].describe())
    # print(mm_2713_pdf['stdevRMSD'].describe())

    # 3. Combine filters:
    dataset_1p0 = mm_2713_pdf.loc[~mask_pidc]
    dataset_1p0 = dataset_1p0.loc[~mask_structural]
    print(dataset_1p0.shape)
    rp_datasets_dir = os.path.join('..', 'data', 'NMR', 'datasets')
    dataset_1p0.to_csv(os.path.join(rp_datasets_dir, 'NMR_1p0.csv'))

