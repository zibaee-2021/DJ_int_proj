### DEFINITIONS:
 
- "IDENTICAL PDB-CHAINS":
  - Both of following must be TRUE:
    - "IDENTICAL SEQUENCES":
      - Both of following must be TRUE:
        - 100 % `pident` 
        - **Either** 100 % `qcov` **or** 100 % `tcov`
    - "IDENTICAL STRUCTURES"
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
    - "EXTREME RMSD":
       - Both of following must be TRUE:
           - mean RMSD > 10
           - stdev RMSD < 0.01
    - "EXTREME TM-SCORE":
       - Both of following must be TRUE:
           - mean TM-score < 0.2
           - stdev TM-score < 0.1 

---

### DATASET 0.9:
| SELECTION CRITERIA                                                 | PDB-CHAINS DROPPED |
|--------------------------------------------------------------------|--------------------|
| INCLUDE ALL HETEROMERIC PDB-CHAINS                                 | 0                  |
| INCLUDE ONLY ONE (RANDOMLY-SELECTED) CHAIN FROM THE HOMOMERIC PDBS | 812                |
| INCLUDES ONLY MULTIMODEL PDB-CHAINS                                | 468                |
v0.9 has 2713 PDB-chains

---

### DATASET 1.0:
SAME AS DATASET 0.9, **PLUS**: 

| SELECTION CRITERIA                                                | PDB-CHAINS DROPPED  |
|-------------------------------------------------------------------|---------------------|
| ANY FURTHER MODIFICATIONS/ADDITIONS REQUESTED BY DJ               | TBC                 |
| EXCLUDE PDB-CHAINS WITH < 3 ALPHA-CARBONS                         | 11                  |
| EXCLUDE PDB-CHAINS WITH "EXTREME RMSD" **AND** "EXTREME TM-SCORE" | 0                   |

v1.0 has **(tbc)** PDB-chains

NOTE: I have not yet (looked for and) removed "IDENTICAL PDB-CHAINS" from v0.9. <br>
Furthermore, the presence of PDB records with identical sequence-structures but under different PDB identifiers 
cannot not be completely ruled out from this dataset.

---

### DATASET 1.1:
SAME AS DATASET 1.0, **PLUS**:

| SELECTION CRITERIA                                                                                                          | PDB-CHAINS DROPPED |
|-----------------------------------------------------------------------------------------------------------------------------|--------------------|
| ANY FURTHER MODIFICATIONS/ADDITIONS REQUESTED BY DJ                                                                         | TBC                |
| EXCLUDE "IDENTICAL PDB-CHAINS"                                                                                              | TBC                |
| EXCLUDE PDB-CHAINS WITH "EXTREME RMSD" **AND** "EXTREME TM-SCORE"<br> WITH MORE STRINGENT THRESHOLDS THAN THAT USED IN v1.0 | TBC                | 

v1.1 has **(tbc)** PDB-chains

NOTE: The method/algorithm for quantifying "IDENTICAL STRUCTURES" when dealing with multi-model PDB-chains has 
not been decided on yet. I may try something more sophisticated than simply mean coordinates though.
<br>(E.g. Gaussian distributions of model coordinates and/or contact maps, followed by calculating similarity via 
KL divergence. I will also look at what's commonly done in Molecular Dynamics field, and/or via MD libraries/tools like MDTraj, MDAnalysis, etc.)

---

### DATASET 1.2:
SAME AS DATASET 1.1, **PLUS**:

| SELECTION CRITERIA                                                                                                      | PDB-CHAINS DROPPED |
|-------------------------------------------------------------------------------------------------------------------------|--------------------|
| ANY FURTHER MODIFICATIONS/ADDITIONS REQUESTED BY DJ                                                                     | TBC                |
| INCLUDE SUBSET OF SINGLE-MODEL PDB-CHAINS WHICH ARE "HOMOLOGOUS" TO OTHER PDB-CHAINS, BUT HAVE NON-IDENTICAL STRUCTURES | TBC                |

v1.2 has **(tbc)** PDB-chains

NOTE: When comparing structures of single-model to multi-model PDB-chains, I will use mean coordinates of the multi-model PDB-chains.