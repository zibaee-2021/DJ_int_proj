- <details>
  <summary><strong>data/</strong></summary>

  - <details>
    <summary><strong>DynDom/</strong></summary>

    - 2 `.html` files & 2 `.csv` files  
      Tables copied from DynDom webpages and scraped to CSV.

    - <details>
      <summary><strong>Hayward_files/</strong></summary>

      - <details>
        <summary><strong>foldseek/</strong></summary>

        - <details>
          <summary><strong>subdirs, 1 per PDB (e.g. 2DN1/)</strong></summary>
          
          - <details>
            <summary><strong>21 Foldseek-related files</strong></summary>

            - <details>
              <summary><strong>14 files output by Foldseek `createdb`</strong></summary>
              
                ```
                 2DN1: Binary database file containing the amino acid sequences (body).
                 2DN1.dbtype: Small file telling FoldSeek what type of file this database is.
                 2DN1.index: Index file mapping internal IDs to byte offsets in 2DN1, so entries can be fetched quickly.
                 2DN1.lookup: Maps internal numeric IDs back to the original identifiers (PDB chain IDs, etc.).
                 2DN1.source: Keeps track of where each entry came from (the input file paths, e.g. nameofpdb.pdb.
                 2DN1_ca: Database of CA-only coordinates extracted from structure. Used for fast structural prefiltering.
                 2DN1_ca.dbtype: Type marker for the CA database.
                 2DN1_ca.index: Index mapping IDs to the positions of the coordinate records inside 2DN1_ca
                 2DN1_h: Database containing headers/identifiers (what would appear after > in a FASTA, i.e. separate from amino acid sequence 'body'.)
                 2DN1_h.dbtype: Type marker: tells Foldseek this DB holds headers, not bodies.
                 2DN1_h.index: Index file for fast lookups in the header database.
                 2DN1_ss: Database containing 3Di alphabet sequence for each residue; symbolic encoding of local 3D geometry.
                 2DN1_ss.dbtype: Type marker for the 3Di database.
                 2DN1_ss.index: Index for fast access to the 3Di database.
                ```

</details>
  - <details>
    <summary><strong>NMR/</strong></summary>
</details>
</details>
</details>
</details>
</details>
</details>
</details>