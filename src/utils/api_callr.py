"""
The two API caller functions below are copy-pasted from a rather obscure, hard-to-find source on RCSB website.
1. Go here: https://www.rcsb.org/docs/search-and-browse/advanced-search/attribute-search
2. Scroll to bottom: "C6. Find NMR structures of homomeric or heteromeric assemblies" in "Some Composite Query examples:"
3. Click on each of the two links:
    [[Run query for homomeric NMR structures]] uses attributes:
        QUERY: Experimental Method (Broader Categories) = "NMR" AND Number of Distinct Protein Entities = 1 AND
               Number of Protein Instances (Chains) per Assembly > 1
    [[Run query for heteromeric NMR structures]] uses attributes
        QUERY: Experimental Method (Broader Categories) = "NMR" AND Number of Distinct Protein Entities > 1
4. This shows the 'Advanced Search Query Builder' GUI above the actual proteins that it returns.
    Click on 'Search API' button in top right (with cogwheel icon). This takes you to the actual API JSON for one of
    the two search queries in 3. Click on 'copy query' button.

These JSONs are what you will find below in:
    `call_rcsb_for_pdbids_of_solution_nmr_homomeric(number)` and
    `call_rcsb_for_pdbids_of_solution_nmr_heteromeric(number)`

HENCE, THEY ARE NOT HAND-CRAFTED.
They are copy-pasted directly from RCSB site.
If you trust the GUI query search tool, then you should also trust these JSONs,
though admittedly they do seem rather long-winded and obscure.
"""
import requests


def call_rcsb_for_cif_or_pdb(pdb_id: str, cif_or_pdb='cif') -> requests.Response:
    """
    Send GET request to https://files.rcsb.org/download/{pdb_id} with PDB identifier of interest to retrieve the
    cif/pdb file of interest.
    :param pdb_id: Alphanumeric 4-character Protein Databank Identifier. e.g. '1OJ6'.
    :param cif_or_pdb: File extension of file desired. Can only be `cif` or `pdb`.
    :return: Response code 200 and text data for given PDB id, or error code such as 404.
    """
    response = None
    pdb_id = pdb_id.upper()  # MUST BE UPPER-CASE
    pdb_id = pdb_id.removesuffix('.cif')
    pdb_id = pdb_id.removesuffix('.pdb')
    url = f'https://files.rcsb.org/download/{pdb_id}.{cif_or_pdb}'
    try:
        print(f'Sending GET request to {url}')
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f'Failed to retrieve data from API: {e}')
    except Exception:
        print(f"Undefined error while trying to fetch {cif_or_pdb} file of '{pdb_id}' from PDB.")
    return response


def call_rcsb_for_pdbids_of_solution_nmr_homomeric(number: int) -> list:
    url = "https://search.rcsb.org/rcsbsearch/v2/query"

    query = {
      "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "exptl.method",
              "operator": "exact_match",
              "negation": False,
              "value": "SOLUTION NMR"
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_entry_info.polymer_entity_count_protein",
              "operator": "equals",
              "negation": False,
              "value": 1
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_assembly_info.polymer_entity_instance_count_protein",
              "operator": "greater",
              "negation": False,
              "value": 1
            }
          }
        ],
        "label": "text"
      },
      "return_type": "entry",
      "request_options": {
        "paginate": {
          "start": 0,
          "rows": number
        },
        "results_content_type": [
          "experimental"
        ],
        "sort": [
          {
            "sort_by": "score",
            "direction": "desc"
          }
        ],
        "scoring_strategy": "combined"
      }
    }

    response = requests.post(url, json=query)
    data = response.json()

    if "result_set" in data:
        pdb_ids = [result['identifier'] for result in data["result_set"]]
        print(f"Retrieved {len(pdb_ids)} PDB IDs")
        return pdb_ids
    else:
        print("No results or invalid query. Full response:")
        print(data)
        return []


def call_rcsb_for_pdbids_of_solution_nmr_heteromeric(number: int) -> list:
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    query = {
      "query": {
        "type": "group",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "exptl.method",
              "operator": "exact_match",
              "value": "SOLUTION NMR",
              "negation": False
            }
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "attribute": "rcsb_entry_info.polymer_entity_count_protein",
              "operator": "greater",
              "negation": False,
              "value": 1
            }
          }
        ],
        "logical_operator": "and",
        "label": "text"
      },
      "return_type": "entry",
      "request_options": {
        "paginate": {
          "start": 0,
          "rows": number
        },
        "results_content_type": [
          "experimental"
        ],
        "sort": [
          {
            "sort_by": "score",
            "direction": "desc"
          }
        ],
        "scoring_strategy": "combined"
      }
    }

    response = requests.post(url, json=query)
    data = response.json()

    if "result_set" in data:
        pdb_ids = [result['identifier'] for result in data["result_set"]]
        print(f"Retrieved {len(pdb_ids)} PDB IDs")
        return pdb_ids
    else:
        print("No results or invalid query. Full response:")
        print(data)
        return []