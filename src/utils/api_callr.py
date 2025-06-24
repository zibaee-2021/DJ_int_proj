
import requests


def call_rcsb_for_cif(pdb_id: str) -> requests.Response:
    """
    Send GET request to https://files.rcsb.org/download/{pdb_id} with PDB identifier of interest.
    :param pdb_id: Alphanumeric 4-character Protein Databank Identifier. e.g. '1OJ6'.
    :return: Response code 200 and text data for given PDB id, or error code such as 404.
    """
    response = None
    pdb_id = pdb_id.upper()  # MUST BE UPPER-CASE
    pdb_id = pdb_id.removesuffix('.cif')
    url = f'https://files.rcsb.org/download/{pdb_id}.cif'
    try:
        print(f'Sending GET request to {url}')
        response = requests.get(url)
        response.raise_for_status()
    except requests.exceptions.RequestException as e:
        print(f'Failed to retrieve data from API: {e}')
    except Exception:
        print(f"Undefined error while trying to fetch '{pdb_id}' from PDB.")
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