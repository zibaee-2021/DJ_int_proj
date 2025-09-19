import os
import re
import requests
from bs4 import BeautifulSoup
import pandas as pd
from collections import OrderedDict

# URL_usercreated = 'https://dyndom.cmp.uea.ac.uk/dyndom/browse.do?type=WEBSITE&id=27d2688ea7053593cd8ca80c002b11'
URL_nonredundant = 'https://dyndom.cmp.uea.ac.uk/dyndom/browse.do?type=NR&id=2dad12c84f009d148c9be288ca8014'
# URL = URL_usercreated
URL = URL_nonredundant


if __name__ == '__main__':

    resp = requests.get(URL, timeout=30)
    resp.raise_for_status()

    xray_dyndom_dir = os.path.join('..', 'data', 'XRAY', 'DynDom')
    os.makedirs(xray_dyndom_dir, exist_ok=True)

    with open(os.path.join(xray_dyndom_dir, 'dyndom_non_redundant.html'), 'w', encoding='utf-8') as f:
    # with open(os.path.join(xray_dyndom_dir, 'dyndom_user_created.html'), 'w', encoding='utf-8') as f:
        f.write(resp.text)

    html = resp.text

    soup = BeautifulSoup(html, 'html.parser')

    # regex extract PDB-like 4-char IDs (start with a digit, then 3 [A-Za-z0-9])
    pdb_id_pattern = re.compile(r'\b[0-9][A-Za-z0-9]{3}\b')

    rows = []

    for tr in soup.find_all('tr'):
        name_td = tr.find('td', class_='ddpagecol1', align='center')
        conf_td = tr.find('td', class_='ddpagecolumn', align='center')
        if not (name_td and conf_td):
            continue  # skip header rows and non-data rows

        protein_name = name_td.get_text(strip=True)
        if not protein_name:
            continue

        conf_text = conf_td.get_text(separator='\n', strip=True)
        ids_found = pdb_id_pattern.findall(conf_text)

        unique_ids = list(OrderedDict.fromkeys(ids_found))

        family_link = tr.select_one('a[href*="Subgroup.do"]')
        family_value = family_link.get_text(strip=True) if family_link else None

        rows.append({
            'Protein name': protein_name,
            'Conformers': unique_ids,  # stored as list, e.g. ['1s3i', '2cfi']
            'Family': family_value,
        })

    user_created_dyndom = pd.DataFrame(rows, columns=['Protein name', 'Conformers', 'Family'])
    # user_created_dyndom.attrs['title'] = 'User-created DynDom'  # optional metadata
    user_created_dyndom.attrs['title'] = 'Non-redundant DynDom'  # optional metadata

    print(user_created_dyndom.head())

    # user_created_dyndom.to_csv(os.path.join(xray_dyndom_dir, 'user_created_dyndom.csv'), index=False)
    user_created_dyndom.to_csv(os.path.join(xray_dyndom_dir, 'non_redundant_dyndom.csv'), index=False)
