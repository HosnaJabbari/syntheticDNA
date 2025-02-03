from selenium import webdriver
from selenium.webdriver.common.by import By
import time
# scraping ecRBPome database for pdb IDs, as there's no output for PDB IDs
driver = webdriver.Chrome()
url = "https://caps.ncbs.res.in/cgi-bin/ecrbpome/uniprot_ID_PDB_link.pl"
driver.get(url)
time.sleep(5) # to load

# click sort twice on "Cross Reference (PDB)" to rows that has IDs come first
pdb_column_header = driver.find_element(By.XPATH, "//th[contains(., 'Cross reference (PDB)')]")
pdb_column_header.click()
time.sleep(3) #to load
pdb_column_header.click()
time.sleep(3)
# scrape
pdb_ids = []
while True:
    # rows from table
    rows = driver.find_elements(By.CSS_SELECTOR, "table#uniprot_ID_PDB_link tbody tr")
    no_table = True
    for row in rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        if len(cells) > 4:  # to make sure row has enough cells
            pdb_links = cells[4].find_elements(By.TAG_NAME, "a")  #5th column PDB
            if pdb_links:
                no_table = False  # non-empty row
                for link in pdb_links:
                    pdb_id = link.get_attribute("href").split('=')[-1]
                    pdb_ids.append(pdb_id)
    if no_table:
        break
    try:#next page
        next_button = driver.find_element(By.LINK_TEXT, "Next")
        if "disabled" in next_button.get_attribute("class"):  # check if disabled
            break
        next_button.click()
        time.sleep(3)  #next page to load
    except Exception as e:
        print("no more pages or error,", e)
        break

with open("ecRBPome_pdb_ids.txt", "w") as file:
    for pdb_id in pdb_ids:
        file.write(f"{pdb_id}\n")
print(f"extracted {len(pdb_ids)} PDB IDs")
driver.quit()