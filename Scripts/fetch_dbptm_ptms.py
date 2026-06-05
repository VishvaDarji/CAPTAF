#!/usr/bin/env python3
"""
Fetch PTM data from dbPTM using the info.php endpoint.
Usage: python fetch_dbptm_ptms.py <fasta_file> <output_file>
"""

import sys
import os
import time
import pandas as pd
import requests
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options
from selenium.common.exceptions import TimeoutException, WebDriverException, NoSuchElementException

def get_dbptm_uniprot_id(fasta_file):
    """Extract dbPTM-compatible UniProt ID from FASTA header."""
    with open(fasta_file) as f:
        header = f.readline().strip()

    print(f"FASTA header: {header}")

    if header.startswith(">"):
        # Example: >sp|Q53ET0|CRTC2_HUMAN CREB-regulated transcription coactivator 2 OS=Homo sapiens...
        # We need CRTC2_HUMAN (the part after the second |)
        parts = header.split("|")
        if len(parts) >= 3:
            dbptm_id = parts[2].split()[0]  # Take first part before space
            print(f"Extracted dbPTM UniProt ID: {dbptm_id}")
            return dbptm_id
        else:
            # Fallback: try to extract from header
            header_clean = header[1:].split()[0]
            print(f"Fallback ID from header: {header_clean}")
            return header_clean

    # Last resort: use filename
    fallback_id = os.path.basename(fasta_file).replace(".fasta", "")
    print(f"Using filename-based ID: {fallback_id}")
    return fallback_id

def setup_chrome_driver():
    """Setup Chrome driver with optimized options."""
    chrome_options = Options()
    chrome_options.add_argument("--headless")
    chrome_options.add_argument("--no-sandbox")
    chrome_options.add_argument("--disable-dev-shm-usage")
    chrome_options.add_argument("--disable-gpu")
    chrome_options.add_argument("--window-size=1920,1080")
    chrome_options.add_argument("--user-agent=Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36")
    chrome_options.add_argument("--disable-extensions")
    chrome_options.add_argument("--disable-plugins")

    # Disable images and CSS for faster loading
    prefs = {
        "profile.managed_default_content_settings.images": 2,
        "profile.default_content_setting_values.stylesheets": 2
    }
    chrome_options.add_experimental_option("prefs", prefs)

    return webdriver.Chrome(options=chrome_options)

def fetch_dbptm_data(uniprot_id, output_file):
    """Fetch PTM data from dbPTM using the info.php endpoint."""
    driver = None
    try:
        print(f"Fetching dbPTM data for {uniprot_id}...")
        driver = setup_chrome_driver()

        # Use the new URL pattern
        url = f"https://biomics.lab.nycu.edu.tw/dbPTM/info.php?id={uniprot_id}"
        print(f"Loading URL: {url}")

        driver.get(url)

        # Wait for page to load
        WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.TAG_NAME, "body"))
        )

        # Wait a bit more for dynamic content
        time.sleep(3)

        # Check if we got an error page or no data
        page_source = driver.page_source.lower()
        if "not found" in page_source or "error" in page_source or "no data" in page_source:
            print(f"No data found for {uniprot_id} or error page returned")
            return False

        # Look for tabs - try to find and click "Experimental PTM sites" tab
        experimental_tab_clicked = False

        # Try different selectors for the experimental PTM tab
        tab_selectors = [
            "//a[contains(text(), 'Experimental PTM sites')]",
            "//a[contains(text(), 'Experimental')]",
            "//li[contains(text(), 'Experimental PTM sites')]",
            "//span[contains(text(), 'Experimental PTM sites')]",
            "//div[contains(text(), 'Experimental PTM sites')]"
        ]

        for selector in tab_selectors:
            try:
                experimental_tab = WebDriverWait(driver, 5).until(
                    EC.element_to_be_clickable((By.XPATH, selector))
                )
                experimental_tab.click()
                print("Successfully clicked 'Experimental PTM sites' tab")
                experimental_tab_clicked = True
                time.sleep(2)  # Wait for tab content to load
                break
            except TimeoutException:
                continue

        if not experimental_tab_clicked:
            print("Could not find or click 'Experimental PTM sites' tab, trying to extract from current page")

        # Wait for the table to load
        time.sleep(3)

        # Look for tables containing PTM data
        tables = driver.find_elements(By.TAG_NAME, "table")
        print(f"Found {len(tables)} tables on the page")

        ptm_data = []

        for i, table in enumerate(tables):
            try:
                # Get all rows
                rows = table.find_elements(By.TAG_NAME, "tr")
                print(f"Table {i+1}: {len(rows)} rows")

                if len(rows) < 2:  # Skip tables with just header or empty tables
                    continue

                # Check if this table has PTM-related headers
                header_row = rows[0]
                header_cells = header_row.find_elements(By.TAG_NAME, "th")
                if not header_cells:
                    header_cells = header_row.find_elements(By.TAG_NAME, "td")

                header_text = [cell.text.strip().lower() for cell in header_cells]
                print(f"Table {i+1} headers: {header_text}")

                # Check if this looks like a PTM table
                ptm_indicators = ['location', 'position', 'modification', 'ptm', 'site']
                if not any(indicator in ' '.join(header_text) for indicator in ptm_indicators):
                    print(f"Table {i+1} doesn't appear to contain PTM data")
                    continue

                print(f"Processing PTM table {i+1}")

                # Extract data from remaining rows
                for row_idx in range(1, len(rows)):
                    row = rows[row_idx]
                    cells = row.find_elements(By.TAG_NAME, "td")

                    if len(cells) >= 3:  # Need at least 3 columns (Location, Modification, Substrate)
                        location = cells[0].text.strip()
                        modification = cells[1].text.strip()
                        substrate = cells[2].text.strip()

                        # Extract position number from location (e.g., "123" from "S123" or "Ser123")
                        position = None
                        for char in location:
                            if char.isdigit():
                                position = ''.join([c for c in location if c.isdigit()])
                                break

                        if position and modification:  # Only add if we have valid data
                            ptm_data.append({
                                'Position': position,
                                'PTM_Type': modification,
                                'Substrate': substrate,
                                'Location_Raw': location
                            })

                # If we found data in this table, we can break
                if ptm_data:
                    print(f"Found {len(ptm_data)} PTM entries in table {i+1}")
                    break

            except Exception as e:
                print(f"Error processing table {i+1}: {e}")
                continue

        # Save the data
        if ptm_data:
            df = pd.DataFrame(ptm_data)
            # Add metadata
            df.insert(0, 'UniProt_ID', uniprot_id)
            df['Source'] = 'dbPTM'

            # Reorder columns to match expected format
            column_order = ['UniProt_ID', 'Position', 'PTM_Type', 'Source', 'Substrate', 'Location_Raw']
            df = df.reindex(columns=[col for col in column_order if col in df.columns])

            df.to_csv(output_file, index=False)
            print(f"Successfully saved {len(ptm_data)} PTM entries to {output_file}")

            # Show preview
            print("Preview of extracted data:")
            print(df.head())
            return True
        else:
            print("No PTM data found in any tables")
            return False

    except Exception as e:
        print(f"Error fetching dbPTM data: {e}")
        return False

    finally:
        if driver:
            try:
                driver.quit()
            except:
                pass

def main():
    if len(sys.argv) != 3:
        print("Usage: python fetch_dbptm_ptms.py <fasta_file> <output_file>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    output_file = sys.argv[2]

    # CHECK IF FILE ALREADY EXISTS
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        print(f"✅ Output file already exists: {output_file}")
        print(f"   Skipping dbPTM fetch (file was created manually or in previous run)")
        # Read and show contents
        try:
            df = pd.read_csv(output_file)
            print(f"   File contains {len(df)} PTM entries")
            if len(df) > 0:
                print(f"   Preview:")
                print(df.head())
            sys.exit(0)  # Success - file exists
        except Exception as e:
            print(f"   Warning: Could not read existing file: {e}")
            print(f"   Will re-fetch from dbPTM")

    # Create output directory
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    # Extract dbPTM-compatible UniProt ID
    uniprot_id = get_dbptm_uniprot_id(fasta_file)

    # Try to fetch data
    try:
        success = fetch_dbptm_data(uniprot_id, output_file)
        if not success:
            # ===== STOP PIPELINE WITH CLEAR ERROR =====
            print(f"\n{'='*60}")
            print(f"❌ ERROR: No PTMs found in dbPTM for target protein!")
            print(f"{'='*60}")
            print(f"Target FASTA: {fasta_file}")
            print(f"UniProt ID: {uniprot_id}")
            print(f"\nPossible reasons:")
            print(f"  1. Protein not in dbPTM database")
            print(f"  2. No experimentally validated PTMs")
            print(f"  3. Scraping/network error")
            print(f"\n🔧 SOLUTION for validation peptides with KNOWN PTM sites:")
            print(f"   Manually create the target PTM CSV file.")
            print(f"\n   Example for RAF1 peptide with phospho-Ser at position 5:")
            print(f"   cat > {output_file} << 'EOF'")
            print(f"   UniProt_ID,Position,PTM_Type,Source")
            print(f"   {uniprot_id},5,Phosphorylation,Literature")
            print(f"   EOF")
            print(f"\n   Then the pipeline will continue processing.")
            print(f"{'='*60}\n")
            sys.exit(1)  # Stop pipeline
            # ==========================================

    except Exception as e:
        print(f"\n{'='*60}")
        print(f"❌ ERROR: Exception while fetching dbPTM data!")
        print(f"{'='*60}")
        print(f"Error: {e}")
        print(f"Target: {fasta_file}")
        print(f"UniProt ID: {uniprot_id}")
        print(f"\nManually create PTM file to continue (see above for example).")
        print(f"{'='*60}\n")
        sys.exit(1)

if __name__ == "__main__":
    main()
