import re
import requests
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
import zipfile
import os, sys

# Create output directory if it does not exist
output_dir = "./afd_output/"
os.makedirs(output_dir, exist_ok=True)

# Define the search terms and precision terms
search_terms = [
    "GFS", "GFS ENSEMBLE", "GEFS", "ECMWF", "ECMWF ENSEMBLE", 
    "EPS", "HRRR", "HREF", "ENSEMBLE", "EFI", "NBM", 
    "NATIONAL BLEND", "CLUSTER"
]

precision_terms = [
    "PERCENTILE", "PROBABILITY", "POSSIBLE", "EXPECTED", 
    "CHANCE", "LIKELY", "CHANCE (", "% CHANCE"
]

# List of WFO sites by region
region_dict ={
    "WR":["BYZ", "BOI", "LKN", "EKA", "FGZ", "GGW", "TFX", "VEF", "LOX", "MFR",
        "MSO", "PDT", "PSR", "PIH", "PQR", "REV", "STO", "SLC", "SGX", "MTR",
        "HNX", "SEW", "OTX", "TWC"],

    "CR":["ABR", "BIS", "CYS", "LOT", "DVN", "BOU", "DMX", "DTX", "DDC", "DLH",
        "FGF", "GLD", "GJT", "GRR", "GRB", "GID", "IND", "JKL", "EAX", "ARX",
        "ILX", "LMK", "MQT", "MKX", "MPX", "LBF", "APX", "IWX", "OAX", "PAH",
        "PUB", "UNR", "RIW", "FSD", "SGF", "LSX", "TOP", "ICT"],

    "ER":["ALY", "LWX", "BGM", "BOX", "BUF", "BTV", "CAR", "CTP", "RLX", "CHS",
        "ILN", "CLE", "CAE", "GSP", "MHX", "OKX", "PHI", "PBZ", "GYX", "RAH",
        "RNK", "AKQ", "ILM"],

    "SR":["ABQ", "AMA", "FFC", "EWX", "BMX", "BRO", "CRP", "EPZ", "FWD", "HGX",
        "HUN", "JAN", "JAX", "KEY", "MRX", "LCH", "LZK", "LUB", "MLB", "MEG",
        "MAF", "MFL", "MOB", "MRX", "OHX", "LIX", "OUN", "SJT", "SHV", "TAE",
        "TBW", "TSA"]
}

# Function to download data for a given WFO and year
def download_data(wfo, year):
    start_date = f"{year}-01-01T00:00Z"
    end_date = f"{year+1}-01-01T00:00Z"
    zip_filename = f"{output_dir}data_{wfo}_{year}.zip"
    extract_dir = f"{output_dir}data_{wfo}_{year}"

    # Check if the extracted directory already exists
    if os.path.exists(extract_dir):
        return extract_dir

    # Check if the zip file is already downloaded
    if os.path.exists(zip_filename):
        with zipfile.ZipFile(zip_filename, 'r') as z:
            z.extractall(extract_dir)
        return extract_dir

    url = f"https://mesonet.agron.iastate.edu/cgi-bin/afos/retrieve.py?sdate={start_date}&edate={end_date}&pil=AFD{wfo}&fmt=zip&limit=99999"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        with open(zip_filename, 'wb') as f:
            f.write(response.content)
        
        with zipfile.ZipFile(zip_filename, 'r') as z:
            z.extractall(extract_dir)
    
    except requests.exceptions.RequestException as e:
        return None

    return extract_dir

# Function to process the downloaded data
def process_downloaded_data(directory):
    search_term_counts = {term: 0 for term in search_terms}
    precision_term_counts = {term: 0 for term in precision_terms}
    file_count = 0
    
    if not os.path.exists(directory):
        return search_term_counts, precision_term_counts, file_count
    
    for filename in os.listdir(directory):
        if filename.endswith(".txt"):
            file_count += 1
            with open(os.path.join(directory, filename), 'r') as file:
                text_content = file.read()

                # Remove any extraneous formatting (e.g., extra spaces, newlines)
                text_content = re.sub(r'\s+', ' ', text_content).strip()

                # Count occurrences of each term using regex with word boundaries
                for term in search_terms:
                    search_term_counts[term] += len(re.findall(r'\b' + re.escape(term) + r'\b', text_content, re.IGNORECASE))

                for term in precision_terms:
                    precision_term_counts[term] += len(re.findall(r'\b' + re.escape(term) + r'\b', text_content, re.IGNORECASE))
    
    return search_term_counts, precision_term_counts, file_count

# Function to aggregate results for multiple years in parallel
def aggregate_results_over_years(wfo, start_year, end_year):
    results = {}
    missing_years = []

    def process_year(year):
        extract_dir = download_data(wfo, year)
        if extract_dir is None:
            return year, None, None, 0
        search_counts, precision_counts, file_count = process_downloaded_data(extract_dir)
        return year, search_counts, precision_counts, file_count
    
    years = range(start_year, end_year + 1)
    with ThreadPoolExecutor(max_workers=30) as executor:
        futures = [executor.submit(process_year, year) for year in years]
        
        with tqdm(total=len(years), desc=f"Processing {wfo}", leave=False) as pbar:
            for future in as_completed(futures):
                year, search_counts, precision_counts, file_count = future.result()
                if search_counts is None or precision_counts is None:
                    missing_years.append(year)
                    pbar.update(1)
                    continue
                
                results[year] = {**search_counts, **precision_counts, 'AFD_Count': file_count}
                pbar.update(1)
    
    return wfo, results, missing_years

# Function to process each WFO in the selected region
def process_region(region):
    all_results = {}
    missing_data = {}

    WFOs = region_dict.get(region, [])

    for wfo in tqdm(WFOs, desc=f"Processing region {region}"):
        csv_file = f"{output_dir}{wfo}_term_counts.csv"
        if os.path.exists(csv_file):
            df_wfo = pd.read_csv(csv_file, index_col='Year')
            all_results[wfo] = df_wfo.to_dict('index')
            continue

        wfo, results, missing_years = aggregate_results_over_years(wfo, 2000, 2025)
        all_results[wfo] = results
        if missing_years:
            missing_data[wfo] = missing_years

        # Save each WFO's output counts to its own CSV file
        df_wfo = pd.DataFrame.from_dict(results, orient='index')
        df_wfo.index.name = 'Year'
        df_wfo.sort_index(inplace=True)
        df_wfo.to_csv(csv_file)

    # Create a combined DataFrame with a multi-index [WFO, Year]
    combined_results = []
    for wfo, yearly_data in all_results.items():
        for year, data in yearly_data.items():
            combined_results.append((wfo, year, *data.values()))

    columns = ['WFO', 'Year'] + list(search_terms) + list(precision_terms) + ['AFD_Count']
    df_combined = pd.DataFrame(combined_results, columns=columns)
    df_combined.set_index(['WFO', 'Year'], inplace=True)
    df_combined.sort_index(inplace=True)

    # Save the combined DataFrame to a CSV file
    df_combined.to_csv(f"{output_dir}combined_term_counts_{region}.csv")

    # Print the combined DataFrame
    print(f"Combined Term Counts for region {region}:")
    print(df_combined)

    # Report missing date ranges
    print(f"\nMissing Data Ranges for region {region}:")
    for wfo, years in missing_data.items():
        print(f"{wfo}: {', '.join(map(str, years))}")


if __name__ == "__main__":
    process_region(sys.argv[1])