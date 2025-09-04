## 📥 Data Sources

### 1. CMNPD

- Website: [https://www.cmnpd.org/](https://www.cmnpd.org/)  
- File format: **TSV (Tab-Separated Values)**  
- Filters applied directly on the CMNPD website before downloading:
  - **Time period**: e.g., compounds reported *after the year 2000*  
  - **Biological source**: e.g., marine bacteria, fungi, animalia, etc.  

⚠️ The CMNPD website provides an interactive browser that allows filtering by **year of discovery/publication** and **source organism**.  
The results can be directly exported as `.tsv` files.

**How to Prepare (CMNPD):**
1. Go to [https://www.cmnpd.org/](https://www.cmnpd.org/)  
2. Use the **search/filter interface**:
   - Select `Year > 2000`  
   - Select desired **organism group(s)** (e.g., Bacteria, Fungi, Animalia)  
3. Export the results as `.tsv` file(s).  
4. Place the downloaded `.tsv` files into the `data/raw/` folder, e.g.: data/raw/CMNPD-animalia-after2000.tsv data/raw/CMNPD-bacteria-after2000.tsv data/raw/CMNPD-fungi-after2000.tsv
5. Use the preprocessing script or notebook to read and merge the `.tsv` files.


### 2. NPAtlas

- Website: [https://www.npatlas.org/](https://www.npatlas.org/)
- File format:  **TSV (Tab-Separated Values)**
- NPAtlas provides curated natural products with compound IDs, SMILES, and organism information.

⚠️ The full dataset can be directly exported as a .tsv file.

**How to Prepare (NPAtlas):**
1. Go to https://www.npatlas.org/
2. Download the latest NPAtlas dataset as .tsv.
3. Place the file into the data/raw/ folder, e.g.: data/raw/NPAtlas_download_2024_09.tsv
4. Use the preprocessing script or notebook to process the `.tsv` file.


## ⚙️Processed Output
- CMNPD processed files → data/processed/data_cmnpd_after2000.csv
- NPAtlas processed file → data/processed/npatlas.csv