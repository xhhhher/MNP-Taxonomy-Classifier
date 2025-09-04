## 📥 Data Sources

The dataset was obtained from **CMNPD (Collection of Marine Natural Products Database)**:

- Website: [https://www.cmnpd.org/](https://www.cmnpd.org/)
- File format: **SDF (Structure Data File)**  
- Filters applied directly on the CMNPD website before downloading:
  - **Time period**: e.g., compounds reported *after the year 2000*  
  - **Biological source**: e.g., marine bacteria, fungi, animalia, etc.  

⚠️ The CMNPD website provides an interactive browser that allows filtering by **year of discovery/publication** and **source organism**. The SDF file can be exported directly after applying these filters.

## ⚙️ How to Prepare

1. Go to [https://www.cmnpd.org/](https://www.cmnpd.org/)  
2. Use the **search/filter interface**:
   - Select `Year > 2000`  
   - Select desired **organism group(s)** (e.g., Bacteria, Fungi, Animalia)  
3. Export the results as `.sdf` file(s).  
4. Place the downloaded `.sdf` files into the `data/raw/` folder, e.g.: data/raw/animalia.sdf (complete set before and after 2000)