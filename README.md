# KEGG-web-scraping

This repository contains a Python-based web scraping tool that extracts Kegg Orthologs from the [KEGG](https://www.kegg.jp/) (Kyoto Encyclopedia of Genes and Genomes) database. Check out the `environment.yml` file to manage dependencies.

---

## 📖 Overview

KEGG is a widely used resource for understanding high-level functions of biological systems. This tool automates the retrieval of protein data from KEGG, which can be used for bioinformatics analysis, research projects, or database integration.

The scraper navigates KEGG's database, extracts detailed protein information, and saves it in a structured format (e.g., CSV or JSON) for downstream analysis.

---

## 🚀 Features

- Scrapes protein information from KEGG pathways and genes.
- Saves output data in a structured format.
- Modular, easy-to-extend codebase.
- Environment management via `conda` for reproducibility.

---

## 🛠️ Installation

Clone the repository:


```bash
git clone https://github.com/basso42/kegg-web-scraping.git
cd kegg-web-scraping
```

## Project structure

```bash
├── environment.yml
├── scrape_kegg_proteins.py
├── utils.py
├── output/
│   └── proteins_data.csv
└── README.md
```
