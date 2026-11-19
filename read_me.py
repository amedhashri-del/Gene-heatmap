# Gene Co-expression Heatmap Generator

## About the Project

The Gene Co-expression Heatmap Generator is a simple Python application that analyzes gene expression data and displays the correlation between different genes using a heatmap.

It uses Pearson correlation to identify how similarly genes behave across different samples. The results are represented using different colours, making gene expression patterns easier to understand.

## Aim

To develop a Python-based application that calculates pairwise Pearson correlation between genes and visualizes the results through a colour-coded heatmap.

## Features

- Load gene expression data from a CSV file.
- Calculate Pearson correlation between genes.
- Generate a colour-coded heatmap.
- Display correlation values inside the heatmap.
- Save the correlation matrix as a CSV file.
- Includes sample gene expression data.

## Technologies Used

- Python
- Tkinter
- CSV
- Math

## How It Works

1. The user loads a gene expression CSV file.
2. The application reads the expression values for each gene.
3. Pearson correlation is calculated between every pair of genes.
4. The correlation values are displayed as a heatmap.
5. The user can save the correlation matrix as a CSV file.

## Understanding the Heatmap

- **Red:** Positive correlation between genes.
- **White:** Correlation close to zero.
- **Blue:** Negative correlation between genes.

Pearson correlation values range from -1 to +1.

## Input Format

The input file should be a CSV containing gene names in the first column and numerical expression values in the remaining columns.

Example:

| Gene | Sample_1 | Sample_2 | Sample_3 |
|---|---:|---:|---:|
| BRCA1 | 12.5 | 15.2 | 14.0 |
| BRCA2 | 11.8 | 14.9 | 13.8 |
| TP53 | 80.1 | 65.4 | 50.2 |

## How to Run

1. Download or clone this repository.
2. Open the project folder in Terminal.
3. Run the following command:

```bash
python3 gene_coexpression_app.py
```

The application will open in a desktop window.

## Project Team

- **A. Medha Sri** — Backend development and correlation analysis
- **Vaishno Kritika** — Input data preparation, testing, documentation and presentation

## Conclusion
