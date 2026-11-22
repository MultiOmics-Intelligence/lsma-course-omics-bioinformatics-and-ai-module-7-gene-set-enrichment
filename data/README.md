# Data

All three files are copied exactly from the data cell and the task 4 table of the lesson. They are **made up for teaching**: no experiment or patient is behind them.

| File | What it is |
|---|---|
| `genes_up.txt` | The 30 genes that went up in a teaching experiment, cells grown in low oxygen, one symbol per line. A made-up list, mostly modelled on genes known to rise when oxygen is short. |
| `pathways.tsv` | Six simplified pathway lists for teaching (glycolysis, TCA cycle, DNA repair, cell cycle, apoptosis, cholesterol synthesis), one row per gene: `pathway`, `gene`. Each list holds only well-established members; real databases hold more. |
| `glucose_links.tsv` | A teaching table of 14 links between consecutive steps, or direct regulation, in glucose use: `gene_a`, `gene_b`. |

The background, the 10,000 genes measured in the experiment, is a number, not a file: `N = 10_000` in the starter script.

Gene symbols are HGNC-approved human gene symbols ([genenames.org](https://www.genenames.org)). The pathway lists are simplified from the kind of content held in Reactome and KEGG. The files are released with the code under the MIT licence.
