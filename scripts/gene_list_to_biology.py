"""Module 7 project: from a gene list to biology.

Run from the top folder of your copy:   python scripts/gene_list_to_biology.py
You can also work in a notebook first (the lesson does), but keep this script
as the final, rerunnable version. Every result goes into results/.

Install once:   pip install pandas scipy networkx
"""
from pathlib import Path

import networkx as nx
import pandas as pd
from scipy.stats import false_discovery_control, hypergeom

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
RESULTS.mkdir(exist_ok=True)

# The data cell of the lesson, read from data/
table = pd.read_csv(DATA / "pathways.tsv", sep="\t")
pathways = {name: list(g["gene"]) for name, g in table.groupby("pathway", sort=False)}
genes_up = (DATA / "genes_up.txt").read_text().split()
links = [tuple(pair) for pair in pd.read_csv(DATA / "glucose_links.tsv", sep="\t").values]
N = 10_000          # background: all genes measured in the experiment
n = len(genes_up)
print(n, "genes up;", len(pathways), "pathways")

# Research question (write it in your own words):
# TODO

# ---------------------------------------------------------------------------
# Task 1 · Count the overlaps (lesson 4)
# ---------------------------------------------------------------------------
# TODO: for each pathway, find which of the 30 genes are in it, and how many
#       (set(members) & set(genes_up)).
# TODO: work out the expected number of hits: list size x pathway size / background.
# TODO: collect one row per pathway: pathway, size, hits, hit_genes, expected
#       (hit_genes: the gene symbols, comma-separated, empty if none).

# ---------------------------------------------------------------------------
# Task 2 · Run the hypergeometric test (lesson 4)
# ---------------------------------------------------------------------------
# TODO: p_value = hypergeom.sf(k - 1, N, K, n) for each pathway; p = 1 if no hits.
# TODO: sort the pathways by p-value.

# ---------------------------------------------------------------------------
# Task 3 · Correct for multiple testing (lesson 4)
# ---------------------------------------------------------------------------
# TODO: bonferroni = p x 6, capped at 1; fdr_bh = false_discovery_control(p, method="bh").
# TODO: save the table as results/enrichment.tsv (tab-separated, index=False) with the
#       columns: pathway, size, hits, hit_genes, expected, p_value, bonferroni, fdr_bh

# ---------------------------------------------------------------------------
# Task 4 · Build a network of the significant genes (lesson 5)
# ---------------------------------------------------------------------------
# TODO: build the full network with nx.Graph(links).
# TODO: keep only the genes in your list (G.subgraph(...)); count nodes, links and
#       pieces (nx.number_connected_components). Which gene is left on its own?
# TODO: put the missing genes back (the full network) and count again.
# TODO: find the gene with the highest degree.
# TODO: save results/network_summary.tsv with two rows, hits_only and
#       with_missing_genes, and the columns:
#       network, nodes, links, pieces, isolated_genes, highest_degree_gene, highest_degree

# ---------------------------------------------------------------------------
# Task 5 · Tell the biological story, with its cautions (lessons 4, 5 and 6)
# ---------------------------------------------------------------------------
# Write results/story.md by hand: a short paragraph that answers the research
# question with your numbers, then a list of at least three cautions, one per line.
