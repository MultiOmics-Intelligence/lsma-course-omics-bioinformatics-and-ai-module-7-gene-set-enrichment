<a href="https://lsma.co.uk/courses"><img src="https://raw.githubusercontent.com/MultiOmics-Intelligence/.github/main/assets/lsma-walk-courses.svg" alt="LSMA, London School of Multiomics and AI: learn by doing, on real data" width="100%"></a>

# Module 7 project: From a gene list to biology

LSMA · Omics, Bioinformatics and AI · Module 7 project

<p>
<img src="https://img.shields.io/badge/LSMA-free%20course-D9481F?style=flat-square" alt="LSMA free course">
<img src="https://img.shields.io/badge/course-Omics%2C%20Bioinformatics%20and%20AI-0F1B33?style=flat-square" alt="Course">
<img src="https://img.shields.io/badge/language-Bash%20%7C%20Python-4EAA25?style=flat-square&logo=gnubash&logoColor=white" alt="Language">
<img src="https://img.shields.io/badge/checks-automatic-2F7D5B?style=flat-square&logo=githubactions&logoColor=white" alt="Automatic checks">
</p>

## Research question

Thirty genes went up when cells were grown in low oxygen. Which biological processes do they point to, how sure can we be, and what story can we honestly tell?

## What you will build

| Files in (`data/`) | Files out (`results/`) |
|---|---|
| `genes_up.txt`, the 30 genes (made up for teaching) | `enrichment.tsv`, one row per pathway: overlaps, expected, p-value, Bonferroni and FDR |
| `pathways.tsv`, six simplified pathway lists | `network_summary.tsv`, nodes, links and pieces of the network, before and after adding the missing genes |
| `glucose_links.tsv`, the teaching table of links in glucose use | `story.md`, the biological story and its cautions |

Your code lives in `scripts/gene_list_to_biology.py`, task by task. It already reads the data into the names the lesson uses: `pathways`, `genes_up`, `links` and `N`.

## Before you start

1. Get your own copy. Use the GitHub Classroom invitation link in the project lesson on [lsma.co.uk](https://lsma.co.uk), or press **Use this template** on this page.
2. Clone your copy to your computer, or open it in GitHub Codespaces.
3. Install the packages once: `pip install pandas scipy networkx` (all already in Colab).
4. Work in `scripts/`. Run everything from the top folder of your copy: `python scripts/gene_list_to_biology.py`. You can explore in a notebook first; the script is what you hand in. (A notebook saved in `scripts/` is accepted too.)
5. Write every output to `results/`, with the exact file names above.

Keep the background fixed. The background is the set of genes that could have appeared in the list: here, the 10,000 genes measured in the experiment. Every p-value depends on it, so state it in your story.

You can run the checks yourself at any time: `python tests/check_results.py`.

## Tasks

### 1. Count the overlaps

Lesson to revisit: Module 7, lesson 4.

For each pathway, find which of the 30 genes are in it, and how many. Work out how many hits you would expect by chance: list size × pathway size ÷ background.

**Checkpoint:** glycolysis has 7 hits of 12 (ALDOA, ENO1, HK2, PFKL, PGK1, PKM, TPI1) against 0.036 expected. TCA cycle (SDHA), cell cycle (CDKN1A) and apoptosis (MCL1) have 1 hit each; DNA repair and cholesterol synthesis have none.

### 2. Run the hypergeometric test

Lesson to revisit: Module 7, lesson 4.

For each pathway, compute the chance of seeing at least that many hits if the 30 genes had been drawn at random from the background. Sort the pathways by p-value.

**Checkpoint:** glycolysis p = 8.1 × 10⁻¹⁶. Three more pathways also fall below 0.05 on a single hit each: TCA cycle 0.027, apoptosis 0.035 and cell cycle 0.047.

### 3. Correct for multiple testing

Lesson to revisit: Module 7, lesson 4.

Apply Bonferroni and Benjamini–Hochberg to the six p-values. List the pathways with FDR below 0.05. Save the full table as `results/enrichment.tsv`.

**Checkpoint:** only glycolysis survives (FDR 4.8 × 10⁻¹⁵). The three single-hit pathways all have FDR 0.070, and Bonferroni values of 0.16 to 0.28.

### 4. Build a network of the significant genes

Lesson to revisit: Module 7, lesson 5.

Use the teaching table of links in `data/glucose_links.tsv`. Keep only the genes in your list, count nodes, links and pieces; then put the missing genes back and count again. Save `results/network_summary.tsv`.

**Checkpoint:** the hits alone form 11 genes, 7 links and 4 pieces; PGK1 is left on its own. Adding GPI, GAPDH and PGAM1, which were measured but not in the list, joins all 14 genes into 1 piece running from SLC2A1 to LDHA. PFKL has the highest degree, 4.

### 5. Tell the biological story, with its cautions

Lesson to revisit: Module 7, lessons 4, 5 and 6.

In `results/story.md`, write a short paragraph answering the research question, with the numbers. Then list at least three cautions, one per line.

**Checkpoint:** your story names glycolysis with 7 of 12 genes (FDR about 5 × 10⁻¹⁵), explains why the network has gaps, says the three single-gene pathways did not survive correction, and proposes a test in another layer, such as measuring lactate.

## Hand in

Commit and push `scripts/` and `results/`. The **Checks** tab (GitHub Actions) then runs. A green tick means every checkpoint passed. You can push as often as you like.

## If a check fails

Open the failed run in the Actions tab and read the log. Each line starts with PASS or FAIL and names the task, for example `FAIL  Task 2 · Run the hypergeometric test: ...`. The message tells you what to look at, not the answer. Fix that task, run your script again and push.

## Getting help

- The project lesson on lsma.co.uk, and the earlier lessons named in each task.
- Helix, the course tutor.
- Email info@lsma.co.uk.

Never post personal data in public, in an issue, a commit or a file.

**Follow LSMA**

<p>
  <a href="https://lsma.co.uk"><img src="https://raw.githubusercontent.com/MultiOmics-Intelligence/.github/main/assets/btn-lsma-website.svg" alt="LSMA: lsma.co.uk" height="54"></a>
  <a href="https://www.linkedin.com/showcase/113252126/"><img src="https://raw.githubusercontent.com/MultiOmics-Intelligence/.github/main/assets/btn-lsma-linkedin.svg" alt="LSMA on LinkedIn" height="54"></a>
  <a href="https://www.instagram.com/lsma.uk/"><img src="https://raw.githubusercontent.com/MultiOmics-Intelligence/.github/main/assets/btn-lsma-instagram.svg" alt="LSMA on Instagram: @lsma.uk" height="54"></a>
  <a href="mailto:info@lsma.co.uk"><img src="https://raw.githubusercontent.com/MultiOmics-Intelligence/.github/main/assets/btn-lsma-email.svg" alt="Email: info@lsma.co.uk" height="54"></a>
</p>

Good luck, and enjoy the project.

**Dr. Anas Alzaeem**<br>
Academic Director, LSMA

[![Dr. Anas Alzaeem, Academic Director, LSMA](https://raw.githubusercontent.com/MultiOmics-Intelligence/.github/main/assets/signature/LSMA_Banner_Anas_Alzaeem.png)](https://lsma.co.uk)

© 2026 MultiOmics Intelligence Ltd. LSMA is the education arm of MultiOmics Intelligence. Code licensed under MIT.
