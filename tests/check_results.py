"""Automatic checks for Module 7 project: From a gene list to biology.

Run from the top folder of your copy:   python tests/check_results.py
It reads only the files you commit in results/ and scripts/, uses the Python
standard library only, and prints PASS or FAIL for each task.
"""
import csv
import hashlib
import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
SCRIPTS = ROOT / "scripts"
DATA = ROOT / "data"
_failed = []


class Fail(Exception):
    pass


def check(label):
    """Run one labelled check and print PASS or FAIL."""
    def wrap(fn):
        try:
            fn()
        except Fail as err:
            print(f"FAIL  {label}: {err}")
            _failed.append(label)
        except Exception as err:  # an unreadable file should guide, not crash
            print(f"FAIL  {label}: could not read your results ({type(err).__name__}: {err})")
            _failed.append(label)
        else:
            print(f"PASS  {label}")
        return fn
    return wrap


def finish():
    print()
    if _failed:
        print(f"{len(_failed)} check(s) failed. Read the messages above, fix that task and push again.")
        sys.exit(1)
    print("Every checkpoint passed. Well done.")


def need(path):
    path = Path(path)
    if not path.exists():
        raise Fail(f"{path.relative_to(ROOT)} not found; create it and commit it")
    return path


def read_rows(path):
    """Read a CSV (or tab-separated) file into a list of dicts with lower-case keys."""
    text = need(path).read_text(encoding="utf-8-sig")
    dialect = "excel-tab" if "\t" in text.splitlines()[0] else "excel"
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    return [{(k or "").strip().lower(): (v or "").strip() for k, v in r.items()} for r in rows]


def column(rows, *names, path=""):
    if not rows:
        raise Fail(f"{path} has no data rows")
    for n in names:
        if n in rows[0]:
            return n
    raise Fail(f"{path} needs a column called {names[0]}; found {', '.join(rows[0])}")


def num(value, what):
    try:
        return float(str(value).replace("%", "").strip())
    except ValueError:
        raise Fail(f"{what}: '{value}' is not a number")


def near(a, b, tol):
    return abs(a - b) <= tol


def _norm(text):
    return re.sub(r"\s+", " ", text.replace("\r", "")).strip()


def sha(path):
    return hashlib.sha256(_norm(Path(path).read_text(encoding="utf-8")).encode()).hexdigest()


def edited(path, starter_sha):
    path = need(path)
    if sha(path) == starter_sha:
        raise Fail(f"{path.relative_to(ROOT)} is still the untouched starter; write your code in it")
    return path


def png_info(path):
    """Width, height and dpi (None if not stored) of a PNG file."""
    data = Path(path).read_bytes()
    if data[:8] != b"\x89PNG\r\n\x1a\n":
        return None
    w, h = struct.unpack(">II", data[16:24])
    dpi = None
    i = 8
    while i + 8 <= len(data):
        length, kind = struct.unpack(">I4s", data[i:i + 8])
        if kind == b"pHYs":
            x, _, unit = struct.unpack(">IIB", data[i + 8:i + 17])
            if unit == 1:
                dpi = x * 0.0254
            break
        if kind == b"IDAT":
            break
        i += 12 + length
    return w, h, dpi


def image(path, min_px=300, min_bytes=5000):
    """Check that a PNG or SVG is a real image of sensible size; return PNG info."""
    path = need(path)
    rel = path.relative_to(ROOT)
    size = path.stat().st_size
    if size < min_bytes:
        raise Fail(f"{rel} is only {size} bytes; it looks empty")
    if path.suffix.lower() == ".svg":
        head = path.read_text(encoding="utf-8", errors="replace")[:5000]
        if "<svg" not in head:
            raise Fail(f"{rel} is not a valid SVG file")
        return None
    info = png_info(path)
    if info is None:
        raise Fail(f"{rel} is not a real PNG file (was it renamed from another format?)")
    w, h, _ = info
    if w < min_px or h < min_px // 3:
        raise Fail(f"{rel} is only {w} x {h} pixels; save it with a higher dpi")
    return info


def text_of(path, min_words=0):
    path = need(path)
    text = path.read_text(encoding="utf-8", errors="replace")
    words = len(re.findall(r"\w+", text))
    if words < min_words:
        raise Fail(f"{path.relative_to(ROOT)} has {words} words; write at least {min_words}")
    return text

import math

N_BG = 10_000


def load_data():
    pathways = {}
    for r in read_rows(DATA / "pathways.tsv"):
        pathways.setdefault(r["pathway"], []).append(r["gene"])
    genes = need(DATA / "genes_up.txt").read_text().split()
    links = [(r["gene_a"], r["gene_b"]) for r in read_rows(DATA / "glucose_links.tsv")]
    return pathways, genes, links


PATHWAYS, GENES, LINKS = load_data()


def hyper_sf(k, N, K, n):
    """P(X >= k) for the hypergeometric distribution, in plain Python."""
    if k <= 0:
        return 1.0
    total = math.comb(N, n)
    return sum(math.comb(K, i) * math.comb(N - K, n - i) for i in range(k, min(K, n) + 1)) / total


def expected():
    rows = []
    for name, members in PATHWAYS.items():
        hits = sorted(set(members) & set(GENES))
        p = hyper_sf(len(hits), N_BG, len(members), len(GENES)) if hits else 1.0
        rows.append({"pathway": name, "size": len(members), "hits": hits,
                     "expected": len(GENES) * len(members) / N_BG, "p": p})
    rows.sort(key=lambda r: r["p"])
    m = len(rows)
    adj, prev = [0.0] * m, 1.0
    for i in range(m - 1, -1, -1):          # Benjamini-Hochberg, step-up
        prev = min(prev, rows[i]["p"] * m / (i + 1))
        adj[i] = prev
    for r, a in zip(rows, adj):
        r["bonf"] = min(1.0, r["p"] * m)
        r["fdr"] = a
    return {r["pathway"].lower(): r for r in rows}


WANT = expected()


def enrichment():
    path = "results/enrichment.tsv"
    rows = read_rows(RESULTS / "enrichment.tsv")
    for c in ("pathway", "size", "hits", "hit_genes", "expected", "p_value", "bonferroni", "fdr_bh"):
        column(rows, c, path=path)
    if len(rows) != len(WANT):
        raise Fail(f"expected {len(WANT)} rows in {path} (one per pathway), found {len(rows)}")
    got = {r["pathway"].lower(): r for r in rows}
    missing = [p for p in WANT if p not in got]
    if missing:
        raise Fail(f"{path} is missing the pathway(s): {', '.join(missing)}; keep the names "
                   "from data/pathways.tsv")
    return rows, got


def rel_near(a, b, rel=0.02, abs_tol=1e-12):
    return abs(a - b) <= max(abs_tol, rel * abs(b))


@check("Setup · Your script")
def _():
    if list(SCRIPTS.glob("*.ipynb")):
        return
    edited(SCRIPTS / "gene_list_to_biology.py", "271908ca1de5477ca17dd7f903678e2b196816590ea8b83bb9940df09364a77e")


@check("Task 1 · Count the overlaps")
def _():
    rows, got = enrichment()
    for name, w in WANT.items():
        r = got[name]
        if int(num(r["size"], "size")) != w["size"]:
            raise Fail(f"the size of {r['pathway']} is not right; count its members")
        if int(num(r["hits"], "hits")) != len(w["hits"]):
            raise Fail(f"the number of hits for {r['pathway']} is not right; use "
                       "set(members) & set(genes_up)")
        listed = sorted(g.strip() for g in re.split(r"[,;\s]+", r["hit_genes"]) if g.strip())
        if listed != w["hits"]:
            raise Fail(f"hit_genes for {r['pathway']} does not list the matching genes")
        if not near(num(r["expected"], "expected"), w["expected"], 0.0006):
            raise Fail(f"the expected number for {r['pathway']} is off; it is list size x "
                       "pathway size / background")


@check("Task 2 · Run the hypergeometric test")
def _():
    rows, got = enrichment()
    for name, w in WANT.items():
        p = num(got[name]["p_value"], "p_value")
        if not rel_near(p, w["p"], 0.01):
            raise Fail(f"the p-value for {got[name]['pathway']} is off; use "
                       "hypergeom.sf(k - 1, N, K, n), and p = 1 when there are no hits")
    ps = [num(r["p_value"], "p_value") for r in rows]
    if ps != sorted(ps):
        raise Fail("results/enrichment.tsv is not sorted by p-value, smallest first")


@check("Task 3 · Correct for multiple testing")
def _():
    rows, got = enrichment()
    for name, w in WANT.items():
        r = got[name]
        if not rel_near(num(r["bonferroni"], "bonferroni"), w["bonf"], 0.01):
            raise Fail(f"the Bonferroni value for {r['pathway']} is off; multiply by the "
                       "number of tests and cap at 1")
        if not rel_near(num(r["fdr_bh"], "fdr_bh"), w["fdr"], 0.01):
            raise Fail(f"the Benjamini-Hochberg value for {r['pathway']} is off; use "
                       "false_discovery_control(p, method='bh') on all six p-values")


def components(nodes, edges):
    adj = {v: set() for v in nodes}
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    seen, pieces = set(), 0
    for v in adj:
        if v not in seen:
            pieces += 1
            stack = [v]
            while stack:
                x = stack.pop()
                if x not in seen:
                    seen.add(x)
                    stack.extend(adj[x] - seen)
    return adj, pieces


@check("Task 4 · Build a network of the significant genes")
def _():
    path = "results/network_summary.tsv"
    rows = read_rows(RESULTS / "network_summary.tsv")
    for c in ("network", "nodes", "links", "pieces", "isolated_genes",
              "highest_degree_gene", "highest_degree"):
        column(rows, c, path=path)
    got = {r["network"].lower(): r for r in rows}
    all_nodes = {g for e in LINKS for g in e}
    hit_nodes = {g for g in all_nodes if g in GENES}
    hit_edges = [e for e in LINKS if e[0] in hit_nodes and e[1] in hit_nodes]
    adj_h, pieces_h = components(hit_nodes, hit_edges)
    adj_f, pieces_f = components(all_nodes, LINKS)
    want = {"hits_only": (len(hit_nodes), len(hit_edges), pieces_h),
            "with_missing_genes": (len(all_nodes), len(LINKS), pieces_f)}
    for name, (nodes, links, pieces) in want.items():
        if name not in got:
            raise Fail(f"{path} needs a row called {name}")
        r = got[name]
        for col, value in (("nodes", nodes), ("links", links), ("pieces", pieces)):
            if int(num(r[col], col)) != value:
                raise Fail(f"{col} for {name} is not right; "
                           + ("keep only the genes in your list with G.subgraph"
                              if name == "hits_only" else "use the full network of all 14 links"))
    alone = sorted(g for g, nb in adj_h.items() if not nb)
    listed = sorted(g.strip() for g in re.split(r"[,;\s]+", got["hits_only"]["isolated_genes"]) if g.strip())
    if listed != alone:
        raise Fail("isolated_genes for hits_only should name the gene(s) left on their own")
    deg = {g: len(nb) for g, nb in adj_f.items()}
    top = max(deg.values())
    r = got["with_missing_genes"]
    if r["highest_degree_gene"].upper() not in {g for g, d in deg.items() if d == top} \
            or int(num(r["highest_degree"], "highest_degree")) != top:
        raise Fail("highest_degree_gene and highest_degree are not right; use G.degree() on "
                   "the full network")


@check("Task 5 · Tell the biological story, with its cautions")
def _():
    text = text_of(RESULTS / "story.md", min_words=150)
    low = text.lower()
    if "glycolysis" not in low:
        raise Fail("results/story.md should name the pathway the list points to")
    if not re.search(r"\b7\b|\bseven\b", low) or not re.search(r"\b12\b|\btwelve\b", low):
        raise Fail("results/story.md should give how many genes of the pathway were hits, "
                   "out of how many (for example 3 of 10)")
    if "fdr" not in low and "false discovery" not in low and "benjamini" not in low:
        raise Fail("results/story.md should quote the corrected value (FDR)")
    if not re.search(r"correct|surviv|did not pass|not significant", low):
        raise Fail("results/story.md should say what happened to the single-gene pathways "
                   "after correction")
    bullets = [l for l in text.splitlines() if re.match(r"\s*(\d+[.)]|[-*])\s+\S", l)]
    if len(bullets) < 3:
        raise Fail(f"results/story.md lists {len(bullets)} caution(s); give at least three, "
                   "one per line, starting with - or a number")
    if not re.search(r"lactate|protein|metabol|hif|another layer|second layer|flux|uptake", low):
        raise Fail("results/story.md should propose a test in another layer")


finish()
