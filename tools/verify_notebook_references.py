"""Check every reference in a notebook's "## References" cell against its PubMed record (found by DOI).

Usage:  python tools/verify_notebook_references.py notebooks/screen.ipynb      (needs nbformat and internet)

Compares each author's surname and initials, year, journal, volume, first page and title. Entries without a DOI (software) or not indexed
in PubMed (for example 1960s physical-chemistry papers) are listed as NOT IN PUBMED and must be checked another way (Crossref). Exit status 1 if any entry fails.
"""
import json, re, sys, time, urllib.error, urllib.parse, urllib.request
import nbformat

nb = nbformat.read(sys.argv[1], as_version=4)
ref_cell = next(c.source for c in nb.cells if c.cell_type == "markdown" and c.source.lstrip().startswith("## References"))
refs = re.findall(r"(?m)^(\d+)\. (.+)$", ref_cell)

def get(url):
    for wait in (0.7, 3, 8, 20):
        time.sleep(wait)
        try:
            return json.loads(urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "sida-docking-course-reference-check"}), timeout=30).read())
        except urllib.error.HTTPError as e:
            if e.code != 429:
                raise
    raise RuntimeError("rate limited")

def clean(t):
    t = re.sub(r"<[^>]+>|\*", "", t)
    return re.sub(r"[^a-z0-9]", "", t.lower())

def initials(s):
    return re.sub(r"[^A-Za-z]", "", s)

summary = {"pubmed_pass": [], "fail": [], "not_in_pubmed": []}
for num, text in refs:
    m_doi = re.search(r"doi:(\S+)$", text)
    if not m_doi:
        summary["not_in_pubmed"].append((num, "no DOI (software / web address)"))
        continue
    doi = m_doi.group(1)
    ids = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&term=" + urllib.parse.quote(doi + "[AID]"))["esearchresult"]["idlist"]
    if not ids:
        summary["not_in_pubmed"].append((num, "DOI %s not found in PubMed" % doi))
        continue
    d = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=%s" % ids[0])["result"][ids[0]]
    problems = []
    # --- split my citation: authors (year) title. journal vol, pages. doi
    m = re.match(r"(.+?) \((\d{4})\) (.+?)\. \*?([^*]+?)\*? (\S+?)(?:, ([^ ]+?))?\. doi:", text)
    if not m:
        problems.append("could not parse the citation text")
    else:
        authors_txt, year, title, journal, vol, pages = m.groups()
        # authors: "Surname, I. J., Surname2, K., and Surname3, L." possibly ending "*et al.*"
        etal = "et al." in authors_txt
        a = re.sub(r"\*et al\.\*,?", "", authors_txt).replace(", and ", ", ").replace(" and ", ", ")
        toks = [t.strip() for t in re.split(r",\s*", a) if t.strip()]
        mine = []
        i = 0
        while i < len(toks):
            sur = toks[i]; ini = toks[i + 1] if i + 1 < len(toks) else ""
            suffix = ""
            if i + 2 < len(toks) and toks[i + 2] in ("Jr.", "III"):
                suffix = toks[i + 2]; i += 1
            mine.append((sur, ini, suffix)); i += 2
        pm = [x["name"] for x in d.get("authors", [])]
        n = len(mine)
        if not etal and n != len(pm):
            problems.append("author count: mine %d, PubMed %d" % (n, len(pm)))
        for k, (sur, ini, suf) in enumerate(mine[:len(pm)]):
            want = pm[k]
            wsur_ini = re.sub(r"\s+(Jr|3rd|2nd|Sr)$", "", want)
            wsur, _, wini = wsur_ini.rpartition(" ")
            if wsur.lower() != sur.lower() or initials(wini).lower() != initials(ini).lower():
                problems.append("author %d: mine '%s, %s' vs PubMed '%s'" % (k + 1, sur, ini, want))
        if year not in d["pubdate"]:
            problems.append("year: mine %s vs PubMed %s" % (year, d["pubdate"]))
        j1, j2 = clean(journal), clean(d["source"].split(".")[0])
        if not (j1 == j2 or j1.startswith(j2) or j2.startswith(j1)):
            problems.append("journal: mine '%s' vs PubMed '%s'" % (journal, d["source"].split(".")[0]))
        if clean(vol) != clean(d.get("volume", "")) and not (not d.get("volume") and True):
            problems.append("volume: mine %s vs PubMed %s" % (vol, d.get("volume")))
        first = re.split(r"[-–]", pages or "")[0]
        if pages and d.get("pages") and clean(first) != clean(re.split(r"[-–]", d["pages"])[0]):
            problems.append("first page: mine %s vs PubMed %s" % (first, d["pages"]))
        if clean(d["title"]).rstrip() not in clean(title) and clean(title) not in clean(d["title"]):
            problems.append("title: mine '%s' vs PubMed '%s'" % (title, d["title"]))
    (summary["fail"] if problems else summary["pubmed_pass"]).append((num, "PMID %s %s" % (ids[0], "; ".join(problems))))
    print("%-3s %s PMID %-9s %s" % (num, "FAIL" if problems else "PASS", ids[0], "; ".join(problems)))
print()
for num, why in summary["not_in_pubmed"]:
    print("%-3s NOT IN PUBMED: %s" % (num, why))
print("\nPubMed PASS: %d | FAIL: %d | not in PubMed: %d" % (len(summary["pubmed_pass"]), len(summary["fail"]), len(summary["not_in_pubmed"])))
sys.exit(1 if summary["fail"] else 0)
