"""Construct check: an explicit residue mapping between two structures.

Author numbering is NEVER changed. The mapping is a table you can read: each row pairs a residue of the
reference with a residue of the query (or marks one side as a gap), and says whether the amino acids are
identical and whether the author numbers agree.
"""
from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import biotite.sequence as bseq
import biotite.sequence.align as balign

AA20 = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class Mapping:
    rows: tuple
    n_ref: int
    n_qry: int
    n_aligned: int
    n_mismatch: int          # aligned pairs whose amino acids differ (point mutations / modified residues)
    n_ref_only: int          # reference residues with no partner (resolved there, not in the query)
    n_qry_only: int
    n_different_number: int  # aligned pairs whose author numbers differ
    offsets: dict            # {qry_num - ref_num: pair count} over aligned pairs whose numbers differ
    n_same_number_conflicts: int  # same author number, different amino acid, and NOT aligned together
    numbering_shared: bool


def _key(r):
    return (r[0], r[1].strip())


def map_residues(ref, qry) -> Mapping:
    """ref/qry: [(author_num:int, icode:str, one_letter_aa:str)] of resolved polymer residues."""
    if not ref or not qry:
        raise ValueError("cannot map an empty residue list (reference %d, query %d)" % (len(ref), len(qry)))
    clean = lambda seq: "".join(a if a in AA20 else "X" for _, _, a in seq)
    s1, s2 = bseq.ProteinSequence(clean(ref).replace("X", "A")), bseq.ProteinSequence(clean(qry).replace("X", "A"))
    aln = balign.align_optimal(s1, s2, balign.SubstitutionMatrix.std_protein_matrix(),
                               gap_penalty=(-10, -1), terminal_penalty=False, local=False, max_number=1)[0]
    rows, aligned = [], []
    for i, j in aln.trace:
        a = ref[i] if i >= 0 else None
        b = qry[j] if j >= 0 else None
        rows.append({
            "ref_num": "%d%s" % (a[0], a[1].strip()) if a else "", "ref_aa": a[2] if a else "",
            "qry_num": "%d%s" % (b[0], b[1].strip()) if b else "", "qry_aa": b[2] if b else "",
            "identical": bool(a and b and a[2] == b[2]),
            "same_number": bool(a and b and _key(a) == _key(b)),
        })
        if a and b:
            aligned.append((a, b))
    n_mis = sum(1 for a, b in aligned if a[2] != b[2])
    diff = [(a, b) for a, b in aligned if _key(a) != _key(b)]
    offsets = dict(Counter(b[0] - a[0] for a, b in diff))
    paired = {(_key(a), _key(b)) for a, b in aligned}
    qry_by_key = {_key(b): b for b in qry}
    conflicts = sum(1 for a in ref if _key(a) in qry_by_key and qry_by_key[_key(a)][2] != a[2]
                    and (_key(a), _key(a)) not in paired)
    return Mapping(
        rows=tuple(rows), n_ref=len(ref), n_qry=len(qry), n_aligned=len(aligned), n_mismatch=n_mis,
        n_ref_only=sum(1 for r in rows if r["ref_aa"] and not r["qry_aa"]),
        n_qry_only=sum(1 for r in rows if r["qry_aa"] and not r["ref_aa"]),
        n_different_number=len(diff), offsets=offsets, n_same_number_conflicts=conflicts,
        numbering_shared=(not diff and conflicts == 0),
    )


def load_polymer_residues(cif_path, chain: str):
    """Resolved polymer residues of one chain in AUTHOR numbering, as map_residues() expects."""
    import gemmi
    st = gemmi.read_structure(str(cif_path))
    out = []
    for r in st[0][chain].get_polymer():
        info = gemmi.find_tabulated_residue(r.name)
        aa = info.one_letter_code.upper() if info and info.one_letter_code != " " else "X"
        out.append((r.seqid.num, r.seqid.icode.strip(), aa))
    if not out:
        raise ValueError("%s chain %s has no polymer residues" % (cif_path, chain))
    return out


def write_mapping_csv(m: Mapping, path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = ["ref_num", "ref_aa", "qry_num", "qry_aa", "identical", "same_number"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(m.rows)
    return path
