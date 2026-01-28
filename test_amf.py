#!/usr/bin/env python3
"""
AMF-style (arbuscular mycorrhizal fungi) *toy* genome simulator
==============================================================

python test_amf.py --total-len 100000 --n-contigs 10 --n-genomes 1000 --outdir assemblies --manifest
Goal
----
Generate small test genomes that *look AMF-like* in structure *without increasing total size*:
- More repeats / TE-like elements (scaled-down proportions)
- Multiple contigs (chromosome fragments) with repeat-enriched "subtelomeric-like" ends
- Occasional low-complexity tracts
- SVs (INS/DEL/DUP/INV/TRA) biased toward repeats/hotspots


This is NOT a biological AMF genome generator; it's a *stress-test generator* for SV/graph tools.

Outputs
-------
- ref.fa                              (multi-contig reference)
- <outdir>/asm_0000.fa ...            (assemblies)
- <outdir>/truth.tsv                  (events applied per assembly)
- <outdir>/manifest.tsv (optional)

Coordinates in truth.tsv
------------------------
- INS/DEL: "pos" is the local position in the evolving contig sequence at time of application.
- DUP/INV/TRA: "start","end" are local positions in the evolving contig at time of application.
- TRA: "target" is the insertion position (local) in the *target* contig at time of application.
- "extra" stores contig ids/names to keep events interpretable.

Default sizes
-------------
Total genome size stays ~100,000 bp by default (split across contigs).
"""

from __future__ import annotations
import argparse
import os
import random
import re
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

DNA = ["A", "C", "G", "T"]

# ----------------------------
# Sequence helpers
# ----------------------------
def random_dna(rng: random.Random, length: int, gc: float = 0.5) -> str:
    """Generate DNA with approximate GC fraction."""
    if length <= 0:
        return ""
    # Probabilities
    pG = gc / 2.0
    pC = gc / 2.0
    pA = (1.0 - gc) / 2.0
    pT = (1.0 - gc) / 2.0
    # CDF sampling
    cdf = [(pA, "A"), (pA + pC, "C"), (pA + pC + pG, "G"), (1.0, "T")]
    out = []
    for _ in range(length):
        x = rng.random()
        for th, b in cdf:
            if x <= th:
                out.append(b)
                break
    return "".join(out)

def revcomp(seq: str) -> str:
    comp = {"A":"T","C":"G","G":"C","T":"A","N":"N",
            "a":"t","c":"g","g":"c","t":"a","n":"n"}
    return "".join(comp.get(b, "N") for b in reversed(seq))

def write_fasta_multi(path: str, contigs: Sequence[Tuple[str,str]], width: int = 60) -> None:
    with open(path, "w") as f:
        for name, seq in contigs:
            f.write(f">{name}\n")
            for i in range(0, len(seq), width):
                f.write(seq[i:i+width] + "\n")

# ----------------------------
# Repeat library + AMF-ish reference
# ----------------------------
def make_repeat_library(rng: random.Random, n: int, min_len: int, max_len: int, gc: float) -> List[str]:
    lib = []
    for _ in range(n):
        L = rng.randint(min_len, max_len)
        # add some internal structure by repeating a short motif
        motif = random_dna(rng, rng.randint(8, 20), gc=gc)
        rep = (motif * (L // len(motif) + 1))[:L]
        # mutate a little
        rep_list = list(rep)
        for i in range(len(rep_list)):
            if rng.random() < 0.02:
                rep_list[i] = rng.choice(DNA)
        lib.append("".join(rep_list))
    return lib

def sprinkle_repeats(rng: random.Random, seq: str, repeat_lib: Sequence[str], density: float) -> str:
    """
    Insert repeat copies with approximate density (fraction of bases contributed by repeats).
    Keeps total length stable by substituting segments rather than inserting.
    """
    if not repeat_lib or density <= 0:
        return seq
    s = list(seq)
    L = len(s)
    # target replaced bases
    target = int(L * density)
    replaced = 0
    while replaced < target:
        rep = rng.choice(repeat_lib)
        rL = len(rep)
        if rL >= L:
            break
        start = rng.randint(0, L - rL)
        s[start:start+rL] = list(rep)
        replaced += rL
    return "".join(s)

def add_low_complexity(rng: random.Random, seq: str, n_tracts: int, tract_len: Tuple[int,int]) -> str:
    s = list(seq)
    L = len(s)
    for _ in range(n_tracts):
        tl = rng.randint(tract_len[0], tract_len[1])
        if tl >= L:
            continue
        start = rng.randint(0, L - tl)
        motif = rng.choice(["A","T","AT","TA","AAA","TTT","AAT","TTA"])
        tract = (motif * (tl // len(motif) + 1))[:tl]
        s[start:start+tl] = list(tract)
    return "".join(s)

def make_amf_like_reference(
    rng: random.Random,
    total_len: int,
    n_contigs: int,
    gc: float,
    repeat_density: float,
    subtel_repeat_density: float,
    subtel_frac: float,
    repeat_lib: Sequence[str],
) -> List[Tuple[str,str]]:
    """
    Produce multiple contigs totaling total_len.
    - Internal region: moderate repeats
    - Ends ("subtelomeric-like"): higher repeat density
    """
    n_contigs = max(1, n_contigs)
    # split lengths (slight variability)
    base = total_len // n_contigs
    lens = [base] * n_contigs
    rem = total_len - base * n_contigs
    for i in range(rem):
        lens[i % n_contigs] += 1
    # jitter while keeping sum fixed
    for _ in range(n_contigs * 2):
        i = rng.randrange(n_contigs)
        j = rng.randrange(n_contigs)
        if i == j:
            continue
        delta = rng.randint(-200, 200)
        if lens[i] + delta < 2000 or lens[j] - delta < 2000:
            continue
        lens[i] += delta
        lens[j] -= delta

    contigs = []
    for ci, clen in enumerate(lens):
        core = random_dna(rng, clen, gc=gc)
        # internal repeats
        core = sprinkle_repeats(rng, core, repeat_lib, density=repeat_density)
        # low complexity
        core = add_low_complexity(rng, core, n_tracts=max(1, clen//25000), tract_len=(60, 250))

        # subtelomeric-like ends
        end_len = int(clen * subtel_frac)
        if end_len >= 200:
            left = core[:end_len]
            mid  = core[end_len:clen-end_len]
            right= core[clen-end_len:]
            left = sprinkle_repeats(rng, left, repeat_lib, density=subtel_repeat_density)
            right= sprinkle_repeats(rng, right, repeat_lib, density=subtel_repeat_density)
            core = left + mid + right

        contigs.append((f"ref_ctg{ci:02d}", core))
    return contigs

# ----------------------------
# SV events & truth
# ----------------------------
@dataclass
class Event:
    asm: str
    eid: int
    kind: str  # INS, DEL, DUP, INV, TRA
    contig: str
    pos: int = 0
    start: int = 0
    end: int = 0
    target_contig: str = ""
    target: int = 0
    length: int = 0
    extra: str = ""

def write_truth_tsv(path: str, truth: Sequence[Event]) -> None:
    header = ["asm","event_id","type","contig","pos","start","end","target_contig","target","length","extra"]
    with open(path, "w") as f:
        f.write("\\t".join(header) + "\\n")
        for ev in truth:
            f.write("\\t".join([
                ev.asm, str(ev.eid), ev.kind, ev.contig,
                str(ev.pos), str(ev.start), str(ev.end),
                ev.target_contig, str(ev.target), str(ev.length),
                ev.extra
            ]) + "\\n")

# ----------------------------
# Hotspot sampling (biased to repeat-enriched ends)
# ----------------------------
def pick_hotspot_pos(rng: random.Random, clen: int, subtel_frac: float) -> int:
    if clen <= 1:
        return 0
    end_len = int(clen * subtel_frac)
    if end_len < 200:
        return rng.randrange(0, clen)
    # 65% chance: choose within ends (repeat-rich)
    if rng.random() < 0.65:
        if rng.random() < 0.5:
            return rng.randrange(0, end_len)
        return rng.randrange(clen - end_len, clen)
    # else choose anywhere
    return rng.randrange(0, clen)

# ----------------------------
# Apply SVs to multi-contig genome
# ----------------------------
def apply_ins(contigs: Dict[str,str], c: str, pos: int, ins: str) -> None:
    s = contigs[c]
    contigs[c] = s[:pos] + ins + s[pos:]

def apply_del(contigs: Dict[str,str], c: str, pos: int, L: int) -> int:
    s = contigs[c]
    end = min(len(s), pos + L)
    contigs[c] = s[:pos] + s[end:]
    return end - pos

def apply_dup(contigs: Dict[str,str], c: str, s0: int, s1: int, target: int) -> None:
    s = contigs[c]
    frag = s[s0:s1]
    contigs[c] = s[:target] + frag + s[target:]

def apply_inv(contigs: Dict[str,str], c: str, s0: int, s1: int) -> None:
    s = contigs[c]
    contigs[c] = s[:s0] + revcomp(s[s0:s1]) + s[s1:]

def apply_tra(contigs: Dict[str,str], c_src: str, s0: int, s1: int, c_tgt: str, target: int) -> None:
    src = contigs[c_src]
    frag = src[s0:s1]
    contigs[c_src] = src[:s0] + src[s1:]
    tgt = contigs[c_tgt]
    target = min(target, len(tgt))
    contigs[c_tgt] = tgt[:target] + frag + tgt[target:]

# ----------------------------
# Simulation driver
# ----------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="AMF-style toy SV simulator (small genomes).")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", type=str, default="assemblies")

    # Size and contigs (keep small for testing)
    ap.add_argument("--total-len", type=int, default=100_000, help="total genome length across contigs (default: 100k)")
    ap.add_argument("--n-contigs", type=int, default=10, help="number of contigs (default: 10)")
    ap.add_argument("--gc", type=float, default=0.35, help="GC fraction (default: 0.35; AMF-ish AT-rich)")

    # AMF-ish repeat profile (scaled-down)
    ap.add_argument("--repeat-density", type=float, default=0.25, help="repeat substitution density in cores (default 0.25)")
    ap.add_argument("--subtel-repeat-density", type=float, default=0.55, help="repeat density at ends (default 0.55)")
    ap.add_argument("--subtel-frac", type=float, default=0.12, help="fraction of each contig treated as end (default 0.12)")

    # Repeat library (toy TE-like)
    ap.add_argument("--repeat-lib-n", type=int, default=40)
    ap.add_argument("--repeat-len-min", type=int, default=80)
    ap.add_argument("--repeat-len-max", type=int, default=600)

    # Assemblies
    ap.add_argument("--n-genomes", type=int, default=1000)
    ap.add_argument("--prefix", type=str, default="asm_")
    ap.add_argument("--digits", type=int, default=4)
    ap.add_argument("--manifest", action="store_true")

    # SV counts per assembly (per-genome totals, distributed across contigs)
    ap.add_argument("--ins", type=int, default=2)
    ap.add_argument("--del", dest="dels", type=int, default=2)
    ap.add_argument("--dup", type=int, default=1)
    ap.add_argument("--inv", type=int, default=1)
    ap.add_argument("--tra", type=int, default=1)

    # SV size ranges (kept small)
    ap.add_argument("--ins-len-min", type=int, default=100)
    ap.add_argument("--ins-len-max", type=int, default=1200)
    ap.add_argument("--del-len-min", type=int, default=100)
    ap.add_argument("--del-len-max", type=int, default=1500)
    ap.add_argument("--seg-len-min", type=int, default=200)
    ap.add_argument("--seg-len-max", type=int, default=5000)

    args = ap.parse_args()
    rng = random.Random(args.seed)

    os.makedirs(args.outdir, exist_ok=True)

    # Build toy TE/repeat library and reference
    repeat_lib = make_repeat_library(
        rng,
        n=args.repeat_lib_n,
        min_len=args.repeat_len_min,
        max_len=args.repeat_len_max,
        gc=args.gc
    )

    ref_contigs = make_amf_like_reference(
        rng,
        total_len=args.total_len,
        n_contigs=args.n_contigs,
        gc=args.gc,
        repeat_density=args.repeat_density,
        subtel_repeat_density=args.subtel_repeat_density,
        subtel_frac=args.subtel_frac,
        repeat_lib=repeat_lib
    )

    write_fasta_multi("ref.fa", ref_contigs)

    # Base reference dict (for per-assembly copy)
    ref_dict: Dict[str,str] = {name: seq for name, seq in ref_contigs}
    contig_names = [name for name, _ in ref_contigs]

    truth: List[Event] = []
    eid = 1

    for i in range(args.n_genomes):
        asm_name = f"{args.prefix}{i:0{args.digits}d}"
        # working copy
        contigs = dict(ref_dict)

        # INS/DEL (biased to ends)
        for _ in range(args.ins):
            c = rng.choice(contig_names)
            pos = pick_hotspot_pos(rng, len(contigs[c]), args.subtel_frac)
            L = rng.randint(args.ins_len_min, args.ins_len_max)
            ins = random_dna(rng, L, gc=args.gc)
            apply_ins(contigs, c, pos, ins)
            truth.append(Event(asm=asm_name, eid=eid, kind="INS", contig=c, pos=pos, length=L,
                               extra="hotspot=subtel" if pos < int(len(contigs[c])*args.subtel_frac) or pos > int(len(contigs[c])*(1-args.subtel_frac)) else "hotspot=any"))
            eid += 1

        for _ in range(args.dels):
            c = rng.choice(contig_names)
            pos = pick_hotspot_pos(rng, len(contigs[c]), args.subtel_frac)
            L = rng.randint(args.del_len_min, args.del_len_max)
            dl = apply_del(contigs, c, pos, L)
            truth.append(Event(asm=asm_name, eid=eid, kind="DEL", contig=c, pos=pos, length=dl,
                               extra="hotspot=subtel" if pos < int(len(contigs[c])*args.subtel_frac) or pos > int(len(contigs[c])*(1-args.subtel_frac)) else "hotspot=any"))
            eid += 1

        # DUP/INV/TRA (operate on evolving contigs)
        for _ in range(args.dup):
            c = rng.choice(contig_names)
            clen = len(contigs[c])
            if clen < args.seg_len_min + 10:
                continue
            L = rng.randint(args.seg_len_min, min(args.seg_len_max, clen - 1))
            s0 = pick_hotspot_pos(rng, clen - L, args.subtel_frac)
            s0 = min(s0, clen - L)
            s1 = s0 + L
            # target biased to nearby (tandem-like) with 60% prob
            if rng.random() < 0.60:
                target = min(clen, s1 + rng.randint(0, 2000))
            else:
                target = rng.randint(0, clen)
            apply_dup(contigs, c, s0, s1, target)
            truth.append(Event(asm=asm_name, eid=eid, kind="DUP", contig=c, start=s0, end=s1, target=target, length=L,
                               extra="mode=tandem" if target >= s1 and target <= s1+2000 else "mode=dispersed"))
            eid += 1

        for _ in range(args.inv):
            c = rng.choice(contig_names)
            clen = len(contigs[c])
            if clen < args.seg_len_min + 10:
                continue
            L = rng.randint(args.seg_len_min, min(args.seg_len_max, clen - 1))
            s0 = pick_hotspot_pos(rng, clen - L, args.subtel_frac)
            s0 = min(s0, clen - L)
            s1 = s0 + L
            apply_inv(contigs, c, s0, s1)
            truth.append(Event(asm=asm_name, eid=eid, kind="INV", contig=c, start=s0, end=s1, length=L,
                               extra=""))
            eid += 1

        for _ in range(args.tra):
            c_src = rng.choice(contig_names)
            c_tgt = rng.choice(contig_names)
            src_len = len(contigs[c_src])
            if src_len < args.seg_len_min + 10:
                continue
            L = rng.randint(args.seg_len_min, min(args.seg_len_max, src_len - 1))
            s0 = pick_hotspot_pos(rng, src_len - L, args.subtel_frac)
            s0 = min(s0, src_len - L)
            s1 = s0 + L
            tgt_len = len(contigs[c_tgt])
            target = pick_hotspot_pos(rng, tgt_len, args.subtel_frac)
            apply_tra(contigs, c_src, s0, s1, c_tgt, target)
            truth.append(Event(asm=asm_name, eid=eid, kind="TRA", contig=c_src, start=s0, end=s1,
                               target_contig=c_tgt, target=target, length=L,
                               extra=""))
            eid += 1

        # write assembly as multi-contig FASTA
        asm_contigs = [(c, contigs[c]) for c in contig_names]
        out_path = os.path.join(args.outdir, f"{asm_name}.fa")
        write_fasta_multi(out_path, asm_contigs)

    # Write per-assembly truth files and a combined truth file
    truth_all_path = os.path.join(args.outdir, "truth_all.tsv")
    write_truth_tsv(truth_all_path, truth)
    per = {}
    for ev in truth:
        per.setdefault(ev.asm, []).append(ev)
    for asm_name, evs in per.items():
        write_truth_tsv(os.path.join(args.outdir, f"{asm_name}.truth.tsv"), evs)

    if args.manifest:
        # Manifest includes per-assembly SV counts + compact positions (joinable to truth.tsv via 'asm').
        # Positions are encoded as:
        #   INS/DEL:  contig:pos(len)
        #   DUP/INV:  contig:start-end
        #   TRA:      contig:start-end->target_contig:target
        # Columns:
        #   asm, fasta_path, ref_path, truth_path,
        #   n_INS, n_DEL, n_DUP, n_INV, n_TRA, n_SV_total,
        #   INS_pos, DEL_pos, DUP_pos, INV_pos, TRA_pos
        mpath = os.path.join(args.outdir, "manifest.tsv")
        outdir_abs = os.path.abspath(args.outdir)
        ref_abs = os.path.abspath("ref.fa")
        truth_abs = os.path.abspath(os.path.join(args.outdir, "truth.tsv"))

        # Aggregate per-assembly
        agg = {}
        for ev in truth:
            a = agg.setdefault(ev.asm, {
                "count": {"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0},
                "pos":   {"INS":[], "DEL":[], "DUP":[], "INV":[], "TRA":[]}
            })
            if ev.kind in a["count"]:
                a["count"][ev.kind] += 1

            if ev.kind == "INS":
                a["pos"]["INS"].append(f"{ev.contig}:{ev.pos}({ev.length})")
            elif ev.kind == "DEL":
                a["pos"]["DEL"].append(f"{ev.contig}:{ev.pos}({ev.length})")
            elif ev.kind == "DUP":
                a["pos"]["DUP"].append(f"{ev.contig}:{ev.start}-{ev.end}")
            elif ev.kind == "INV":
                a["pos"]["INV"].append(f"{ev.contig}:{ev.start}-{ev.end}")
            elif ev.kind == "TRA":
                a["pos"]["TRA"].append(f"{ev.contig}:{ev.start}-{ev.end}->{ev.target_contig}:{ev.target}")

        with open(mpath, "w") as mf:
            mf.write(
                "asm\tfasta_path\tref_path\ttruth_path"
                "\tn_INS\tn_DEL\tn_DUP\tn_INV\tn_TRA\tn_SV_total"
                "\tINS_pos\tDEL_pos\tDUP_pos\tINV_pos\tTRA_pos\n"
            )

            for i in range(args.n_genomes):
                asm_name = f"{args.prefix}{i:0{args.digits}d}"
                fasta_abs = os.path.join(outdir_abs, asm_name + ".fa")

                a = agg.get(asm_name, {"count":{"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0},
                                       "pos":{"INS":[],"DEL":[],"DUP":[],"INV":[],"TRA":[]}})
                c = a["count"]
                total = c["INS"] + c["DEL"] + c["DUP"] + c["INV"] + c["TRA"]

                # Compact, stable ordering (sort by contig then numeric position)
                def sort_key(tok: str):
                    # token starts with contig:...
                    contig, rest = tok.split(":", 1)
                    # extract first number in rest
                    num = 0
                    m = re.search(r"(\d+)", rest)
                    if m: num = int(m.group(1))
                    return (contig, num, rest)

                def join_list(xs):
                    xs2 = sorted(xs, key=sort_key)
                    return ";".join(xs2)

                mf.write(
                    f"{asm_name}\t{fasta_abs}\t{ref_abs}\t{truth_abs}"
                    f"\t{c['INS']}\t{c['DEL']}\t{c['DUP']}\t{c['INV']}\t{c['TRA']}\t{total}"
                    f"\t{join_list(a['pos']['INS'])}"
                    f"\t{join_list(a['pos']['DEL'])}"
                    f"\t{join_list(a['pos']['DUP'])}"
                    f"\t{join_list(a['pos']['INV'])}"
                    f"\t{join_list(a['pos']['TRA'])}\n"
                )
        # Count SVs per assembly from the in-memory truth list
        counts = {}
        for ev in truth:
            d = counts.setdefault(ev.asm, {"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0})
            if ev.kind in d:
                d[ev.kind] += 1

        with open(mpath, "w") as mf:
            mf.write("asm\tfasta_path\tref_path\ttruth_path\tn_INS\tn_DEL\tn_DUP\tn_INV\tn_TRA\tn_SV_total\n")
            for i in range(args.n_genomes):
                asm_name = f"{args.prefix}{i:0{args.digits}d}"
                fasta_abs = os.path.join(outdir_abs, asm_name + ".fa")
                d = counts.get(asm_name, {"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0})
                total = d["INS"] + d["DEL"] + d["DUP"] + d["INV"] + d["TRA"]
                mf.write(
                    f"{asm_name}\t{fasta_abs}\t{ref_abs}\t{truth_abs}"
                    f"\t{d['INS']}\t{d['DEL']}\t{d['DUP']}\t{d['INV']}\t{d['TRA']}\t{total}\n"
                )
    print("AMF-style simulation complete:")
    print("  ref.fa (multi-contig)")
    first = f"{args.prefix}{0:0{args.digits}d}"
    last = f"{args.prefix}{(args.n_genomes-1):0{args.digits}d}"
    print(f"  {os.path.join(args.outdir, first + '.fa')} ... {os.path.join(args.outdir, last + '.fa')}  (n={args.n_genomes})")
    print(f"  {os.path.join(args.outdir, 'truth.tsv')}")
    if args.manifest:
        # Manifest includes per-assembly SV counts + compact positions + per-assembly truth path.
        # Per-assembly truth is written to: <outdir>/<asm>.truth.tsv
        # Positions are encoded as:
        #   INS/DEL:  contig:pos(len)
        #   DUP/INV:  contig:start-end
        #   TRA:      contig:start-end->target_contig:target
        mpath = os.path.join(args.outdir, "manifest.tsv")
        outdir_abs = os.path.abspath(args.outdir)
        ref_abs = os.path.abspath("ref.fa")
        truth_all_abs = os.path.abspath(os.path.join(args.outdir, "truth_all.tsv"))

        # Aggregate per-assembly from the combined truth list
        agg = {}
        for ev in truth:
            a = agg.setdefault(ev.asm, {
                "count": {"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0},
                "pos":   {"INS":[], "DEL":[], "DUP":[], "INV":[], "TRA":[]}
            })
            if ev.kind in a["count"]:
                a["count"][ev.kind] += 1
            if ev.kind == "INS":
                a["pos"]["INS"].append(f"{ev.contig}:{ev.pos}({ev.length})")
            elif ev.kind == "DEL":
                a["pos"]["DEL"].append(f"{ev.contig}:{ev.pos}({ev.length})")
            elif ev.kind == "DUP":
                a["pos"]["DUP"].append(f"{ev.contig}:{ev.start}-{ev.end}")
            elif ev.kind == "INV":
                a["pos"]["INV"].append(f"{ev.contig}:{ev.start}-{ev.end}")
            elif ev.kind == "TRA":
                a["pos"]["TRA"].append(f"{ev.contig}:{ev.start}-{ev.end}->{ev.target_contig}:{ev.target}")

        def sort_key(tok: str):
            contig, rest = tok.split(":", 1)
            num = 0
            m = re.search(r"(\d+)", rest)
            if m: num = int(m.group(1))
            return (contig, num, rest)

        def join_list(xs):
            xs2 = sorted(xs, key=sort_key)
            return ";".join(xs2)

        with open(mpath, "w") as mf:
            mf.write(
                "asm\tfasta_path\tref_path\ttruth_path\ttruth_all_path"
                "\tn_INS\tn_DEL\tn_DUP\tn_INV\tn_TRA\tn_SV_total"
                "\tINS_pos\tDEL_pos\tDUP_pos\tINV_pos\tTRA_pos\n"
            )

            for i in range(args.n_genomes):
                asm_name = f"{args.prefix}{i:0{args.digits}d}"
                fasta_abs = os.path.join(outdir_abs, asm_name + ".fa")
                truth_abs = os.path.join(outdir_abs, asm_name + ".truth.tsv")

                a = agg.get(asm_name, {"count":{"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0},
                                       "pos":{"INS":[],"DEL":[],"DUP":[],"INV":[],"TRA":[]}})
                c = a["count"]
                total = c["INS"] + c["DEL"] + c["DUP"] + c["INV"] + c["TRA"]

                mf.write(
                    f"{asm_name}\t{fasta_abs}\t{ref_abs}\t{truth_abs}\t{truth_all_abs}"
                    f"\t{c['INS']}\t{c['DEL']}\t{c['DUP']}\t{c['INV']}\t{c['TRA']}\t{total}"
                    f"\t{join_list(a['pos']['INS'])}"
                    f"\t{join_list(a['pos']['DEL'])}"
                    f"\t{join_list(a['pos']['DUP'])}"
                    f"\t{join_list(a['pos']['INV'])}"
                    f"\t{join_list(a['pos']['TRA'])}\n"
                )
if __name__ == "__main__":
    main()
