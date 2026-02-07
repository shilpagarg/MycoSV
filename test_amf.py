#!/usr/bin/env python3
"""
AMF-style (arbuscular mycorrhizal fungi) *toy* genome simulator
==============================================================

python test_amf.py --total-len 100000 --n-contigs 10 --n-genomes 4 --outdir assemblies --manifest

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
import concurrent.futures
import threading

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
        f.write("\t".join(header) + "\n")
        for ev in truth:
            f.write("\t".join([
                ev.asm, str(ev.eid), ev.kind, ev.contig,
                str(ev.pos), str(ev.start), str(ev.end),
                ev.target_contig, str(ev.target), str(ev.length),
                ev.extra
            ]) + "\n")

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
# Simulation driver (streaming; scales to thousands of genomes)
# ----------------------------
def _write_truth_header(fh) -> None:
    fh.write("\t".join(["asm","event_id","type","contig","pos","start","end","target_contig","target","length","extra"]) + "\n")

def _write_truth_event(fh, ev: Event) -> None:
    fh.write("\t".join([
        ev.asm, str(ev.eid), ev.kind, ev.contig,
        str(ev.pos), str(ev.start), str(ev.end),
        ev.target_contig, str(ev.target), str(ev.length),
        ev.extra
    ]) + "\n")

def main() -> None:
    ap = argparse.ArgumentParser(description="AMF-style toy SV simulator (supports thousands of genomes efficiently).")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outdir", type=str, default="assemblies")

    # Size and contigs
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
    ap.add_argument("--start-idx", type=int, default=0, help="start index for assembly names (useful for batching)")
    ap.add_argument("--prefix", type=str, default="asm_")
    ap.add_argument("--digits", type=int, default=4)
    ap.add_argument("--threads", type=int, default=1, help="parallelism for per-assembly simulation (default 1)")

    # Truth/manifest outputs (streaming)
    ap.add_argument("--truth-all-mode", choices=["write","append","none"], default="write",
                    help="combined truth file mode for <outdir>/truth_all.tsv (default: write). Use 'none' for fastest runs.")
    ap.add_argument("--per-asm-truth", action="store_true", default=True,
                    help="write <outdir>/<asm>.truth.tsv for each genome (default: on)")
    ap.add_argument("--no-per-asm-truth", dest="per_asm_truth", action="store_false",
                    help="disable per-assembly truth files")
    ap.add_argument("--manifest", action="store_true", help="write <outdir>/manifest.tsv (streamed)")
    ap.add_argument("--manifest-positions", action="store_true",
                    help="include compact position strings in manifest (can be large; off by default)")
    ap.add_argument("--ref-out", type=str, default="ref.fa",
                    help="where to write the reference FASTA (default: ref.fa in current working directory)")

    # SV counts per assembly
    ap.add_argument("--ins", type=int, default=2)
    ap.add_argument("--del", dest="dels", type=int, default=2)
    ap.add_argument("--dup", type=int, default=1)
    ap.add_argument("--inv", type=int, default=1)
    ap.add_argument("--tra", type=int, default=1)

    # SV size ranges
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

    # Write reference (explicit destination)
    write_fasta_multi(args.ref_out, ref_contigs)

    # Base reference dict (for per-assembly copy)
    ref_dict: Dict[str,str] = {name: seq for name, seq in ref_contigs}
    contig_names = [name for name, _ in ref_contigs]

    # Open combined truth file (optional), streaming
    truth_all_fh = None
    truth_all_path = os.path.join(args.outdir, "truth_all.tsv")
    if args.truth_all_mode != "none":
        mode = "w" if args.truth_all_mode == "write" else "a"
        first_write = (mode == "w") or (not os.path.exists(truth_all_path)) or (os.path.getsize(truth_all_path) == 0)
        truth_all_fh = open(truth_all_path, mode)
        if first_write:
            _write_truth_header(truth_all_fh)

    # Open manifest (optional), streaming
    manifest_fh = None
    if args.manifest:
        mpath = os.path.join(args.outdir, "manifest.tsv")
        manifest_fh = open(mpath, "w")
        # Keep manifest compact by default for large N
        cols = [
            "asm","fasta_path","ref_path","truth_path","truth_all_path",
            "n_INS","n_DEL","n_DUP","n_INV","n_TRA","n_SV_total"
        ]
        if args.manifest_positions:
            cols += ["INS_pos","DEL_pos","DUP_pos","INV_pos","TRA_pos"]
        manifest_fh.write("\t".join(cols) + "\n")

    outdir_abs = os.path.abspath(args.outdir)
    ref_abs = os.path.abspath(args.ref_out)
    truth_all_abs = os.path.abspath(truth_all_path) if args.truth_all_mode != "none" else ""

    io_lock = threading.Lock()
    eid_lock = threading.Lock()
    eid = 1

    def next_eid() -> int:
        nonlocal eid
        with eid_lock:
            x = eid
            eid += 1
            return x

    def sort_key(tok: str):
        contig, rest = tok.split(":", 1)
        num = 0
        m = re.search(r"(\d+)", rest)
        if m: num = int(m.group(1))
        return (contig, num, rest)

    def _run_one(i: int) -> None:
        # Use a per-assembly RNG for determinism + thread safety
        lrng = random.Random(args.seed + i * 1000003)

        asm_name = f"{args.prefix}{i:0{args.digits}d}"

        # Working copy per genome
        contigs = dict(ref_dict)

        # Per-assembly truth file (optional)
        per_truth_path = os.path.join(args.outdir, f"{asm_name}.truth.tsv")
        per_truth_fh = None
        if args.per_asm_truth:
            per_truth_fh = open(per_truth_path, "w")
            _write_truth_header(per_truth_fh)

        # Per-genome manifest stats (streamed)
        counts = {"INS":0,"DEL":0,"DUP":0,"INV":0,"TRA":0}
        posbags = {"INS":[], "DEL":[], "DUP":[], "INV":[], "TRA":[]}

        def emit(ev: Event):
            # Write truth (truth_all.tsv is shared across threads; per-assembly truth is private)
            if per_truth_fh is not None:
                _write_truth_event(per_truth_fh, ev)
            if truth_all_fh is not None:
                with io_lock:
                    _write_truth_event(truth_all_fh, ev)

            if ev.kind in counts:
                counts[ev.kind] += 1
                if args.manifest_positions:
                    if ev.kind == "INS":
                        posbags["INS"].append(f"{ev.contig}:{ev.pos}({ev.length})")
                    elif ev.kind == "DEL":
                        posbags["DEL"].append(f"{ev.contig}:{ev.pos}({ev.length})")
                    elif ev.kind == "DUP":
                        posbags["DUP"].append(f"{ev.contig}:{ev.start}-{ev.end}")
                    elif ev.kind == "INV":
                        posbags["INV"].append(f"{ev.contig}:{ev.start}-{ev.end}")
                    elif ev.kind == "TRA":
                        posbags["TRA"].append(f"{ev.contig}:{ev.start}-{ev.end}->{ev.target_contig}:{ev.target}")

        # Apply DUP/INV/TRA first (avoid nested composition)
        for _ in range(args.dup):
            c = lrng.choice(contig_names)
            clen = len(contigs[c])
            if clen < args.seg_len_min + 10:
                continue
            L = lrng.randint(args.seg_len_min, min(args.seg_len_max, clen - 1))
            s0 = pick_hotspot_pos(lrng, clen - L, args.subtel_frac)
            s0 = min(s0, clen - L)
            s1 = s0 + L
            if lrng.random() < 0.60:
                target = min(clen, s1 + lrng.randint(0, 2000))
            else:
                target = lrng.randint(0, clen)
            apply_dup(contigs, c, s0, s1, target)
            emit(Event(asm=asm_name, eid=next_eid(), kind="DUP", contig=c, start=s0, end=s1, target=target, length=L,
                       extra="mode=tandem" if target >= s1 and target <= s1+2000 else "mode=dispersed"))

        for _ in range(args.inv):
            c = lrng.choice(contig_names)
            clen = len(contigs[c])
            if clen < args.seg_len_min + 10:
                continue
            L = lrng.randint(args.seg_len_min, min(args.seg_len_max, clen - 1))
            s0 = pick_hotspot_pos(lrng, clen - L, args.subtel_frac)
            s0 = min(s0, clen - L)
            s1 = s0 + L
            apply_inv(contigs, c, s0, s1)
            emit(Event(asm=asm_name, eid=next_eid(), kind="INV", contig=c, start=s0, end=s1, length=L, extra=""))

        for _ in range(args.tra):
            c_src = lrng.choice(contig_names)
            c_tgt = lrng.choice(contig_names)
            src_len = len(contigs[c_src])
            if src_len < args.seg_len_min + 10:
                continue
            L = lrng.randint(args.seg_len_min, min(args.seg_len_max, src_len - 1))
            s0 = pick_hotspot_pos(lrng, src_len - L, args.subtel_frac)
            s0 = min(s0, src_len - L)
            s1 = s0 + L
            tgt_len = len(contigs[c_tgt])
            target = pick_hotspot_pos(lrng, tgt_len, args.subtel_frac)
            apply_tra(contigs, c_src, s0, s1, c_tgt, target)
            emit(Event(asm=asm_name, eid=next_eid(), kind="TRA", contig=c_src, start=s0, end=s1,
                       target_contig=c_tgt, target=target, length=L, extra=""))

        # INS/DEL after
        for _ in range(args.ins):
            c = lrng.choice(contig_names)
            pos = pick_hotspot_pos(lrng, len(contigs[c]), args.subtel_frac)
            L = lrng.randint(args.ins_len_min, args.ins_len_max)
            ins = random_dna(lrng, L, gc=args.gc)
            apply_ins(contigs, c, pos, ins)
            emit(Event(asm=asm_name, eid=next_eid(), kind="INS", contig=c, pos=pos, length=L,
                       extra="hotspot=subtel" if pos < int(len(contigs[c])*args.subtel_frac) or pos > int(len(contigs[c])*(1-args.subtel_frac)) else "hotspot=any"))

        for _ in range(args.dels):
            c = lrng.choice(contig_names)
            pos = pick_hotspot_pos(lrng, len(contigs[c]), args.subtel_frac)
            L = lrng.randint(args.del_len_min, args.del_len_max)
            dl = apply_del(contigs, c, pos, L)
            emit(Event(asm=asm_name, eid=next_eid(), kind="DEL", contig=c, pos=pos, length=dl,
                       extra="hotspot=subtel" if pos < int(len(contigs[c])*args.subtel_frac) or pos > int(len(contigs[c])*(1-args.subtel_frac)) else "hotspot=any"))

        # Write assembly as multi-contig FASTA
        asm_contigs = [(c, contigs[c]) for c in contig_names]
        out_path = os.path.join(args.outdir, f"{asm_name}.fa")
        write_fasta_multi(out_path, asm_contigs)

        if per_truth_fh is not None:
            per_truth_fh.close()

        # Stream manifest line (if enabled)
        if manifest_fh is not None:
            with io_lock:
                cols = [
                    asm_name,
                    os.path.join(outdir_abs, f"{asm_name}.fa"),
                    ref_abs,
                    os.path.join(outdir_abs, f"{asm_name}.truth.tsv") if args.per_asm_truth else "",
                    truth_all_abs,
                    str(counts["INS"]), str(counts["DEL"]), str(counts["DUP"]), str(counts["INV"]), str(counts["TRA"]),
                    str(sum(counts.values()))
                ]
                if args.manifest_positions:
                    for k in ["INS","DEL","DUP","INV","TRA"]:
                        posbags[k].sort(key=sort_key)
                    cols += [
                        ",".join(posbags["INS"]),
                        ",".join(posbags["DEL"]),
                        ",".join(posbags["DUP"]),
                        ",".join(posbags["INV"]),
                        ",".join(posbags["TRA"]),
                    ]
                manifest_fh.write("\t".join(cols) + "\n")
    

    idxs = list(range(args.start_idx, args.start_idx + args.n_genomes))
    if args.threads <= 1:
        for i in idxs:
            _run_one(i)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as ex:
            list(ex.map(_run_one, idxs))

    if truth_all_fh is not None:
        truth_all_fh.close()
    if manifest_fh is not None:
        manifest_fh.close()

    # Final message
    first = f"{args.prefix}{args.start_idx:0{args.digits}d}"
    last = f"{args.prefix}{(args.start_idx + args.n_genomes - 1):0{args.digits}d}"
    print("AMF-style simulation complete:")
    print(f"  {args.ref_out} (multi-contig)")
    print(f"  {os.path.join(args.outdir, first + '.fa')} ... {os.path.join(args.outdir, last + '.fa')}  (n={args.n_genomes})")
    if args.truth_all_mode != "none":
        print(f"  {os.path.join(args.outdir, 'truth_all.tsv')}")
    if args.per_asm_truth:
        print(f"  {os.path.join(args.outdir, first + '.truth.tsv')} ...")
    if args.manifest:
        print(f"  {os.path.join(args.outdir, 'manifest.tsv')}")

if __name__ == '__main__':
    main()
