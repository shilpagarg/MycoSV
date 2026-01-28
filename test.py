#!/usr/bin/env python3
"""
Simulate two assemblies from a reference genome with:
  - Insertions (INS) and deletions (DEL) applied relative to the reference coordinates
  - Duplications (DUP), inversions (INV), and translocations (TRA) applied as rearrangements
    on the evolving sequence.

Outputs:
  ref.fa
  assemblies/asm_a.fa
  assemblies/asm_b.fa
  assemblies/truth.tsv   (events used to generate each assembly)

Notes:
- Coordinates in truth.tsv are 0-based half-open on the sequence *at the time the event is applied*.
- Indels are applied first (reference-based positions), then rearrangements.
"""
from __future__ import annotations

import argparse
import os
import random
from dataclasses import dataclass
from typing import List, Sequence, Tuple

DNA = ["A", "C", "G", "T"]


# ----------------------------
# Utility functions
# ----------------------------
def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choices(DNA, k=length))


def revcomp(seq: str) -> str:
    comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(comp.get(b, "N") for b in reversed(seq))


def write_fasta(path: str, name: str, seq: str) -> None:
    with open(path, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")


def clamp(x: int, lo: int, hi: int) -> int:
    return lo if x < lo else hi if x > hi else x


# ----------------------------
# Event representation
# ----------------------------
@dataclass
class Event:
    asm: str
    eid: int
    kind: str  # INS, DEL, DUP, INV, TRA
    pos: int = 0
    length: int = 0
    start: int = 0
    end: int = 0
    target: int = 0
    extra: str = ""


# ----------------------------
# SV operations
# ----------------------------
def apply_ins_del_ref_based(
    ref: str,
    indels: Sequence[Tuple[str, int, int]],
    rng: random.Random,
    asm_name: str,
    eid_start: int,
    truth: List[Event],
) -> Tuple[str, int]:
    """
    Apply ("INS", ref_pos, len) and ("DEL", ref_pos, len) relative to the original reference.
    We apply them in ascending ref_pos with an offset to keep them deterministic.

    Returns (sequence, next_eid)
    """
    seq = ref
    offset = 0
    eid = eid_start

    for kind, ref_pos, length in sorted(indels, key=lambda x: x[1]):
        pos = ref_pos + offset
        pos = clamp(pos, 0, len(seq))

        if kind == "INS":
            ins_seq = random_dna(rng, length)
            seq = seq[:pos] + ins_seq + seq[pos:]
            truth.append(Event(asm=asm_name, eid=eid, kind="INS", pos=pos, length=length, extra=f"ref_pos={ref_pos}"))
            eid += 1
            offset += length

        elif kind == "DEL":
            del_len = clamp(length, 0, len(seq) - pos)
            seq = seq[:pos] + seq[pos + del_len :]
            truth.append(Event(asm=asm_name, eid=eid, kind="DEL", pos=pos, length=del_len, extra=f"ref_pos={ref_pos}"))
            eid += 1
            offset -= del_len

        else:
            raise ValueError(f"Unknown indel kind: {kind}")

    return seq, eid


def apply_dup(seq: str, s: int, e: int, t: int) -> str:
    frag = seq[s:e]
    return seq[:t] + frag + seq[t:]


def apply_inv(seq: str, s: int, e: int) -> str:
    frag = revcomp(seq[s:e])
    return seq[:s] + frag + seq[e:]


def apply_tra(seq: str, s: int, e: int, t: int) -> str:
    frag = seq[s:e]
    seq2 = seq[:s] + seq[e:]  # excise
    if t > s:
        t -= (e - s)  # adjust if target was after cut
    t = clamp(t, 0, len(seq2))
    return seq2[:t] + frag + seq2[t:]


def apply_rearrangements_evolving(
    seq: str,
    rearr_ops: Sequence[Tuple],
    asm_name: str,
    eid_start: int,
    truth: List[Event],
) -> Tuple[str, int]:
    """
    Apply DUP/INV/TRA on the *current* sequence.
    ops may be:
      ("DUP", start, end, target_pos)
      ("INV", start, end)
      ("TRA", start, end, target_pos)

    Coordinates are interpreted on the current sequence at time of application.
    Returns (sequence, next_eid)
    """
    eid = eid_start
    for op in rearr_ops:
        kind = op[0]
        L = len(seq)

        if kind == "DUP":
            _, s, e, t = op
            s = clamp(s, 0, L)
            e = clamp(e, s, L)
            t = clamp(t, 0, L)
            if e - s <= 0:
                continue
            seq = apply_dup(seq, s, e, t)
            truth.append(Event(asm=asm_name, eid=eid, kind="DUP", start=s, end=e, target=t, length=(e - s)))
            eid += 1

        elif kind == "INV":
            _, s, e = op
            s = clamp(s, 0, L)
            e = clamp(e, s, L)
            if e - s <= 1:
                continue
            seq = apply_inv(seq, s, e)
            truth.append(Event(asm=asm_name, eid=eid, kind="INV", start=s, end=e, length=(e - s)))
            eid += 1

        elif kind == "TRA":
            _, s, e, t = op
            s = clamp(s, 0, L)
            e = clamp(e, s, L)
            if e - s <= 0:
                continue
            t = clamp(t, 0, L)
            if t >= s and t <= e:
                t = clamp(e + 1, 0, L)  # avoid trivial no-op
            seq = apply_tra(seq, s, e, t)
            truth.append(Event(asm=asm_name, eid=eid, kind="TRA", start=s, end=e, target=t, length=(e - s)))
            eid += 1

        else:
            raise ValueError(f"Unknown rearrangement kind: {kind}")

    return seq, eid


# ----------------------------
# Random event generators
# ----------------------------
def gen_indels(
    rng: random.Random,
    ref_len: int,
    n_ins: int,
    n_del: int,
    ins_len: Tuple[int, int],
    del_len: Tuple[int, int],
) -> List[Tuple[str, int, int]]:
    ops: List[Tuple[str, int, int]] = []
    for _ in range(n_ins):
        pos = rng.randrange(0, ref_len + 1)
        length = rng.randrange(ins_len[0], ins_len[1] + 1)
        ops.append(("INS", pos, length))
    for _ in range(n_del):
        pos = rng.randrange(0, ref_len + 1)
        length = rng.randrange(del_len[0], del_len[1] + 1)
        ops.append(("DEL", pos, length))
    return ops


def gen_rearrangements(
    rng: random.Random,
    seq_len: int,
    n_dup: int,
    n_inv: int,
    n_tra: int,
    seg_len: Tuple[int, int],
) -> List[Tuple]:
    """
    Generate rearrangement ops with coordinates on an *approximate* post-indel sequence length.
    (They are clamped safely at application time.)
    """
    ops: List[Tuple] = []

    def rand_interval(L: int) -> Tuple[int, int]:
        ln = rng.randrange(seg_len[0], min(seg_len[1], max(seg_len[0], L)) + 1)
        if ln >= L:
            return 0, L
        s = rng.randrange(0, L - ln + 1)
        return s, s + ln

    for _ in range(n_dup):
        s, e = rand_interval(seq_len)
        t = e if rng.random() < 0.5 else rng.randrange(0, seq_len + 1)  # 50% tandem, 50% dispersed
        ops.append(("DUP", s, e, t))

    for _ in range(n_inv):
        s, e = rand_interval(seq_len)
        ops.append(("INV", s, e))

    for _ in range(n_tra):
        s, e = rand_interval(seq_len)
        t = rng.randrange(0, seq_len + 1)
        ops.append(("TRA", s, e, t))

    rng.shuffle(ops)
    return ops


# ----------------------------
# Main
# ----------------------------
def write_truth_tsv(path: str, truth: Sequence[Event]) -> None:
    header = ["asm", "event_id", "type", "pos", "start", "end", "target", "length", "extra"]
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for ev in truth:
            f.write(
                "\t".join(
                    [
                        ev.asm,
                        str(ev.eid),
                        ev.kind,
                        str(ev.pos),
                        str(ev.start),
                        str(ev.end),
                        str(ev.target),
                        str(ev.length),
                        ev.extra,
                    ]
                )
                + "\n"
            )


def simulate_one(
    ref: str,
    rng: random.Random,
    asm_name: str,
    indels: List[Tuple[str, int, int]],
    rearr: List[Tuple],
    eid_start: int,
    truth: List[Event],
) -> Tuple[str, int]:
    seq, eid = apply_ins_del_ref_based(ref, indels, rng, asm_name, eid_start, truth)
    seq, eid = apply_rearrangements_evolving(seq, rearr, asm_name, eid, truth)
    return seq, eid


def main() -> None:
    ap = argparse.ArgumentParser(description="Simulate assemblies with INS/DEL/DUP/INV/TRA.")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--ref-len", type=int, default=100_000)
    ap.add_argument("--outdir", type=str, default="assemblies")

    # Indels per assembly
    ap.add_argument("--ins", type=int, default=2, help="number of insertions per assembly")
    ap.add_argument("--del", dest="dels", type=int, default=2, help="number of deletions per assembly")
    ap.add_argument("--ins-len-min", type=int, default=200)
    ap.add_argument("--ins-len-max", type=int, default=1200)
    ap.add_argument("--del-len-min", type=int, default=200)
    ap.add_argument("--del-len-max", type=int, default=1500)

    # Rearrangements per assembly
    ap.add_argument("--dup", type=int, default=1)
    ap.add_argument("--inv", type=int, default=1)
    ap.add_argument("--tra", type=int, default=1)
    ap.add_argument("--seg-len-min", type=int, default=500)
    ap.add_argument("--seg-len-max", type=int, default=5000)

    args = ap.parse_args()
    rng = random.Random(args.seed)

    # Reference
    ref = random_dna(rng, args.ref_len)
    os.makedirs(args.outdir, exist_ok=True)
    write_fasta("ref.fa", "ref", ref)

    truth: List[Event] = []
    eid = 1

    for asm_name in ("asm_a", "asm_b"):
        indels = gen_indels(
            rng,
            ref_len=len(ref),
            n_ins=args.ins,
            n_del=args.dels,
            ins_len=(args.ins_len_min, args.ins_len_max),
            del_len=(args.del_len_min, args.del_len_max),
        )

        approx_len = len(ref) + sum(l for k, _, l in indels if k == "INS") - sum(l for k, _, l in indels if k == "DEL")
        approx_len = max(1, approx_len)

        rearr = gen_rearrangements(
            rng,
            seq_len=approx_len,
            n_dup=args.dup,
            n_inv=args.inv,
            n_tra=args.tra,
            seg_len=(args.seg_len_min, args.seg_len_max),
        )

        seq, eid = simulate_one(ref, rng, asm_name, indels, rearr, eid, truth)
        write_fasta(os.path.join(args.outdir, f"{asm_name}.fa"), asm_name, seq)

    write_truth_tsv(os.path.join(args.outdir, "truth.tsv"), truth)

    print("Simulation complete:")
    print("  ref.fa")
    print(f"  {args.outdir}/asm_a.fa")
    print(f"  {args.outdir}/asm_b.fa")
    print(f"  {args.outdir}/truth.tsv")


if __name__ == "__main__":
    main()

