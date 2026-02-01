#!/usr/bin/env bash
set -euo pipefail

# run_phylum_refcoords_bench_updated.sh
#
# Simulate fungal genomes per phylum with test_amf.py, run main.cpp SV caller,
# lift truth into reference-coordinate VCF (replay liftover), then compare truth vs calls.
#
# Usage:
#   PHYLUM_FILTER=Basidiomycota \
#   bash run_phylum_refcoords_bench_updated.sh /path/test_amf.py /path/main.cpp results_dir

#	chmod +x /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/run_phylum_refcoords_bench.sh

#	PHYLUM_FILTER=Basidiomycota \
#	bash /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/run_phylum_refcoords_bench.sh \
#	  /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/test_amf.py \
#	  /mnt/bmh01-rds/Shilpa_Group/2024/projects/fungi/AMF/main.cpp \
#	  results_basidio_candidate_split
#
# Inputs:
#   $1 = test_amf.py
#   $2 = main.cpp (or main_candidate_split.cpp)
#   $3 = output root directory
#
# Environment variables (optional):
#   PHYLUM_FILTER   : if set, only runs that phylum (exact match)
#   N_GENOMES       : override n-genomes per phylum (default per table)
#   TOTAL_LEN       : override total-len per phylum (default per table)
#   TOL_BP          : breakpoint tolerance for evaluation (default 5000)
#
# Outputs:
#   results_dir/<phylum>/
#     ref.fa
#     assemblies/
#       asm_*.fa
#       truth_all.tsv
#       manifest.tsv (if produced)
#     truth.refcoords.vcf
#     out.vcf
#     out.gfa
#     stdout.log
#     stderr.log
#     metrics.txt
#   results_dir/metrics.tsv

TEST_AMF="${1:-/mnt/data/test_amf.py}"
MAIN_CPP="${2:-/mnt/data/main.cpp}"
ROOT_OUT="${3:-results_phyla_refcoords}"

if [[ ! -f "$TEST_AMF" ]]; then
  echo "[error] missing test_amf.py at: $TEST_AMF" >&2
  exit 1
fi
if [[ ! -f "$MAIN_CPP" ]]; then
  echo "[error] missing main.cpp at: $MAIN_CPP" >&2
  exit 1
fi

# Make ROOT_OUT absolute to avoid nested path duplication after pushd
ROOT_OUT="$(realpath -m "$ROOT_OUT")"
mkdir -p "$ROOT_OUT"

BIN="$ROOT_OUT/fungi_pangenome"
echo "[info] Compiling $MAIN_CPP -> $BIN"
g++ -O3 -std=c++17 -pthread "$MAIN_CPP" -o "$BIN"

# -----------------------------------------------------------------------------
# truth_all.tsv -> truth.refcoords.vcf (reference coordinate space) via replay liftover
# -----------------------------------------------------------------------------
TRUTH_LIFT="$ROOT_OUT/truth_to_refcoords_vcf.py"
cat > "$TRUTH_LIFT" <<'PY'
#!/usr/bin/env python3
import sys
from dataclasses import dataclass
from typing import List, Optional, Dict, Tuple
from collections import defaultdict

@dataclass
class Block:
    cur0: int
    cur1: int
    ref_ctg: str
    ref0: Optional[int]
    ref1: Optional[int]
    strand: int  # +1 or -1

def load_fasta_lengths(ref_fa: str) -> Dict[str,int]:
    lens={}
    name=None
    cur=0
    with open(ref_fa) as f:
        for line in f:
            line=line.strip()
            if not line: continue
            if line.startswith(">"):
                if name is not None:
                    lens[name]=cur
                name=line[1:].split()[0]
                cur=0
            else:
                cur += len(line)
    if name is not None:
        lens[name]=cur
    return lens

def parse_truth_tsv(path: str):
    with open(path) as f:
        hdr=f.readline().rstrip("\n").split("\t")
        idx={h:i for i,h in enumerate(hdr)}
        def g(parts,k,default=""):
            return parts[idx[k]] if k in idx and idx[k] < len(parts) else default
        for line in f:
            line=line.rstrip("\n")
            if not line: continue
            p=line.split("\t")
            yield {
                "asm": g(p,"asm"),
                "event_id": int(g(p,"event_id","0") or "0"),
                "type": (g(p,"type") or g(p,"kind") or "").upper(),
                "contig": g(p,"contig") or "",
                "pos": int(g(p,"pos","0") or "0"),
                "start": int(g(p,"start","0") or "0"),
                "end": int(g(p,"end","0") or "0"),
                "target_contig": g(p,"target_contig") or g(p,"to_contig") or "",
                "target": int(g(p,"target","0") or "0"),
                "length": int(g(p,"length","0") or "0"),
            }

def blocks_total_len(blks: List[Block]) -> int:
    return blks[-1].cur1 if blks else 0

def split_at(blks: List[Block], x: int) -> None:
    i=0
    while i < len(blks):
        b=blks[i]
        if x <= b.cur0:
            return
        if b.cur0 < x < b.cur1:
            left_len = x - b.cur0
            if b.ref0 is None:
                left = Block(b.cur0, x, b.ref_ctg, None, None, b.strand)
                right = Block(x, b.cur1, b.ref_ctg, None, None, b.strand)
            else:
                if b.strand == +1:
                    left_ref0 = b.ref0
                    left_ref1 = b.ref0 + left_len
                    right_ref0 = left_ref1
                    right_ref1 = b.ref1
                else:
                    left_ref1 = b.ref1
                    left_ref0 = b.ref1 - left_len
                    right_ref1 = left_ref0
                    right_ref0 = b.ref0
                left = Block(b.cur0, x, b.ref_ctg, left_ref0, left_ref1, b.strand)
                right = Block(x, b.cur1, b.ref_ctg, right_ref0, right_ref1, b.strand)
            blks[i:i+1] = [left,right]
            return
        i += 1

def shift_blocks(blks: List[Block], start_idx: int, delta: int) -> None:
    if delta == 0: return
    for j in range(start_idx, len(blks)):
        blks[j].cur0 += delta
        blks[j].cur1 += delta

def map_point_to_ref(blks: List[Block], x: int) -> Tuple[str, Optional[int]]:
    if not blks:
        return ("", None)
    if x < 0: x = 0
    L = blocks_total_len(blks)
    if x > L: x = L
    probe = x-1 if x > 0 else x
    for b in blks:
        if b.cur0 <= probe < b.cur1:
            if b.ref0 is None:
                return (b.ref_ctg, None)
            off = probe - b.cur0
            if b.strand == +1:
                return (b.ref_ctg, b.ref0 + off)
            else:
                return (b.ref_ctg, b.ref1 - 1 - off)
    return (blks[0].ref_ctg, None)

def extract_range(blks: List[Block], s: int, e: int) -> List[Block]:
    frag=[]
    for b in blks:
        if b.cur1 <= s: continue
        if b.cur0 >= e: break
        if b.cur0 >= s and b.cur1 <= e:
            frag.append(b)
    out=[]
    cur=0
    for b in frag:
        L = b.cur1 - b.cur0
        out.append(Block(cur, cur+L, b.ref_ctg, b.ref0, b.ref1, b.strand))
        cur += L
    return out

def delete_range(blks: List[Block], s: int, e: int) -> None:
    i=0
    while i < len(blks):
        b=blks[i]
        if b.cur1 <= s:
            i += 1
            continue
        if b.cur0 >= e:
            break
        if b.cur0 >= s and b.cur1 <= e:
            del blks[i]
            continue
        i += 1
    shift = e - s
    j=0
    while j < len(blks) and blks[j].cur0 < e:
        j += 1
    for k in range(j, len(blks)):
        blks[k].cur0 -= shift
        blks[k].cur1 -= shift

def insert_blocks(blks: List[Block], x: int, ins: List[Block]) -> None:
    if not ins: return
    idx=0
    while idx < len(blks) and blks[idx].cur0 < x:
        idx += 1
    ins_len = ins[-1].cur1
    shift_blocks(blks, idx, ins_len)
    reb=[]
    for b in ins:
        reb.append(Block(x + b.cur0, x + b.cur1, b.ref_ctg, b.ref0, b.ref1, b.strand))
    blks[idx:idx] = reb

def invert_range(blks: List[Block], s: int, e: int) -> None:
    idx_s=None
    idx_e=None
    for i,b in enumerate(blks):
        if b.cur0 == s: idx_s=i
        if b.cur0 == e: idx_e=i; break
    if idx_s is None: return
    if idx_e is None: idx_e=len(blks)
    frag=blks[idx_s:idx_e]
    new=[]
    cur=0
    for b in reversed(frag):
        L = b.cur1 - b.cur0
        new.append(Block(cur, cur+L, b.ref_ctg, b.ref0, b.ref1, -b.strand))
        cur += L
    for i,b in enumerate(new):
        new[i] = Block(s + b.cur0, s + b.cur1, b.ref_ctg, b.ref0, b.ref1, b.strand)
    blks[idx_s:idx_e] = new

def make_inserted_block(ref_ctg: str, L: int) -> List[Block]:
    return [Block(0, L, ref_ctg, None, None, +1)] if L > 0 else []

def main():
    if len(sys.argv) != 4:
        print("usage: truth_to_refcoords_vcf.py ref.fa truth_all.tsv out.vcf", file=sys.stderr)
        sys.exit(1)

    ref_fa, truth_tsv, out_vcf = sys.argv[1:]
    ref_lens = load_fasta_lengths(ref_fa)

    evs_by_asm=defaultdict(list)
    for ev in parse_truth_tsv(truth_tsv):
        if not ev["asm"] or not ev["type"] or not ev["contig"]:
            continue
        evs_by_asm[ev["asm"]].append(ev)

    samples = sorted(evs_by_asm.keys())

    maps: Dict[str, Dict[str, List[Block]]] = {}
    for asm in samples:
        maps[asm] = {}
        for ctg,L in ref_lens.items():
            maps[asm][ctg] = [Block(0, L, ctg, 0, L, +1)]

    records=[]

    for asm in samples:
        for ev in evs_by_asm[asm]:
            t = ev["type"]
            if t not in ("INS","DEL","DUP","INV","TRA"):
                continue

            if t == "INS":
                c = ev["contig"]
                blks = maps[asm][c]
                p = ev["pos"]
                L = ev["length"]
                split_at(blks, p)
                ref_chr, refp0 = map_point_to_ref(blks, p)
                pos1 = (refp0 + 1) if refp0 is not None else 1
                info = f"SVTYPE=INS;SVLEN={L};END={pos1}"
                rid = f"{asm}:truth:INS:{ref_chr}:{pos1}:{ev['event_id']}"
                records.append((ref_chr, pos1, rid, info, asm))
                insert_blocks(blks, p, make_inserted_block(ref_chr, L))

            elif t == "DEL":
                c = ev["contig"]
                blks = maps[asm][c]
                p = ev["pos"]
                L = ev["length"]
                s = p
                e = p + L
                split_at(blks, s); split_at(blks, e)
                ref_chr, ref_s0 = map_point_to_ref(blks, s)
                _, ref_e0 = map_point_to_ref(blks, e)
                pos1 = (ref_s0 + 1) if ref_s0 is not None else 1
                end1 = (ref_e0 + 1) if ref_e0 is not None else pos1
                svlen = -abs(end1 - pos1)
                info = f"SVTYPE=DEL;SVLEN={svlen};END={end1}"
                rid = f"{asm}:truth:DEL:{ref_chr}:{pos1}:{ev['event_id']}"
                records.append((ref_chr, pos1, rid, info, asm))
                delete_range(blks, s, e)

            elif t == "DUP":
                c = ev["contig"]
                blks = maps[asm][c]
                s0 = ev["start"]
                s1 = ev["end"]
                target = ev["target"]
                split_at(blks, s0); split_at(blks, s1); split_at(blks, target)
                ref_chr, ref_s0 = map_point_to_ref(blks, s0)
                _, ref_s1 = map_point_to_ref(blks, s1)
                pos1 = (ref_s0 + 1) if ref_s0 is not None else 1
                end1 = (ref_s1 + 1) if ref_s1 is not None else pos1
                svlen = abs(end1 - pos1)
                info = f"SVTYPE=DUP;SVLEN={svlen};END={end1}"
                rid = f"{asm}:truth:DUP:{ref_chr}:{pos1}:{ev['event_id']}"
                records.append((ref_chr, pos1, rid, info, asm))
                frag = extract_range(blks, s0, s1)
                insert_blocks(blks, target, frag)

            elif t == "INV":
                c = ev["contig"]
                blks = maps[asm][c]
                s0 = ev["start"]
                s1 = ev["end"]
                split_at(blks, s0); split_at(blks, s1)
                ref_chr, ref_s0 = map_point_to_ref(blks, s0)
                _, ref_s1 = map_point_to_ref(blks, s1)
                pos1 = (ref_s0 + 1) if ref_s0 is not None else 1
                end1 = (ref_s1 + 1) if ref_s1 is not None else pos1
                svlen = abs(end1 - pos1)
                info = f"SVTYPE=INV;SVLEN={svlen};END={end1}"
                rid = f"{asm}:truth:INV:{ref_chr}:{pos1}:{ev['event_id']}"
                records.append((ref_chr, pos1, rid, info, asm))
                invert_range(blks, s0, s1)

            elif t == "TRA":
                c_src = ev["contig"]
                c_tgt = ev["target_contig"]
                s0 = ev["start"]
                s1 = ev["end"]
                target = ev["target"]

                blks_src = maps[asm][c_src]
                blks_tgt = maps[asm][c_tgt]
                split_at(blks_src, s0); split_at(blks_src, s1)
                split_at(blks_tgt, target)

                ref_chr1, ref_p10 = map_point_to_ref(blks_src, s0)
                ref_chr2, ref_p20 = map_point_to_ref(blks_tgt, target)

                pos1 = (ref_p10 + 1) if ref_p10 is not None else 1
                pos2 = (ref_p20 + 1) if ref_p20 is not None else 1
                info = f"SVTYPE=TRA;SVLEN=0;END={pos1};CHR2={ref_chr2};POS2={pos2}"
                rid = f"{asm}:truth:TRA:{ref_chr1}:{pos1}:{ev['event_id']}"
                records.append((ref_chr1, pos1, rid, info, asm))

                frag = extract_range(blks_src, s0, s1)
                delete_range(blks_src, s0, s1)
                insert_blocks(blks_tgt, target, frag)

    records.sort(key=lambda x: (x[0], x[1], x[2]))

    with open(out_vcf, "w") as o:
        o.write("##fileformat=VCFv4.2\n")
        o.write("##source=test_amf_truth_refcoords\n")
        o.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">\n")
        o.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">\n")
        o.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"SV length\">\n")
        o.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second chromosome for TRA\">\n")
        o.write("##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Second position for TRA\">\n")
        o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
        for s in samples:
            o.write("\t" + s)
        o.write("\n")

        for chrom,pos1,rid,info,asm in records:
            svt = info.split("SVTYPE=")[1].split(";")[0]
            alt = f"<{svt}>"
            o.write(f"{chrom}\t{pos1}\t{rid}\tN\t{alt}\t.\tPASS\t{info}\tGT")
            for s in samples:
                o.write("\t1" if s == asm else "\t0")
            o.write("\n")

if __name__ == "__main__":
    main()
PY
chmod +x "$TRUTH_LIFT"

# -----------------------------------------------------------------------------
# Evaluate: fuzzy match truth.refcoords.vcf vs out.vcf (both reference coords)
# -----------------------------------------------------------------------------
EVAL="$ROOT_OUT/eval_vcfs.py"
cat > "$EVAL" <<'PY'
#!/usr/bin/env python3
import sys
from collections import defaultdict

def parse_info(info):
    d={}
    for p in info.split(";"):
        if "=" in p:
            k,v=p.split("=",1)
            d[k]=v
    return d

def read_vcf(path):
    samples=[]
    recs=[]
    with open(path) as f:
        for line in f:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples=line.rstrip("\n").split("\t")[9:]
                continue
            if not line.strip() or line.startswith("#"): continue
            P=line.rstrip("\n").split("\t")
            chrom=P[0]; pos=int(P[1])
            info=parse_info(P[7])
            svt=info.get("SVTYPE","")
            end=int(info.get("END",str(pos)))
            chr2=info.get("CHR2","")
            pos2=int(info.get("POS2","0")) if "POS2" in info else 0
            present=set()
            for s,gt in zip(samples, P[9:]):
                gt=gt.strip()
                if gt not in ("0","0/0","./."):
                    present.add(s)
            recs.append((chrom,pos,end,svt,chr2,pos2,present))
    return recs

def bucket_key(r, binw):
    chrom,pos,end,svt,chr2,pos2,_=r
    if svt=="TRA":
        # Canonicalize TRA so that (CHROM,POS,CHR2,POS2) and swapped representations
        # hash to the same bucket.
        if (chrom, chr2) <= (chr2, chrom):
            c1,p1,c2,p2 = chrom,pos,chr2,pos2
        else:
            c1,p1,c2,p2 = chr2,pos2,chrom,pos
        return (c1,svt,c2,p1//binw,p2//binw)
    return (chrom,svt,chr2,0,pos//binw)

def is_match(t,c,tol):
    ctg_t,p_t,e_t,svt_t,chr2_t,pos2_t,_=t
    ctg_c,p_c,e_c,svt_c,chr2_c,pos2_c,_=c
    if svt_t!=svt_c:
        return False
    if svt_t=="TRA":
        # Compare canonicalized endpoints (unordered TRA).
        if (ctg_t, chr2_t) <= (chr2_t, ctg_t):
            t1, tp1, t2, tp2 = ctg_t, p_t, chr2_t, pos2_t
        else:
            t1, tp1, t2, tp2 = chr2_t, pos2_t, ctg_t, p_t
        if (ctg_c, chr2_c) <= (chr2_c, ctg_c):
            c1, cp1, c2, cp2 = ctg_c, p_c, chr2_c, pos2_c
        else:
            c1, cp1, c2, cp2 = chr2_c, pos2_c, ctg_c, p_c
        if t1!=c1 or t2!=c2:
            return False
        return abs(tp1-cp1)<=tol and abs(tp2-cp2)<=tol
    if ctg_t!=ctg_c:
        return False
    if ctg_t!=ctg_c:
        return False
    return abs(p_t-p_c)<=tol and abs(e_t-e_c)<=tol

def main():
    if len(sys.argv) < 3:
        print("usage: eval_vcfs.py truth.vcf out.vcf [tol=1000]", file=sys.stderr)
        sys.exit(1)
    truth_vcf=sys.argv[1]
    out_vcf=sys.argv[2]
    tol=int(sys.argv[3]) if len(sys.argv)>3 else 1000
    binw=max(200,tol)

    truth=read_vcf(truth_vcf)
    calls=read_vcf(out_vcf)

    idx=defaultdict(list)
    for t in truth:
        idx[bucket_key(t,binw)].append(t)

    used=set()
    TP=FP=FN=0

    for c in calls:
        for s in c[6]:
            found=False
            chrom,svt,chr2,pos2b,b = bucket_key(c,binw)
            for db in (-1,0,1):
                cand = idx.get((chrom,svt,chr2,pos2b,b+db), [])
                for t in cand:
                    if s not in t[6]:
                        continue
                    tid=(id(t),s)
                    if tid in used:
                        continue
                    if is_match(t,c,tol):
                        used.add(tid)
                        TP += 1
                        found=True
                        break
                if found:
                    break
            if not found:
                FP += 1

    for t in truth:
        for s in t[6]:
            if (id(t),s) not in used:
                FN += 1

    prec = TP/(TP+FP) if (TP+FP) else 0.0
    rec  = TP/(TP+FN) if (TP+FN) else 0.0

    print(f"TP\t{TP}")
    print(f"FP\t{FP}")
    print(f"FN\t{FN}")
    print(f"precision\t{prec:.6f}")
    print(f"recall\t{rec:.6f}")

if __name__ == "__main__":
    main()
PY
chmod +x "$EVAL"

# -----------------------------------------------------------------------------
# Phylum presets (heuristics)
# Columns: PHYLUM GC REP SUBTEL_REP SUBTEL_FRAC TOTAL_LEN N_CONTIGS N_GENOMES
# -----------------------------------------------------------------------------
PHYLUMS=(
  "Basidiomycota        0.50 0.18 0.35 0.10 200000 12 5"
  "Ascomycota           0.48 0.12 0.25 0.10 200000 12 5"
  "Mucoromycota         0.40 0.20 0.40 0.12 200000 12 5"
  "Zoopagomycota        0.42 0.22 0.45 0.12 200000 12 5"
  "Chytridiomycota      0.45 0.15 0.30 0.10 200000 12 5"
  "Blastocladiomycota   0.46 0.16 0.32 0.10 200000 12 5"
  "Microsporidia        0.25 0.05 0.10 0.06 200000 12 5"
  "Cryptomycota         0.38 0.20 0.40 0.12 200000 12 5"
)

# Allow overrides
PHYLUM_FILTER="${PHYLUM_FILTER:-}"
TOL_BP="${TOL_BP:-5000}"
OVR_NGEN="${N_GENOMES:-}"
OVR_TLEN="${TOTAL_LEN:-}"

# SV profile per genome
INS=2
DEL=2
DUP=1
INV=1
TRA=1

# Event length profile (bp)
# The caller is a prototype; using mid/large events makes the benchmark stable and
# ensures we can meaningfully score DUP/INV/TRA rather than getting swamped by repeat-driven noise.
INS_LEN_MIN=800
INS_LEN_MAX=1200
DEL_LEN_MIN=800
DEL_LEN_MAX=1500
SEG_LEN_MIN=2000   # used for DUP/INV/TRA segments
SEG_LEN_MAX=5000

# Caller parameters
MIN_SEEDS=5
TOP_CONTIGS=1000
SV_MIN=500          # minimum SV size for INV/DUP (and general heuristics)
SV_MIN_INDEL=500    # minimum INS/DEL size extracted from alignment

# Candidate/split settings (for main_candidate_split.cpp; harmless if ignored)
CANDIDATES=3
MAPQ_RATIO=1.15
SPLIT_MAP=1

# Extra runtime-stabilization knobs (if supported; ignored otherwise)
VCF_CHECKPOINT=1

echo -e "phylum\tTP\tFP\tFN\tprecision\trecall" > "$ROOT_OUT/metrics.tsv"

for row in "${PHYLUMS[@]}"; do
  # shellcheck disable=SC2206
  F=($row)
  PHYLUM="${F[0]}"

  if [[ -n "$PHYLUM_FILTER" && "$PHYLUM" != "$PHYLUM_FILTER" ]]; then
    continue
  fi

  GC="${F[1]}"
  REP="${F[2]}"
  SUBREP="${F[3]}"
  SUBFRAC="${F[4]}"
  TLEN="${F[5]}"
  NCONTIG="${F[6]}"
  NGEN="${F[7]}"

  if [[ -n "$OVR_TLEN" ]]; then TLEN="$OVR_TLEN"; fi
  if [[ -n "$OVR_NGEN" ]]; then NGEN="$OVR_NGEN"; fi

  OUT="$ROOT_OUT/$PHYLUM"
  mkdir -p "$OUT"
  pushd "$OUT" >/dev/null

  ASM="assemblies"
  mkdir -p "$ASM"

  echo "[info] ---- $PHYLUM ----"
  echo "[info] Simulating: total-len=$TLEN n-contigs=$NCONTIG n-genomes=$NGEN"

  python3 "$TEST_AMF" \
    --seed 42 \
    --outdir "$ASM" \
    --total-len "$TLEN" \
    --n-contigs "$NCONTIG" \
    --gc "$GC" \
    --repeat-density "$REP" \
    --subtel-repeat-density "$SUBREP" \
    --subtel-frac "$SUBFRAC" \
    --n-genomes "$NGEN" \
    --ins "$INS" --del "$DEL" --dup "$DUP" --inv "$INV" --tra "$TRA" \
    --ins-len-min "$INS_LEN_MIN" --ins-len-max "$INS_LEN_MAX" \
    --del-len-min "$DEL_LEN_MIN" --del-len-max "$DEL_LEN_MAX" \
    --seg-len-min "$SEG_LEN_MIN" --seg-len-max "$SEG_LEN_MAX" \
    --manifest

  if [[ ! -f "ref.fa" ]]; then
    echo "[error] ref.fa missing for $PHYLUM" >&2
    exit 1
  fi
  if [[ ! -f "$ASM/truth_all.tsv" ]]; then
    echo "[error] truth_all.tsv missing at $ASM/truth_all.tsv for $PHYLUM" >&2
    exit 1
  fi

  echo "[info] Converting truth_all.tsv -> truth.refcoords.vcf"
  python3 "$TRUTH_LIFT" "ref.fa" "$ASM/truth_all.tsv" "$OUT/truth.refcoords.vcf"

  echo "[info] Running caller -> out.vcf"
  "$BIN" \
    --ref "$OUT/ref.fa" \
    --asm-dir "$ASM" \
    --out "$OUT/out.gfa" \
    --vcf "$OUT/out.vcf" \
    --sv "$SV_MIN" \
    --sv-indel "$SV_MIN_INDEL" \
    --min-seeds "$MIN_SEEDS" \
    --top-contigs "$TOP_CONTIGS" \
    --min-contig 1 \
    --candidates "$CANDIDATES" \
    --mapq-ratio "$MAPQ_RATIO" \
    --split-map "$SPLIT_MAP" \
    --vcf-checkpoint "$VCF_CHECKPOINT" \
    > "$OUT/stdout.log" 2> "$OUT/stderr.log" || true

  if [[ ! -f "$OUT/out.vcf" ]]; then
    echo "[warn] out.vcf missing for $PHYLUM (see stderr.log)" >&2
    popd >/dev/null
    continue
  fi

  echo "[info] Evaluating out.vcf vs truth.refcoords.vcf (tol=${TOL_BP})"
  python3 "$EVAL" "$OUT/truth.refcoords.vcf" "$OUT/out.vcf" "$TOL_BP" | tee "$OUT/metrics.txt"

  TP=$(awk '$1=="TP"{print $2}' "$OUT/metrics.txt")
  FP=$(awk '$1=="FP"{print $2}' "$OUT/metrics.txt")
  FN=$(awk '$1=="FN"{print $2}' "$OUT/metrics.txt")
  PREC=$(awk '$1=="precision"{print $2}' "$OUT/metrics.txt")
  REC=$(awk '$1=="recall"{print $2}' "$OUT/metrics.txt")

  echo -e "${PHYLUM}\t${TP}\t${FP}\t${FN}\t${PREC}\t${REC}" | tee -a "$ROOT_OUT/metrics.tsv"

  popd >/dev/null
done

echo "[done] Summary written to: $ROOT_OUT/metrics.tsv"

