// main.cpp
// Build: g++ -O3 -std=c++17 -pthread -o fungi_pangenome main.cpp
//
// What this prototype does:
//  - Builds a reference backbone graph (GFA) by segmenting the reference FASTA (all contigs).
//  - For each assembly contig, uses syncmer seeding + simple chaining to infer a mapped locus.
//  - Runs a banded GLOBAL affine-gap alignment (correct state-specific traceback).
//  - Calls large INDEL SVs (>= --sv) from the resulting CIGAR.
//  - Inserts SV alleles as detours into the graph + edits a sample path.
//
// Key update for large fungal genomes (e.g. Rhizophagus irregularis C2 / DAOM-197198):
//  - If the reference/contig is too long for full alignment, we align a mapped WINDOW instead
//    (computed from the seed chain + padding), then lift SV coordinates back to global ref coords.
//  - This preserves the same SV->graph behaviour as the small "test" example, but scales.
//
// NOTE: Still a prototype. For production pangenomes you'd want better mapping, split alignment,
//       robust SV normalization, and proper allele graph construction.

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

// ---------------- Utils ----------------
static inline uint8_t base_to_bits(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return 4;
    }
}

static inline std::string to_upper_acgt(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        char u = std::toupper(static_cast<unsigned char>(c));
        if (u == 'A' || u == 'C' || u == 'G' || u == 'T') out.push_back(u);
        else out.push_back('N');
    }
    return out;
}

static inline std::string to_lower_ascii(std::string s) {
    for (char& c : s) c = (char)std::tolower((unsigned char)c);
    return s;
}

static inline bool icontains(const std::string& hay, const std::string& needle) {
    return to_lower_ascii(hay).find(to_lower_ascii(needle)) != std::string::npos;
}

// ---------------- FASTA ----------------
struct FastaRecord {
    std::string name;
    std::string seq;
};

static std::vector<FastaRecord> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Failed to open FASTA: " + path);

    std::vector<FastaRecord> recs;
    std::string line;
    FastaRecord cur;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!cur.name.empty()) {
                cur.seq = to_upper_acgt(cur.seq);
                recs.push_back(std::move(cur));
                cur = FastaRecord{};
            }
            cur.name = line.substr(1);
            auto pos = cur.name.find_first_of(" \t");
            if (pos != std::string::npos) cur.name.resize(pos);
        } else {
            cur.seq += line;
        }
    }
    if (!cur.name.empty()) {
        cur.seq = to_upper_acgt(cur.seq);
        recs.push_back(std::move(cur));
    }
    return recs;
}

// ---------------- Hashing ----------------
static inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline uint64_t hash_kmer(const std::string& s, size_t pos, int k) {
    uint64_t x = 0;
    for (int i = 0; i < k; i++) {
        uint8_t b = base_to_bits(s[pos + (size_t)i]);
        if (b >= 4) return 0ULL;
        x = (x << 2) | b;
    }
    return splitmix64(x);
}

static inline std::pair<uint64_t, int> min_smer_hash_in_kmer(
    const std::string& s, size_t pos, int k, int smer)
{
    uint64_t best = std::numeric_limits<uint64_t>::max();
    int best_i = -1;

    for (int i = 0; i <= k - smer; i++) {
        bool ok = true;
        uint64_t x = 0;
        for (int j = 0; j < smer; j++) {
            uint8_t b = base_to_bits(s[pos + (size_t)i + (size_t)j]);
            if (b >= 4) { ok = false; break; }
            x = (x << 2) | b;
        }
        if (!ok) continue;
        uint64_t h = splitmix64(x);
        if (h < best) { best = h; best_i = i; }
    }
    return {best, best_i};
}

// ---------------- Syncmers ----------------
struct Syncmer { uint64_t h; uint32_t pos; };

static std::vector<Syncmer> compute_syncmers(
    const std::string& seq, int k, int smer, int t)
{
    std::vector<Syncmer> out;
    if ((int)seq.size() < k || smer > k) return out;

    for (size_t i = 0; i + (size_t)k <= seq.size(); i++) {
        auto [minh, argmin] = min_smer_hash_in_kmer(seq, i, k, smer);
        (void)minh;
        if (argmin == t) {
            uint64_t hk = hash_kmer(seq, i, k);
            if (hk != 0ULL) out.push_back({hk, (uint32_t)i});
        }
    }
    return out;
}

// ---------------- Reference index (bucketed, on concatenated ref) ----------------
struct RefIndex {
    int interval_w = 2000;
    // bucket -> hash -> list of global positions in concatenated ref
    std::unordered_map<uint32_t, std::unordered_map<uint64_t, std::vector<uint32_t>>> buckets;
};

static RefIndex build_ref_index(
    const std::string& refseq_concat, int k, int smer, int t, int interval_w)
{
    RefIndex idx;
    idx.interval_w = interval_w;
    auto syncs = compute_syncmers(refseq_concat, k, smer, t);
    idx.buckets.reserve(syncs.size() / 8 + 1);

    for (const auto& sm : syncs) {
        uint32_t b = sm.pos / (uint32_t)interval_w;
        idx.buckets[b][sm.h].push_back(sm.pos);
    }
    return idx;
}

// ---------------- Seeding ----------------
struct SeedPair { uint32_t qpos, rpos; }; // rpos is global position in concatenated ref

static std::vector<SeedPair> find_seed_pairs(
    const std::string& qry,
    const RefIndex& ridx,
    int k, int smer, int t,
    int interval_w,
    int max_buckets_scan = 2)
{
    std::vector<SeedPair> seeds;
    auto qsyncs = compute_syncmers(qry, k, smer, t);
    seeds.reserve(qsyncs.size() * 2);

    for (const auto& qs : qsyncs) {
        uint32_t qb = qs.pos / (uint32_t)interval_w;
        for (int delta = -max_buckets_scan; delta <= max_buckets_scan; delta++) {
            int64_t b = (int64_t)qb + delta;
            if (b < 0) continue;
            auto itb = ridx.buckets.find((uint32_t)b);
            if (itb == ridx.buckets.end()) continue;

            auto ith = itb->second.find(qs.h);
            if (ith == itb->second.end()) continue;

            for (uint32_t rp : ith->second) seeds.push_back({qs.pos, rp});
        }
    }

    std::sort(seeds.begin(), seeds.end(),
              [](const SeedPair& a, const SeedPair& b) {
                  return (a.qpos < b.qpos) || (a.qpos == b.qpos && a.rpos < b.rpos);
              });
    return seeds;
}

// ---------------- Chaining (simple LIS with diagonal constraint) ----------------
struct Chain { std::vector<SeedPair> points; int32_t score = 0; };

static Chain chain_seeds(const std::vector<SeedPair>& seeds, int max_diag_gap = 2000) {
    Chain best;
    if (seeds.empty()) return best;

    const int n = (int)seeds.size();
    std::vector<int> dp(n, 1), prev(n, -1);

    auto diag = [&](int i) -> int64_t { return (int64_t)seeds[i].rpos - (int64_t)seeds[i].qpos; };

    int best_i = 0;
    for (int i = 0; i < n; i++) {
        int best_here = 1;
        int best_prev = -1;
        int64_t di = diag(i);

        for (int j = 0; j < i; j++) {
            if (seeds[j].qpos < seeds[i].qpos && seeds[j].rpos < seeds[i].rpos) {
                int64_t dj = diag(j);
                if (std::llabs(di - dj) <= max_diag_gap) {
                    int cand = dp[j] + 1;
                    if (cand > best_here) { best_here = cand; best_prev = j; }
                }
            }
        }
        dp[i] = best_here;
        prev[i] = best_prev;
        if (dp[i] > dp[best_i]) best_i = i;
    }

    std::vector<SeedPair> chain;
    for (int i = best_i; i != -1; i = prev[i]) chain.push_back(seeds[i]);
    std::reverse(chain.begin(), chain.end());

    best.points = std::move(chain);
    best.score = dp[best_i];
    return best;
}

// ---------------- Auto band ----------------
static int compute_band_full_alignment(
    const Chain& chain,
    int ref_len,
    int qry_len,
    int min_band,
    int band_cap)
{
    int64_t len_diff = std::llabs((int64_t)ref_len - (int64_t)qry_len);
    int band = (int)len_diff + 256;

    if (chain.points.size() >= 2) {
        int64_t dmin = std::numeric_limits<int64_t>::max();
        int64_t dmax = std::numeric_limits<int64_t>::min();
        for (const auto& p : chain.points) {
            int64_t d = (int64_t)p.rpos - (int64_t)p.qpos;
            dmin = std::min(dmin, d);
            dmax = std::max(dmax, d);
        }
        band = (int)std::max<int64_t>(band, (dmax - dmin) + len_diff + 256);
    }

    band = std::max(band, min_band);
    band = std::min(band, band_cap);
    return band;
}

// ---------------- Affine-gap banded global alignment with correct backtrace ----------------
struct AlnOp { char op; int len; };
struct AlignmentResult { int score = 0; int edit = 0; std::vector<AlnOp> cigar; };

// Safety guard: banded DP with traceback is O(n * band) memory/time.
// For very long sequences/windows, we fall back to anchor-splitting (below).
static constexpr uint64_t MAX_BT_CELLS = 250000000ULL; // ~250M cells => ~250MB for 1 byte/cell BT

static AlignmentResult banded_global_affine_bt(
    const std::string& ref,
    const std::string& qry,
    int band,
    int mismatch_cost,
    int gap_open,
    int gap_ext)
{
    const int n = (int)ref.size();
    const int m = (int)qry.size();
    const int INF = 1000000000;

    if (std::abs(n - m) > band) return AlignmentResult{INF, INF, {}};

    const int W = 2 * band + 1;
    auto idx = [&](int i, int j) -> int { return (j - i) + band; };

    std::vector<int> Mp(W, INF), Ip(W, INF), Dp(W, INF);
    std::vector<int> Mc(W, INF), Ic(W, INF), Dc(W, INF);

    // predecessor state for each state at each cell, packed into 1 byte:
    // bits 0-1: pm, bits 2-3: pi, bits 4-5: pd. state: 0=M, 1=I, 2=D
    const uint64_t bt_cells = (uint64_t)(n + 1) * (uint64_t)W;
    if (bt_cells > MAX_BT_CELLS) return AlignmentResult{INF, INF, {}};
    std::vector<uint8_t> BT(bt_cells, 0);

    auto bt_at = [&](int i, int k) -> uint8_t& { return BT[(uint64_t)i * (uint64_t)W + (uint64_t)k]; };
    auto set_pm = [&](int i, int k, uint8_t which) {
        uint8_t& x = bt_at(i, k);
        x = (uint8_t)((x & 0b11111100u) | (which & 0b11u));
    };
    auto set_pi = [&](int i, int k, uint8_t which) {
        uint8_t& x = bt_at(i, k);
        x = (uint8_t)((x & 0b11110011u) | ((which & 0b11u) << 2));
    };
    auto set_pd = [&](int i, int k, uint8_t which) {
        uint8_t& x = bt_at(i, k);
        x = (uint8_t)((x & 0b11001111u) | ((which & 0b11u) << 4));
    };
    auto get_pm = [&](int i, int k) -> uint8_t { return (uint8_t)(bt_at(i, k) & 0b11u); };
    auto get_pi = [&](int i, int k) -> uint8_t { return (uint8_t)((bt_at(i, k) >> 2) & 0b11u); };
    auto get_pd = [&](int i, int k) -> uint8_t { return (uint8_t)((bt_at(i, k) >> 4) & 0b11u); };

    auto argmin3 = [&](int a, int b, int c, int& best, uint8_t& which) {
        best = a; which = 0;
        if (b < best) { best = b; which = 1; }
        if (c < best) { best = c; which = 2; }
    };

    // (0,0)
    Mp[idx(0, 0)] = 0;

    // row 0 => insertions
    for (int j = 1; j <= m; j++) {
        if (std::abs(j) > band) continue;
        int k = idx(0, j);
        Ip[k] = gap_open + gap_ext * (j - 1);
        set_pi(0, k, 1); // from I
    }

    for (int i = 1; i <= n; i++) {
        std::fill(Mc.begin(), Mc.end(), INF);
        std::fill(Ic.begin(), Ic.end(), INF);
        std::fill(Dc.begin(), Dc.end(), INF);

        int jlo = std::max(0, i - band);
        int jhi = std::min(m, i + band);

        for (int j = jlo; j <= jhi; j++) {
            int k = idx(i, j);

            // col 0 => deletions
            if (j == 0) {
                Dc[k] = gap_open + gap_ext * (i - 1);
                set_pd(i, k, 2); // from D
                continue;
            }

            // D(i,j) from (i-1,j)
            if (std::abs((i - 1) - j) <= band) {
                int ku = idx(i - 1, j);
                int fromM = Mp[ku] + gap_open;
                int fromI = Ip[ku] + gap_open;
                int fromD = Dp[ku] + gap_ext;
                int best; uint8_t which;
                argmin3(fromM, fromI, fromD, best, which);
                Dc[k] = best;
                set_pd(i, k, which);
            }

            // I(i,j) from (i,j-1)
            if (std::abs(i - (j - 1)) <= band) {
                int kl = idx(i, j - 1);
                int fromM = Mc[kl] + gap_open;
                int fromI = Ic[kl] + gap_ext;
                int fromD = Dc[kl] + gap_open;
                int best; uint8_t which;
                argmin3(fromM, fromI, fromD, best, which);
                Ic[k] = best;
                set_pi(i, k, which);
            }

            // M(i,j) from (i-1,j-1)
            if (std::abs((i - 1) - (j - 1)) <= band) {
                int kd = idx(i - 1, j - 1);
                int sub = (ref[i - 1] == qry[j - 1]) ? 0 : mismatch_cost;
                int fromM = Mp[kd] + sub;
                int fromI = Ip[kd] + sub;
                int fromD = Dp[kd] + sub;
                int best; uint8_t which;
                argmin3(fromM, fromI, fromD, best, which);
                Mc[k] = best;
                set_pm(i, k, which);
            }
        }

        std::swap(Mp, Mc);
        std::swap(Ip, Ic);
        std::swap(Dp, Dc);
    }

    int kend = idx(n, m);
    int best; uint8_t state;
    argmin3(Mp[kend], Ip[kend], Dp[kend], best, state);

    AlignmentResult res;
    res.score = best;
    res.edit = best;

    // traceback
    int i = n, j = m;
    std::vector<char> ops;
    ops.reserve((size_t)n + (size_t)m);

    while (i > 0 || j > 0) {
        int k = idx(i, j);
        if (state == 0) { // M
            ops.push_back('M');
            uint8_t prev = get_pm(i, k);
            i--; j--;
            state = prev;
        } else if (state == 1) { // I
            ops.push_back('I');
            uint8_t prev = get_pi(i, k);
            j--;
            state = prev;
        } else { // D
            ops.push_back('D');
            uint8_t prev = get_pd(i, k);
            i--;
            state = prev;
        }
    }

    std::reverse(ops.begin(), ops.end());
    for (char op : ops) {
        if (res.cigar.empty() || res.cigar.back().op != op) res.cigar.push_back({op, 1});
        else res.cigar.back().len++;
    }
    return res;
}

// ---------------- Anchor-splitting alignment (scales to chromosomes) ----------------
// Align between chained seed anchors, then stitch CIGAR. Each DP problem stays small.
// If there are very large anchor gaps, increase max_chunk and/or make seeding denser.
static void cigar_push(std::vector<AlnOp>& cigar, char op, int len) {
    if (len <= 0) return;
    if (!cigar.empty() && cigar.back().op == op) cigar.back().len += len;
    else cigar.push_back({op, len});
}

static void cigar_append(std::vector<AlnOp>& cigar, const std::vector<AlnOp>& add) {
    for (const auto& c : add) cigar_push(cigar, c.op, c.len);
}

static Chain shift_and_clip_chain(const Chain& chain, uint32_t r0, uint32_t r1, uint32_t q0, uint32_t q1) {
    Chain out;
    out.score = chain.score;
    out.points.reserve(chain.points.size());
    for (const auto& p : chain.points) {
        if (p.rpos < r0 || p.rpos >= r1) continue;
        if (p.qpos < q0 || p.qpos >= q1) continue;
        out.points.push_back(SeedPair{(uint32_t)(p.qpos - q0), (uint32_t)(p.rpos - r0)});
    }
    return out;
}

// --- Robust anchor-splitting alignment (scales to chromosomes) ---
// Fixes "Alignment failed even after fallback" by:
//  1) Equalizing large |rlen-qlen| with a big I/D before DP.
//  2) Recursively splitting a chunk if DP still fails.
//  3) Printing diagnostics for failing chunks.

static bool align_chunk_recursive(
    AlignmentResult& out,
    const std::string& rseg_in,
    const std::string& qseg_in,
    int min_band,
    int band_cap,
    int mismatch,
    int gap_open,
    int gap_ext,
    int depth)
{
    const int INF = 1000000000;

    // Hard stop to avoid infinite recursion
    if (depth > 20) return false;

    std::string rseg = rseg_in;
    std::string qseg = qseg_in;

    int rlen = (int)rseg.size();
    int qlen = (int)qseg.size();
    int diff = rlen - qlen; // >0 means ref longer => deletion

    // If net length diff is too large for banded DP, "equalize" first.
    // This is SV-scale and prevents immediate INF due to |n-m| > band.
    if (std::abs(diff) > (band_cap - 256)) {
        if (diff > 0) {
            cigar_push(out.cigar, 'D', diff);
            int cost = gap_open + gap_ext * std::max(0, diff - 1);
            out.score += cost; out.edit += cost;
            rseg = rseg.substr((size_t)diff);
        } else {
            int ins = -diff;
            cigar_push(out.cigar, 'I', ins);
            int cost = gap_open + gap_ext * std::max(0, ins - 1);
            out.score += cost; out.edit += cost;
            qseg = qseg.substr((size_t)ins);
        }
    }

    rlen = (int)rseg.size();
    qlen = (int)qseg.size();
    if (rlen == 0 && qlen == 0) return true;

    int band = std::max(min_band, std::abs(rlen - qlen) + 256);
    band = std::min(band, band_cap);

    auto part = banded_global_affine_bt(rseg, qseg, band, mismatch, gap_open, gap_ext);
    if (part.edit < INF) {
        cigar_append(out.cigar, part.cigar);
        out.score += part.score;
        out.edit  += part.edit;
        return true;
    }

    // DP failed: split the chunk and try again (this is the main robustness booster).
    // Split around midpoints, align left then right.
    if (rlen < 2000 && qlen < 2000) {
        std::cerr << "    [debug] DP failed on small chunk even after equalize: "
                  << "rlen=" << rlen << " qlen=" << qlen
                  << " band=" << band << " cap=" << band_cap << "\n";
        return false;
    }

    int rm = rlen / 2;
    int qm = qlen / 2;

    std::string rL = rseg.substr(0, (size_t)rm);
    std::string rR = rseg.substr((size_t)rm);
    std::string qL = qseg.substr(0, (size_t)qm);
    std::string qR = qseg.substr((size_t)qm);

    // Try aligning halves. If that fails, try an alternate split (1/3-2/3).
    if (align_chunk_recursive(out, rL, qL, min_band, band_cap, mismatch, gap_open, gap_ext, depth + 1) &&
        align_chunk_recursive(out, rR, qR, min_band, band_cap, mismatch, gap_open, gap_ext, depth + 1)) {
        return true;
    }

    int r1 = rlen / 3;
    int q1 = qlen / 3;
    rL = rseg.substr(0, (size_t)r1);
    rR = rseg.substr((size_t)r1);
    qL = qseg.substr(0, (size_t)q1);
    qR = qseg.substr((size_t)q1);

    if (align_chunk_recursive(out, rL, qL, min_band, band_cap, mismatch, gap_open, gap_ext, depth + 1) &&
        align_chunk_recursive(out, rR, qR, min_band, band_cap, mismatch, gap_open, gap_ext, depth + 1)) {
        return true;
    }

    std::cerr << "    [debug] DP failed chunk (after equalize + split attempts): "
              << "rlen=" << rlen << " qlen=" << qlen
              << " |diff|=" << std::abs(rlen - qlen)
              << " band_cap=" << band_cap << "\n";
    return false;
}

static AlignmentResult align_by_anchors(
    const std::string& ref_sub,
    const std::string& qry_sub,
    const Chain& chain_local,
    int k,
    int min_band,
    int band_cap,
    int mismatch,
    int gap_open,
    int gap_ext,
    int max_chunk)
{
    const int INF = 1000000000;
    AlignmentResult out;

    if (chain_local.points.size() < 2) return AlignmentResult{INF, INF, {}};

    uint32_t ri = 0, qi = 0;

    for (const auto& a : chain_local.points) {
        uint32_t ar = a.rpos;
        uint32_t aq = a.qpos;

        if (ar < ri || aq < qi) continue;

        uint32_t rlen = ar - ri;
        uint32_t qlen = aq - qi;

        if (rlen > 0 || qlen > 0) {
            if ((int)rlen > max_chunk || (int)qlen > max_chunk) return AlignmentResult{INF, INF, {}};

            std::string rseg = ref_sub.substr(ri, rlen);
            std::string qseg = qry_sub.substr(qi, qlen);

            if (!align_chunk_recursive(out, rseg, qseg,
                                       min_band, band_cap, mismatch, gap_open, gap_ext, 0)) {
                return AlignmentResult{INF, INF, {}};
            }
        }

        // Anchor itself: do NOT force "exact match" blindly.
        // Still record as M for stitching, but keep it short (k) and clipped.
        uint32_t mk = (uint32_t)k;
        if (ar + mk > (uint32_t)ref_sub.size()) mk = (uint32_t)ref_sub.size() - ar;
        if (aq + mk > (uint32_t)qry_sub.size()) mk = std::min<uint32_t>(mk, (uint32_t)qry_sub.size() - aq);
        cigar_push(out.cigar, 'M', (int)mk);

        ri = ar + mk;
        qi = aq + mk;
    }

    // Tail
    if (ri <= (uint32_t)ref_sub.size() && qi <= (uint32_t)qry_sub.size()) {
        uint32_t rlen = (uint32_t)ref_sub.size() - ri;
        uint32_t qlen = (uint32_t)qry_sub.size() - qi;

        if (rlen > 0 || qlen > 0) {
            if ((int)rlen > max_chunk || (int)qlen > max_chunk) return AlignmentResult{INF, INF, {}};

            std::string rseg = ref_sub.substr(ri, rlen);
            std::string qseg = qry_sub.substr(qi, qlen);

            if (!align_chunk_recursive(out, rseg, qseg,
                                       min_band, band_cap, mismatch, gap_open, gap_ext, 0)) {
                return AlignmentResult{INF, INF, {}};
            }
        }
    }

    return out;
}

// ---------------- SV extraction ----------------
struct SVEvent {
    // Local position (within the reference contig) of the breakpoint / start of event.
    uint32_t ref_pos = 0;

    // Types:
    //   INS / DEL  : extracted from CIGAR on the aligned window
    //   DUP / INV  : inferred from k-mer block mapping (allele stored in allele_seq; INV stored as reverse-complement)
    //   TRA        : inferred from k-mer block mapping (source contig -> target contig jump)
    std::string type;       // INS / DEL / DUP / INV / TRA

    // For INS/DUP/INV we store the alternative allele sequence (for INV it is stored as reverse-complement).
    // For DEL/TRA this is "*".
    std::string allele_seq = "*";

    // Length of allele or deleted span in bp (TRA uses 0).
    uint32_t len = 0;

    // Optional contig names (used for k-mer inferred events).
    std::string ref_contig1;
    std::string ref_contig2; // for TRA
    uint32_t ref_pos2 = 0;   // breakpoint on target contig for TRA
    char orient1 = '+';      // '+' or '-'
    char orient2 = '+';
};

static std::vector<SVEvent> extract_svs_from_cigar(
    const std::string& ref,
    const std::string& qry,
    const std::vector<AlnOp>& cigar,
    uint32_t sv_min_len)
{
    std::vector<SVEvent> svs;
    uint32_t r = 0, q = 0;

    for (const auto& c : cigar) {
        if (c.op == 'M') { r += (uint32_t)c.len; q += (uint32_t)c.len; }
        else if (c.op == 'D') {
            if ((uint32_t)c.len >= sv_min_len) {
                { SVEvent ev; ev.ref_pos = r; ev.type = "DEL"; ev.allele_seq = "*"; ev.len = (uint32_t)c.len; svs.push_back(ev); }
            }
            r += (uint32_t)c.len;
        } else if (c.op == 'I') {
            if ((uint32_t)c.len >= sv_min_len) {
                std::string ins = qry.substr(q, (size_t)c.len);
                { SVEvent ev; ev.ref_pos = r; ev.type = "INS"; ev.allele_seq = ins; ev.len = (uint32_t)c.len; svs.push_back(ev); }
            }
            q += (uint32_t)c.len;
        }
    }
    return svs;
}

// ---------------- Simple GFA graph ----------------
struct Graph {
    std::vector<std::string> seg_name;
    std::vector<std::string> seg_seq;

    struct Edge { int a; char ao; int b; char bo; };
    std::vector<Edge> edges;

    struct PathStep { int seg; char orient; };
    std::unordered_map<std::string, std::vector<PathStep>> paths;

    int add_segment(const std::string& name, const std::string& seq) {
        seg_name.push_back(name);
        seg_seq.push_back(seq);
        return (int)seg_name.size() - 1;
    }

    void add_edge(int a, char ao, int b, char bo) { edges.push_back({a, ao, b, bo}); }

    void add_path_step(const std::string& pname, int seg, char orient) {
        paths[pname].push_back({seg, orient});
    }

    void write_gfa(const std::string& out) const {
        std::ofstream o(out);
        if (!o) throw std::runtime_error("Failed to write: " + out);

        o << "H\tVN:Z:1.0\n";
        for (size_t i = 0; i < seg_name.size(); i++) {
            o << "S\t" << seg_name[i] << "\t" << seg_seq[i] << "\n";
        }
        for (const auto& e : edges) {
            o << "L\t" << seg_name[e.a] << "\t" << e.ao
              << "\t" << seg_name[e.b] << "\t" << e.bo
              << "\t0M\n";
        }
        for (const auto& kv : paths) {
            o << "P\t" << kv.first << "\t";
            const auto& steps = kv.second;
            for (size_t i = 0; i < steps.size(); i++) {
                if (i) o << ",";
                o << seg_name[steps[i].seg] << steps[i].orient;
            }
            o << "\t*\n";
        }
    }
};

static int find_segment_id(const Graph& g, const std::string& sname) {
    for (int i = 0; i < (int)g.seg_name.size(); i++) {
        if (g.seg_name[i] == sname) return i;
    }
    return -1;
}

static int ref_segment_index(uint32_t ref_pos_local, int segment_size) {
    return (int)(ref_pos_local / (uint32_t)segment_size);
}

static int find_ref_segment_id(const Graph& g, const std::string& ref_contig_name, int seg_idx) {
    return find_segment_id(g, ref_contig_name + "_seg" + std::to_string(seg_idx));
}

static void ensure_sample_path(Graph& g, const std::string& sample_name) {
    std::string pname = "sample:" + sample_name;
    if (g.paths.find(pname) == g.paths.end()) {
        // copy reference path
        auto it = g.paths.find("ref");
        if (it != g.paths.end()) g.paths[pname] = it->second;
        else g.paths[pname] = {};
    }
}

static void add_sv_to_graph(Graph& g,
                            const std::string& sample_name,
                            const std::string& ref_contig_name,
                            const std::string& asm_contig_name,
                            uint32_t ref_pos_local,
                            const SVEvent& sv_local,
                            int segment_size)
{
    const bool is_ins  = (sv_local.type == "INS");
    const bool is_del  = (sv_local.type == "DEL");
    const bool is_dup  = (sv_local.type == "DUP");
    const bool is_inv  = (sv_local.type == "INV");
    const bool is_tra  = (sv_local.type == "TRA");
    const bool is_ins_like = (is_ins || is_dup || is_inv);

    int left_seg = ref_segment_index(ref_pos_local, segment_size);

    // By default, connect to the next segment after the breakpoint.
    int right_seg = left_seg + 1;

    if (is_del) {
        // Jump to the segment just after the deleted interval.
        uint32_t del_end = ref_pos_local + sv_local.len;
        right_seg = ref_segment_index(del_end, segment_size) + 1;
        if (right_seg <= left_seg) right_seg = left_seg + 1;
    }

    std::ostringstream nm;
    nm << sample_name << "|" << asm_contig_name << "|" << sv_local.type
       << "|ref:" << ref_contig_name << ":" << ref_pos_local << "|len" << sv_local.len;

    if (is_tra) {
        nm << "|to:" << sv_local.ref_contig2 << ":" << sv_local.ref_pos2
           << "|o:" << sv_local.orient1 << ">" << sv_local.orient2;
    }

    std::string allele_seq = is_ins_like ? sv_local.allele_seq : "*";
    int alt_id = g.add_segment(nm.str(), allele_seq);

    int left_id  = find_ref_segment_id(g, ref_contig_name, left_seg);

    int right_id = -1;
    if (is_tra) {
        int tgt_seg = ref_segment_index(sv_local.ref_pos2, segment_size);
        right_id = find_ref_segment_id(g, sv_local.ref_contig2, tgt_seg);
    } else {
        right_id = find_ref_segment_id(g, ref_contig_name, right_seg);
    }

    // Add detour edges
    if (left_id  >= 0) g.add_edge(left_id,  '+', alt_id, '+');
    if (right_id >= 0) g.add_edge(alt_id, '+', right_id, '+');

    // Update sample path
    ensure_sample_path(g, sample_name);
    std::string pname = "sample:" + sample_name;
    auto& steps = g.paths[pname];

    // Find the first occurrence of left_id in the current steps
    int insert_pos = -1;
    for (int i = 0; i < (int)steps.size(); i++) {
        if (steps[i].seg == left_id) { insert_pos = i + 1; break; }
    }
    if (insert_pos < 0) {
        steps.push_back({alt_id, '+'});
        return;
    }

    steps.insert(steps.begin() + insert_pos, Graph::PathStep{alt_id, '+'});

    if (is_del && right_id >= 0) {
        // Remove deleted reference segments between left_seg and right_seg.
        std::unordered_set<int> to_remove;
        for (int si = left_seg + 1; si <= right_seg - 1; si++) {
            int sid = find_ref_segment_id(g, ref_contig_name, si);
            if (sid >= 0) to_remove.insert(sid);
        }

        int cur = insert_pos + 1;
        while (cur < (int)steps.size()) {
            int sid = steps[cur].seg;
            if (sid == right_id) break;
            if (to_remove.find(sid) != to_remove.end()) {
                steps.erase(steps.begin() + cur);
                continue;
            }
            break;
        }
    }

    if (is_tra && right_id >= 0) {
        // For translocations, splice the sample path to jump to the target segment.
        int j = -1;
        for (int ii = insert_pos + 1; ii < (int)steps.size(); ii++) {
            if (steps[ii].seg == right_id) { j = ii; break; }
        }
        if (j >= 0) {
            steps.erase(steps.begin() + (insert_pos + 1), steps.begin() + j);
        }
    }
}

// ---------------- Reference contig bookkeeping ----------------
struct RefContigInfo {
    std::string name;
    uint32_t start_global = 0; // start position in concatenated reference
    uint32_t len = 0;
};

static bool global_to_contig(const std::vector<RefContigInfo>& infos,
                             uint32_t pos_global,
                             size_t& contig_id,
                             uint32_t& pos_local)
{
    for (size_t i = 0; i < infos.size(); i++) {
        uint32_t st = infos[i].start_global;
        uint32_t en = st + infos[i].len;
        if (pos_global >= st && pos_global < en) {
            contig_id = i;
            pos_local = pos_global - st;
            return true;
        }
    }
    return false;
}


// ---------------- K-mer based complex SV discovery (INV/DUP/TRA) ----------------
static inline bool encode_kmer_2bit(const std::string& s, uint32_t pos, uint32_t k, uint64_t& out) {
    out = 0;
    if (pos + k > (uint32_t)s.size()) return false;
    for (uint32_t i = 0; i < k; i++) {
        uint8_t b = base_to_bits(s[pos + i]);
        if (b > 3) return false; // skip N/other
        out = (out << 2) | (uint64_t)b;
    }
    return true;
}

static inline uint64_t revcomp_2bit(uint64_t x, uint32_t k) {
    // Reverse complement for 2-bit encoding A=0,C=1,G=2,T=3 (comp = 3-b)
    uint64_t y = 0;
    for (uint32_t i = 0; i < k; i++) {
        uint64_t b = (x & 3ULL);
        b = 3ULL - b;
        y = (y << 2) | b;
        x >>= 2;
    }
    return y;
}

static inline std::string rev_comp_str(const std::string& s) {
    std::string out;
    out.resize(s.size());
    for (size_t i = 0; i < s.size(); i++) {
        char c = s[s.size() - 1 - i];
        switch (std::toupper(static_cast<unsigned char>(c))) {
            case 'A': out[i] = 'T'; break;
            case 'C': out[i] = 'G'; break;
            case 'G': out[i] = 'C'; break;
            case 'T': out[i] = 'A'; break;
            default:  out[i] = 'N'; break;
        }
    }
    return out;
}

static std::unordered_map<uint64_t, uint32_t> build_unique_kmer_index_2bit(const std::string& ref, uint32_t k) {
    std::unordered_map<uint64_t, uint32_t> pos;
    std::unordered_map<uint64_t, uint32_t> cnt;
    pos.reserve(ref.size() / 2);
    cnt.reserve(ref.size() / 2);

    for (uint32_t i = 0; i + k <= (uint32_t)ref.size(); i++) {
        uint64_t key = 0;
        if (!encode_kmer_2bit(ref, i, k, key)) continue;
        uint32_t c = ++cnt[key];
        if (c == 1) pos[key] = i;
        else pos.erase(key); // keep only unique
    }
    return pos;
}

struct KmerHit {
    uint32_t q = 0;
    uint32_t r = 0;   // global ref coordinate (forward)
    char orient = '+'; // '+' or '-'
};

struct MapBlock {
    uint32_t q0 = 0, q1 = 0;   // half-open query interval
    uint32_t r0 = 0, r1 = 0;   // half-open ref interval on forward strand (global)
    char orient = '+';
    int hits = 0;
};


static inline bool kmer_intervals_overlap(uint32_t a0, uint32_t a1, uint32_t b0, uint32_t b1, uint32_t tol = 0) {
    uint32_t x0 = std::max(a0, b0);
    uint32_t x1 = std::min(a1, b1);
    if (x1 <= x0) return false;
    return (x1 - x0) > tol;
}


// Build collinear blocks from hits. For '+' blocks, r should increase with q.
// For '-' blocks, r should decrease with q (since r is in forward coordinates).
static std::vector<MapBlock> build_blocks_from_hits(std::vector<KmerHit> hits,
                                                   uint32_t k,
                                                   uint32_t max_q_gap,
                                                   uint32_t max_diag_dev,
                                                   uint32_t min_span,
                                                   int min_hits)
{
    std::vector<MapBlock> out;
    if (hits.empty()) return out;
    std::sort(hits.begin(), hits.end(), [](const KmerHit& a, const KmerHit& b){ return a.q < b.q; });

    auto start_block = [&](const KmerHit& h) {
        MapBlock b;
        b.q0 = h.q;
        b.q1 = h.q + k;
        b.r0 = h.r;
        b.r1 = h.r + k;
        b.orient = h.orient;
        b.hits = 1;
        return b;
    };

    MapBlock cur = start_block(hits[0]);
    // Diagonal measure: '+' uses (r-q), '-' uses (r+q)
    int64_t cur_diag = (cur.orient == '+') ? ((int64_t)hits[0].r - (int64_t)hits[0].q)
                                          : ((int64_t)hits[0].r + (int64_t)hits[0].q);
    uint32_t prev_q = hits[0].q;
    uint32_t prev_r = hits[0].r;

    for (size_t i = 1; i < hits.size(); i++) {
        const auto& h = hits[i];
        uint32_t q = h.q;
        uint32_t r = h.r;

        uint32_t q_gap = (q > prev_q) ? (q - prev_q) : (prev_q - q);
        bool mono_ok = true;
        if (cur.orient == '+') mono_ok = (r >= prev_r);
        else mono_ok = (r <= prev_r);

        int64_t diag = (cur.orient == '+') ? ((int64_t)r - (int64_t)q) : ((int64_t)r + (int64_t)q);
        bool diag_ok = (std::llabs(diag - cur_diag) <= (int64_t)max_diag_dev);

        if (q_gap <= max_q_gap && mono_ok && diag_ok) {
            cur.q1 = std::max(cur.q1, q + k);
            cur.r0 = std::min(cur.r0, r);
            cur.r1 = std::max(cur.r1, r + k);
            cur.hits++;
            // slowly adapt diag
            cur_diag = (cur_diag * 7 + diag) / 8;
        } else {
            uint32_t span = (cur.q1 > cur.q0) ? (cur.q1 - cur.q0) : 0;
            if (span >= min_span && cur.hits >= min_hits) out.push_back(cur);
            cur = start_block(h);
            cur_diag = diag;
        }

        prev_q = q;
        prev_r = r;
    }

    uint32_t span = (cur.q1 > cur.q0) ? (cur.q1 - cur.q0) : 0;
    if (span >= min_span && cur.hits >= min_hits) out.push_back(cur);
    return out;
}

// Pick up to N non-overlapping blocks (in query coordinates), preferring more hits.
static std::vector<MapBlock> pick_blocks_nonoverlap(std::vector<MapBlock> blocks, int N, uint32_t min_q_sep) {
    std::vector<MapBlock> picked;
    if (blocks.empty()) return picked;

    std::sort(blocks.begin(), blocks.end(), [](const MapBlock& a, const MapBlock& b){
        if (a.hits != b.hits) return a.hits > b.hits;
        return (a.q1 - a.q0) > (b.q1 - b.q0);
    });

    for (const auto& b : blocks) {
        bool ok = true;
        for (const auto& p : picked) {
            if (kmer_intervals_overlap(p.q0, p.q1, b.q0, b.q1, min_q_sep)) { ok = false; break; }
        }
        if (ok) {
            picked.push_back(b);
            if ((int)picked.size() >= N) break;
        }
    }

    std::sort(picked.begin(), picked.end(), [](const MapBlock& a, const MapBlock& b){ return a.q0 < b.q0; });
    return picked;
}

static std::vector<SVEvent> infer_dup_inv_tra_from_blocks(const std::string& qry,
                                                         const std::vector<RefContigInfo>& ref_infos,
                                                         const std::vector<MapBlock>& blocks,
                                                         uint32_t k,
                                                         uint32_t tra_min_dist,
                                                         uint32_t min_dup_ovl)
{
    std::vector<SVEvent> out;
    if (blocks.empty()) return out;

    // Map global ref coords to contigs/locals for each block (skip blocks landing in separator Ns).
    struct BlkX { MapBlock b; std::string contig; uint32_t r0l=0, r1l=0; };
    std::vector<BlkX> bx;
    bx.reserve(blocks.size());
    for (const auto& b : blocks) {
        size_t cid0=0, cid1=0;
        uint32_t r0l=0, r1l=0;
        if (!global_to_contig(ref_infos, b.r0, cid0, r0l)) continue;
        // b.r1 is half-open; map r1-1
        if (b.r1 == 0) continue;
        if (!global_to_contig(ref_infos, b.r1 - 1, cid1, r1l)) continue;
        if (cid0 != cid1) {
            // block straddles contigs (rare with unique kmers); skip
            continue;
        }
        BlkX x; x.b = b; x.contig = ref_infos[cid0].name; x.r0l = r0l; x.r1l = r1l + 1;
        bx.push_back(std::move(x));
    }
    if (bx.empty()) return out;

    // INV: merge '-' blocks into maximal inversion intervals (prevents multi-counting one inversion).
{
    struct InvBlk { std::string contig; uint32_t q0,q1,r0,r1; int hits; };
    std::vector<InvBlk> invs;
    invs.reserve(bx.size());
    for (const auto& x : bx) {
        if (x.b.orient != '-') continue;
        invs.push_back({x.contig, x.b.q0, x.b.q1, x.b.r0, x.b.r1, x.b.hits});
    }
    if (!invs.empty()) {
        std::sort(invs.begin(), invs.end(), [](const InvBlk& a, const InvBlk& b){
            if (a.contig != b.contig) return a.contig < b.contig;
            if (a.r0 != b.r0) return a.r0 < b.r0;
            return a.q0 < b.q0;
        });

        const uint32_t MERGE_GAP = 200;   // bp tolerance to merge adjacent blocks
        // Minimum inversion length. This removes short spurious "-" blocks that typically
        // come from local repeats/noise and inflate the SV count.
        const uint32_t MIN_SPAN = 800;

        uint32_t mq0=invs[0].q0, mq1=invs[0].q1;
        uint32_t mr0=invs[0].r0, mr1=invs[0].r1;
        int mhits=invs[0].hits;
        std::string mcontig=invs[0].contig;

        auto flush = [&](const std::string& contig, uint32_t q0, uint32_t q1, uint32_t r0, uint32_t r1, int hits){
            if (q1 <= q0) return;
            uint32_t span = q1 - q0;
            if (span < MIN_SPAN) return;
            SVEvent ev;
            ev.type = "INV";
            ev.ref_contig1 = contig;
            // Map global r0 to local
            size_t cid=0;
            uint32_t r0l=0;
            if (!global_to_contig(ref_infos, r0, cid, r0l)) return;
            ev.ref_pos = r0l;
            ev.len = span;
            std::string seg = qry.substr(q0, std::min<size_t>((size_t)span, qry.size() - q0));
            ev.allele_seq = rev_comp_str(seg);
            out.push_back(std::move(ev));
        };

        for (size_t i=1;i<invs.size();i++){
            const auto& b = invs[i];
            bool same = (b.contig == mcontig);
            bool ref_close = (b.r0 <= mr1 + MERGE_GAP);
            bool qry_close = (b.q0 <= mq1 + MERGE_GAP);
            if (same && ref_close && qry_close) {
                mq0 = std::min(mq0, b.q0); mq1 = std::max(mq1, b.q1);
                mr0 = std::min(mr0, b.r0); mr1 = std::max(mr1, b.r1);
                mhits += b.hits;
            } else {
                flush(mcontig, mq0, mq1, mr0, mr1, mhits);
                mcontig=b.contig; mq0=b.q0; mq1=b.q1; mr0=b.r0; mr1=b.r1; mhits=b.hits;
            }
        }
        flush(mcontig, mq0, mq1, mr0, mr1, mhits);
    }
}
    // DUP/TRA: analyze adjacency between '+' blocks by query order
    // Track visited ref intervals (on same contig) to infer duplication
    struct Interval { std::string contig; uint32_t a0, a1; };
    std::vector<Interval> visited;

    auto overlap_len = [](uint32_t a0, uint32_t a1, uint32_t b0, uint32_t b1)->uint32_t{
        uint32_t x0 = std::max(a0,b0), x1 = std::min(a1,b1);
        return (x1 > x0) ? (x1 - x0) : 0;
    };

    // collect '+' blocks only
    std::vector<BlkX> plus;
    for (auto& x : bx) if (x.b.orient == '+') plus.push_back(x);
    if (plus.size() < 2) return out;

    for (size_t i = 0; i < plus.size(); i++) {
        // mark visited
        visited.push_back({plus[i].contig, plus[i].r0l, plus[i].r1l});

        if (i + 1 >= plus.size()) break;
        const auto& A = plus[i];
        const auto& B = plus[i+1];

        // Breakpoint at end of A on source contig
        uint32_t bp = A.r1l;

        // Determine if DUP-like: B overlaps any earlier visited interval.
        // Important: a transposed duplication usually involves a big reference jump; we
        // must prefer DUP over TRA when overlap evidence is present.
        bool is_dup = false;
        for (size_t v = 0; v < visited.size(); v++) {
            if (visited[v].contig != B.contig) continue;
            uint32_t ov = overlap_len(visited[v].a0, visited[v].a1, B.r0l, B.r1l);
            if (ov >= min_dup_ovl) { is_dup = true; break; }
        }

        // Determine if TRA-like jump (only if not a duplication)
        bool is_tra = false;
        if (!is_dup) {
            if (A.contig != B.contig) {
                is_tra = true;
            } else {
                uint32_t jump = (B.r0l > A.r1l) ? (B.r0l - A.r1l) : (A.r1l - B.r0l);
                if (jump >= tra_min_dist) is_tra = true;
            }
        }

        if (is_dup) {
            SVEvent ev;
            ev.type = "DUP";
            ev.ref_contig1 = A.contig;
            ev.ref_pos = bp;
            // Use the *reference* span of the duplicated block as DUP length. Using query span can
            // substantially overestimate DUP size when the block covers multiple nearby events.
            ev.len = (B.r1l > B.r0l) ? (B.r1l - B.r0l) : 0;
            size_t take = std::min<size_t>(ev.len, (B.b.q1 > B.b.q0) ? (B.b.q1 - B.b.q0) : 0);
            if (take > 0 && B.b.q0 < qry.size()) {
                ev.allele_seq = qry.substr(B.b.q0, std::min(take, qry.size() - (size_t)B.b.q0));
            }
            out.push_back(std::move(ev));
        } else if (is_tra) {
            SVEvent ev;
            ev.type = "TRA";
            ev.ref_contig1 = A.contig;
            ev.ref_pos = bp;
            ev.ref_contig2 = B.contig;
            ev.ref_pos2 = B.r0l;
            ev.len = 0;
            ev.allele_seq = "*";
            ev.orient1 = '+';
            ev.orient2 = '+';
            out.push_back(std::move(ev));
        }
    }


return out;
}





// Deduplicate SV events that are nearly identical (prevents over-counting from fragmented evidence).
static std::vector<SVEvent> dedup_svs(std::vector<SVEvent> evs,
                                     uint32_t pos_bin = 25,
                                     uint32_t len_bin = 25)
{
    auto bin = [](uint32_t x, uint32_t b)->uint32_t { return b ? (x / b) : x; };

    std::vector<SVEvent> out;
    out.reserve(evs.size());
    std::unordered_set<std::string> seen;
    seen.reserve(evs.size() * 2);

    for (auto& e : evs) {
        std::ostringstream key;
        // For intra-contig TRA, treat the two breakpoints as an undirected pair.
        // This merges A->B and B->A representations of the same transposition.
        if (e.type == "TRA" && e.ref_contig1 == e.ref_contig2) {
            uint32_t a = bin(e.ref_pos, pos_bin);
            uint32_t b = bin(e.ref_pos2, pos_bin);
            if (b < a) std::swap(a, b);
            key << "TRA|" << e.ref_contig1 << ":" << a << "-" << b;
        } else {
            key << e.type << "|"
                << e.ref_contig1 << ":" << bin(e.ref_pos, pos_bin) << "|"
                << "l" << bin(e.len, len_bin);
            if (e.type == "TRA") {
                key << "|to:" << e.ref_contig2 << ":" << bin(e.ref_pos2, pos_bin);
            }
        }
        std::string k = key.str();
        if (seen.insert(k).second) out.push_back(std::move(e));
    }
    return out;
}

// ---------------- Build reference backbone graph for ALL contigs ----------------
static Graph build_ref_backbone_graph_all(const std::vector<FastaRecord>& ref_recs,
                                         int segment_size,
                                         std::vector<RefContigInfo>& out_infos,
                                         std::string& out_concat,
                                         int concat_sep_len = 100)
{
    Graph g;
    out_infos.clear();
    out_concat.clear();

    uint32_t cursor = 0;
    bool first_contig = true;

    for (const auto& rr : ref_recs) {
        if (rr.seq.empty()) continue;

        RefContigInfo info;
        info.name = rr.name.empty() ? ("ref_contig" + std::to_string(out_infos.size())) : rr.name;
        info.start_global = cursor;
        info.len = (uint32_t)rr.seq.size();
        out_infos.push_back(info);

        // add segments for this contig
        int nseg = (int)((rr.seq.size() + (size_t)segment_size - 1) / (size_t)segment_size);
        std::vector<int> ids(nseg);

        for (int i = 0; i < nseg; i++) {
            size_t st = (size_t)i * (size_t)segment_size;
            size_t en = std::min(rr.seq.size(), st + (size_t)segment_size);

            std::ostringstream nm;
            nm << info.name << "_seg" << i;
            ids[i] = g.add_segment(nm.str(), rr.seq.substr(st, en - st));
            g.add_path_step("ref", ids[i], '+');
            if (i > 0) g.add_edge(ids[i - 1], '+', ids[i], '+');
        }

        // concatenate ref sequence for indexing/alignment
        if (!first_contig) {
            out_concat.append((size_t)concat_sep_len, 'N');
            cursor += (uint32_t)concat_sep_len;
        }
        first_contig = false;

        out_concat += rr.seq;
        cursor += info.len;
    }

    return g;
}

// ---------------- File discovery ----------------
static bool is_fasta_ext(std::string ext) {
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return (ext == ".fa" || ext == ".fasta" || ext == ".fna" || ext == ".fas");
}

static std::vector<std::string> list_fasta_files_in_dirs(
    const std::vector<std::string>& dirs,
    bool recursive)
{
    std::vector<std::string> files;
    std::unordered_set<std::string> seen;

    auto add_file = [&](const fs::path& p) {
        std::error_code ec;
        fs::path canon = fs::weakly_canonical(p, ec);
        std::string key = ec ? p.lexically_normal().string() : canon.string();
        if (seen.insert(key).second) files.push_back(key);
    };

    for (const auto& d : dirs) {
        fs::path root(d);
        if (!fs::exists(root)) throw std::runtime_error("Assembly dir not found: " + d);
        if (!fs::is_directory(root)) throw std::runtime_error("Not a directory: " + d);

        if (recursive) {
            for (const auto& entry : fs::recursive_directory_iterator(root)) {
                if (!entry.is_regular_file()) continue;
                auto p = entry.path();
                if (is_fasta_ext(p.extension().string())) add_file(p);
            }
        } else {
            for (const auto& entry : fs::directory_iterator(root)) {
                if (!entry.is_regular_file()) continue;
                auto p = entry.path();
                if (is_fasta_ext(p.extension().string())) add_file(p);
            }
        }
    }

    std::sort(files.begin(), files.end());
    return files;
}

// ---------------- Sample naming (Rhizophagus C2 / DAOM-197198) ----------------
static std::string normalize_sample_name_from_path(const std::string& path) {
    std::string base = fs::path(path).filename().string();
    if (auto dot = base.find_last_of('.'); dot != std::string::npos) base = base.substr(0, dot);

    // Rhizophagus irregularis strain heuristics
    if (icontains(base, "c2") && (icontains(base, "rhizophagus") || icontains(base, "irregularis") || icontains(base, "rhi"))) {
        return "rhizophagus_irregularis_C2";
    }
    if ((icontains(base, "daom") || icontains(base, "197198")) &&
        (icontains(base, "rhizophagus") || icontains(base, "irregularis") || icontains(base, "rhi"))) {
        return "rhizophagus_irregularis_DAOM-197198";
    }
    // also accept directory naming
    std::string parent = fs::path(path).parent_path().filename().string();
    if (icontains(parent, "c2") && (icontains(parent, "rhizophagus") || icontains(parent, "irregularis") || icontains(parent, "rhi"))) {
        return "rhizophagus_irregularis_C2";
    }
    if ((icontains(parent, "daom") || icontains(parent, "197198")) &&
        (icontains(parent, "rhizophagus") || icontains(parent, "irregularis") || icontains(parent, "rhi"))) {
        return "rhizophagus_irregularis_DAOM-197198";
    }

    return base;
}

// ---------------- Mapping window selection ----------------
struct Window { uint32_t r0=0,r1=0,q0=0,q1=0; };

static bool compute_mapped_window(
    const Chain& chain,
    int k,
    uint32_t ref_len,
    uint32_t qry_len,
    uint32_t pad,
    uint32_t max_ref_window,
    Window& out)
{
    if (chain.points.size() < 2) return false;

    uint32_t qmin = std::numeric_limits<uint32_t>::max();
    uint32_t qmax = 0;
    uint32_t rmin = std::numeric_limits<uint32_t>::max();
    uint32_t rmax = 0;

    for (const auto& p : chain.points) {
        qmin = std::min(qmin, p.qpos);
        qmax = std::max(qmax, (uint32_t)(p.qpos + (uint32_t)k));
        rmin = std::min(rmin, p.rpos);
        rmax = std::max(rmax, (uint32_t)(p.rpos + (uint32_t)k));
    }

    if (qmin >= qry_len || rmin >= ref_len) return false;

    uint32_t q0 = (qmin > pad) ? (qmin - pad) : 0;
    uint32_t q1 = std::min<uint32_t>(qry_len, qmax + pad);

    uint32_t r0 = (rmin > pad) ? (rmin - pad) : 0;
    uint32_t r1 = std::min<uint32_t>(ref_len, rmax + pad);

    // Cap ref window size (keep centered-ish around chain)
    if (max_ref_window > 0 && r1 > r0 && (r1 - r0) > max_ref_window) {
        uint64_t mid = (uint64_t)r0 + (uint64_t)(r1 - r0) / 2;
        uint64_t half = (uint64_t)max_ref_window / 2;
        uint64_t nr0 = (mid > half) ? (mid - half) : 0;
        uint64_t nr1 = std::min<uint64_t>((uint64_t)ref_len, nr0 + (uint64_t)max_ref_window);
        // try to keep length exactly max_ref_window
        if (nr1 - nr0 < max_ref_window && nr0 > 0) {
            uint64_t need = (uint64_t)max_ref_window - (nr1 - nr0);
            if (need > nr0) need = nr0;
            nr0 -= need;
        }
        r0 = (uint32_t)nr0;
        r1 = (uint32_t)nr1;
    }

    // Ensure non-empty
    if (q1 <= q0 || r1 <= r0) return false;

    out = Window{r0,r1,q0,q1};
    return true;
}

// ---------------- CLI ----------------
struct Args {
    std::string ref_path;
    std::vector<std::string> asm_dirs;
    std::string out_gfa = "pangenome.gfa";

    int k = 15;
    int s = 9;
    int t = 5;
    int interval_w = 2000;

    int band = 256;        // min band
    int band_cap = 8192;   // cap for band

    int mismatch = 3;
    int gap_open = 5;
    int gap_ext  = 1;

    int sv_min = 50;
    int ref_segment = 1000;

    bool recursive = true;

    int top_contigs = 50;
    int min_contig_len = 1;

    bool debug = true;

    // For large genomes: align windows if sequences exceed this.
    int max_full_align = 200000;

    // Window extraction for big genomes
    int window_pad = 20000;       // padding around chained seeds
    int max_ref_window = 6000000; // cap mapped ref window length (0=unlimited)
};

static void usage() {
    std::cerr
        << "fungi_pangenome --ref ref.fa --asm-dir DIR [--asm-dir DIR2 ...] --out out.gfa [options]\n"
        << "Options:\n"
        << "  --k INT              syncmer k (default 15)\n"
        << "  --s INT              syncmer s (default 9)\n"
        << "  --t INT              open-syncmer offset t (default 5)\n"
        << "  --w INT              bucket size (default 2000)\n"
        << "  --band INT           minimum band (default 256)\n"
        << "  --band-cap INT       maximum band cap (default 8192)\n"
        << "  --mismatch INT       mismatch cost (default 3)\n"
        << "  --gap-open INT       gap open cost (default 5)\n"
        << "  --gap-ext INT        gap extend cost (default 1)\n"
        << "  --sv INT             SV threshold (default 50)\n"
        << "  --seg INT            ref segment size (default 1000)\n"
        << "  --max-full-align INT full align only if ref+contig <= this (default 200000)\n"
        << "  --window-pad INT     pad around mapped region (default 20000)\n"
        << "  --max-ref-window INT max ref window length (default 6000000; 0=unlimited)\n"
        << "  --no-recursive       do not scan subdirectories\n"
        << "  --top-contigs INT    process top N contigs per assembly file (default 50)\n"
        << "  --min-contig INT     skip contigs shorter than this (default 1)\n";
}

static Args parse_args(int argc, char** argv) {
    Args a;
    for (int i = 1; i < argc; i++) {
        std::string x = argv[i];
        auto need = [&](const std::string& opt) {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + opt);
            return std::string(argv[++i]);
        };

        if (x == "--ref") a.ref_path = need(x);
        else if (x == "--asm-dir") a.asm_dirs.push_back(need(x));
        else if (x == "--out") a.out_gfa = need(x);

        else if (x == "--k") a.k = std::stoi(need(x));
        else if (x == "--s") a.s = std::stoi(need(x));
        else if (x == "--t") a.t = std::stoi(need(x));
        else if (x == "--w") a.interval_w = std::stoi(need(x));

        else if (x == "--band") a.band = std::stoi(need(x));
        else if (x == "--band-cap") a.band_cap = std::stoi(need(x));

        else if (x == "--mismatch") a.mismatch = std::stoi(need(x));
        else if (x == "--gap-open") a.gap_open = std::stoi(need(x));
        else if (x == "--gap-ext") a.gap_ext = std::stoi(need(x));

        else if (x == "--sv") a.sv_min = std::stoi(need(x));
        else if (x == "--seg") a.ref_segment = std::stoi(need(x));
        else if (x == "--max-full-align") a.max_full_align = std::stoi(need(x));
        else if (x == "--window-pad") a.window_pad = std::stoi(need(x));
        else if (x == "--max-ref-window") a.max_ref_window = std::stoi(need(x));

        else if (x == "--top-contigs") a.top_contigs = std::stoi(need(x));
        else if (x == "--min-contig") a.min_contig_len = std::stoi(need(x));

        else if (x == "--no-recursive") a.recursive = false;
        else if (x == "--help" || x == "-h") { usage(); std::exit(0); }
        else throw std::runtime_error("Unknown option: " + x);
    }

    if (a.ref_path.empty() || a.asm_dirs.empty()) {
        usage();
        throw std::runtime_error("Need --ref and at least one --asm-dir");
    }
    if (a.s > a.k) throw std::runtime_error("--s must be <= --k");
    if (a.t < 0 || a.t > a.k - a.s) throw std::runtime_error("--t must be in [0, k-s]");
    a.band = std::max(16, a.band);
    a.band_cap = std::max(a.band, a.band_cap);
    a.ref_segment = std::max(50, a.ref_segment);
    a.top_contigs = std::max(1, a.top_contigs);
    a.min_contig_len = std::max(0, a.min_contig_len);
    a.max_full_align = std::max(1000, a.max_full_align);
    a.window_pad = std::max(0, a.window_pad);
    a.max_ref_window = std::max(0, a.max_ref_window);
    return a;
}

// ---------------- Main ----------------
int main(int argc, char** argv) {
    try {
        Args args = parse_args(argc, argv);

        auto ref_recs = read_fasta(args.ref_path);
        if (ref_recs.empty()) throw std::runtime_error("Empty reference FASTA");

        std::cerr << "[info] Reference contigs: " << ref_recs.size() << "\n";

        std::vector<RefContigInfo> ref_infos;
        std::string ref_concat;
        Graph g = build_ref_backbone_graph_all(ref_recs, args.ref_segment, ref_infos, ref_concat, 100);

        std::cerr << "[info] Concatenated reference length: " << ref_concat.size() << "\n";

        // Build unique 31-mer indices for complex SV discovery (INV/DUP/TRA).
        const uint32_t k_complex = 31;
        auto ref_kmer_pos_fwd = build_unique_kmer_index_2bit(ref_concat, k_complex);
        std::string ref_rc = rev_comp_str(ref_concat);
        auto ref_kmer_pos_rc  = build_unique_kmer_index_2bit(ref_rc, k_complex);
        std::cerr << "[info] Unique " << k_complex << "-mers (fwd)=" << ref_kmer_pos_fwd.size()
                  << " (rc)=" << ref_kmer_pos_rc.size() << "\n";


        std::cerr << "[info] Building reference syncmer index...\n";
        RefIndex ridx = build_ref_index(ref_concat, args.k, args.s, args.t, args.interval_w);
        std::cerr << "[info] Buckets: " << ridx.buckets.size() << "\n";

        auto asm_files = list_fasta_files_in_dirs(args.asm_dirs, args.recursive);
        if (asm_files.empty()) throw std::runtime_error("No FASTA assemblies found in --asm-dir");
        std::cerr << "[info] Found " << asm_files.size() << " assembly FASTA files\n";

        for (const auto& ap : asm_files) {
            auto arecs = read_fasta(ap);
            if (arecs.empty()) { std::cerr << "[warn] Empty assembly FASTA: " << ap << "\n"; continue; }

            std::string sample = normalize_sample_name_from_path(ap);

            std::sort(arecs.begin(), arecs.end(),
                      [](const FastaRecord& a, const FastaRecord& b) { return a.seq.size() > b.seq.size(); });

            int N = std::min<int>(args.top_contigs, (int)arecs.size());
            std::cerr << "[info] Sample " << sample << ": contigs=" << arecs.size()
                      << " processing top " << N << "\n";

            int total_svs_added = 0;

            for (int ci = 0; ci < N; ci++) {
                const auto& contig = arecs[ci];
                if ((int)contig.seq.size() < args.min_contig_len) continue;

                const std::string& qry_seq_full = contig.seq;
                std::cerr << "  [info] Contig " << ci << " " << contig.name
                          << " len=" << qry_seq_full.size() << "\n";

                // --- Seeding/chaining to estimate mapped region ---
                auto seeds = find_seed_pairs(qry_seq_full, ridx, args.k, args.s, args.t, args.interval_w);
                auto chain = chain_seeds(seeds);

                std::cerr << "    [info] Seeds=" << seeds.size()
                          << " chain_points=" << chain.points.size()
                          << " chain_score=" << chain.score << "\n";

                // Choose full vs windowed alignment
                std::string ref_sub;
                std::string qry_sub;
                uint32_t ref_sub_global0 = 0;
                uint32_t qry_sub0 = 0;

                Chain chain_local;

                // Never try to align the entire concatenated genome backbone for huge refs/contigs.
                // If too large, align a mapped WINDOW and lift coordinates.
                bool use_full = ((int)qry_seq_full.size() <= args.max_full_align) &&
                                ((args.max_ref_window == 0) || ((int)ref_concat.size() <= args.max_ref_window)) &&
                                ((int)ref_concat.size() <= args.max_full_align);

                if (use_full) {
                    ref_sub = ref_concat;
                    qry_sub = qry_seq_full;
                    ref_sub_global0 = 0;
                    qry_sub0 = 0;
                    chain_local = chain; // already in same coordinate system
                } else {
                    Window win;
                    bool ok = compute_mapped_window(chain, args.k,
                                                    (uint32_t)ref_concat.size(), (uint32_t)qry_seq_full.size(),
                                                    (uint32_t)args.window_pad, (uint32_t)args.max_ref_window, win);
                    if (!ok) {
                        std::cerr << "    [warn] No reliable chain/window (or too few seeds). Skipping contig.\n";
                        continue;
                    }

                    ref_sub = ref_concat.substr(win.r0, (size_t)(win.r1 - win.r0));
                    qry_sub = qry_seq_full.substr(win.q0, (size_t)(win.q1 - win.q0));
                    ref_sub_global0 = win.r0;
                    qry_sub0 = win.q0;

                    chain_local = shift_and_clip_chain(chain, win.r0, win.r1, win.q0, win.q1);

                    std::cerr << "    [info] Windowed alignment: ref[" << win.r0 << "," << win.r1 << ")"
                              << " qry[" << win.q0 << "," << win.q1 << ")\n";
                }

                // Band on the region we actually align (windowed or full).
                int band = compute_band_full_alignment(chain_local, (int)ref_sub.size(), (int)qry_sub.size(),
                                                       args.band, args.band_cap);

                // Try full banded affine DP first (on the window). If DP is too big, the aligner returns INF.
                auto aln = banded_global_affine_bt(ref_sub, qry_sub, band,
                                                   args.mismatch, args.gap_open, args.gap_ext);

                // Retry with a larger band (cheap) before fallback.
                if (aln.edit >= 1000000000) {
                    int band2 = std::min(args.band_cap, band * 2);
                    if (band2 > band) {
                        std::cerr << "    [warn] Alignment failed at band=" << band
                                  << "; retry band=" << band2 << "\n";
                        aln = banded_global_affine_bt(ref_sub, qry_sub, band2,
                                                      args.mismatch, args.gap_open, args.gap_ext);
                        band = band2;
                    }
                }

                // Debug: show max anchor gaps in the aligned coordinate system (helps choose max_chunk/seed density).
                if (chain_local.points.size() >= 2) {
                    uint32_t max_r_gap = 0, max_q_gap = 0;
                    for (size_t ii = 1; ii < chain_local.points.size(); ii++) {
                        max_r_gap = std::max(max_r_gap, chain_local.points[ii].rpos - chain_local.points[ii - 1].rpos);
                        max_q_gap = std::max(max_q_gap, chain_local.points[ii].qpos - chain_local.points[ii - 1].qpos);
                    }
                    std::cerr << "    [debug] max anchor gap: ref=" << max_r_gap << " qry=" << max_q_gap << "\n";
                }

                // If DP is still impossible (too big or band miss), fall back to anchor-splitting.
                if (aln.edit >= 1000000000) {
                aln = align_by_anchors(ref_sub, qry_sub, chain_local,
                       args.k,
                       args.band, args.band_cap,
                       args.mismatch, args.gap_open, args.gap_ext,
                       /*max_chunk=*/2000000);
		}

                if (aln.edit >= 1000000000) {
                    std::cerr << "    [warn] Alignment failed even after fallback.\n";
                    continue;
                }

                if (args.debug) {
                    int maxI = 0, maxD = 0, bigI = 0, bigD = 0;
                    for (auto& c : aln.cigar) {
                        if (c.op == 'I') { maxI = std::max(maxI, c.len); if (c.len >= args.sv_min) bigI++; }
                        if (c.op == 'D') { maxD = std::max(maxD, c.len); if (c.len >= args.sv_min) bigD++; }
                    }
                    std::cerr << "    [debug] band=" << band << " score=" << aln.score
                              << " maxI=" << maxI << " maxD=" << maxD
                              << " bigI=" << bigI << " bigD=" << bigD << "\n";
                }

                auto svs = extract_svs_from_cigar(ref_sub, qry_sub, aln.cigar, (uint32_t)args.sv_min);
                std::cerr << "    [info] SVs >= " << args.sv_min << "bp: " << svs.size() << "\n";

                // Collect INDEL SVs (do NOT add to graph yet). We'll filter out INDELs that are
                // likely artifacts of larger rearrangements (INV/DUP/TRA) to avoid double-counting.
                std::vector<SVEvent> indel_evs;
                indel_evs.reserve(svs.size());
                for (const auto& sv : svs) {
                    uint32_t sv_ref_global = ref_sub_global0 + sv.ref_pos;

                    // map SV start to a reference contig
                    size_t ref_cid = 0;
                    uint32_t ref_local = 0;
                    if (!global_to_contig(ref_infos, sv_ref_global, ref_cid, ref_local)) {
                        // likely within separator Ns; skip
                        continue;
                    }
                    const std::string& ref_contig_name = ref_infos[ref_cid].name;

                    SVEvent sv_adj = sv;
                    sv_adj.ref_contig1 = ref_contig_name;
                    sv_adj.ref_pos = ref_local;

                    // If deletion crosses contig boundary, clamp len to contig end.
                    if (sv_adj.type == "DEL") {
                        uint32_t contig_len = ref_infos[ref_cid].len;
                        uint32_t end_local = ref_local + sv_adj.len;
                        if (end_local > contig_len) {
                            sv_adj.len = contig_len - ref_local;
                        }
                    }
                    indel_evs.push_back(std::move(sv_adj));
                }
                // --- K-mer block based complex SV discovery (INV/DUP/TRA) ---
                // Tuned for ~100 kb contigs and ~100–3000 bp events.
                {
                    const uint32_t min_split_piece = (uint32_t)std::max<uint32_t>(150, (uint32_t)args.sv_min);
                    const uint32_t tra_min_dist = 5000; // intra-contig "translocation-like" jump threshold
                    const uint32_t min_dup_ovl = (uint32_t)std::max<uint32_t>(50, min_split_piece / 5);

                    std::vector<KmerHit> hits_plus;
                    std::vector<KmerHit> hits_minus;
                    hits_plus.reserve(qry_seq_full.size() / 4);
                    hits_minus.reserve(qry_seq_full.size() / 4);

                    const uint32_t Rlen = (uint32_t)ref_concat.size();
                    const uint32_t Qlen = (uint32_t)qry_seq_full.size();

                    for (uint32_t q = 0; q + k_complex <= Qlen; q++) {
                        uint64_t key = 0;
                        if (!encode_kmer_2bit(qry_seq_full, q, k_complex, key)) continue;

                        auto itp = ref_kmer_pos_fwd.find(key);
                        if (itp != ref_kmer_pos_fwd.end()) {
                            hits_plus.push_back({q, itp->second, '+'});
                        }

                        auto itn = ref_kmer_pos_rc.find(key);
                        if (itn != ref_kmer_pos_rc.end()) {
                            uint32_t prc = itn->second;
                            if (prc + k_complex <= Rlen) {
                                uint32_t r_fwd = Rlen - (prc + k_complex);
                                hits_minus.push_back({q, r_fwd, '-'});
                            }
                        }
                    }

                    auto blocks_p = build_blocks_from_hits(hits_plus,  k_complex, /*max_q_gap=*/300, /*max_diag_dev=*/250,
                                                           min_split_piece, /*min_hits=*/12);
                    auto blocks_m = build_blocks_from_hits(hits_minus, k_complex, /*max_q_gap=*/300, /*max_diag_dev=*/250,
                                                           min_split_piece, /*min_hits=*/12);
// Pick blocks per orientation to ensure we keep inversion evidence even when '+' blocks dominate.
// IMPORTANT for DUP:
//   Duplications often create multiple '+' blocks that are *disjoint in query* but map to overlapping
//   reference intervals. If we keep only a tiny number of '+' blocks (e.g. 4), we can easily drop the
//   duplicated copy and never call DUP. With ~100 kb contigs, keeping up to ~20 non-overlapping '+'
//   blocks is still cheap and makes DUP detection reliable.
auto picked_p = pick_blocks_nonoverlap(blocks_p, /*N=*/20, /*min_q_sep=*/25);
auto picked_m = pick_blocks_nonoverlap(blocks_m, /*N=*/4,  /*min_q_sep=*/50);

std::vector<MapBlock> picked;
picked.reserve(picked_p.size() + picked_m.size());
picked.insert(picked.end(), picked_p.begin(), picked_p.end());
picked.insert(picked.end(), picked_m.begin(), picked_m.end());

// As a backstop, if we have any '-' blocks but none were picked (e.g., overlap),
// force-include the single best '-' block by hits.
if (!blocks_m.empty() && picked_m.empty()) {
    auto best_m = *std::max_element(blocks_m.begin(), blocks_m.end(),
        [](const MapBlock& a, const MapBlock& b){ return a.hits < b.hits; });
    picked.push_back(best_m);
}

std::sort(picked.begin(), picked.end(),
          [](const MapBlock& a, const MapBlock& b){ return a.q0 < b.q0; });

                    auto complex_evs = infer_dup_inv_tra_from_blocks(qry_seq_full, ref_infos, picked,
                                                                     k_complex, tra_min_dist, min_dup_ovl);

                    // Avoid counting the same event multiple times from fragmented blocks.
                    // Coarser bins reduce multi-counting of the same breakpoint produced by fragmented blocks.
                    // Values tuned for 100 kb contigs and ~100–3000 bp SVs.
                    complex_evs = dedup_svs(std::move(complex_evs), /*pos_bin=*/500, /*len_bin=*/500);

                    // Filter INDEL calls that are likely just manifestations of larger rearrangements.
                    // Heuristic: if an INDEL lies within +/-500 bp of a complex breakpoint on the same contig,
                    // drop it. This prevents double-counting INV/DUP/TRA as multiple small INDELs.
                    if (!complex_evs.empty() && !indel_evs.empty()) {
                        const uint32_t PAD = 200;
                        std::vector<SVEvent> kept;
                        kept.reserve(indel_evs.size());
                        for (const auto& e : indel_evs) {
                            bool near = false;
                            for (const auto& c : complex_evs) {
                                if (c.ref_contig1 != e.ref_contig1) continue;
                                uint32_t c0 = (c.ref_pos > PAD) ? (c.ref_pos - PAD) : 0;
                                uint32_t c1 = c.ref_pos + std::max<uint32_t>(c.len, 1u) + PAD;
                                if (e.ref_pos >= c0 && e.ref_pos <= c1) { near = true; break; }
                            }
                            if (!near) kept.push_back(e);
                        }
                        indel_evs.swap(kept);
                    }

                    // Add all SV alleles to graph (INDEL first, then complex) after filtering.
                    for (auto& e : indel_evs) {
                        add_sv_to_graph(g, sample, e.ref_contig1, contig.name, e.ref_pos, e, args.ref_segment);
                        total_svs_added++;
                    }
                    if (!complex_evs.empty()) {
                        std::cerr << "    [info] Complex SVs (kmer) : " << complex_evs.size() << "\n";
                        for (auto& ev : complex_evs) {
                            if (ev.ref_contig1.empty()) continue;
                            add_sv_to_graph(g, sample, ev.ref_contig1, contig.name, ev.ref_pos, ev, args.ref_segment);
                            total_svs_added++;
                        }
                    }
                }

            }

            std::cerr << "[info] Sample " << sample << ": total SV alleles added=" << total_svs_added << "\n";
        }

        g.write_gfa(args.out_gfa);
        std::cerr << "[done] Wrote " << args.out_gfa << "\n";
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "[error] " << e.what() << "\n";
        return 1;
    }
}

