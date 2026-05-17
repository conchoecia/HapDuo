#!/usr/bin/env python3
"""
Find inversions by detecting strand switches in a chained PAF.

Pipeline:
  1. Chain the PAF with tight parameters so each chain represents a single
     local alignment block (no merging across gaps that look like SVs).
  2. Drop chains shorter than --min-anchor on the shorter axis. These are
     the "anchors" used to call strand context.
  3. For each matched chr<N> ↔ chr<N> pair, sort anchors by target start.
  4. Walk the sorted list and identify runs of consecutive "-" anchors
     bounded by "+" anchors. Each such run is a candidate inversion.
  5. The breakpoints are the strand-switch positions, located at the
     boundary between the flanking "+" anchor and the first/last "-" anchor.

This intentionally does NOT merge "-" blocks separated by a gap that has
a "+" alignment between them — that's two separate inversions, not one.

Output TSV: same schema as find_inversions_chained.py.
"""
from __future__ import annotations
import argparse
import os
import re
import sys
from collections import defaultdict

from .paf_chain import parse_paf, chain_paf

CHR_RE = re.compile(r"_chr(\d+)$")


def chr_label(name):
    m = CHR_RE.search(name)
    return m.group(1) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paf")
    ap.add_argument("--min-anchor", type=int, default=50_000,
                    help="ignore chained blocks shorter than this on the "
                         "shorter axis (these are the strand-context anchors)")
    ap.add_argument("--min-inv", type=int, default=20_000,
                    help="minimum inversion length on the shorter axis")
    ap.add_argument("--chain-max-gap", type=int, default=50_000,
                    help="conservative max gap for chaining (tighter than "
                         "the plotting chainer so anchors stay local)")
    ap.add_argument("--chain-rel-skew", type=float, default=0.2)
    ap.add_argument("--min-input-block", type=int, default=1000)
    ap.add_argument("--length-ratio", type=float, default=0.5)
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--max-intrusion-gap", type=int, default=1_500_000,
                    help="t-axis gap tolerated between successive '-' anchors "
                         "in the same inversion. Short '+' intrusions inside "
                         "an inversion are absorbed if both flanking '-' "
                         "anchors are within this distance.")
    args = ap.parse_args()

    blocks = parse_paf(args.paf)
    print(f"# loaded {len(blocks)} PAF records", file=sys.stderr)

    chained = chain_paf(blocks,
                        max_gap=args.chain_max_gap,
                        rel_skew=args.chain_rel_skew,
                        min_block_len=args.min_input_block)
    print(f"# {len(chained)} tight-chained blocks", file=sys.stderr)

    # filter to anchors above min_anchor on the shorter axis
    anchors = []
    for b in chained:
        if chr_label(b.tname) is None or chr_label(b.tname) != chr_label(b.qname):
            continue
        t_len = b.t_end - b.t_start
        q_len = b.q_end - b.q_start
        if min(t_len, q_len) < args.min_anchor:
            continue
        anchors.append(b)
    print(f"# {len(anchors)} anchors >= {args.min_anchor:,} bp on matched chr pairs",
          file=sys.stderr)

    # group by chr pair
    pairs = defaultdict(list)
    for b in anchors:
        pairs[(b.tname, b.qname)].append(b)

    def keep_antiparallel_chain(minus_run):
        """Within a run of '-' anchors, keep the longest subsequence whose
        q-coords are monotonically decreasing as t increases (true
        antiparallel synteny). This trims rogue anchors that point to a
        distant hapB region (e.g. tandem inverted repeats elsewhere)."""
        if not minus_run:
            return []
        run = sorted(minus_run, key=lambda b: b.t_start)
        # LIS on negated q_start — equivalently LDS on q_start
        n = len(run)
        # Use bp-weighted longest chain: each block contributes its shorter-axis
        # length to a "score" so we pick the chain that covers the most sequence
        weight = [min(b.t_end - b.t_start, b.q_end - b.q_start) for b in run]
        best = [weight[k] for k in range(n)]
        prev = [-1] * n
        for k in range(n):
            for m in range(k):
                # m must come before k in t and AFTER k in q (antiparallel)
                if run[m].q_start > run[k].q_start:
                    if best[m] + weight[k] > best[k]:
                        best[k] = best[m] + weight[k]
                        prev[k] = m
        end = max(range(n), key=lambda k: best[k])
        chain_idx = []
        while end != -1:
            chain_idx.append(end)
            end = prev[end]
        chain_idx.reverse()
        return [run[k] for k in chain_idx]

    candidates = []
    for (tname, qname), group in pairs.items():
        group.sort(key=lambda b: b.t_start)
        n = len(group)
        consumed = [False] * n
        i = 0
        while i < n:
            if consumed[i] or group[i].strand != "-":
                i += 1
                continue
            # Greedy extension: gather all unconsumed '-' anchors whose
            # t_start lies within --max-intrusion-gap of the running '-' run's
            # right edge. '+' anchors are skipped (not added, not used as
            # terminator) as long as a '-' anchor reappears before the gap
            # closes — this handles inversions like ----+----+---- where the
            # inner '+' anchors are repeat artifacts inside one true inversion.
            minus_idx = [i]
            run_t_right = group[i].t_end
            k = i + 1
            while k < n:
                b = group[k]
                if b.t_start - run_t_right > args.max_intrusion_gap:
                    break
                if b.strand == "-" and not consumed[k]:
                    minus_idx.append(k)
                    run_t_right = max(run_t_right, b.t_end)
                k += 1
            minus_run = keep_antiparallel_chain([group[m] for m in minus_idx])
            if not minus_run:
                consumed[i] = True
                i += 1
                continue
            t_left  = min(b.t_start for b in minus_run)
            t_right = max(b.t_end   for b in minus_run)
            q_left  = min(b.q_start for b in minus_run)
            q_right = max(b.q_end   for b in minus_run)
            inv_len = min(t_right - t_left, q_right - q_left)
            longer  = max(t_right - t_left, q_right - q_left)
            if longer == 0 or inv_len / longer < args.length_ratio:
                consumed[i] = True
                i += 1
                continue
            if inv_len < args.min_inv:
                consumed[i] = True
                i += 1
                continue
            # Refine boundaries using the nearest '+' anchor that lies wholly
            # outside [t_left, t_right] — i.e. flanking the inversion, not
            # buried inside it. Only refine if the flanking '+' anchor is
            # within --max-intrusion-gap of the chain edge; otherwise the
            # nearest '+' is in unrelated territory and the raw chain bound is
            # already the safer estimate.
            refine_cap = args.max_intrusion_gap
            t_left_break = t_left
            for prev in range(i - 1, -1, -1):
                if group[prev].strand != "+":
                    continue
                if group[prev].t_end > t_left:
                    continue
                if t_left - group[prev].t_end > refine_cap:
                    break
                t_left_break = (group[prev].t_end + t_left) // 2
                break
            t_right_break = t_right
            for nxt in range(k, n):
                if group[nxt].strand != "+":
                    continue
                if group[nxt].t_start < t_right:
                    continue
                if group[nxt].t_start - t_right > refine_cap:
                    break
                t_right_break = (t_right + group[nxt].t_start) // 2
                break

            # On query axis: use the raw bounds of the trimmed antiparallel
            # chain. (An earlier attempt to interpolate with the flanking '+'
            # anchor's q_end collapsed both q breakpoints toward the inversion
            # midpoint — wrong.)
            q_right_break = q_right
            q_left_break  = q_left

            # Dedup: mark every "-" anchor that participated in this chain
            # (or that overlaps the inversion's t-extent) as consumed, so the
            # next outer iteration doesn't re-emit the same inversion starting
            # from an interior anchor.
            chain_ids = set(id(b) for b in minus_run)
            for m in range(n):
                b = group[m]
                if b.strand != "-":
                    continue
                if id(b) in chain_ids:
                    consumed[m] = True
                    continue
                if b.t_end > t_left and b.t_start < t_right:
                    consumed[m] = True

            candidates.append({
                "tname": tname, "qname": qname,
                "t_left": t_left_break, "t_right": t_right_break,
                "q_left": q_left_break, "q_right": q_right_break,
                "inv_len": inv_len,
                "n_anchors": len(minus_run),
                "n_total": k - i,
            })
            i += 1

    candidates.sort(key=lambda c: c["inv_len"], reverse=True)

    # Post-hoc dedup: drop any candidate whose t-extent substantially overlaps
    # a higher-ranked candidate on the same chromosome pair. "Substantial"
    # means reciprocal overlap >= 0.30 (Jaccard-like). Two candidates that
    # share <30% of their t-bases are kept as distinct calls.
    def t_overlap_frac(a, b):
        if a["tname"] != b["tname"] or a["qname"] != b["qname"]:
            return 0.0
        lo = max(a["t_left"], b["t_left"])
        hi = min(a["t_right"], b["t_right"])
        inter = max(0, hi - lo)
        union = max(a["t_right"], b["t_right"]) - min(a["t_left"], b["t_left"])
        return inter / union if union > 0 else 0.0

    kept = []
    for c in candidates:
        if any(t_overlap_frac(c, k) >= 0.30 for k in kept):
            continue
        kept.append(c)
    top = kept[:args.top]

    print(f"# top {len(top)} inversions:", file=sys.stderr)
    for i, c in enumerate(top, 1):
        print(f"#  {i}: {c['tname']}:{c['t_left']}-{c['t_right']}  vs  "
              f"{c['qname']}:{c['q_left']}-{c['q_right']}  "
              f"len={c['inv_len']:,}  n_anchors={c['n_anchors']}",
              file=sys.stderr)

    print("rank\thapA_chr\thapA_pos\thapB_chr\thapB_pos\tend\tinv_len\tn_anchors")
    for rank, c in enumerate(top, 1):
        # left end: hapA = t_left, hapB = q_right (antiparallel pairing)
        print(f"{rank}\t{c['tname']}\t{c['t_left']}\t{c['qname']}\t{c['q_right']}\t"
              f"left\t{c['inv_len']}\t{c['n_anchors']}")
        # right end: hapA = t_right, hapB = q_left
        print(f"{rank}\t{c['tname']}\t{c['t_right']}\t{c['qname']}\t{c['q_left']}\t"
              f"right\t{c['inv_len']}\t{c['n_anchors']}")


if __name__ == "__main__":
    main()
