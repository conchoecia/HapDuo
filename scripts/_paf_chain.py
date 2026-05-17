"""
Merge adjacent collinear PAF records into syntenic blocks.

minimap2 emits the same chain as many smaller alignments because indels,
repeats, and low-similarity patches break alignments apart. For the
ribbon plot we want one polygon per real synteny block, not 50.

Algorithm:
  1. Bucket records by (target_scaf, query_scaf, strand).
  2. Sort each bucket by target_start.
  3. Greedy-merge consecutive records when:
       - target_gap = next.t_start - prev.t_end  is in [-overlap, max_gap]
       - query_gap  = corresponding gap on the query (sign depends on strand)
         is in [-overlap, max_gap]
       - |t_gap - q_gap| / max(|t_gap|, |q_gap|, 1) <= rel_skew
     The relative-skew check rejects junctions where the query jumps a
     much larger or smaller distance than the target — that's not a single
     chain anymore.

For minus-strand records, "next on the target axis" means "previous on the
query axis," so the query-gap is prev.q_start - next.q_end (which is
positive when properly ordered).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from collections import defaultdict


@dataclass
class Block:
    tname: str
    qname: str
    strand: str
    t_start: int
    t_end: int
    q_start: int
    q_end: int
    n_records: int = 1
    blocklen: int = 0
    matched: int = 0

    def merge(self, other: "Block"):
        self.t_start = min(self.t_start, other.t_start)
        self.t_end   = max(self.t_end,   other.t_end)
        self.q_start = min(self.q_start, other.q_start)
        self.q_end   = max(self.q_end,   other.q_end)
        self.n_records += other.n_records
        self.blocklen += other.blocklen
        self.matched  += other.matched


def parse_paf(path: str):
    out = []
    with open(path) as f:
        for line in f:
            if not line.strip():
                continue
            p = line.rstrip("\n").split("\t")
            out.append(Block(
                tname=p[5], qname=p[0], strand=p[4],
                t_start=int(p[7]), t_end=int(p[8]),
                q_start=int(p[2]), q_end=int(p[3]),
                blocklen=int(p[10]),
                matched=int(p[9]),
            ))
    return out


def chain_paf(blocks, max_gap=200_000, overlap_tol=20_000, rel_skew=0.5,
              min_block_len=1000):
    """Merge adjacent collinear blocks.

    max_gap     -- maximum allowed t-gap and q-gap between consecutive
                   alignments to still be considered part of the same chain
    overlap_tol -- allow this much overlap on either axis (negative gap)
    rel_skew    -- |t_gap - q_gap| / max_gap must be <= this fraction
    min_block_len -- discard input records shorter than this before chaining
    """
    blocks = [b for b in blocks if b.blocklen >= min_block_len]
    buckets = defaultdict(list)
    for b in blocks:
        buckets[(b.tname, b.qname, b.strand)].append(b)
    merged = []
    for key, group in buckets.items():
        group.sort(key=lambda b: b.t_start)
        active = None
        for b in group:
            if active is None:
                active = Block(**b.__dict__)
                continue
            t_gap = b.t_start - active.t_end
            if b.strand == "+":
                q_gap = b.q_start - active.q_end
            else:
                # minus strand: as t increases, q decreases
                # active spans [active.q_start, active.q_end]
                # next   spans [b.q_start,      b.q_end]
                # gap on q axis = active.q_start - b.q_end
                q_gap = active.q_start - b.q_end
            if (-overlap_tol <= t_gap <= max_gap
                    and -overlap_tol <= q_gap <= max_gap
                    and abs(t_gap - q_gap) <=
                    max(max_gap * rel_skew, overlap_tol)):
                active.merge(b)
            else:
                merged.append(active)
                active = Block(**b.__dict__)
        if active is not None:
            merged.append(active)
    return merged
