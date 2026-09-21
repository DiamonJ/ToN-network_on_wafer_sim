#!/usr/bin/env python3
"""
LAMMPS MPI Communication Feature Analysis Tool
==============================================
Parses CCDG files from DUMPI traces and generates:
  1. MPI operation type distribution (pie/bar charts)
  2. Rank-to-rank communication heatmap
  3. Compute vs communication ratio
  4. Cross-dimensional comparison (atom count vs rank count)
"""

import json
import os
import sys
import re
import glob
from collections import Counter, defaultdict

import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    plt.rcParams['font.size'] = 12
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not installed — will only output JSON stats")


# ─── Node type mapping ──────────────────────────────────────────────────────
COMM_TYPES = {
    'SEND': 'SEND', 'RECV': 'RECV',
    'ISEND': 'ISEND', 'IRECV': 'IRECV',
    'WAIT': 'WAIT', 'WAITALL': 'WAITALL',
    'ALLREDUCE': 'ALLREDUCE', 'BARRIER': 'BARRIER',
    'BCAST': 'BCAST', 'GATHER': 'GATHER',
    'ALLGATHER': 'ALLGATHER', 'SCATTER': 'SCATTER',
    'ALLTOALL': 'ALLTOALL', 'REDUCE': 'REDUCE',
    'OTHER': 'OTHER',
}
COMM_TYPE_ORDER = [
    'ALLREDUCE', 'BARRIER', 'BCAST', 'SEND', 'RECV',
    'ISEND', 'IRECV', 'WAIT', 'WAITALL',
    'GATHER', 'ALLGATHER', 'SCATTER', 'ALLTOALL', 'REDUCE',
    'OTHER'
]

# Colors for MPI operation types
COMM_COLORS = {
    'ALLREDUCE': '#E74C3C', 'BARRIER': '#8E44AD', 'BCAST': '#3498DB',
    'SEND': '#2ECC71', 'RECV': '#1ABC9C',
    'ISEND': '#27AE60', 'IRECV': '#16A085',
    'WAIT': '#F39C12', 'WAITALL': '#E67E22',
    'GATHER': '#2980B9', 'ALLGATHER': '#2471A3',
    'SCATTER': '#5DADE2', 'ALLTOALL': '#A569BD',
    'REDUCE': '#EC7063', 'OTHER': '#95A5A6',
}


# ─── CCDG parsing ────────────────────────────────────────────────────────────

def parse_ccdg(filepath):
    """Parse a CCDG file, extract nodes and cross_rank_edges.
    Returns dict with num_ranks, nodes, cross_rank_edges, and stats header.
    """
    if not os.path.exists(filepath):
        return None

    with open(filepath, 'r') as f:
        content = f.read()

    # Find JSON start
    json_start = content.find('{')
    if json_start == -1:
        return None

    # Header stats (text before JSON)
    header = content[:json_start]
    stats = {
        'total_nodes': _extract_int(header, r'Total nodes:\s+(\d+)'),
        'compute_nodes': _extract_int(header, r'Compute nodes:\s+(\d+)'),
        'comm_nodes': _extract_int(header, r'Communication nodes:\s+(\d+)'),
        'total_compute_time': _extract_float(header, r'Total compute time:\s+([\d.]+)'),
        'cross_edges': _extract_int(header, r'Cross-rank edges:\s+(\d+)'),
    }

    try:
        data = json.loads(content[json_start:])
    except json.JSONDecodeError as e:
        print(f"  WARNING: JSON parse error: {e}")
        return None

    return {
        'stats': stats,
        'num_ranks': data.get('num_ranks', 0),
        'nodes': data.get('nodes', []),
        'cross_rank_edges': data.get('cross_rank_edges', []),
    }


def _extract_int(text, pattern):
    m = re.search(pattern, text)
    return int(m.group(1)) if m else None


def _extract_float(text, pattern):
    m = re.search(pattern, text)
    return float(m.group(1)) if m else None


# ─── Communication feature extraction ────────────────────────────────────────

def analyze_ccdg(ccdg_data):
    """Extract communication features from parsed CCDG data."""
    if ccdg_data is None:
        return None

    nodes = ccdg_data['nodes']
    num_ranks = ccdg_data['num_ranks']
    stats = ccdg_data['stats']

    # 1. MPI operation type distribution
    type_counter = Counter()
    for node in nodes:
        ntype = node.get('type', 'UNKNOWN')
        if ntype in COMM_TYPES:
            type_counter[ntype] += 1
        elif ntype == 'COMPUTE':
            pass  # Not a comm type
        else:
            type_counter['OTHER'] += 1

    total_comm = sum(type_counter.values())

    # 2. Rank-to-rank communication matrix (from SEND/ISEND with comm_bytes)
    comm_matrix = np.zeros((num_ranks, num_ranks), dtype=np.float64)
    comm_count_matrix = np.zeros((num_ranks, num_ranks), dtype=np.int64)
    for node in nodes:
        ntype = node.get('type', '')
        if ntype in ('SEND', 'ISEND'):
            src = node.get('comm_src')
            dst = node.get('comm_dst')
            bytes_val = node.get('comm_bytes', 0)
            if src is not None and dst is not None and bytes_val:
                comm_matrix[src, dst] += bytes_val
                comm_count_matrix[src, dst] += 1

    # 3. Compute vs communication ratio
    total_compute_nodes = stats.get('compute_nodes', 0)
    total_comm_nodes = stats.get('comm_nodes', 0)
    total_compute_time = stats.get('total_compute_time', 0)

    # Per-rank node counts
    rank_node_count = Counter()
    rank_compute_count = Counter()
    rank_comm_count = Counter()
    for node in nodes:
        rank = node.get('rank')
        ntype = node.get('type', '')
        rank_node_count[rank] += 1
        if ntype == 'COMPUTE':
            rank_compute_count[rank] += 1
        elif ntype in COMM_TYPES:
            rank_comm_count[rank] += 1

    # 4. Communication volume per rank
    rank_send_bytes = np.zeros(num_ranks, dtype=np.float64)
    rank_recv_bytes = np.zeros(num_ranks, dtype=np.float64)
    for node in nodes:
        ntype = node.get('type', '')
        if ntype == 'SEND' or ntype == 'ISEND':
            src = node.get('comm_src')
            dst = node.get('comm_dst')
            bytes_val = node.get('comm_bytes', 0)
            if bytes_val and src is not None:
                rank_send_bytes[src] += bytes_val
            if bytes_val and dst is not None:
                rank_recv_bytes[dst] += bytes_val

    # 5. Collective operation breakdown
    collective_counter = Counter()
    for ntype in ('ALLREDUCE', 'BARRIER', 'BCAST', 'GATHER', 'ALLGATHER',
                  'SCATTER', 'ALLTOALL', 'REDUCE'):
        collective_counter[ntype] = type_counter.get(ntype, 0)

    # 6. P2P vs collective
    p2p_types = {'SEND', 'RECV', 'ISEND', 'IRECV', 'WAIT', 'WAITALL'}
    coll_types = {'ALLREDUCE', 'BARRIER', 'BCAST', 'GATHER', 'ALLGATHER',
                  'SCATTER', 'ALLTOALL', 'REDUCE'}
    p2p_count = sum(type_counter.get(t, 0) for t in p2p_types)
    coll_count = sum(type_counter.get(t, 0) for t in coll_types)

    return {
        'num_ranks': num_ranks,
        'type_distribution': dict(type_counter),
        'total_comm_nodes': int(total_comm_nodes),
        'total_compute_nodes': int(total_compute_nodes),
        'total_comm_ops': int(total_comm),
        'compute_comm_ratio': total_compute_nodes / max(total_comm_nodes, 1),
        'total_compute_time_sec': total_compute_time,
        'comm_matrix': comm_matrix.tolist(),
        'comm_count_matrix': comm_count_matrix.tolist(),
        'rank_send_bytes': rank_send_bytes.tolist(),
        'rank_recv_bytes': rank_recv_bytes.tolist(),
        'collective_breakdown': dict(collective_counter),
        'p2p_count': int(p2p_count),
        'coll_count': int(coll_count),
        'rank_node_counts': {
            int(r): {'total': rank_node_count[r],
                     'compute': rank_compute_count[r],
                     'comm': rank_comm_count[r]}
            for r in sorted(rank_node_count.keys())
        },
    }


# ─── Visualization ──────────────────────────────────────────────────────────

def plot_operation_distribution(results, output_dir):
    """Plot MPI operation type distribution as pie + bar chart."""
    type_dist = results['type_distribution']
    if not type_dist:
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Pie chart (top 8 types + others)
    sorted_types = sorted(type_dist.items(), key=lambda x: -x[1])
    labels = []
    sizes = []
    colors = []
    other = 0
    for t, c in sorted_types:
        if len(labels) < 8:
            labels.append(t)
            sizes.append(c)
            colors.append(COMM_COLORS.get(t, '#95A5A6'))
        else:
            other += c
    if other > 0:
        labels.append('OTHER')
        sizes.append(other)
        colors.append('#95A5A6')

    axes[0].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%',
                startangle=90, textprops={'fontsize': 10})
    axes[0].set_title(f'MPI Operation Distribution\n({results["total_comm_ops"]} total ops)',
                      fontsize=13, fontweight='bold')

    # Bar chart
    all_types = sorted(type_dist.items(), key=lambda x: -x[1])
    type_names = [t for t, _ in all_types]
    type_counts = [c for _, c in all_types]
    bar_colors = [COMM_COLORS.get(t, '#95A5A6') for t in type_names]

    bars = axes[1].barh(range(len(type_names)), type_counts,
                        color=bar_colors, edgecolor='white')
    axes[1].set_yticks(range(len(type_names)))
    axes[1].set_yticklabels(type_names, fontsize=10)
    axes[1].set_xlabel('Count')
    axes[1].set_title('MPI Operation Counts', fontsize=13, fontweight='bold')

    for bar, count in zip(bars, type_counts):
        axes[1].text(bar.get_width() + max(type_counts) * 0.01,
                     bar.get_y() + bar.get_height() / 2,
                     str(count), va='center', fontsize=9)

    plt.tight_layout()
    path = os.path.join(output_dir, '01_operation_distribution.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_heatmap(results, output_dir):
    """Plot rank-to-rank communication heatmap."""
    comm_matrix = np.array(results['comm_matrix'])
    comm_count_matrix = np.array(results['comm_count_matrix'])
    num_ranks = results['num_ranks']

    if num_ranks == 0 or comm_matrix.size == 0:
        return

    has_comm = comm_matrix.max() > 0

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Bytes heatmap
    total_bytes = comm_matrix.sum()
    if has_comm:
        min_b = comm_matrix[comm_matrix > 0].min()
        vmin = max(min_b, 1)
        im1 = axes[0].imshow(comm_matrix, cmap='YlOrRd',
                             norm=LogNorm(vmin=vmin,
                                          vmax=max(comm_matrix.max(), 1)),
                             aspect='auto')
    else:
        im1 = axes[0].imshow(comm_matrix, cmap='YlOrRd',
                             aspect='auto')
    axes[0].set_title(f'Comm Volume (bytes)\nTotal: {total_bytes/1024:.1f} KB',
                      fontsize=13, fontweight='bold')
    axes[0].set_xlabel('Destination Rank')
    axes[0].set_ylabel('Source Rank')
    plt.colorbar(im1, ax=axes[0], label='Bytes (log scale)')

    # Add text annotations (only for non-zero)
    for i in range(num_ranks):
        for j in range(num_ranks):
            if comm_matrix[i, j] > 0:
                val_kb = comm_matrix[i, j] / 1024
                text = f'{val_kb:.0f}K' if val_kb >= 1 else f'{comm_matrix[i, j]:.0f}'
                axes[0].text(j, i, text, ha='center', va='center',
                            fontsize=7, color='black' if comm_matrix[i, j] < comm_matrix.max()/2 else 'white')

    # Count heatmap
    total_packets = comm_count_matrix.sum()
    if has_comm:
        vmin_c = max(comm_count_matrix[comm_count_matrix > 0].min(), 1)
        im2 = axes[1].imshow(comm_count_matrix, cmap='Blues',
                             norm=LogNorm(vmin=vmin_c,
                                          vmax=max(comm_count_matrix.max(), 1)),
                             aspect='auto')
    else:
        im2 = axes[1].imshow(comm_count_matrix, cmap='Blues',
                             aspect='auto')
    axes[1].set_title(f'Comm Count (messages)\nTotal: {int(total_packets)} pkts',
                      fontsize=13, fontweight='bold')
    axes[1].set_xlabel('Destination Rank')
    axes[1].set_ylabel('Source Rank')
    plt.colorbar(im2, ax=axes[1], label='Count (log scale)')

    for i in range(num_ranks):
        for j in range(num_ranks):
            if comm_count_matrix[i, j] > 0:
                axes[1].text(j, i, str(int(comm_count_matrix[i, j])),
                            ha='center', va='center',
                            fontsize=7, color='black' if comm_count_matrix[i, j] < comm_count_matrix.max()/2 else 'white')

    plt.tight_layout()
    path = os.path.join(output_dir, '02_rank_comm_heatmap.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_compute_comm_ratio(results, output_dir):
    """Plot compute vs communication node count comparison."""
    compute = results['total_compute_nodes']
    comm = results['total_comm_nodes']
    ratio = results['compute_comm_ratio']

    fig, ax = plt.subplots(figsize=(8, 5))
    categories = ['Compute', 'Communication']
    values = [compute, comm]
    colors = ['#2ECC71', '#E74C3C']

    bars = ax.bar(categories, values, color=colors, edgecolor='white', width=0.5)
    ax.set_ylabel('Node Count')
    ax.set_title(f'Compute vs Communication Nodes\nRatio={ratio:.2f}',
                 fontsize=14, fontweight='bold')

    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(values) * 0.02,
                str(val), ha='center', va='bottom', fontsize=12)

    plt.tight_layout()
    path = os.path.join(output_dir, '03_compute_comm_ratio.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_rank_volume(results, output_dir):
    """Plot per-rank send/recv volume."""
    send = np.array(results['rank_send_bytes'])
    recv = np.array(results['rank_recv_bytes'])
    ranks = np.arange(len(send))

    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(ranks))
    width = 0.35

    ax.bar(x - width / 2, send / 1024, width, label='Send', color='#3498DB',
           edgecolor='white')
    ax.bar(x + width / 2, recv / 1024, width, label='Recv', color='#E74C3C',
           edgecolor='white')

    ax.set_xlabel('Rank')
    ax.set_ylabel('Volume (KB)')
    ax.set_title('Per-Rank Communication Volume', fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels([str(r) for r in ranks])
    ax.legend()

    plt.tight_layout()
    path = os.path.join(output_dir, '04_rank_comm_volume.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_collective_breakdown(results, output_dir):
    """Bar chart of collective operation breakdown."""
    collective = results.get('collective_breakdown', {})
    if not collective:
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    types = []
    counts = []
    colors = []
    for t in COMM_TYPE_ORDER:
        if t in collective and collective[t] > 0:
            types.append(t)
            counts.append(collective[t])
            colors.append(COMM_COLORS.get(t, '#95A5A6'))

    bars = ax.bar(range(len(types)), counts, color=colors, edgecolor='white')
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(types, rotation=30, ha='right', fontsize=10)
    ax.set_ylabel('Count')
    ax.set_title('Collective Operation Breakdown', fontsize=14, fontweight='bold')

    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(counts) * 0.01,
                str(c), ha='center', fontsize=10)

    plt.tight_layout()
    path = os.path.join(output_dir, '05_collective_breakdown.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def generate_all_plots(results, output_dir):
    """Generate all visualization plots for one experiment."""
    if not HAS_MPL:
        print("  SKIP: matplotlib not available")
        return

    os.makedirs(output_dir, exist_ok=True)
    plot_operation_distribution(results, output_dir)
    plot_heatmap(results, output_dir)
    plot_compute_comm_ratio(results, output_dir)
    plot_rank_volume(results, output_dir)
    plot_collective_breakdown(results, output_dir)


# ─── Comparative analysis ───────────────────────────────────────────────────

def plot_comparative_scaling(all_results, output_dir):
    """Generate comparative plots across experiments.
    all_results: list of dicts with keys: case, atoms, ranks, results
    """
    if not HAS_MPL or len(all_results) < 2:
        return

    # Group by case
    cases_data = defaultdict(list)
    for entry in all_results:
        case_key = entry['case']
        cases_data[case_key].append(entry)

    # 1. Comm ratio vs ranks (lines, one per case)
    fig, ax = plt.subplots(figsize=(10, 6))
    markers = ['o', 's', '^', 'D']
    for idx, (case, entries) in enumerate(sorted(cases_data.items())):
        entries_sorted = sorted(entries, key=lambda x: int(x['ranks']))
        ranks = [int(e['ranks']) for e in entries_sorted]
        ratios = [
            e['results']['compute_comm_ratio']
            if e['results'] else 0 for e in entries_sorted
        ]
        label_map = {'lj_0.1k': 'LJ 0.1k', 'lj_1k': 'LJ 1k',
                     'lj_10k': 'LJ 10k', 'lj_mix_1k': 'LJ-mix 1k'}
        ax.plot(ranks, ratios, marker=markers[idx % len(markers)],
                label=label_map.get(case, case), linewidth=2, markersize=8)

    ax.set_xlabel('Number of Ranks', fontsize=13)
    ax.set_ylabel('Compute / Communication Ratio', fontsize=13)
    ax.set_title('Compute-Comm Ratio vs Ranks (Strong Scaling)',
                 fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log', base=2)
    ax.set_xticks([4, 8, 16, 32])
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

    plt.tight_layout()
    path = os.path.join(output_dir, '10_scaling_comm_ratio.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

    # 2. Comm operation count vs ranks
    fig, ax = plt.subplots(figsize=(10, 6))
    for idx, (case, entries) in enumerate(sorted(cases_data.items())):
        entries_sorted = sorted(entries, key=lambda x: int(x['ranks']))
        ranks = [int(e['ranks']) for e in entries_sorted]
        comm_counts = [
            e['results']['total_comm_ops']
            if e['results'] else 0 for e in entries_sorted
        ]
        label_map = {'lj_0.1k': 'LJ 0.1k', 'lj_1k': 'LJ 1k',
                     'lj_10k': 'LJ 10k', 'lj_mix_1k': 'LJ-mix 1k'}
        ax.plot(ranks, comm_counts, marker=markers[idx % len(markers)],
                label=label_map.get(case, case), linewidth=2, markersize=8)

    ax.set_xlabel('Number of Ranks', fontsize=13)
    ax.set_ylabel('Total Communication Operations', fontsize=13)
    ax.set_title('Comm Operations vs Ranks', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log', base=2)
    ax.set_xticks([4, 8, 16, 32])
    ax.get_xaxis().set_major_formatter(plt.ScalarFormatter())

    plt.tight_layout()
    path = os.path.join(output_dir, '11_scaling_comm_ops.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")

    # 3. P2P vs Collective breakdown (stacked bar, one per experiment)
    fig, ax = plt.subplots(figsize=(14, 6))
    labels = []
    p2p_vals = []
    coll_vals = []
    coll_colors_list = []

    # Group by case then sort by ranks
    for case in ['lj_0.1k', 'lj_1k', 'lj_10k', 'lj_mix_1k']:
        if case not in cases_data:
            continue
        entries_sorted = sorted(cases_data[case], key=lambda x: int(x['ranks']))
        label_map = {'lj_0.1k': 'LJ 0.1k', 'lj_1k': 'LJ 1k',
                     'lj_10k': 'LJ 10k', 'lj_mix_1k': 'LJ-mix 1k'}
        for e in entries_sorted:
            r = e['results']
            if r:
                labels.append(f"{label_map.get(case, case)}\n{e['ranks']}r")
                p2p_vals.append(r['p2p_count'])
                coll_vals.append(r['coll_count'])

    x = np.arange(len(labels))
    width = 0.6

    bars1 = ax.bar(x, p2p_vals, width, label='P2P (Send/Recv/Wait)',
                   color='#3498DB', edgecolor='white')
    bars2 = ax.bar(x, coll_vals, width, bottom=p2p_vals,
                   label='Collective (AR/Barrier/BCast/...)',
                   color='#E74C3C', edgecolor='white')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9, rotation=45, ha='right')
    ax.set_ylabel('Count')
    ax.set_title('P2P vs Collective Communication', fontsize=14, fontweight='bold')
    ax.legend(fontsize=11)

    plt.tight_layout()
    path = os.path.join(output_dir, '12_p2p_vs_coll.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ─── Rank heatmap grid (comparative, across rank counts) ─────────────────────

def plot_rank_heatmap_grid(all_results, output_dir):
    """Generate a 2×2 grid of rank-to-rank heatmaps for 4/8/16/32 ranks.
    Uses a representative case (e.g., lj_1k) to show the scaling trend.
    If lj_1k isn't available, uses the first case with data for all 4 rank counts.
    """
    if not HAS_MPL or len(all_results) < 4:
        return

    # Find the best case (lj_1k preferred, otherwise first complete case)
    case_groups = defaultdict(list)
    for entry in all_results:
        case_groups[entry['case']].append(entry)

    target_case = None
    for preferred in ['lj_1k', 'lj_0.1k', 'lj_10k', 'lj_mix_1k']:
        if preferred in case_groups:
            entries = case_groups[preferred]
            ranks_present = sorted(int(e['ranks']) for e in entries if e['results'])
            if set([4, 8, 16, 32]).issubset(set(ranks_present)):
                target_case = preferred
                break

    if target_case is None:
        # Just use the first case with the most rank entries
        target_case = max(case_groups.keys(),
                          key=lambda c: len(case_groups[c]))
        entries = [e for e in all_results if e['case'] == target_case
                   and e['results']]
    else:
        entries = [e for e in all_results if e['case'] == target_case
                   and e['results']]

    entries_sorted = sorted(entries, key=lambda x: int(x['ranks']))

    label_map = {'lj_0.1k': 'LJ 0.1k', 'lj_1k': 'LJ 1k',
                 'lj_10k': 'LJ 10k', 'lj_mix_1k': 'LJ-mix 1k'}
    case_label = label_map.get(target_case, target_case)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes_flat = axes.flatten()

    # Determine global color range across all rank counts
    all_mats = []
    for e in entries_sorted:
        mat = np.array(e['results']['comm_matrix'])
        if mat.max() > 0:
            all_mats.append(mat)

    if not all_mats:
        plt.close()
        return

    global_min = min(m[m > 0].min() for m in all_mats if m.max() > 0)
    global_max = max(m.max() for m in all_mats)

    for idx, e in enumerate(entries_sorted):
        if idx >= 4:
            break
        ax = axes_flat[idx]
        mat = np.array(e['results']['comm_matrix'])
        num_ranks = e['results']['num_ranks']

        # Show bytes in KB
        mat_kb = mat / 1024.0

        im = ax.imshow(mat_kb, cmap='YlOrRd',
                       norm=LogNorm(vmin=max(global_min / 1024, 0.01),
                                    vmax=global_max / 1024),
                       aspect='auto')
        ax.set_title(f'{case_label} · {num_ranks} Ranks',
                     fontsize=14, fontweight='bold')
        ax.set_xlabel('Destination Rank', fontsize=11)
        ax.set_ylabel('Source Rank', fontsize=11)

        # Add text annotations
        for i in range(num_ranks):
            for j in range(num_ranks):
                val_kb = mat_kb[i, j]
                if val_kb > 0:
                    text = f'{val_kb:.1f}K' if val_kb >= 1 else f'{mat[i,j]:.0f}B'
                    color = 'white' if mat[i, j] > global_max * 0.5 else 'black'
                    ax.text(j, i, text, ha='center', va='center',
                            fontsize=8, color=color)

        ax.set_xticks(range(num_ranks))
        ax.set_yticks(range(num_ranks))

    # Share colorbar
    cbar = fig.colorbar(im, ax=axes_flat, fraction=0.02, pad=0.02)
    cbar.set_label('Communication Volume (KB, log scale)', fontsize=11)

    plt.suptitle('Rank-to-Rank Communication Volume\n'
                 f'(Case: {case_label}, CPU Freq: 2.8 GHz)',
                 fontsize=16, fontweight='bold', y=1.01)

    plt.tight_layout()
    path = os.path.join(output_dir, '20_rank_heatmap_grid.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def plot_rank_heatmap_grid_bytes(all_results, output_dir):
    """Same as above but uses comm_count_matrix (message count) instead of bytes."""
    if not HAS_MPL or len(all_results) < 4:
        return

    case_groups = defaultdict(list)
    for entry in all_results:
        case_groups[entry['case']].append(entry)

    target_case = None
    for preferred in ['lj_1k', 'lj_0.1k', 'lj_10k', 'lj_mix_1k']:
        if preferred in case_groups:
            entries = case_groups[preferred]
            ranks_present = sorted(int(e['ranks']) for e in entries if e['results'])
            if set([4, 8, 16, 32]).issubset(set(ranks_present)):
                target_case = preferred
                break

    if target_case is None:
        target_case = max(case_groups.keys(),
                          key=lambda c: len(case_groups[c]))
        entries = [e for e in all_results if e['case'] == target_case
                   and e['results']]
    else:
        entries = [e for e in all_results if e['case'] == target_case
                   and e['results']]

    entries_sorted = sorted(entries, key=lambda x: int(x['ranks']))
    label_map = {'lj_0.1k': 'LJ 0.1k', 'lj_1k': 'LJ 1k',
                 'lj_10k': 'LJ 10k', 'lj_mix_1k': 'LJ-mix 1k'}
    case_label = label_map.get(target_case, target_case)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes_flat = axes.flatten()

    all_mats = []
    for e in entries_sorted:
        mat = np.array(e['results']['comm_count_matrix'])
        if mat.max() > 0:
            all_mats.append(mat)

    if not all_mats:
        plt.close()
        return

    global_min = min(m[m > 0].min() for m in all_mats if m.max() > 0)
    global_max = max(m.max() for m in all_mats)

    for idx, e in enumerate(entries_sorted):
        if idx >= 4:
            break
        ax = axes_flat[idx]
        mat = np.array(e['results']['comm_count_matrix'])
        num_ranks = e['results']['num_ranks']

        im = ax.imshow(mat, cmap='Blues',
                       norm=LogNorm(vmin=max(global_min, 1),
                                    vmax=global_max),
                       aspect='auto')
        ax.set_title(f'{case_label} · {num_ranks} Ranks',
                     fontsize=14, fontweight='bold')
        ax.set_xlabel('Destination Rank', fontsize=11)
        ax.set_ylabel('Source Rank', fontsize=11)

        for i in range(num_ranks):
            for j in range(num_ranks):
                val = mat[i, j]
                if val > 0:
                    color = 'white' if val > global_max * 0.5 else 'black'
                    ax.text(j, i, str(int(val)), ha='center', va='center',
                            fontsize=8, color=color)

        ax.set_xticks(range(num_ranks))
        ax.set_yticks(range(num_ranks))

    cbar = fig.colorbar(im, ax=axes_flat, fraction=0.02, pad=0.02)
    cbar.set_label('Message Count (log scale)', fontsize=11)

    plt.suptitle('Rank-to-Rank Message Count\n'
                 f'(Case: {case_label}, CPU Freq: 2.8 GHz)',
                 fontsize=16, fontweight='bold', y=1.01)

    plt.tight_layout()
    path = os.path.join(output_dir, '21_rank_msgcount_grid.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ─── Main entry point ────────────────────────────────────────────────────────

def process_single_experiment(ccdg_file, output_dir):
    """Process a single CCDG file and generate all outputs."""
    print(f"  Parsing: {ccdg_file}")
    data = parse_ccdg(ccdg_file)
    if data is None:
        print("  FAILED: Could not parse CCDG")
        return None

    results = analyze_ccdg(data)
    if results is None:
        print("  FAILED: Could not analyze CCDG")
        return None

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Save JSON stats
    stats_json = {k: v for k, v in results.items()
                  if k not in ('comm_matrix', 'comm_count_matrix',
                               'rank_send_bytes', 'rank_recv_bytes')}
    stats_path = os.path.join(output_dir, 'stats.json')
    with open(stats_path, 'w') as f:
        json.dump(stats_json, f, indent=2)
    print(f"  Saved: {stats_path}")

    # Generate plots
    generate_all_plots(results, output_dir)

    return results


def scan_runs(runs_dir):
    """Scan runs directory and find all CCDG files."""
    results_db = []
    run_dirs = sorted(glob.glob(os.path.join(runs_dir, '*')))
    for run_dir in run_dirs:
        if not os.path.isdir(run_dir):
            continue
        dirname = os.path.basename(run_dir)
        if dirname == 'results':
            continue

        # Find CCDG file
        ccdg_files = glob.glob(os.path.join(run_dir, '*_global.ccdg'))
        if not ccdg_files:
            ccdg_files = glob.glob(os.path.join(run_dir, '*.ccdg'))
        if not ccdg_files:
            continue

        ccdg_file = ccdg_files[0]

        # Extract case name and rank count from dirname
        # Format: {case}_{rank}r_{timestamp}
        match = re.match(r'(.+)_(\d+)r_\d+', dirname)
        if match:
            case_name = match.group(1)
            ranks = int(match.group(2))
        else:
            case_name = dirname
            ranks = 0

        results_db.append({
            'case': case_name,
            'ranks': ranks,
            'dirname': dirname,
            'ccdg_file': ccdg_file,
            'run_dir': run_dir,
        })

    return results_db


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description='LAMMPS MPI Communication Feature Analysis')
    parser.add_argument('--runs-dir', default=None,
                        help='Path to runs directory (default: lammps_trace/runs)')
    parser.add_argument('--output-dir', default=None,
                        help='Output directory for analysis results')
    parser.add_argument('--single', default=None,
                        help='Process a single CCDG file')
    parser.add_argument('--name', default='experiment',
                        help='Experiment name (for single mode)')
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    runs_dir = args.runs_dir or os.path.join(base_dir, '..', 'runs')
    runs_dir = os.path.abspath(runs_dir)
    output_dir = args.output_dir or os.path.join(base_dir, 'results')

    if args.single:
        # Single experiment mode
        ccdg_file = args.single
        exp_output = os.path.join(output_dir, args.name)
        os.makedirs(exp_output, exist_ok=True)
        print(f"\n{'='*60}")
        print(f"Processing: {ccdg_file}")
        print(f"{'='*60}")
        results = process_single_experiment(ccdg_file, exp_output)
        if results:
            print(f"\n  Summary:")
            print(f"    Ranks:              {results['num_ranks']}")
            print(f"    Total comm ops:     {results['total_comm_ops']}")
            print(f"    Compute/Comm ratio: {results['compute_comm_ratio']:.2f}")
            print(f"    P2P ops:            {results['p2p_count']}")
            print(f"    Collective ops:     {results['coll_count']}")
        return

    # Batch mode: scan and process all runs
    print(f"\n{'='*60}")
    print(f"LAMMPS Communication Feature Analysis (Batch Mode)")
    print(f"  Runs directory: {runs_dir}")
    print(f"  Output directory: {output_dir}")
    print(f"{'='*60}")

    experiments = scan_runs(runs_dir)
    if not experiments:
        print("No CCDG files found in runs directory.")
        print("Run batch_trace_capture.sh first to generate traces.")
        return

    print(f"Found {len(experiments)} experiments")

    all_results = []
    for exp in experiments:
        print(f"\n{'='*60}")
        print(f"Case: {exp['case']} | Ranks: {exp['ranks']}")
        print(f"  Dir: {exp['dirname']}")
        print(f"{'='*60}")

        exp_output = os.path.join(output_dir, f"{exp['case']}_{exp['ranks']}r")
        results = process_single_experiment(exp['ccdg_file'], exp_output)
        if results:
            all_results.append({
                'case': exp['case'],
                'ranks': exp['ranks'],
                'results': results,
                'dirname': exp['dirname'],
            })
            print(f"\n  Summary:")
            print(f"    Total comm ops:     {results['total_comm_ops']}")
            print(f"    Compute/Comm ratio: {results['compute_comm_ratio']:.2f}")
            print(f"    P2P ops:            {results['p2p_count']}")
            print(f"    Collective ops:     {results['coll_count']}")

    # Comparative analysis
    if len(all_results) >= 4:
        print(f"\n{'='*60}")
        print("Generating comparative analysis...")
        print(f"{'='*60}")
        os.makedirs(output_dir, exist_ok=True)
        plot_comparative_scaling(all_results, output_dir)
        plot_rank_heatmap_grid(all_results, output_dir)
        plot_rank_heatmap_grid_bytes(all_results, output_dir)

    # Save master summary
    summary_path = os.path.join(output_dir, 'all_results.json')
    summary_data = []
    for entry in all_results:
        r = entry['results']
        summary_data.append({
            'case': entry['case'],
            'ranks': entry['ranks'],
            'num_ranks': r['num_ranks'],
            'total_comm_ops': r['total_comm_ops'],
            'total_compute_nodes': r['total_compute_nodes'],
            'total_comm_nodes': r['total_comm_nodes'],
            'compute_comm_ratio': r['compute_comm_ratio'],
            'total_compute_time_sec': r['total_compute_time_sec'],
            'p2p_count': r['p2p_count'],
            'coll_count': r['coll_count'],
            'collective_breakdown': r['collective_breakdown'],
            'type_distribution': r['type_distribution'],
        })

    with open(summary_path, 'w') as f:
        json.dump(summary_data, f, indent=2)
    print(f"\nSaved master summary: {summary_path}")

    # Print markdown summary table
    print(f"\n{'='*60}")
    print("Summary Table")
    print(f"{'='*60}")
    print(f"{'Case':<12} {'Ranks':<6} {'Comm Ops':<10} {'C/C Ratio':<10} "
          f"{'P2P':<8} {'Coll':<8} {'Compute(s)':<10}")
    print("-" * 64)
    for entry in all_results:
        r = entry['results']
        print(f"{entry['case']:<12} {r['num_ranks']:<6} {r['total_comm_ops']:<10} "
              f"{r['compute_comm_ratio']:<10.2f} {r['p2p_count']:<8} "
              f"{r['coll_count']:<8} {r['total_compute_time_sec']:<10.4f}")

    print(f"\nDone! Results in: {output_dir}")


if __name__ == '__main__':
    main()
