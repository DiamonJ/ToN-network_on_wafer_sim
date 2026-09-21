#!/usr/bin/env python3
"""
ccdg_unroll.py — 把单步(迭代段)CCDG 按程序序展开 N 步，用于 BookSim 多步仿真。

用法: ccdg_unroll.py <输入.ccdg> <步数N> <输出.ccdg>

规则:
  - 节点复制 N 份，copy k 的节点 id 偏移 = k × 单步节点总数
  - 每 rank 按程序序把 copy k 的末节点链接为 copy k+1 首节点的 predecessor
  - cross_rank_edges 按相同偏移逐份复制（步 k 的 SEND 配步 k 的 WAIT，
    模式与单步一致）
  - 首 COMPUTE 处理: 新版 Loop-time 窗口裁剪（2026-08-21 起）下首 COMPUTE
    已是裁剪后的稳态计算量，展开时原样复制（默认）。仅对旧版裁剪（首
    COMPUTE = "setup 尾部→首迭代"一次性间隙）的 CCDG，设环境变量
    CCDG_UNROLL_ZERO_FIRST=1 将其置 0，避免 N 倍放大。
  - 输出前校验 id 唯一性与 predecessor 引用完整性
前提: 输入 CCDG 已经过 setup 裁剪（CCDG_TRIM_SETUP），仅含迭代段。
"""
import json
import os
import sys


def main():
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} <输入.ccdg> <步数N> <输出.ccdg>", file=sys.stderr)
        sys.exit(2)
    src, dst = sys.argv[1], sys.argv[3]
    try:
        n_steps = int(sys.argv[2])
    except ValueError:
        print("ERROR: 步数须为整数", file=sys.stderr)
        sys.exit(2)
    if n_steps < 1:
        print("ERROR: 步数须 >= 1", file=sys.stderr)
        sys.exit(2)

    with open(src) as f:
        d = json.load(f)
    nodes = d["nodes"]
    n0 = len(nodes)
    num_ranks = d["num_ranks"]
    edges = d.get("cross_rank_edges", [])

    # id -> 位置索引（校验 id 连续性，兼容任意起始编号）
    id2idx = {n["id"]: i for i, n in enumerate(nodes)}
    if len(id2idx) != n0:
        print("ERROR: 输入 CCDG 存在重复节点 id", file=sys.stderr)
        sys.exit(2)
    for n in nodes:
        for p in n.get("predecessors", []):
            if p not in id2idx:
                print(f"ERROR: predecessor {p} 引用不存在的节点", file=sys.stderr)
                sys.exit(2)
    for e in edges:
        if e["src_node"] not in id2idx or e["dst_node"] not in id2idx:
            print(f"ERROR: cross edge 引用不存在的节点: {e}", file=sys.stderr)
            sys.exit(2)

    # 每 rank 的首/末节点位置（nodes 按 rank 分组有序排列）
    rank_first, rank_last = {}, {}
    for i, n in enumerate(nodes):
        r = n["rank"]
        if r not in rank_first:
            rank_first[r] = i
        rank_last[r] = i
    if len(rank_first) != num_ranks:
        print(f"WARNING: 节点覆盖 {len(rank_first)} 个 rank, num_ranks={num_ranks}",
              file=sys.stderr)

    bytes_single = sum(n.get("comm_bytes", 0) for n in nodes)

    # 首 COMPUTE 清零（仅旧版裁剪 CCDG，环境变量 CCDG_UNROLL_ZERO_FIRST=1）:
    # 新版 Loop-time 窗口裁剪后首 COMPUTE 是稳态计算量，不应清零
    zero_first = bool(os.environ.get("CCDG_UNROLL_ZERO_FIRST"))
    rank_first_comp = {}
    if zero_first:
        for i, n in enumerate(nodes):
            if n["type"] == "COMPUTE" and n["rank"] not in rank_first_comp:
                rank_first_comp[n["rank"]] = i
    gap_cycles = sum(nodes[idx].get("compute_cycles", 0)
                     for idx in rank_first_comp.values())

    out_nodes = []
    for k in range(n_steps):
        off = k * n0
        for i, n in enumerate(nodes):
            m = dict(n)
            m["id"] = n["id"] + off
            if n.get("predecessors"):
                m["predecessors"] = [p + off for p in n["predecessors"]]
            else:
                m.pop("predecessors", None)
            if n["type"] == "COMPUTE" and n["rank"] in rank_first_comp \
                    and i == rank_first_comp[n["rank"]]:
                m["compute_cycles"] = 0
                m["compute_time_sec"] = 0
            out_nodes.append(m)
    # 相邻步链接: copy k 的每 rank 末节点 -> copy k+1 的首节点
    for k in range(n_steps - 1):
        for r in rank_first:
            first_idx = rank_first[r]
            last_id = nodes[rank_last[r]]["id"] + k * n0
            tgt = out_nodes[(k + 1) * n0 + first_idx]
            tgt.setdefault("predecessors", []).append(last_id)

    out_edges = []
    for k in range(n_steps):
        off = k * n0
        for e in edges:
            out_edges.append({"src_node": e["src_node"] + off,
                              "dst_node": e["dst_node"] + off})

    out = {"num_ranks": num_ranks, "nodes": out_nodes,
           "cross_rank_edges": out_edges}
    with open(dst, "w") as f:
        json.dump(out, f)

    # 自检
    ids = [n["id"] for n in out_nodes]
    assert len(set(ids)) == len(ids), "展开后出现重复 id"
    idset = set(ids)
    for n in out_nodes:
        for p in n.get("predecessors", []):
            assert p in idset, f"展开后 predecessor {p} 悬空"
    print(f"unroll: {src} x{n_steps} -> {dst}")
    print(f"  nodes: {n0} -> {len(out_nodes)}, cross edges: {len(edges)} -> {len(out_edges)}")
    print(f"  comm bytes: {bytes_single} -> {bytes_single * n_steps} (x{n_steps})")
    cyc_in = sum(n.get("compute_cycles", 0) for n in nodes if n["type"] == "COMPUTE")
    cyc_out = sum(n.get("compute_cycles", 0) for n in out_nodes if n["type"] == "COMPUTE")
    print(f"  compute cycles: {cyc_in:.0f} -> {cyc_out:.0f}"
          + (f" (首COMPUTE间隙 {gap_cycles:.0f} 已清零)" if zero_first
             else " (首COMPUTE 为稳态量, 原样复制)"))


if __name__ == "__main__":
    main()
