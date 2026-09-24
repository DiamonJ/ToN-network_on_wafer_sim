#!/usr/bin/env python3
"""
gen_marching_ccdg.py — 合成 marching multicast 风格的单步 CCDG，
用于验证 WSE wavelet 门控在"链路争用真实存在"场景下的收益。

复刻 SC24 论文 III-B 的争用前提（b > 1 时多个源共享同一 mesh 链路）:
  - KxK 方形 mesh，rank r 位于 (r%K, r//K)
  - 每个 rank 将自己的数据多播到 +x 方向的 b 个邻居（条带内东向推进）
  - 边界裁剪（东侧 rank 少发、西侧 rank 少收，对应论文 wafer 边缘条带）
  - 行内东向链路 (c→c+1) 被 x∈[c+1-b, c] 共 b 个源共享 → 真实链路争用

节点结构（每 rank，程序序即注入序）:
  BARRIER → COMPUTE → [SEND×b | IRECV×b] → [IRECV×b | SEND×b] → WAITALL → BARRIER

  - order=send_first（默认）: 先发后收，所有 rank 同时注入 → 同步波前，
    争用最大化（对应论文 b>1 的原始担忧）；依赖边挂到 WAITALL（IRECV 非阻塞）
  - order=recv_first: 先收后发，依赖边挂到各自 IRECV（阻塞收），流量自西向东
    逐跳串行化 → 模拟 wavelet 逐条带推进（marching）的依赖链

用法:
  python3 gen_marching_ccdg.py --k 4 --b 2 --bytes 16384 --compute 1000 \
      --order send_first -o runs/synthetic/mm_16r_b2_16k.ccdg
  python3 ccdg_unroll.py <单步.ccdg> 100 <输出.ccdg>
  cd booksim2 && ./run_ccdg_mesh.sh <展开文件> 3600 100
"""
import argparse
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=4, help="mesh 每边节点数（总 rank = k*k）")
    ap.add_argument("--b", type=int, default=2, help="多播宽度 b（>1 才有链路共享）")
    ap.add_argument("--bytes", type=int, default=16384, help="每条消息字节数")
    ap.add_argument("--compute", type=int, default=1000,
                    help="每 rank 每步计算周期（各 rank 相同，保持同步波前）")
    ap.add_argument("--order", choices=["send_first", "recv_first"], default="send_first")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args()

    K, B = args.k, args.b
    num_ranks = K * K
    nodes = []
    edges = []
    nid = 0

    def add(rank, ntype, **kw):
        nonlocal nid
        n = {"id": nid, "rank": rank, "type": ntype}
        n.update(kw)
        nodes.append(n)
        nid += 1
        return n["id"]

    # pass 1: 为每 rank 预生成结构（记录各节点的通信伙伴，便于 pass 2 连边）
    # rank r = (x, y)；东向多播: r -> y*K + (x+d), d=1..b（x+d < K 才发）
    # 西向接收:   来源 y*K + (x-d), d=1..b（x-d >= 0 才收）
    send_list = {r: [] for r in range(num_ranks)}   # (dst, tag)
    recv_list = {r: [] for r in range(num_ranks)}   # (src, tag)
    for r in range(num_ranks):
        x, y = r % K, r // K
        for d in range(1, B + 1):
            if x + d < K:
                send_list[r].append((y * K + (x + d), d))
    for r in range(num_ranks):
        x, y = r % K, r // K
        for d in range(1, B + 1):
            if x - d >= 0:
                recv_list[r].append((y * K + (x - d), d))

    # pass 2: 按 rank 分组生成节点（ccdg_unroll.py 要求 nodes 按 rank 分组有序）
    for r in range(num_ranks):
        pred = None

        def chain(node_id):
            nonlocal pred
            if pred is not None:
                for n in nodes:
                    if n["id"] == node_id:
                        n["predecessors"] = [pred]
                        break
            pred = node_id

        bar1 = add(r, "BARRIER")
        chain(bar1)
        comp = add(r, "COMPUTE", compute_cycles=args.compute)
        chain(comp)

        send_ids = {}    # (dst, tag) -> SEND node id
        recv_ids = {}    # (src, tag) -> IRECV node id

        def emit_sends():
            for dst, tag in send_list[r]:
                sid = add(r, "SEND", comm_src=r, comm_dst=dst, comm_tag=tag,
                          comm_bytes=args.bytes, comm_count=max(1, args.bytes // 8))
                chain(sid)
                send_ids[(dst, tag)] = sid

        def emit_recvs():
            for src, tag in recv_list[r]:
                rid = add(r, "IRECV", comm_src=src, comm_dst=r, comm_tag=tag,
                          comm_bytes=args.bytes, comm_count=max(1, args.bytes // 8))
                chain(rid)
                recv_ids[(src, tag)] = rid

        if args.order == "send_first":
            emit_sends()
            emit_recvs()
        else:
            emit_recvs()
            emit_sends()

        wa = add(r, "WAITALL")
        chain(wa)
        bar2 = add(r, "BARRIER")
        chain(bar2)

        # 依赖边: SEND(src rank) -> 目标 rank 的收侧节点
        #   send_first: 挂 WAITALL（IRECV 非阻塞，等待集中在 WAITALL）
        #   recv_first: 挂对应 IRECV（阻塞收，形成自西向东的 marching 依赖链）
        # 本 rank 的收侧目标节点 id 在 pass 3 统一连接（此时对端 WAITALL 尚未建好）
        nodes[-1]["_meta"] = {"send_ids": dict(send_ids),
                              "recv_ids": dict(recv_ids),
                              "waitall": wa}

    # pass 3: 连 cross_rank_edges
    meta = {n["id"]: n.pop("_meta") for n in nodes if "_meta" in n}
    for r in range(num_ranks):
        m = meta[  # 该 rank 的 meta 挂在最后一个 BARRIER 上
            next(n["id"] for n in reversed(nodes) if n["rank"] == r)
        ]
        for (dst, tag), sid in m["send_ids"].items():
            tm = meta[
                next(n["id"] for n in reversed(nodes) if n["rank"] == dst)
            ]
            if args.order == "send_first":
                target = tm["waitall"]
            else:
                target = tm["recv_ids"][(r, tag)]
            edges.append({"src_node": sid, "dst_node": target})

    out = {"num_ranks": num_ranks, "nodes": nodes, "cross_rank_edges": edges}
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(out, f)

    # 自检: SEND/IRECV 数量匹配（Σsend == Σrecv）、边两端类型合法
    ns = sum(len(send_list[r]) for r in range(num_ranks))
    nr = sum(len(recv_list[r]) for r in range(num_ranks))
    assert ns == nr, f"send/recv 数不匹配: {ns} vs {nr}"
    assert len(edges) == ns, f"边数 {len(edges)} != SEND 数 {ns}"
    type_by_id = {n["id"]: n["type"] for n in nodes}
    for e in edges:
        assert type_by_id[e["src_node"]] == "SEND"
        assert type_by_id[e["dst_node"]] in ("IRECV", "WAITALL")

    # 链路负载估算（XY 路由, 东向单跳）: 每条东向链路被 b 个源共享
    total_bytes = ns * args.bytes
    per_link = {}
    for r in range(num_ranks):
        for (dst, _) in send_list[r]:
            link = (r % K, r // K)  # 东向相邻链路 (x→x+1, y)
            per_link[link] = per_link.get(link, 0) + args.bytes
    hottest = max(per_link.values()) if per_link else 0

    print(f"marching multicast CCDG -> {args.output}")
    print(f"  mesh {K}x{K}, ranks {num_ranks}, b={B}, "
          f"msg {args.bytes}B, compute {args.compute} cyc, order={args.order}")
    print(f"  SENDs/rank avg {ns/num_ranks:.2f}, total msgs {ns}, "
          f"total bytes {total_bytes/1024:.0f} KB")
    print(f"  hottest eastbound link: {hottest/1024:.0f} KB "
          f"({hottest/args.bytes:.0f} msgs share it) "
          f"~ serial time {hottest} cyc @1 flit/cyc (flit_size=1B)")
    print(f"  compute dominates: {args.compute} cyc vs link serial {hottest} cyc "
          f"-> {'CONTENTION' if hottest > args.compute else 'compute-bound'} regime")


if __name__ == "__main__":
    main()
