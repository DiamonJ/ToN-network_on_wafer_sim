#!/usr/bin/env python3
"""生成 rank 间 2D mesh 的 SimGrid platform XML。

rank r -> host host_r, 坐标 (r % KX, r // KX), 水平/垂直双向链路,
全显式路由（routing="Full"）。host speed 1.42 Gf, 链路 6.8 GB/s / 10 ns。

用法: python3 gen_mesh_platform.py KX KY [outdir]
"""
import sys

BW = "6800MBps"        # 6.8 GB/s
LAT = "0.00000001s"    # 10 ns
HOST_SPEED = "1.42Gf"


def gen(kx, ky):
    n = kx * ky
    lines = []
    lines.append("<?xml version='1.0'?>")
    lines.append("<!DOCTYPE platform SYSTEM \"https://simgrid.org/simgrid.dtd\">")
    lines.append("<platform version=\"4.1\">")
    lines.append(f"  <zone id=\"mesh{kx}x{ky}\" routing=\"Full\">")
    for r in range(n):
        lines.append(f"    <host id=\"host{r}\" speed=\"{HOST_SPEED}\"/>")
    # loopback
    lines.append(f"    <link id=\"loopback\" bandwidth=\"100000MBps\" latency=\"{LAT}\"/>")
    # mesh links: 水平 (x 方向) 与垂直 (y 方向)
    for y in range(ky):
        for x in range(kx):
            r = y * kx + x
            if x + 1 < kx:
                lines.append(f"    <link id=\"h{r}_{r+1}\" bandwidth=\"{BW}\" latency=\"{LAT}\"/>")
            if y + 1 < ky:
                rb = r + kx
                lines.append(f"    <link id=\"v{r}_{rb}\" bandwidth=\"{BW}\" latency=\"{LAT}\"/>")
    # 路由: BFS 最短路（先 x 后 y）；FullZone 路由默认对称，只需 a<=b
    def path(a, b):
        if a == b:
            return ["loopback"]
        ax, ay = a % kx, a // kx
        bx, by = b % kx, b // kx
        links = []
        x, y = ax, ay
        while x != bx:
            r = y * kx + x
            if x < bx:
                links.append(f"h{r}_{r+1}")
                x += 1
            else:
                x -= 1
                r = y * kx + x
                links.append(f"h{r}_{r+1}")
        while y != by:
            r = y * kx + x
            if y < by:
                links.append(f"v{r}_{r+kx}")
                y += 1
            else:
                y -= 1
                r = y * kx + x
                links.append(f"v{r}_{r+kx}")
        return links
    for a in range(n):
        for b in range(a, n):
            ctns = "".join(f"<link_ctn id=\"{l}\"/>" for l in path(a, b))
            lines.append(f"    <route src=\"host{a}\" dst=\"host{b}\">{ctns}</route>")
    lines.append("  </zone>")
    lines.append("</platform>")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    kx, ky = int(sys.argv[1]), int(sys.argv[2])
    outdir = sys.argv[3] if len(sys.argv) > 3 else "."
    fn = f"{outdir}/platform_mesh_{kx}x{ky}.xml"
    with open(fn, "w") as f:
        f.write(gen(kx, ky))
    print(f"wrote {fn} ({kx}x{ky} = {kx*ky} hosts)")
