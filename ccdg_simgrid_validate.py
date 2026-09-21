#!/usr/bin/env python3
"""
Single-run CCDG validation driver: SimGrid DAG 仿真 vs 真实 LAMMPS 耗时.
用法: ccdg_simgrid_validate.py <run_dir> <num_ranks> [ccdg文件]
  不指定 ccdg文件时自动取 run_dir 下按名排序的第一个 *_global.ccdg
环境变量:
  SIMGRID_CPU_FREQ=2.80  → 固定频率模式(跨 rank 对比研究用，避免通信开销被校准吸收)
  未设置                → 自动迭代校准频率(单次验证用，验证 CCDG 结构正确性)
退出码: 0=PASS(误差≤5%)  1=FAIL  2=流水线错误
复用 validate_simgrid.py 的 DAG 生成/编译/运行/校准函数(全互联 backbone 平台)。
"""
import json
import os
import re
import sys

import validate_simgrid as vs


def die(msg):
    print(f"ERROR: {msg}")
    sys.exit(2)


def main():
    if len(sys.argv) not in (3, 4):
        print("usage: ccdg_simgrid_validate.py <run_dir> <num_ranks> [ccdg文件]")
        sys.exit(2)
    run_dir, num_ranks = sys.argv[1], int(sys.argv[2])
    env_freq = os.environ.get("SIMGRID_CPU_FREQ")

    # T_real
    t_real = vs.get_t_real(os.path.join(run_dir, "lammps.log"))
    if t_real is None:
        t_real = vs.get_t_real(os.path.join(run_dir, "log.lammps"))
    if t_real is None:
        die("cannot find 'Loop time' (T_real) in log")
    print(f"T_real     = {t_real:.6f} s")

    # CCDG（可显式指定，默认取 run_dir 下排序第一个 *_global.ccdg）
    if len(sys.argv) >= 4:
        ccdg_path = sys.argv[3]
        if not os.path.isfile(ccdg_path):
            die(f"ccdg file not found: {ccdg_path}")
        print(f"CCDG       = {ccdg_path}")
    else:
        ccdgs = [f for f in sorted(os.listdir(run_dir)) if f.endswith("_global.ccdg")]
        if not ccdgs:
            die(f"no *_global.ccdg in {run_dir}")
        ccdg_path = os.path.join(run_dir, ccdgs[0])
        print(f"CCDG       = {ccdgs[0]}")

    cpu_freq = float(env_freq) if env_freq else None
    calibrated = cpu_freq is None

    if calibrated:
        # 自动迭代校准频率（单次验证口径）
        cpu_freq, iter_err, t_sim_iter, iter_stats = \
            vs.compute_cpu_freq_iterative(run_dir, num_ranks, t_real, ccdg_path)
        print(f"CPU freq   = {cpu_freq:.4f} GHz (iterative-calibrated, iter_error={iter_err:.4f}%)")

        # 物理上限检查：校准频率不得超过机器物理 CPU 频率上限，
        # 超限说明 CCDG 计算量建模超出真实迭代时长（如 setup 残留未裁净）。
        freq_cap = float(os.environ.get("SIMGRID_FREQ_CAP_GHZ", "4.0"))
        if cpu_freq > freq_cap:
            print(f"E_freq     = FAIL ✗ 校准频率 {cpu_freq:.2f} GHz 超过物理上限 "
                  f"{freq_cap:.1f} GHz（CCDG 计算量建模异常）")
            sys.exit(1)

        # 用校准频率生成正式产物 + 最终结果
        dag_dir = os.path.join(run_dir, "simgrid_dag")
        os.makedirs(dag_dir, exist_ok=True)
        platform_xml = os.path.join(dag_dir, "platform.xml")
        vs.generate_platform_xml(num_ranks, cpu_freq, platform_xml)
        cpp_path = os.path.join(dag_dir, "ccdg_dag_sim.cpp")
        stats = vs.ccdg_to_dag_cpp_fixed(ccdg_path, cpp_path, cpu_freq)
        print(f"Activities = {stats['num_acts']} (Exec {stats['exec_count']}, Comm {stats['comm_count']})")
        bin_path = os.path.join(dag_dir, "ccdg_dag_sim")
        if not vs.compile_simgrid_dag(cpp_path, bin_path):
            die("SimGrid DAG 编译失败")
        res = vs.run_simgrid_dag(bin_path, platform_xml)
        if res.returncode != 0:
            die(f"SimGrid DAG 运行失败 rc={res.returncode}: {res.stderr[:300]}")
    else:
        # 固定频率模式
        print(f"CPU freq   = {cpu_freq} GHz (fixed, SIMGRID_CPU_FREQ)")
        dag_dir = os.path.join(run_dir, "simgrid_dag")
        os.makedirs(dag_dir, exist_ok=True)
        platform_xml = os.path.join(dag_dir, "platform.xml")
        vs.generate_platform_xml(num_ranks, cpu_freq, platform_xml)
        cpp_path = os.path.join(dag_dir, "ccdg_dag_sim.cpp")
        stats = vs.ccdg_to_dag_cpp_fixed(ccdg_path, cpp_path, cpu_freq)
        print(f"Activities = {stats['num_acts']} (Exec {stats['exec_count']}, Comm {stats['comm_count']})")
        bin_path = os.path.join(dag_dir, "ccdg_dag_sim")
        if not vs.compile_simgrid_dag(cpp_path, bin_path):
            die("SimGrid DAG 编译失败")
        res = vs.run_simgrid_dag(bin_path, platform_xml)
        if res.returncode != 0:
            die(f"SimGrid DAG 运行失败 rc={res.returncode}: {res.stderr[:300]}")

    t_sim = None
    for line in res.stdout.strip().splitlines():
        if "T_simgrid_sec" in line:
            try:
                t_sim = json.loads(line)["T_simgrid_sec"]
            except json.JSONDecodeError:
                m = re.search(r"([\d.]+)", line)
                if m:
                    t_sim = float(m.group(1))
            break
    if t_sim is None:
        die(f"cannot parse T_simgrid from output: {res.stdout[:200]}")

    e_total = abs(t_sim - t_real) / t_real * 100.0
    passed = e_total <= 5.0
    print(f"T_simgrid  = {t_sim:.6f} s")
    print(f"E_total    = {e_total:.4f}%  {'PASS ✓ (≤5%)' if passed else 'FAIL ✗ (>5%)'}")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()