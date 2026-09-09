#!/usr/bin/env bash
# smoke_test.sh — 复现包一键冒烟验证
# 用法: bash smoke_test.sh [--skip-lammps-build] [--skip-pipeline]
# 步骤: ① 构建 sst-dumpi / dumpi2ccdg / BookSim（已有产物则跳过）
#       ② 样例 run 的 CCDG 直接注入 BookSim（不需要 MPI/LAMMPS）
#       ③ 完整流水线: LAMMPS 16r 捕获 → 五档编译/注入 → evaluation
set -uo pipefail
PKG="$(cd "$(dirname "$0")" && pwd)"
NPROC=$(nproc)
FAIL=0
log() { echo "[smoke $(date +%H:%M:%S)] $*"; }

SKIP_LMP_BUILD=0; SKIP_PIPELINE=0
for a in "$@"; do
  case $a in
    --skip-lammps-build) SKIP_LMP_BUILD=1 ;;
    --skip-pipeline) SKIP_PIPELINE=1 ;;
    *) echo "未知参数 $a"; exit 2 ;;
  esac
done

command -v g++ >/dev/null || { echo "ERROR: 缺少 g++"; exit 2; }
command -v python3 >/dev/null || { echo "ERROR: 缺少 python3"; exit 2; }
command -v mpirun >/dev/null || { echo "ERROR: 缺少 mpirun (OpenMPI)"; exit 2; }

# ── ① 构建 ────────────────────────────────────────────────────────────
if [ ! -x "$PKG/install/lib/libundumpi.so" ] && [ ! -f "$PKG/install/lib/libundumpi.so" ]; then
  log "构建 sst-dumpi → install/ ..."
  ( cd "$PKG/sst-dumpi" && ./configure --prefix="$PKG/install" --enable-libdumpi CC=mpicc CXX=mpic++ >/dev/null \
    && make -j"$NPROC" >/dev/null && make install >/dev/null ) || { echo "ERROR: sst-dumpi 构建失败"; exit 2; }
else
  log "sst-dumpi 已构建，跳过"
fi

if [ ! -x "$PKG/dumpi2ccdg/dumpi2ccdg" ]; then
  log "构建 dumpi2ccdg ..."
  ( cd "$PKG/dumpi2ccdg" && make ) || { echo "ERROR: dumpi2ccdg 构建失败"; exit 2; }
else
  log "dumpi2ccdg 已构建，跳过"
fi

if [ ! -x "$PKG/booksim2/src/booksim" ]; then
  log "构建 BookSim ..."
  ( cd "$PKG/booksim2/src" && make -j"$NPROC" ) || { echo "ERROR: BookSim 构建失败"; exit 2; }
else
  log "BookSim 已构建，跳过"
fi

LMP_NEEDS_BUILD=0
[ ! -x "$PKG/install/bin/lmp" ] && LMP_NEEDS_BUILD=1
[ "$PKG/lammps-src/src/comm_brick.cpp" -nt "$PKG/install/bin/lmp" ] && LMP_NEEDS_BUILD=1
[ "$PKG/lammps-src/src/comm_brick.h" -nt "$PKG/install/bin/lmp" ] && LMP_NEEDS_BUILD=1
if [ "$SKIP_LMP_BUILD" != 1 ] && [ "$LMP_NEEDS_BUILD" = 1 ]; then
  # PATH 里可能有 vendor 自带的坏 cmake；选第一个能真正运行的
  CMAKE_BIN=""
  for c in /usr/bin/cmake /usr/local/bin/cmake "$(command -v cmake)"; do
    [ -n "$c" ] && [ -x "$c" ] && "$c" --version >/dev/null 2>&1 && CMAKE_BIN="$c" && break
  done
  [ -n "$CMAKE_BIN" ] || { echo "ERROR: 找不到可用的 cmake"; exit 2; }
  log "构建 LAMMPS (KSPACE+MANYBODY, ~10 分钟, cmake=$CMAKE_BIN) ..."
  ( cd "$PKG/lammps-src" && mkdir -p build && cd build \
    && "$CMAKE_BIN" ../cmake -D PKG_KSPACE=yes -D PKG_MANYBODY=yes \
         -D BUILD_MPI=yes -D WITH_PNG=no -D WITH_JPEG=no -D WITH_FFMPEG=no -D WITH_GZIP=no \
         -D CMAKE_BUILD_TYPE=Release -D CMAKE_INSTALL_PREFIX="$PKG/install" >/dev/null \
    && make -j"$NPROC" >/dev/null && make install >/dev/null ) \
    || { echo "ERROR: LAMMPS 构建失败"; exit 2; }
else
  log "LAMMPS 已构建（或 --skip-lammps-build），跳过"
fi

# ── ② 样例 CCDG 直接注入（不需要 MPI）────────────────────────────────
SAMPLE="$PKG/runs/pipeline/short_lialocl_2688a_16r_20260831_195028/trimonly_16ranks_global.ccdg"
if [ -f "$SAMPLE" ]; then
  log "样例 CCDG 注入 BookSim (free, 1 步) ..."
  ( cd "$PKG/booksim2" && CCDG_COMPUTE_CAP=2.5e10 ./run_ccdg_mesh.sh "$SAMPLE" 600 1 \
      > /tmp/smoke_inject.log 2>&1 )
  CYCLES=$(grep -oP '完成: cycles=\K[0-9]+' /tmp/smoke_inject.log | tail -1)
  UNRES=$(grep -oP 'unresolved=\K[0-9]+' /tmp/smoke_inject.log | tail -1)
  if [ -n "$CYCLES" ] && [ "$UNRES" = "0" ]; then
    log "② PASS: cycles=$CYCLES unresolved=0 (参考 174,000±20%)"
  else
    log "② FAIL: cycles=$CYCLES unresolved=$UNRES"; FAIL=1
  fi
else
  log "② 跳过: 样例 CCDG 缺失"
fi

# ── ③ 完整流水线 ──────────────────────────────────────────────────────
if [ "$SKIP_PIPELINE" != 1 ]; then
  log "完整流水线: run_noc_pipeline.sh short 16 lialocl 2688 both (约 1 分钟) ..."
  cd "$PKG" && bash run_noc_pipeline.sh short 16 lialocl 2688 both > /tmp/smoke_pipeline.log 2>&1
  RC=$?
  OVERALL=$(grep -oP 'OVERALL: \K\w+' /tmp/smoke_pipeline.log | tail -1)
  EVAL=$(grep -oP '运行目录: \K\S+' /tmp/smoke_pipeline.log | tail -1)
  if [ "$OVERALL" = "PASS" ] && [ $RC -eq 0 ]; then
    log "③ PASS: $EVAL/evaluation.txt"
  else
    log "③ FAIL: rc=$RC OVERALL=$OVERALL (日志: /tmp/smoke_pipeline.log)"; FAIL=1
  fi
fi

[ $FAIL -eq 0 ] && log "SMOKE: ALL PASS ✓" || log "SMOKE: FAILED ✗"
exit $FAIL
