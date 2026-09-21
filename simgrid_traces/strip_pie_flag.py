#!/usr/bin/env python3
"""清除 ELF 的 DF_1_PIE 标志，使 PIE 可执行文件可被 glibc dlopen（SMPI 需要）。"""
import struct, sys

path = sys.argv[1]
data = bytearray(open(path, "rb").read())
e_phoff   = struct.unpack_from("<Q", data, 0x20)[0]
e_phentsz = struct.unpack_from("<H", data, 0x36)[0]
e_phnum   = struct.unpack_from("<H", data, 0x38)[0]
dyn_off = dyn_sz = None
for i in range(e_phnum):
    off = e_phoff + i * e_phentsz
    if struct.unpack_from("<I", data, off)[0] == 2:  # PT_DYNAMIC
        dyn_off = struct.unpack_from("<Q", data, off + 0x08)[0]
        dyn_sz  = struct.unpack_from("<Q", data, off + 0x20)[0]
        break
assert dyn_off, "no PT_DYNAMIC"
DT_FLAGS_1, DF_1_PIE = 0x6FFFFFFB, 0x08000000
for i in range(dyn_sz // 16):
    ent = dyn_off + i * 16
    tag = struct.unpack_from("<q", data, ent)[0]
    if tag == 0:
        break
    if tag == DT_FLAGS_1:
        val = struct.unpack_from("<Q", data, ent + 8)[0]
        struct.pack_into("<Q", data, ent + 8, val & ~DF_1_PIE)
        open(path, "wb").write(data)
        print(f"{path}: FLAGS_1 {val:#x} -> {val & ~DF_1_PIE:#x}")
        sys.exit(0)
sys.exit("no DT_FLAGS_1")
