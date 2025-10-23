# CMM_Common

Shared **common library** used by the [FlatnessScan](https://github.com/LucianoRistori/FlatnessScan) and [CompareScan](https://github.com/LucianoRistori/CompareScan) projects.

This repository contains reusable C++ classes and utilities for coordinate-measurement-machine (CMM) data analysis, grid detection, and point-cloud handling.  
It is designed to promote **code reuse**, **consistency**, and **simplified maintenance** across the CMM toolchain.

---

## 📦 Contents

| Folder / File | Description |
|----------------|-------------|
| `Points.h`, `Points.cpp` | Core data structure for 3D measurement points (`X, Y, Z`, residuals, etc.). |
| `GridCheckResult.h`, `GridCheckResult.cpp` | Struct and methods for grid-finding results (missing points, grid size, success flag, etc.). |
| *(future)* Additional shared headers for plane fits, utilities, or geometry helpers. |

---

## 🧰 Build Integration

Both **FlatnessScan** and **CompareScan** include the `common` folder as a sub-directory in their Makefiles.  
Typical compilation line (macOS + ROOT example):

```bash
clang++ -std=c++17 -O2 -Wall \
    ../common/Points.cpp ../common/GridCheckResult.cpp \
    FlatnessScan.cpp -o FlatnessScan \
    -I../common $(root-config --cflags --libs)
