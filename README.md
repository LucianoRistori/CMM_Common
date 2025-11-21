# common

Shared C++ utilities for all CMM-related analysis programs (FlatnessScan, CompareScan, AddDisplacedPoints, MakeHexGrid, PlotPoints, OptimizePath, etc.).  
This repository provides a stable, centralized implementation of the **Point** data structure and the **readPoints()** function used across all geometry and measurement tools.

---

## Contents

### Points.h / Points.cpp

**Points.h** defines the `Point` struct:
- `std::string label` — point identifier (e.g. C1, H7, P12)
- `double coords[3]` — X, Y, Z coordinates in millimeters

**Points.cpp** implements:
- `std::vector<Point> readPoints(const std::string& filename, int nExpected)`

`readPoints()`:
- Accepts both CSV and space-separated files  
- Skips empty lines  
- Reads optional labels  
- Reads the first **three numeric fields** (X, Y, Z)  
- Ignores extra fields  
- Warns on malformed lines  
- Warns if the number of points differs from `nExpected`  
- Returns a clean `std::vector<Point>`

---

## Purpose

This repository is the **single source of truth** for point-reading and point-representation code used by all downstream analysis tools.  
Centralizing this logic ensures:
- stable behavior across all programs  
- no duplicated code  
- consistent handling of labels and coordinates  
- easy maintenance  

---

## Usage

Other repositories include this library via:

```cpp
#include "../common/Points.h"
```

Link with:

```bash
../common/Points.o
```

Example:

```cpp
std::vector<Point> pts = readPoints("myfile.csv", 3);
```

---

## File Format

Each line should contain:

```
LABEL   X   Y   Z   [optional extra fields...]
```

Examples:

```
C14, 123.45, 67.89, -0.12
P7   100.0  -50.0   0.00  extra data
42   0   0   0
```

Commas or spaces accepted.  
Extra fields ignored.

---

## Versioning

Repository uses annotated tags:
- v1.0 — initial release  
- v1.2.x — updated readPoints()  
- v2.x — structural improvements and documentation  

---

## Development Notes

- Only the first 3 numeric fields are used (X, Y, Z).  
- Labels are optional and not validated.  
- Extra fields allow forward-compatible extensions.  
- No external dependencies: pure C++17, lightweight, and portable.  

---

## Future Enhancements

Optional future improvements:
- support for quoted CSV  
- support for comment lines (# or //)  
- configurable column positions  
- automatic label generation  
- storing additional metadata  

---

## License

Internal use — part of the CMM analysis ecosystem maintained by Luciano Ristori.
