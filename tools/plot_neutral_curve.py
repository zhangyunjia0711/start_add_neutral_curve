#!/usr/bin/env python3
"""Plot START neutral curves and optional growth fields from Tecplot text data."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
import re
import sys
import warnings

import numpy as np

COLUMNS = (
    "outer", "inner", "iloc", "x", "Re", "alpha_r", "alpha_i",
    "beta_r", "beta_i", "omega_r", "omega_i", "growth", "branch", "status",
)
INDEX = {name.casefold(): index for index, name in enumerate(COLUMNS)}


@dataclass
class TecplotData:
    values: np.ndarray
    row_numbers: np.ndarray
    sample_numbers: np.ndarray
    zones: tuple[str, ...]

    def column(self, name: str) -> np.ndarray:
        return self.values[:, INDEX[name.casefold()]]


def read_tecplot(path: str | Path) -> TecplotData:
    """Read POINT-format tables, retaining empty zones and original scan indices.

    A row_N zone identifies the original scan row. Failed samples are retained:
    removing them before triangulation could draw a false field across a failure.
    Empty files and header-only files return a correctly shaped empty table.
    """
    path = Path(path)
    records, row_numbers, sample_numbers, zones = [], [], [], []
    variable_lines = []
    in_variables = False
    current_row, sample = 1, 0
    saw_zone = False
    for line_number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "!")):
            continue
        upper = line.upper()
        if upper.startswith("VARIABLES"):
            in_variables = True
            variable_lines.append(line.split("=", 1)[-1])
            continue
        if upper.startswith("ZONE"):
            in_variables = False
            title_match = re.search(r'\bT\s*=\s*"([^"]*)"', line, re.I)
            title = title_match.group(1) if title_match else f"zone_{len(zones) + 1}"
            row_match = re.fullmatch(r"row_(\d+)", title, re.I)
            current_row = int(row_match.group(1)) if row_match else len(zones) + 1
            sample = 0
            zones.append(title)
            saw_zone = True
            if re.search(r"\b(?:DATAPACKING\s*=\s*BLOCK|F\s*=\s*BLOCK)\b", line, re.I):
                raise ValueError(f"{path}:{line_number}: BLOCK data are unsupported; use POINT data")
            continue
        if upper.startswith(("TITLE", "TEXT", "AUXDATA", "DATASETAUXDATA")):
            continue
        if in_variables and not re.match(r"^[+\-.\d]", line):
            variable_lines.append(line)
            continue
        in_variables = False
        tokens = line.replace(",", " ").split()
        if len(tokens) != len(COLUMNS):
            raise ValueError(
                f"{path}:{line_number}: expected {len(COLUMNS)} columns, found {len(tokens)}"
            )
        try:
            record = [float(token.replace("D", "E").replace("d", "e")) for token in tokens]
        except ValueError as error:
            raise ValueError(f"{path}:{line_number}: invalid numeric record") from error
        if not saw_zone:
            zones.append("row_1")
            saw_zone = True
        records.append(record)
        row_numbers.append(current_row)
        sample_numbers.append(sample)
        sample += 1
    if variable_lines:
        header = " ".join(variable_lines)
        names = re.findall(r'"([^"]+)"', header)
        if not names:
            names = header.replace(",", " ").split()
        if [name.casefold() for name in names] != [name.casefold() for name in COLUMNS]:
            raise ValueError(
                f"{path}: VARIABLES must be in this order: {', '.join(COLUMNS)}"
            )
    return TecplotData(
        np.asarray(records, dtype=float).reshape((-1, len(COLUMNS))),
        np.asarray(row_numbers, dtype=int),
        np.asarray(sample_numbers, dtype=int),
        tuple(zones),
    )


def successful(data: TecplotData) -> np.ndarray:
    """Only status=0 samples with finite coordinates and growth may be drawn."""
    return (
        (data.column("status") == 0)
        & np.isfinite(data.column("outer"))
        & np.isfinite(data.column("inner"))
        & np.isfinite(data.column("growth"))
    )


def root_segments(data: TecplotData) -> list[np.ndarray]:
    """Split by branch, failed points, and missing original rows; preserve order."""
    if not len(data.values):
        return []
    good = successful(data)
    branches = data.column("branch")
    segments = []
    for branch in np.unique(branches[np.isfinite(branches)]):
        current, previous_row = [], None
        for index in np.flatnonzero(branches == branch):
            row = data.row_numbers[index]
            contiguous = previous_row is None or 0 <= row - previous_row <= 1
            if not good[index] or not contiguous:
                if current:
                    segments.append(np.asarray(current, dtype=int))
                current = []
            if good[index]:
                current.append(index)
            previous_row = row
        if current:
            segments.append(np.asarray(current, dtype=int))
    return segments


def rectilinear_field(data: TecplotData):
    """Return a masked rectilinear grid, or None for an irregular sample cloud."""
    coordinates = np.isfinite(data.column("outer")) & np.isfinite(data.column("inner"))
    indices = np.flatnonzero(coordinates)
    if not len(indices):
        return None
    xs = np.unique(data.column("outer")[indices])
    ys = np.unique(data.column("inner")[indices])
    if len(xs) < 2 or len(ys) < 2:
        return None
    # A field zone must have constant outer coordinate.
    for row in np.unique(data.row_numbers[indices]):
        same_row = indices[data.row_numbers[indices] == row]
        if len(np.unique(data.column("outer")[same_row])) != 1:
            return None
        # Local root corrections add inner coordinates that differ by row.
        # Such clouds must not become a sparse rectilinear grid merely because
        # repeated evaluations increase the number of input records.
        if not np.array_equal(np.unique(data.column("inner")[same_row]), ys):
            return None
    grid = np.full((len(ys), len(xs)), np.nan)
    rows_for_x = {}
    good = successful(data)
    unknown_rows = np.unique(data.row_numbers[~coordinates & ~good])
    seen = {}
    for index in indices:
        x_index = int(np.searchsorted(xs, data.column("outer")[index]))
        y_index = int(np.searchsorted(ys, data.column("inner")[index]))
        key = (y_index, x_index)
        rows_for_x[x_index] = data.row_numbers[index]
        # Duplicates containing a failure remain masked.
        previous_good = seen.get(key, True)
        seen[key] = previous_good and bool(good[index]) and data.row_numbers[index] not in unknown_rows
        grid[key] = data.column("growth")[index] if seen[key] else np.nan
    # An empty row_N zone still represents a gap. Insert a masked column so
    # contourf cannot connect rows on either side of it.
    expanded_x, expanded_grid = [], []
    typical_x = typical_spacing(np.diff(xs))
    for x_index, value in enumerate(xs):
        if x_index and (abs(rows_for_x[x_index] - rows_for_x[x_index - 1]) > 1
                        or value - xs[x_index - 1] > 2.5 * typical_x):
            expanded_x.append(0.5 * (xs[x_index - 1] + value))
            expanded_grid.append(np.full(len(ys), np.nan))
        expanded_x.append(value)
        expanded_grid.append(grid[:, x_index])
    expanded_grid = np.asarray(expanded_grid).T
    expanded_y, final_grid = [], []
    typical_y = typical_spacing(np.diff(ys))
    for y_index, value in enumerate(ys):
        if y_index and value - ys[y_index - 1] > 2.5 * typical_y:
            expanded_y.append(0.5 * (ys[y_index - 1] + value))
            final_grid.append(np.full(len(expanded_x), np.nan))
        expanded_y.append(value)
        final_grid.append(expanded_grid[y_index])
    return np.asarray(expanded_x), np.asarray(expanded_y), np.ma.masked_invalid(final_grid)


def typical_spacing(steps: np.ndarray) -> float:
    """Estimate the scan spacing while limiting bias from local root refinement."""
    steps = np.asarray(steps, dtype=float)
    steps = steps[np.isfinite(steps) & (steps > 0)]
    if not len(steps):
        return np.inf
    # Root corrections add many geometrically shrinking intervals. The densest
    # logarithmic spacing band usually retains the original scan resolution.
    bands = np.floor(np.log10(steps) / 0.25).astype(int)
    unique, counts = np.unique(bands, return_counts=True)
    chosen = unique[np.flatnonzero(counts == np.max(counts))[-1]]
    return float(np.median(steps[bands == chosen]))


def field_triangulation(data: TecplotData):
    """Triangulate irregular scans without skipping failed/missing scan samples."""
    import matplotlib.tri as mtri

    finite = np.isfinite(data.column("outer")) & np.isfinite(data.column("inner"))
    selected = np.flatnonzero(finite)
    if len(selected) < 3:
        return None
    coordinates = np.column_stack((data.column("outer")[selected], data.column("inner")[selected]))
    unique, first, inverse = np.unique(coordinates, axis=0, return_index=True, return_inverse=True)
    if len(unique) < 3:
        return None
    original = selected[first]
    all_good = successful(data)
    good = all_good[original].copy()
    for offset, unique_index in enumerate(inverse):
        good[unique_index] &= all_good[selected[offset]]
    # A failure without finite coordinates cannot be located inside its row.
    # Conservatively hide that row instead of interpolating through it.
    unknown_rows = np.unique(data.row_numbers[~finite & ~all_good])
    good[np.isin(data.row_numbers[original], unknown_rows)] = False
    try:
        mesh = mtri.Triangulation(unique[:, 0], unique[:, 1])
    except (RuntimeError, ValueError):
        return None
    triangles = mesh.triangles
    masked = ~np.all(good[triangles], axis=1)
    row = data.row_numbers[original]
    rank = np.zeros(len(unique), dtype=int)
    x_steps = np.diff(np.unique(unique[:, 0]))
    typical_x = typical_spacing(x_steps)
    y_steps = []
    for row_number in np.unique(row):
        row_indices = np.flatnonzero(row == row_number)
        ordered = row_indices[np.argsort(unique[row_indices, 1])]
        rank[ordered] = np.arange(len(ordered))
        differences = np.diff(unique[ordered, 1])
        y_steps.extend(differences[differences > 0])
    typical_y = typical_spacing(np.asarray(y_steps))
    for left, right in ((0, 1), (1, 2), (2, 0)):
        a, b = triangles[:, left], triangles[:, right]
        masked |= np.abs(row[a] - row[b]) > 1
        same_row = row[a] == row[b]
        masked |= same_row & (np.abs(rank[a] - rank[b]) > 1)
        masked |= np.abs(unique[a, 0] - unique[b, 0]) > 2.5 * typical_x
        masked |= np.abs(unique[a, 1] - unique[b, 1]) > 2.5 * typical_y
    mesh.set_mask(masked)
    if np.all(masked):
        return None
    growth = data.column("growth")[original].copy()
    growth[~np.isfinite(growth)] = 0.0  # all associated triangles are masked
    return mesh, growth


def plot_data(roots: TecplotData, field: TecplotData | None = None, *,
              x_label: str = "outer", y_label: str = "inner"):
    import matplotlib.pyplot as plt

    if not np.any(successful(roots)) and (field is None or not np.any(successful(field))):
        raise ValueError("No successful neutral roots or field samples to plot")
    figure, axes = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    if field is not None and np.any(successful(field)):
        values = field.column("growth")[successful(field)]
        bound = max(float(np.max(np.abs(values))), np.finfo(float).eps)
        levels = np.linspace(-bound, bound, 25)
        image = None
        grid = rectilinear_field(field)
        if grid is not None and min(grid[2].shape) >= 2:
            # corner_mask=False masks each whole cell touching a failed point.
            image = axes.contourf(*grid, levels=levels, cmap="RdBu_r", corner_mask=False)
        else:
            triangulated = field_triangulation(field)
            if triangulated is not None:
                image = axes.tricontourf(*triangulated, levels=levels, cmap="RdBu_r")
        if image is not None:
            figure.colorbar(image, ax=axes, label="growth")
        else:
            warnings.warn("Insufficient adjacent successful samples for a field background")
    branches_seen = set()
    for indices in root_segments(roots):
        branch = roots.column("branch")[indices[0]]
        label = f"branch {branch:g}" if branch not in branches_seen else "_nolegend_"
        axes.plot(roots.column("outer")[indices], roots.column("inner")[indices],
                  ".-", color="black", linewidth=1.3, markersize=4, label=label)
        branches_seen.add(branch)
    if branches_seen:
        axes.legend(loc="best", fontsize="small")
    else:
        axes.text(0.02, 0.98, "No successful neutral roots", transform=axes.transAxes,
                  ha="left", va="top", bbox={"facecolor": "white", "alpha": 0.8})
    axes.set(xlabel=x_label, ylabel=y_label, title="Neutral curve")
    axes.grid(alpha=0.2)
    return figure, axes


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", type=Path, help="prefix_roots.plt")
    parser.add_argument("--field", type=Path, help="optional prefix_field.plt growth background")
    parser.add_argument("--output", type=Path, help="save a figure instead of opening a window")
    parser.add_argument("--x-label", default="outer")
    parser.add_argument("--y-label", default="inner")
    parser.add_argument("--dpi", type=int, default=180)
    arguments = parser.parse_args(argv)
    if arguments.output:
        import matplotlib
        matplotlib.use("Agg")
    try:
        roots = read_tecplot(arguments.roots)
        field = read_tecplot(arguments.field) if arguments.field else None
        figure, _ = plot_data(roots, field, x_label=arguments.x_label, y_label=arguments.y_label)
        if arguments.output:
            figure.savefig(arguments.output, dpi=arguments.dpi)
            print(f"Saved {arguments.output}")
        else:
            import matplotlib.pyplot as plt
            plt.show()
        return 0
    except (OSError, ValueError) as error:
        print(f"plot_neutral_curve: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
