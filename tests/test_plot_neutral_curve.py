#!/usr/bin/env python3
"""Regression tests for parsing and safe drawing of neutral-curve output."""
import contextlib
import io
from pathlib import Path
import tempfile
import unittest

import matplotlib
matplotlib.use("Agg")
import numpy as np

from plot_neutral_curve import (
    COLUMNS, TecplotData, field_triangulation, main, read_tecplot,
    rectilinear_field, root_segments, successful,
)


def record(outer, inner, *, growth=0.0, branch=1, status=0):
    return [outer, inner, 270, 0.3, 79100, 3.5, 0, -3.2059, 0, -0.03933479,
            growth, growth, branch, status]


def write_table(path, zones):
    header = 'TITLE = "neutral test"\nVARIABLES = "outer" "inner" "iloc" "x" "Re"\n'
    header += ' "alpha_r" "alpha_i" "beta_r" "beta_i" "omega_r" "omega_i"\n'
    header += ' "growth" "branch" "status"\n'
    lines = [header]
    for row_number, samples in zones:
        lines.append(f'ZONE T="row_{row_number}", DATAPACKING=POINT\n')
        for sample in samples:
            lines.append(" ".join(f"{value:.10E}".replace("E", "D") for value in sample) + "\n")
    path.write_text("".join(lines), encoding="utf-8")


class PlotTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def test_multiline_variables_fortran_exponent_and_branch_gaps(self):
        path = self.directory / "roots.plt"
        write_table(path, [
            (1, [record(1, 0.1)]),
            (2, [record(2, 0.2), record(2, 0.8, branch=2)]),
            (3, []),
            (4, [record(4, 0.4), record(4, 0.6, branch=2)]),
            (5, [record(5, 0.5, status=4)]),
            (6, [record(6, 0.6)]),
        ])
        data = read_tecplot(path)
        self.assertEqual(data.values.shape, (7, 14))
        self.assertEqual(data.zones[2], "row_3")
        np.testing.assert_allclose(data.column("Re"), 79100)
        self.assertEqual(sorted(map(len, root_segments(data))), [1, 1, 1, 1, 2])
        self.assertEqual(int(np.sum(successful(data))), 6)
        for segment in root_segments(data):
            self.assertTrue(np.all(np.diff(data.row_numbers[segment]) <= 1))
            self.assertTrue(np.all(successful(data)[segment]))

    def test_empty_and_header_only_files_have_clear_cli_failure(self):
        for content in ("", 'VARIABLES = ' + " ".join(f'"{name}"' for name in COLUMNS)):
            path = self.directory / "empty.plt"
            path.write_text(content, encoding="utf-8")
            data = read_tecplot(path)
            self.assertEqual(data.values.shape, (0, 14))
            with contextlib.redirect_stderr(io.StringIO()) as error:
                code = main([str(path), "--output", str(self.directory / "empty.png")])
            self.assertEqual(code, 2)
            self.assertIn("No successful neutral roots", error.getvalue())
            self.assertFalse((self.directory / "empty.png").exists())

    def test_bad_variables_and_incomplete_record_are_rejected(self):
        path = self.directory / "bad.plt"
        path.write_text('VARIABLES = "inner" "outer"\n', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "VARIABLES must be"):
            read_tecplot(path)
        path.write_text('ZONE T="row_1"\n1 2 3\n', encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "expected 14 columns"):
            read_tecplot(path)

    def test_regular_field_masks_failed_points_and_empty_rows(self):
        path = self.directory / "field.plt"
        zones = []
        for row in (1, 2, 4, 5):
            samples = [record(row, inner, growth=inner - 2, branch=0,
                              status=5 if (row, inner) == (2, 2) else 0)
                       for inner in range(5)]
            zones.append((row, samples))
        zones.insert(2, (3, []))
        write_table(path, zones)
        data = read_tecplot(path)
        xs, ys, growth = rectilinear_field(data)
        np.testing.assert_allclose(xs, [1, 2, 3, 4, 5])
        self.assertTrue(np.all(np.ma.getmaskarray(growth)[:, 2]))
        self.assertTrue(np.ma.getmaskarray(growth)[2, 1])
        self.assertEqual(int(growth.count()), 19)

    def test_irregular_field_never_skips_failed_sample_or_missing_row(self):
        path = self.directory / "irregular.plt"
        zones = []
        for row in range(1, 5):
            samples = [record(row, inner + row * 0.01, growth=inner - 2, branch=0,
                              status=2 if (row, inner) == (2, 2) else 0)
                       for inner in range(5)]
            zones.append((row, samples))
        write_table(path, zones)
        data = read_tecplot(path)
        self.assertIsNone(rectilinear_field(data))
        mesh, growth = field_triangulation(data)
        visible = mesh.triangles[~mesh.mask]
        self.assertGreater(len(visible), 0)
        for triangle in visible:
            vertices = np.column_stack((mesh.x[triangle], mesh.y[triangle]))
            self.assertFalse(any(np.allclose(vertex, [2, 2.02]) for vertex in vertices))
            self.assertLessEqual(np.ptp(mesh.x[triangle]), 1)
        write_table(path, [zones[0], (2, []), zones[2]])
        self.assertIsNone(field_triangulation(read_tecplot(path)))

    def test_irregular_field_masks_excessive_coordinate_gaps(self):
        path = self.directory / "large_gap.plt"
        zones = []
        for row, outer in enumerate((0, 1, 100, 101), 1):
            samples = [record(outer, inner + row * 0.01, growth=inner - 1, branch=0)
                       for inner in range(4)]
            zones.append((row, samples))
        write_table(path, zones)
        mesh, _ = field_triangulation(read_tecplot(path))
        visible = mesh.triangles[~mesh.mask]
        self.assertGreater(len(visible), 0)
        self.assertTrue(np.all(np.ptp(mesh.x[visible], axis=1) <= 1))

    def test_refinement_order_and_duplicate_failure_do_not_make_false_field(self):
        path = self.directory / "refinement.plt"
        zones = []
        for row in range(1, 4):
            samples = [record(row, inner, growth=inner - 2, branch=0)
                       for inner in range(5)]
            samples.extend(record(row, inner, growth=inner - 2, branch=0)
                           for inner in (2.01, 1.995, 2.005, 2.0))
            if row == 2:
                samples.append(record(row, 2.0, branch=0, status=2))
            zones.append((row, samples))
        write_table(path, zones)
        mesh, _ = field_triangulation(read_tecplot(path))
        visible = mesh.triangles[~mesh.mask]
        self.assertGreater(len(visible), 0)
        for triangle in visible:
            vertices = np.column_stack((mesh.x[triangle], mesh.y[triangle]))
            self.assertFalse(any(np.allclose(vertex, [2, 2]) for vertex in vertices))

    def test_arclength_zone_preserves_turning_point_sequence(self):
        path = self.directory / "arc.plt"
        samples = [record(outer, inner) for outer, inner in
                   ((0, 0), (1, 1), (0, 2), (-1, 1), (0, 0))]
        write_table(path, [(1, samples)])
        path.write_text(path.read_text().replace('T="row_1"', 'T="arclength"'))
        data = read_tecplot(path)
        segments = root_segments(data)
        self.assertEqual(len(segments), 1)
        np.testing.assert_allclose(data.column("outer")[segments[0]], [0, 1, 0, -1, 0])

    def test_duplicate_evaluations_do_not_force_irregular_cloud_into_grid(self):
        path = self.directory / "duplicates.plt"
        first = [record(1, inner, growth=inner - 1, branch=0) for inner in (0, 1, 1.5, 2)]
        first.extend([record(1, 1.5, growth=0.5, branch=0)] * 100)
        second = [record(2, inner, growth=inner - 1, branch=0) for inner in (0, 1, 2)]
        write_table(path, [(1, first), (2, second)])
        self.assertIsNone(rectilinear_field(read_tecplot(path)))
        self.assertIsNotNone(field_triangulation(read_tecplot(path)))

    def test_rectilinear_field_leaves_large_inner_gap_blank(self):
        path = self.directory / "large_inner_gap.plt"
        write_table(path, [(row, [record(row, inner, growth=inner, branch=0)
                                  for inner in (0, 1, 100, 101)]) for row in (1, 2, 3)])
        xs, ys, growth = rectilinear_field(read_tecplot(path))
        self.assertEqual(len(ys), 5)
        self.assertTrue(np.all(np.ma.getmaskarray(growth)[2]))

    def test_cli_exports_valid_image_with_field_and_root_gaps(self):
        roots, field = self.directory / "roots.plt", self.directory / "field.plt"
        write_table(roots, [(row, [record(row, 2 + row * 0.05)]) for row in (1, 2, 4, 5)])
        write_table(field, [(row, [record(row, inner, growth=inner - 2,
                                         branch=0, status=1 if (row, inner) == (3, 2) else 0)
                                   for inner in range(5)]) for row in range(1, 6)])
        output = self.directory / "test.png"
        with contextlib.redirect_stdout(io.StringIO()):
            code = main([str(roots), "--field", str(field), "--output", str(output),
                         "--x-label", "beta", "--y-label", "alpha"])
        self.assertEqual(code, 0)
        self.assertEqual(output.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
        self.assertGreater(output.stat().st_size, 10000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
