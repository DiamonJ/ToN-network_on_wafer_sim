import json
import tempfile
import unittest
from pathlib import Path

from estimate_lammps_cost import estimate, parse_system
from summarize_compute_profile import summarize
from validate_static_cost import compare


class StaticBoundsTest(unittest.TestCase):
    def test_estimator_emits_uncalibrated_bounds_without_midpoint(self):
        root = Path(__file__).resolve().parents[1]
        result = estimate(parse_system(root / "cases/c1_validation_cu_108/in.lammps"))
        self.assertEqual(result["schema_version"], 2)
        self.assertEqual(result["assumptions"]["calibration"], "none")
        self.assertNotIn("C1_steady_ops_midpoint", result["summary"])
        for row in result["ranks"]:
            self.assertLessEqual(
                row["C1_steady_ops_lower"], row["C1_steady_ops_upper"]
            )
            self.assertEqual(
                row["C1_steady_ops_budget"], row["C1_steady_ops_upper"]
            )
            self.assertLessEqual(
                row["T1_send_bytes_lower"], row["T1_send_bytes_upper"]
            )
        assumptions = json.dumps(result["assumptions"])
        self.assertNotIn("160 bytes", assumptions)


class ComputeProfileTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write_perf(self, name, cycles, instructions, duration_ns):
        (self.root / name).write_text(
            "# synthetic perf output\n"
            f"{cycles};;cycles:u;1;100.00;;\n"
            f"{instructions};;instructions:u;1;100.00;;\n"
            f"{duration_ns};ns;duration_time:u;1;100.00;;\n",
            encoding="utf-8",
        )

    def test_wall_time_and_frequency_are_not_fitted(self):
        self.write_perf("baseline.repeat0.rank0000.csv", 100, 200, 1_000_000)
        self.write_perf("run.repeat0.rank0000.csv", 1100, 2200, 3_000_000)
        result = summarize(
            self.root, 1, 10, 1, 2.5, 2.0, 3.0, "test", "core",
            ["cycles", "instructions", "duration_time"],
        )
        rank = result["ranks"][0]
        self.assertEqual(result["calibration"], "none")
        self.assertEqual(rank["retired_instructions_per_step"], 200)
        self.assertEqual(rank["hardware_cycles_per_step"], 100)
        self.assertAlmostEqual(rank["wall_seconds_per_step"], 0.0002)
        self.assertEqual(rank["wall_cycles_lower_per_step"], 400_000)
        self.assertEqual(rank["wall_cycles_nominal_per_step"], 500_000)
        self.assertEqual(rank["wall_cycles_upper_per_step"], 600_000)
        self.assertNotIn("dp_ops_per_step", result)
        self.assertEqual(
            result["C1_measured"]["total_retired_instructions_per_step"],
            200,
        )
        self.assertEqual(
            result["C1_measured"]["critical_rank_wall_cycle_estimate_per_step"],
            {"lower": 400_000, "nominal": 500_000, "upper": 600_000},
        )


class BoundsValidationTest(unittest.TestCase):
    def fixtures(self):
        estimate_data = {
            "input_file": "/tmp/in.lammps",
            "input": {"num_ranks": 1, "atoms": 4},
            "ranks": [{
                "rank": 0,
                "T1_steady_send_bytes": 100,
                "T1_rebuild_send_bytes": 150,
                "C1_steady_ops_min": 10,
                "C1_steady_ops_max": 20,
                "C1_steady_ops_lower": 10,
                "C1_steady_ops_upper": 20,
            }],
        }
        plan = {
            "num_ranks": 1,
            "records": [
                {"kind": "message", "scope": "run", "rank": 0, "bytes": 100},
                {
                    "kind": "message",
                    "scope": "setup",
                    "phase": "borders",
                    "rank": 0,
                    "bytes": 150,
                },
            ],
        }
        profile = {
            "input_file": "/tmp/in.lammps",
            "num_ranks": 1,
            "method": "synthetic",
            "event_encoding": "core",
            "steps": 10,
            "ranks": [{
                "rank": 0,
                "retired_instructions_per_step": 1000,
                "hardware_cycles_per_step": 500,
                "wall_seconds_per_step": 1e-6,
                "wall_cycles_lower_per_step": 400,
                "wall_cycles_nominal_per_step": 500,
                "wall_cycles_upper_per_step": 600,
            }],
        }
        return estimate_data, plan, profile

    def test_validator_checks_bounds_without_cross_unit_fit(self):
        estimate_data, plan, profile = self.fixtures()
        result = compare(estimate_data, plan, profile)
        self.assertTrue(result["passed"])
        self.assertEqual(result["calibration"], "none")
        self.assertNotIn("global_relative_error", result["c1"])
        self.assertTrue(
            result["c1"]["hardware_cycles_inside_wall_frequency_bounds"]
        )
        self.assertEqual(
            result["t1_iteration_send_bytes"]["observed_upper_sum"], 150
        )

    def test_validator_rejects_mismatched_input(self):
        estimate_data, plan, profile = self.fixtures()
        profile["input_file"] = "/tmp/other.lammps"
        with self.assertRaisesRegex(ValueError, "different input files"):
            compare(estimate_data, plan, profile)


if __name__ == "__main__":
    unittest.main()
