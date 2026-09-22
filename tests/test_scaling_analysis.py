import unittest

from analyze_scaling import measured_compute, render_markdown


class ScalingAnalysisTest(unittest.TestCase):
    def test_measured_compute_keeps_total_and_critical_rank_metrics(self):
        profile = {
            "ranks": [
                {
                    "retired_instructions_per_step": 10,
                    "hardware_cycles_per_step": 5,
                    "wall_seconds_per_step": 1e-6,
                    "wall_cycles_nominal_per_step": 2500,
                },
                {
                    "retired_instructions_per_step": 30,
                    "hardware_cycles_per_step": 15,
                    "wall_seconds_per_step": 2e-6,
                    "wall_cycles_nominal_per_step": 5000,
                },
            ]
        }
        result = measured_compute(profile)
        self.assertEqual(result["measured_retired_instructions_per_step"], 40)
        self.assertEqual(
            result["measured_retired_instructions_per_rank_average"], 20
        )
        self.assertEqual(
            result["measured_retired_instructions_per_rank_maximum"], 30
        )
        self.assertEqual(result["measured_hardware_cycles_per_step_sum"], 20)
        self.assertEqual(result["measured_critical_wall_seconds_per_step"], 2e-6)

    def test_summary_table_uses_critical_path_compute_to_communication_ratio(self):
        row = {
            "scenario_scale": "2688 atoms / 16 ranks (4x4)",
            "ranks": 16,
            "compute_to_communication_ratio": 2.0,
            "expected_cycles": 300,
            "theoretical_critical_rank_ops_lower": 10,
            "theoretical_critical_rank_ops_upper": 20,
            "theoretical_compute_ops_lower": 100,
            "theoretical_compute_ops_upper": 200,
            "measured_dp_ops_per_step_optional": 150,
            "measured_retired_instructions_per_rank_maximum": 30,
            "measured_retired_instructions_per_step": 40,
            "expected_compute_cycles": 200,
            "expected_communication_cycles": 100,
            "booksim_measured_cycles": 301,
            "booksim_relative_error": 1 / 300,
        }
        table = render_markdown([row])
        self.assertIn("| 场景规模 | 计算通信比 | 预期 cycles |", table)
        self.assertIn("| 2688 atoms / 16 ranks (4x4) | 2.0000 | 300 |", table)
        self.assertIn("200/100", table)


if __name__ == "__main__":
    unittest.main()
