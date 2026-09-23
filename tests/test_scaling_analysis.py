import unittest

from analyze_scaling import measured_compute, render_markdown, render_result_document


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
            "atoms": 2688,
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
            "booksim_average_packet_queue_cycles": 1.5,
            "booksim_average_flit_queue_cycles": 3.5,
            "booksim_average_injection_rate": 0.2,
            "booksim_injection_saturation_ratio": 0.25,
            "booksim_saturated_injection_rate": 0.8,
            "booksim_aggregate_communication_to_compute_ratio": 0.5,
        }
        table = render_markdown([row])
        self.assertIn("| 场景规模 | 计算通信比 | 预期 cycles |", table)
        self.assertIn("| 2688 atoms / 16 ranks (4x4) | 2.0000 | 300 |", table)
        self.assertIn("200/100", table)
        report = render_result_document(
            [row], {"compute_capability_ops_s": 2.5e10}
        )
        self.assertIn("平均 packet 排队", report)
        self.assertIn("注入饱和占比", report)
        self.assertIn("饱和时注入率", report)
        self.assertIn("BookSim 通信/计算比", report)
        self.assertIn("## 1. 计算量改造怎么实现", report)
        self.assertIn("## 2. 实验结果", report)


if __name__ == "__main__":
    unittest.main()
