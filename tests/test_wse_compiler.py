import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from wse_compiler import build_replay_ccdg, compile_program  # noqa: E402


class WseCompilerTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.plan_path = root / "plan.json"
        self.cost_path = root / "cost.json"
        self.cfg_path = root / "booksim.cfg"
        self.cfg_path.write_text(
            """
            sim_type = ccdg;
            topology = mesh;
            k = 2;
            n = 2;
            flit_size_bytes = 8;
            noc_frequency_ghz = 2.0;
            ccdg_compute_capability = 2.5e9;
            routing_delay = 1;
            vc_alloc_delay = 1;
            sw_alloc_delay = 1;
            """,
            encoding="utf-8",
        )
        metadata = []
        for rank, (x, y) in enumerate(((0, 0), (0, 1), (1, 0), (1, 1))):
            metadata.append({"rank": rank, "myloc": [x, y, 0]})
        records = []
        seq = 0
        for phase, direction in (("forward", 1), ("reverse", -1)):
            for rank in range(4):
                x, y = metadata[rank]["myloc"][:2]
                dst = next(
                    item["rank"]
                    for item in metadata
                    if item["myloc"][:2] == [1 - x, y]
                )
                records.append(
                    {
                        "kind": "message",
                        "scope": "run",
                        "timestep": 1,
                        "phase": phase,
                        "dimension": 0,
                        "direction": direction,
                        "round": 0,
                        "seq": seq,
                        "src": rank,
                        "dst": dst,
                        "bytes": 16 + rank,
                    }
                )
                seq += 1
        self.plan = {
            "num_ranks": 4,
            "procgrid": [2, 2, 1],
            "rank_metadata": metadata,
            "records": records,
        }
        self.cost = {
            "ranks": [
                {"rank": rank, "C1_steady_ops_midpoint": 100 + rank}
                for rank in range(4)
            ]
        }
        self.plan_path.write_text(json.dumps(self.plan), encoding="utf-8")
        self.cost_path.write_text(json.dumps(self.cost), encoding="utf-8")

    def tearDown(self):
        self.tmp.cleanup()

    def test_short_plan_compiles_with_conservation_and_reduction(self):
        program, est, report = compile_program(
            self.plan_path, self.cost_path, self.cfg_path
        )
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["conservation"]["input_messages"], 8)
        self.assertEqual(
            report["conservation"]["input_bytes"],
            report["conservation"]["emitted_logical_bytes"],
        )
        self.assertTrue(report["assertions"]["link_schedule_conflict_free"])
        self.assertTrue(program["hardware"]["fold_pbc"])
        self.assertEqual(len({tuple(item["physical"]) for item in program["placement"]}), 4)
        stages = program["timesteps"][0]["stages"]
        self.assertEqual([item["mode"] for item in stages], ["multicast", "reduction"])
        self.assertEqual([item["phase_count"] for item in stages], [2, 2])
        self.assertEqual(est, sorted(est, key=lambda item: (item[1], item[0])))
        replay, schedule, stats = build_replay_ccdg(program)
        self.assertEqual(replay["num_ranks"], 4)
        self.assertEqual(stats["commands"], 4)
        self.assertEqual(stats["branches"], 8)
        self.assertEqual(stats["packets"], stats["commands"] + stats["branches"])
        self.assertEqual(
            stats["branch_link_bytes"], report["conservation"]["input_bytes"]
        )
        self.assertEqual(len(schedule), len(replay["nodes"]))

    def test_long_range_phase_is_rejected(self):
        self.plan["records"][0]["phase"] = "fft_remap"
        self.plan_path.write_text(json.dumps(self.plan), encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "cannot lower phases: fft_remap"):
            compile_program(self.plan_path, self.cost_path, self.cfg_path)

    def test_fast_profile_matches_supported_iq_minimum(self):
        program, _, _ = compile_program(
            self.plan_path,
            self.cost_path,
            self.cfg_path,
            wse_fast_profile=True,
        )
        self.assertEqual(program["hardware"]["profile"], "wse_fast")
        self.assertEqual(program["hardware"]["flit_size_bytes"], 4)
        self.assertEqual(program["hardware"]["hop_stride_cycles"], 3)

    def test_fold_can_be_disabled(self):
        program, _, report = compile_program(
            self.plan_path,
            self.cost_path,
            self.cfg_path,
            fold_pbc=False,
        )
        self.assertFalse(program["hardware"]["fold_pbc"])
        self.assertTrue(report["assertions"]["placement_bijective"])
        self.assertTrue(
            all(item["logical"] == item["physical"] for item in program["placement"])
        )

    def test_compute_capability_can_match_experiment_baseline(self):
        program, _, _ = compile_program(
            self.plan_path,
            self.cost_path,
            self.cfg_path,
            compute_capability=2.5e10,
        )
        self.assertEqual(
            program["hardware"]["compute_capability_ops_s"], 2.5e10
        )

    def test_new_cost_schema_uses_uncalibrated_upper_bound(self):
        self.cost["ranks"] = [
            {
                "rank": rank,
                "C1_steady_ops_lower": 10 + rank,
                "C1_steady_ops_upper": 20 + rank,
            }
            for rank in range(4)
        ]
        self.cost_path.write_text(json.dumps(self.cost), encoding="utf-8")
        program, _, _ = compile_program(
            self.plan_path, self.cost_path, self.cfg_path
        )
        self.assertEqual(
            [row["ops"] for row in program["compute_model"]["blocks"]],
            [20, 21, 22, 23],
        )


if __name__ == "__main__":
    unittest.main()
