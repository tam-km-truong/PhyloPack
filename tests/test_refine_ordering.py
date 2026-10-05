import argparse
import os
import shutil
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from phylopack.preorder.refine_ordering import (
    add_refine_args,
    add_refine_parser,
    run_refine,
)


class TestRefineOrdering(unittest.TestCase):

    def setUp(self):
        self.test_dir = tempfile.mkdtemp()
        self.addCleanup(lambda: shutil.rmtree(self.test_dir, ignore_errors=True))

        # Create dummy input file with 10 lines
        self.input_file = os.path.join(self.test_dir, "input_genomes.txt")
        self.sample_genomes = [f"/path/to/genome_{i}.fa" for i in range(10)]
        with open(self.input_file, "w") as f:
            for g in self.sample_genomes:
                f.write(g + "\n")

    def test_window_validation(self):
        """Window size smaller than 4 must raise ValueError."""
        output_file = os.path.join(self.test_dir, "out.txt")
        args = argparse.Namespace(
            input_genomes=self.input_file,
            output=output_file,
            window=3,
            min_size=10,
            verbose=False,
            debug=False,
        )
        with self.assertRaises(ValueError) as ctx:
            run_refine(args)
        self.assertIn("--window must be at least 4", str(ctx.exception))

    def test_nonexistent_input_file(self):
        """Nonexistent input path must raise FileNotFoundError."""
        output_file = os.path.join(self.test_dir, "out.txt")
        args = argparse.Namespace(
            input_genomes=os.path.join(self.test_dir, "does_not_exist.txt"),
            output=output_file,
            window=5,
            min_size=10,
            verbose=False,
            debug=False,
        )
        with self.assertRaises(FileNotFoundError):
            run_refine(args)

    def test_fallback_when_below_min_size(self):
        """When total genomes < min_size, copy input directly to output, creating directories as needed."""
        # Output is in a new nested directory
        nested_out = os.path.join(self.test_dir, "nested", "dir", "out.txt")
        args = argparse.Namespace(
            input_genomes=self.input_file,
            output=nested_out,
            window=5,
            min_size=20,  # 10 < 20
            verbose=False,
            debug=False,
        )
        run_refine(args)

        self.assertTrue(os.path.exists(nested_out))
        with open(nested_out) as f:
            lines = [line.strip() for line in f if line.strip()]
        self.assertEqual(lines, self.sample_genomes)

    @patch("phylopack.preorder.refine_ordering.attotree")
    def test_run_refine_window_processing_mocked(self, mock_attotree):
        """Verify windowing and concatenation logic when min_size is met."""
        output_file = os.path.join(self.test_dir, "out_refined.txt")

        # Mock attotree to simulate writing leaf_order output file
        def fake_attotree(args):
            # Write lines in reverse to simulate reordering
            with open(args.input_genomes) as in_f:
                lines = [l.strip() for l in in_f if l.strip()]
            with open(args.leaf_order, "w") as out_f:
                for line in reversed(lines):
                    out_f.write(line + "\n")

        mock_attotree.side_effect = fake_attotree

        args = argparse.Namespace(
            input_genomes=self.input_file,
            output=output_file,
            window=5,
            min_size=5,  # 10 genomes, window 5 -> 2 batches of 5
            k=21,
            s=10000,
            t=10,
            m="nj",
            verbose=False,
            debug=False,
        )

        run_refine(args)

        self.assertTrue(os.path.exists(output_file))
        self.assertEqual(mock_attotree.call_count, 2)

        with open(output_file) as f:
            refined_lines = [l.strip() for l in f if l.strip()]

        # Each batch of 5 was reversed: [4..0] + [9..5]
        expected = list(reversed(self.sample_genomes[:5])) + list(reversed(self.sample_genomes[5:]))
        self.assertEqual(refined_lines, expected)

    @patch("phylopack.preorder.refine_ordering.attotree")
    def test_debug_mode_preserves_temp_dir(self, mock_attotree):
        """When --debug is set, the phylopack_tmp directory must be preserved."""
        mock_attotree.side_effect = lambda a: open(a.leaf_order, "w").close()

        output_file = os.path.join(self.test_dir, "out.txt")
        args = argparse.Namespace(
            input_genomes=self.input_file,
            output=output_file,
            window=5,
            min_size=5,
            k=21,
            s=10000,
            t=10,
            m="nj",
            verbose=False,
            debug=True,
        )

        run_refine(args)
        tmp_dir = os.path.join(self.test_dir, "phylopack_tmp")
        self.assertTrue(os.path.exists(tmp_dir))

    def test_cli_parser_defaults(self):
        """Verify argparse options and defaults configured for refine command."""
        subparsers = argparse.ArgumentParser().add_subparsers(dest="command")
        add_refine_parser(subparsers)

        parser = subparsers.choices["refine"]
        parsed = parser.parse_args(["in.txt", "out.txt"])

        self.assertEqual(parsed.input_genomes, "in.txt")
        self.assertEqual(parsed.output, "out.txt")
        self.assertEqual(parsed.window, 5000)
        self.assertEqual(parsed.min_size, 10000)
        self.assertEqual(parsed.k, 21)
        self.assertEqual(parsed.s, 10000)
        self.assertEqual(parsed.t, 10)
        self.assertEqual(parsed.m, "nj")
        self.assertFalse(parsed.verbose)
        self.assertFalse(parsed.debug)

    def test_real_attotree_integration(self):
        """Integration test with real attotree binary if installed on the system."""
        if shutil.which("attotree") is None:
            self.skipTest("attotree binary not found in PATH")

        real_data = "tests/data/genomes.txt"
        if not os.path.exists(real_data):
            self.skipTest("tests/data/genomes.txt not found")

        output_file = os.path.join(self.test_dir, "real_refined.txt")
        args = argparse.Namespace(
            input_genomes=real_data,
            output=output_file,
            window=5,
            min_size=4,
            k=21,
            s=10000,
            t=2,
            m="nj",
            verbose=False,
            debug=False,
        )

        run_refine(args)
        self.assertTrue(os.path.exists(output_file))
        with open(output_file) as f:
            lines = [l.strip() for l in f if l.strip()]
        self.assertEqual(len(lines), 10)


if __name__ == "__main__":
    unittest.main()
