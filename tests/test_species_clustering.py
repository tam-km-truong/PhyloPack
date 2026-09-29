import os
import tempfile
import unittest
from argparse import Namespace
from pathlib import Path

from phylopack.macro.species_clustering import (
    normalize_delimiter,
    resolve_columns,
    count_and_rank_species,
    create_sorted_chunks,
    chunk_line_generator,
    merge_and_emit,
    run_species_clustering,
)


class TestSpeciesClustering(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.tmp_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_normalize_delimiter(self):
        self.assertEqual(normalize_delimiter("\\t"), "\t")
        self.assertEqual(normalize_delimiter("tab"), "\t")
        self.assertEqual(normalize_delimiter("\\s"), " ")
        self.assertEqual(normalize_delimiter("space"), " ")
        self.assertEqual(normalize_delimiter(","), ",")
        self.assertEqual(normalize_delimiter("\t"), "\t")

    def test_resolve_columns(self):
        header = "acc\tspecies\tpath\textra"
        # Test named resolution
        acc_idx, sp_idx, path_idx = resolve_columns(header, "\t", "acc", "species", "path")
        self.assertEqual((acc_idx, sp_idx, path_idx), (0, 1, 2))

        # Test optional path
        acc_idx, sp_idx, path_idx = resolve_columns(header, "\t", "acc", "species", None)
        self.assertEqual((acc_idx, sp_idx, path_idx), (0, 1, None))

        # Test numeric indices
        acc_idx, sp_idx, path_idx = resolve_columns(header, "\t", "3", "1", "0")
        self.assertEqual((acc_idx, sp_idx, path_idx), (3, 1, 0))

        # Test invalid column name
        with self.assertRaises(ValueError):
            resolve_columns(header, "\t", "nonexistent", "species", None)

    def test_count_and_rank_species(self):
        tsv_path = self.tmp_path / "metadata.tsv"
        with open(tsv_path, "w") as f:
            f.write("accession\tspecies_id\tgenome_path\n")
            f.write("acc1\tsp_large\t/path1\n")
            f.write("acc2\tsp_large\t/path2\n")
            f.write("acc3\tsp_large\t/path3\n")
            f.write("acc4\tsp_medium\t/path4\n")
            f.write("acc5\tsp_medium\t/path5\n")
            f.write("acc6\tsp_small\t/path6\n")

        # Test descending order (largest first)
        counts, sp_to_rank, sp_order = count_and_rank_species(str(tsv_path), "\t", 1, "desc")
        self.assertEqual(counts["sp_large"], 3)
        self.assertEqual(counts["sp_medium"], 2)
        self.assertEqual(counts["sp_small"], 1)
        self.assertEqual(sp_order, ["sp_large", "sp_medium", "sp_small"])
        self.assertEqual(sp_to_rank["sp_large"], 0)
        self.assertEqual(sp_to_rank["sp_small"], 2)

        # Test ascending order (smallest first)
        _, _, sp_order_asc = count_and_rank_species(str(tsv_path), "\t", 1, "asc")
        self.assertEqual(sp_order_asc, ["sp_small", "sp_medium", "sp_large"])

    def test_create_sorted_chunks_and_fallback(self):
        tsv_path = self.tmp_path / "metadata_fallback.tsv"
        with open(tsv_path, "w") as f:
            f.write("accession\tspecies_id\tgenome_path\n")
            f.write("accB\tsp1\t\n")  # Empty path -> fallback to accession
            f.write("accA\tsp1\trelative/path.fna\n")

        col_indices = (0, 1, 2)
        species_to_rank = {"sp1": 0}

        chunk_dir = self.tmp_path / "chunks"
        chunk_dir.mkdir()

        chunk_files = create_sorted_chunks(
            metadata_path=str(tsv_path),
            delim="\t",
            col_indices=col_indices,
            species_to_rank=species_to_rank,
            temp_dir=str(chunk_dir),
            chunk_size=10,
            path_prefix="/data/genomes"
        )

        self.assertEqual(len(chunk_files), 1)
        with open(chunk_files[0]) as f:
            lines = [line.strip().split("\t") for line in f]

        # Sorted by accession: accA then accB
        self.assertEqual(lines[0][2], "accA")
        self.assertEqual(lines[0][3], "/data/genomes/relative/path.fna")
        self.assertEqual(lines[1][2], "accB")
        self.assertEqual(lines[1][3], "accB")  # Fallback to accession

    def test_merge_and_emit(self):
        chunk_dir = self.tmp_path / "merge_chunks"
        chunk_dir.mkdir()

        # Create 2 simulated chunk files
        chunk1 = chunk_dir / "chunk_0000.tsv"
        with open(chunk1, "w") as f:
            # rank, species_id, accession, path
            f.write("0\tspA\taccA1\t/path/A1\n")
            f.write("1\tspB\taccB2\t/path/B2\n")

        chunk2 = chunk_dir / "chunk_0001.tsv"
        with open(chunk2, "w") as f:
            f.write("0\tspA\taccA2\t/path/A2\n")
            f.write("0\tspA\taccA3\t/path/A3\n")
            f.write("1\tspB\taccB1\t/path/B1\n")

        output_path = self.tmp_path / "ordered_genomes.txt"
        sep_path = self.tmp_path / "separators.tsv"
        species_counts = {"spA": 3, "spB": 2}

        separators = merge_and_emit(
            chunk_paths=[str(chunk1), str(chunk2)],
            output_path=str(output_path),
            species_counts=species_counts,
            boundaries_path=str(sep_path)
        )

        # Check ordered genome paths
        with open(output_path) as f:
            emitted_paths = [l.strip() for l in f]

        expected_paths = [
            "/path/A1",
            "/path/A2",
            "/path/A3",
            "/path/B1",
            "/path/B2",
        ]
        self.assertEqual(emitted_paths, expected_paths)

        # Check separators
        # spA (count 3, start 0): median index is 0 + (3 // 2) = 1 -> /path/A2
        # spB (count 2, start 3): median index is 3 + (2 // 2) = 4 -> /path/B2
        self.assertEqual(len(separators), 2)
        self.assertEqual(separators[0], ("spA", 0, 3, "/path/A2"))
        self.assertEqual(separators[1], ("spB", 3, 2, "/path/B2"))

        # Verify separators file written to disk
        self.assertTrue(sep_path.exists())
        with open(sep_path) as f:
            sep_lines = [l.strip() for l in f]
        self.assertEqual(sep_lines[0], "species_id\tstart_idx\tcount\trepresentative")
        self.assertEqual(sep_lines[1], "spA\t0\t3\t/path/A2")
        self.assertEqual(sep_lines[2], "spB\t3\t2\t/path/B2")

    def test_run_species_clustering_end_to_end(self):
        metadata_file = self.tmp_path / "input_metadata.tsv"
        with open(metadata_file, "w") as f:
            f.write("acc_id\tspecies\n")
            f.write("G003\tEscherichia_coli\n")
            f.write("G001\tEscherichia_coli\n")
            f.write("G002\tEscherichia_coli\n")
            f.write("G004\tSalmonella_enterica\n")
            f.write("G005\tSalmonella_enterica\n")
            f.write("G006\tBacillus_subtilis\n")

        out_genomes = self.tmp_path / "sub_dir" / "out_genomes.txt"
        out_separators = self.tmp_path / "sub_dir" / "separators.tsv"
        out_species_order = self.tmp_path / "sub_dir" / "species_order.txt"

        args = Namespace(
            metadata=str(metadata_file),
            output=str(out_genomes),
            delimiter="\t",
            separators_output=str(out_separators),
            species_order_output=str(out_species_order),
            accession_col="acc_id",
            species_col="species",
            path_col=None,
            path_prefix="",
            sort_species_by="desc",
            verbose=False
        )

        separators, species_order = run_species_clustering(args)

        # E. coli (3) > S. enterica (2) > B. subtilis (1)
        self.assertEqual(species_order, ["Escherichia_coli", "Salmonella_enterica", "Bacillus_subtilis"])

        # Check emitted genomes: no path given, so accessions emitted in order
        with open(out_genomes) as f:
            ordered_accessions = [l.strip() for l in f]

        self.assertEqual(
            ordered_accessions,
            ["G001", "G002", "G003", "G004", "G005", "G006"]
        )

        # Check species order file
        with open(out_species_order) as f:
            saved_species_order = [l.strip() for l in f]
        self.assertEqual(saved_species_order, species_order)

        # Check separators
        self.assertEqual(len(separators), 3)
        self.assertEqual(separators[0], ("Escherichia_coli", 0, 3, "G002"))
        self.assertEqual(separators[1], ("Salmonella_enterica", 3, 2, "G005"))
        self.assertEqual(separators[2], ("Bacillus_subtilis", 5, 1, "G006"))


if __name__ == "__main__":
    unittest.main()
