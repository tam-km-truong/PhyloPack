import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from tests.generate_mock_data import (
    generate_mock_dataset,
    mutate_sequence,
    generate_random_sequence,
    partition_species_counts,
)


class TestMockDataGenerator(unittest.TestCase):

    def test_random_sequence_generation(self):
        seq = generate_random_sequence(1000, gc_content=0.5)
        self.assertEqual(len(seq), 1000)
        self.assertTrue(set(seq).issubset({'A', 'C', 'G', 'T'}))

    def test_sequence_mutation(self):
        seq = "A" * 1000
        mutated = mutate_sequence(seq, mutation_rate=0.05)
        self.assertEqual(len(mutated), 1000)
        diff_count = sum(1 for a, b in zip(seq, mutated) if a != b)
        self.assertEqual(diff_count, 50)

    def test_partition_species_counts(self):
        # Equal partition
        counts_equal = partition_species_counts(num_species=4, total_genomes=40, unequal=False)
        self.assertEqual(counts_equal, [10, 10, 10, 10])

        # Unequal partition preserves exact total
        counts_unequal = partition_species_counts(num_species=4, total_genomes=40, unequal=True)
        self.assertEqual(sum(counts_unequal), 40)
        self.assertTrue(len(set(counts_unequal)) > 1)
        self.assertTrue(all(c >= 1 for c in counts_unequal))

    def test_generate_mock_dataset_neutral_accessions_and_tsv(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir) / "test_synth"
            num_species = 3
            strains_per_species = 4
            result = generate_mock_dataset(
                output_dir=output_dir,
                num_species=num_species,
                strains_per_species=strains_per_species,
                genome_length=5000,
                seed=123,
                accession_prefix="GENOME",
                relative_manifest=False
            )

            manifest_path = result.manifest_path
            tsv_path = result.tsv_path

            self.assertTrue(manifest_path.exists())
            self.assertTrue(tsv_path.exists())

            with open(manifest_path) as f:
                lines = [line.strip() for line in f if line.strip()]

            expected_count = num_species * strains_per_species
            self.assertEqual(len(lines), expected_count)

            # Check that each file exists, is valid gzip, and has neutral accession
            for line in lines:
                p = Path(line)
                self.assertTrue(p.exists())
                self.assertTrue(p.name.endswith(".fa.gz"))

                accession = p.name.split(".")[0]
                self.assertTrue(accession.startswith("GENOME"))
                self.assertNotIn("species", accession.lower())
                self.assertNotIn("sp0", accession.lower())

                with gzip.open(p, "rt") as f:
                    content = f.read()
                    self.assertTrue(content.startswith(f">{accession}"))
                    self.assertNotIn("species", content.splitlines()[0].lower())
                    seq = "".join(l.strip() for l in content.splitlines() if not l.startswith(">"))
                    self.assertEqual(len(seq), 5000)

            # Check TSV mapping
            with open(tsv_path) as f:
                reader = csv.DictReader(f, delimiter="\t")
                self.assertEqual(reader.fieldnames, ["accession", "species_id", "genome_path"])
                rows = list(reader)

            self.assertEqual(len(rows), expected_count)

            species_counts = {}
            for row in rows:
                sp = row["species_id"]
                species_counts[sp] = species_counts.get(sp, 0) + 1
                self.assertTrue(row["accession"].startswith("GENOME"))
                self.assertTrue(Path(row["genome_path"]).exists())

            self.assertEqual(len(species_counts), num_species)
            for count in species_counts.values():
                self.assertEqual(count, strains_per_species)

    def test_generate_mock_dataset_unequal_species(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir) / "test_synth_unequal"
            total_genomes = 40
            num_species = 4
            result = generate_mock_dataset(
                output_dir=output_dir,
                num_species=num_species,
                total_genomes=total_genomes,
                unequal_species=True,
                genome_length=2000,
                seed=42,
                relative_manifest=False
            )

            counts = list(result.species_counts.values())
            self.assertEqual(sum(counts), total_genomes)
            self.assertEqual(len(counts), num_species)
            # Verify counts are not all identical
            self.assertTrue(len(set(counts)) > 1)

            # Verify TSV matches counts
            with open(result.tsv_path) as f:
                rows = list(csv.DictReader(f, delimiter="\t"))
            self.assertEqual(len(rows), total_genomes)

            tsv_counts = {}
            for row in rows:
                tsv_counts[row["species_id"]] = tsv_counts.get(row["species_id"], 0) + 1
            self.assertEqual(tsv_counts, result.species_counts)

    def test_generate_mock_dataset_explicit_counts(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            output_dir = Path(tmp_dir) / "test_synth_explicit"
            explicit = [12, 7, 3]
            result = generate_mock_dataset(
                output_dir=output_dir,
                species_counts=explicit,
                genome_length=2000,
                relative_manifest=False
            )

            self.assertEqual(list(result.species_counts.values()), explicit)
            self.assertEqual(sum(result.species_counts.values()), sum(explicit))


if __name__ == "__main__":
    unittest.main()
