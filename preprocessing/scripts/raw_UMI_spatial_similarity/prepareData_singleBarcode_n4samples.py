#!/usr/bin/env python3

import gzip
from pathlib import Path

from tqdm import tqdm


# ==============================================================================
# Ustawienia
# ==============================================================================

project_dir = Path("/home/mateusz/projects/ippas-kunevicius-spatial")
raw_dir = project_dir / "raw"

output_dir = (
    project_dir
    / "results"
    / "fastq_barcode_comparison"
    / "selected_fragments_n4samples"
)

barcode = "CATGAACCTCTTATCA"

# Maksymalna liczba pierwszych par odczytów analizowana dla każdej próbki.
# Ustaw None, aby przeskanować cały FASTQ.
max_read_pairs_to_scan = 5_000_000

# Maksymalna liczba pasujących par zapisana dla każdej próbki.
# Ustaw None, aby zapisać wszystkie znalezione odczyty.
max_matching_pairs = 100_000


samples = {
    "20_1F": {
        "R1": raw_dir / "20_1F" / "20_1F_S5_L001_R1_001.fastq.gz",
        "R2": raw_dir / "20_1F" / "20_1F_S5_L001_R2_001.fastq.gz",
    },
    "15_1M": {
        "R1": raw_dir / "15_1M" / "15_1M_S5_L001_R1_001.fastq.gz",
        "R2": raw_dir / "15_1M" / "15_1M_S5_L001_R2_001.fastq.gz",
    },
    "12_3F": {
        "R1": raw_dir / "12_3F" / "12_3F_S10_L002_R1_001.fastq.gz",
        "R2": raw_dir / "12_3F" / "12_3F_S10_L002_R2_001.fastq.gz",
    },
    "20_3M": {
        "R1": raw_dir / "20_3M" / "20_3M_S10_L002_R1_001.fastq.gz",
        "R2": raw_dir / "20_3M" / "20_3M_S10_L002_R2_001.fastq.gz",
    },
}


# ==============================================================================
# Przygotowanie
# ==============================================================================

output_dir.mkdir(parents=True, exist_ok=True)

summary_file = output_dir / "summary.tsv"

with summary_file.open("w") as summary:
    summary.write(
        "sample_ID\tprocessed_pairs\tmatched_pairs\tselected_R1\tselected_R2\n"
    )

    for sample_id, paths in samples.items():

        r1_file = paths["R1"]
        r2_file = paths["R2"]

        if not r1_file.exists():
            raise FileNotFoundError(f"Brak pliku R1: {r1_file}")

        if not r2_file.exists():
            raise FileNotFoundError(f"Brak pliku R2: {r2_file}")

        sample_output_dir = output_dir / sample_id
        sample_output_dir.mkdir(parents=True, exist_ok=True)

        selected_r1_file = (
            sample_output_dir / f"{sample_id}_selected_R1.fastq.gz"
        )

        selected_r2_file = (
            sample_output_dir / f"{sample_id}_selected_R2.fastq.gz"
        )

        processed_pairs = 0
        matched_pairs = 0

        print(f"\nPrzetwarzanie próbki: {sample_id}")

        with (
            gzip.open(r1_file, "rt") as r1_handle,
            gzip.open(r2_file, "rt") as r2_handle,
            gzip.open(
                selected_r1_file,
                "wt",
                compresslevel=1,
            ) as output_r1,
            gzip.open(
                selected_r2_file,
                "wt",
                compresslevel=1,
            ) as output_r2,
            tqdm(
                total=max_read_pairs_to_scan,
                desc=sample_id,
                unit=" par",
                unit_scale=True,
            ) as progress,
        ):

            while True:

                r1_header = r1_handle.readline()
                r2_header = r2_handle.readline()

                if not r1_header and not r2_header:
                    break

                if not r1_header or not r2_header:
                    raise RuntimeError(
                        f"{sample_id}: R1 i R2 mają różną liczbę rekordów."
                    )

                r1_sequence = r1_handle.readline()
                r1_plus = r1_handle.readline()
                r1_quality = r1_handle.readline()

                r2_sequence = r2_handle.readline()
                r2_plus = r2_handle.readline()
                r2_quality = r2_handle.readline()

                if not all([
                    r1_sequence,
                    r1_plus,
                    r1_quality,
                    r2_sequence,
                    r2_plus,
                    r2_quality,
                ]):
                    raise RuntimeError(
                        f"{sample_id}: znaleziono niepełny rekord FASTQ."
                    )

                processed_pairs += 1
                progress.update(1)

                if r1_sequence.startswith(barcode):

                    output_r1.write(r1_header)
                    output_r1.write(r1_sequence)
                    output_r1.write(r1_plus)
                    output_r1.write(r1_quality)

                    output_r2.write(r2_header)
                    output_r2.write(r2_sequence)
                    output_r2.write(r2_plus)
                    output_r2.write(r2_quality)

                    matched_pairs += 1

                    if (
                        max_matching_pairs is not None
                        and matched_pairs >= max_matching_pairs
                    ):
                        break

                if (
                    max_read_pairs_to_scan is not None
                    and processed_pairs >= max_read_pairs_to_scan
                ):
                    break

        summary.write(
            f"{sample_id}\t"
            f"{processed_pairs}\t"
            f"{matched_pairs}\t"
            f"{selected_r1_file}\t"
            f"{selected_r2_file}\n"
        )

        print(f"Przeskanowane pary: {processed_pairs:,}")
        print(f"Pasujące pary:      {matched_pairs:,}")
        print(f"R1: {selected_r1_file}")
        print(f"R2: {selected_r2_file}")


print(f"\nGotowe. Podsumowanie:\n{summary_file}")