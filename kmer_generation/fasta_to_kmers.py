import gzip
import shutil
import subprocess
from pathlib import Path


# Default absolute path to the dicey executable.
DICEY_PATH = "/g/korbel/shared/software/dicey/bin/dicey"


def _run_command(command):
    # Execute an external command and fail immediately if it returns a non-zero exit code.
    subprocess.run(command, check=True)


def _sort_and_uniq(input_path, output_path, add_counts, gzipped_output):
    # Build the uniq command in case-insensitive mode and optionally include counts.
    uniq_command = ["uniq", "-i"]
    if add_counts:
        uniq_command.append("-c")

    output_path = Path(output_path)
    if gzipped_output:
        # Write plain output first, then gzip in a separate step.
        if output_path.suffix == ".gz":
            plain_output_path = output_path.with_suffix("")
        else:
            plain_output_path = Path(str(output_path) + ".tmp")
    else:
        plain_output_path = output_path

    # Start sort and stream its output directly into uniq, writing to disk.
    sort_process = subprocess.Popen(["sort", str(input_path)], stdout=subprocess.PIPE)
    if sort_process.stdout is None:
        raise RuntimeError("Failed to capture sort stdout")

    with open(plain_output_path, "wb") as plain_output_file:
        uniq_process = subprocess.Popen(uniq_command, stdin=sort_process.stdout, stdout=plain_output_file)
        sort_process.stdout.close()

        # Wait for both processes and surface pipeline failures.
        uniq_return_code = uniq_process.wait()
        sort_return_code = sort_process.wait()

    if sort_return_code != 0 or uniq_return_code != 0:
        raise subprocess.CalledProcessError(
            uniq_return_code if uniq_return_code != 0 else sort_return_code,
            "sort|uniq",
        )

    # Optionally gzip the already materialized plain-text output.
    if gzipped_output:
        with open(plain_output_path, "rb") as input_file, gzip.open(output_path, "wb") as output_file:
            shutil.copyfileobj(input_file, output_file)

        if plain_output_path.exists():
            plain_output_path.unlink()


def fasta_to_kmers(
    fasta_path,
    sample,
    workdir,
    functionality,
    gzipped_output=False,
    dicey_path=DICEY_PATH,
    kmer_size=61,
):
    """
    Generate kmers from a FASTA file.

    functionality values:
    - "unique": unique forward+reverse-complement kmers
    - "unique_counts": counted unique forward+reverse-complement kmers
    - "forward_unique_counts": counted unique forward-only kmers
    """

    ### Check for valid functionality
    valid_functionalities = {"unique", "unique_counts", "forward_unique_counts"}
    if functionality not in valid_functionalities:
        raise ValueError(f"functionality must be one of {valid_functionalities}")

    ### Prepare paths
    workdir_path = Path(workdir)
    raw_list_path = workdir_path / f"{sample}.dicey_raw.list"
    normal_prefix = workdir_path / f"{sample}.dicey_raw_normal"

    ### Run first dicey chop for normal forward kmers, save in {normal_prefix}.fq.gz
    _run_command(
        [
            dicey_path,
            "chop",
            "--length",
            str(kmer_size),
            "--se",
            "-f",
            str(normal_prefix),
            str(fasta_path),
        ]
    )
    
    ### Collect the forward fastq.gz output and decide whether reverse-complement kmers are needed.
    fastq_gz_paths = [Path(f"{normal_prefix}.fq.gz")]
    include_revcomp = functionality in {"unique", "unique_counts"}

    ### Optionally generate reverse-complement kmers (using --revcomp for dicey), save this in {revcomp_prefix}.fq.gz and append their output file.
    if include_revcomp:
        revcomp_prefix = workdir_path / f"{sample}.dicey_raw_revcomp"
        _run_command(
            [
                dicey_path,
                "chop",
                "--length",
                str(kmer_size),
                "--se",
                "--revcomp",
                "-f",
                str(revcomp_prefix),
                str(fasta_path),
            ]
        )
        fastq_gz_paths.append(Path(f"{revcomp_prefix}.fq.gz"))

    ### Parse FASTQ records and write only sequence lines (line index 1 in each 4-line record).
    with open(raw_list_path, "w") as kmer_list_file:
        for fastq_gz_path in fastq_gz_paths:
            with gzip.open(fastq_gz_path, "rt") as fastq_file:
                for line_number, line in enumerate(fastq_file):
                    if (line_number % 4) == 1:
                        kmer_list_file.write(line)

    ### Remove intermediate dicey FASTQ files once kmer sequences are extracted.
    for fastq_gz_path in fastq_gz_paths:
        if fastq_gz_path.exists():
            fastq_gz_path.unlink()

    ### Decide whether output should include counts and build output filename.
    add_counts = functionality in {"unique_counts", "forward_unique_counts"}
    if add_counts:
        output_name = f"{sample}.dicey_unique.counts.list"
    else:
        output_name = f"{sample}.dicey_unique.list"
    if gzipped_output:
        output_name += ".gz"

    # Sort and collapse kmers to unique entries (with or without counts).
    output_path = workdir_path / output_name
    _sort_and_uniq(raw_list_path, output_path, add_counts=add_counts, gzipped_output=gzipped_output)

    ### Remove temporary raw kmer list file.
    if raw_list_path.exists():
        raw_list_path.unlink()