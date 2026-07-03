#!/usr/bin/env python3

from pathlib import Path


OUTDIR = Path(__file__).resolve().parent
N = 520


def wrap(seq, width=80):
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def write_fasta(path, prefix, count, seq):
    with open(path, "w") as out:
        for index in range(1, count + 1):
            out.write(f">{prefix}{index:04d}\n")
            out.write(wrap(seq(index)) + "\n")


def coding_seq(index):
    body = ("GCCGAACTG" * 18)[:162]
    return "ATG" + body + "TAA"


def noncoding_seq(index):
    return (("ATATGCGT" if index % 2 else "TGTACACA") * 22)[:170]


def write_reference_and_novel():
    genome = ["N"] * 400000
    gtf_lines = []
    pos = 1

    for index in range(1, N + 1):
        tid = f"mrna{index:04d}"
        gid = f"gene{index:04d}"
        seq = coding_seq(index)
        start = pos
        end = pos + len(seq) - 1
        genome[start - 1 : end] = seq
        attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";'
        gtf_lines.extend(
            [
                f"chrTest\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attrs}",
                f"chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attrs}",
                f"chrTest\ttest\tCDS\t{start}\t{end - 3}\t.\t+\t0\t{attrs}",
            ]
        )
        pos = end + 30

    for index in range(1, N + 1):
        tid = f"lncrna{index:04d}"
        gid = f"lncgene{index:04d}"
        seq = noncoding_seq(index)
        start = pos
        end = pos + len(seq) - 1
        genome[start - 1 : end] = seq
        attrs = f'gene_id "{gid}"; transcript_id "{tid}"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";'
        gtf_lines.extend(
            [
                f"chrTest\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attrs}",
                f"chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attrs}",
            ]
        )
        pos = end + 30

    with open(OUTDIR / "novel_transcripts.gtf", "w") as gtf, open(
        OUTDIR / "novel_transcripts.tmap", "w"
    ) as tmap:
        tmap.write("ref_gene_id\tref_id\tclass_code\tqry_gene_id\tqry_id\n")
        pos += 1000
        for index in range(1, N + 1):
            tid = f"novel{index:04d}"
            gid = f"novel_gene{index:04d}"
            seq = coding_seq(index)
            start = pos
            end = pos + len(seq) - 1
            genome[start - 1 : end] = seq
            attrs = f'gene_id "{gid}"; transcript_id "{tid}";'
            gtf.write(f"chrTest\ttest\ttranscript\t{start}\t{end}\t.\t+\t.\t{attrs}\n")
            gtf.write(f"chrTest\ttest\texon\t{start}\t{end}\t.\t+\t.\t{attrs}\n")
            tmap.write(f"-\t-\tu\t{gid}\t{tid}\n")
            pos = end + 80

    with open(OUTDIR / "genome.fa", "w") as out:
        out.write(">chrTest\n")
        out.write(wrap("".join(genome)) + "\n")

    with open(OUTDIR / "reference.gtf", "w") as out:
        out.write("\n".join(gtf_lines) + "\n")


def main():
    write_reference_and_novel()
    write_fasta(OUTDIR / "cds.fa", "cds", N, coding_seq)
    write_fasta(OUTDIR / "mrna.fa", "mrna", N, coding_seq)
    write_fasta(OUTDIR / "noncoding.fa", "noncoding", N, noncoding_seq)
    write_fasta(OUTDIR / "novel_transcripts.fa", "novel", N, coding_seq)


if __name__ == "__main__":
    main()
