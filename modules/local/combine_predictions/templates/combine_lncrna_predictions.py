#!/usr/bin/env python3

from Bio import SeqIO
import platform


CPAT_CUTOFF = float("${task.ext.cpat_cutoff}" if "${task.ext.cpat_cutoff}" != "null" else 0.364)


def read_cpat(path):
    data = {}
    with open(path) as handle:
        header = handle.readline()
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 11:
                continue
            score = float(fields[10])
            data[fields[0]] = ["Non-coding" if score < CPAT_CUTOFF else "Coding", score]
    return data


def read_simple_prediction(path, noncoding_labels):
    data = {}
    with open(path) as handle:
        header = handle.readline()
        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if len(fields) < 3:
                continue
            prediction = "Non-coding" if fields[1] in noncoding_labels else "Coding"
            data[fields[0]] = [prediction, fields[2]]
    return data


def read_tmap(path):
    data = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\\n").lstrip("#").split("\\t")
        columns = {name: index for index, name in enumerate(header)}

        for line in handle:
            fields = line.rstrip("\\n").split("\\t")
            if not fields:
                continue
            query_id = fields[columns.get("qry_id", 4)]
            data[query_id] = {
                "class_code": fields[columns.get("class_code", 2)] if len(fields) > columns.get("class_code", 2) else "NA",
                "ref_gene_id": fields[columns.get("ref_gene_id", 0)] if len(fields) > columns.get("ref_gene_id", 0) else "NA",
                "ref_id": fields[columns.get("ref_id", 1)] if len(fields) > columns.get("ref_id", 1) else "NA",
            }
    return data


def consensus(votes, mode):
    coding = votes.count("Coding")
    total = len(votes)
    if total == 0:
        return "No Prediction"
    if mode == "strict":
        return "Protein-coding" if coding == total else "Non-coding"
    if mode == "lenient":
        return "Protein-coding" if coding >= 1 else "Non-coding"
    return "Protein-coding" if coding >= 2 else "Non-coding"


def write_summary(path, transcript_ids, predictions, tmap_data):
    header = [
        "transcript_id", "class_code", "ref_gene_id", "ref_id", "consensus",
        "votes", "cpat", "feelnc", "plek", "cpat_score", "feelnc_score", "plek_score"
    ]
    with open(path, "w") as out:
        out.write("\\t".join(header) + "\\n")
        for tid in sorted(transcript_ids):
            row = predictions[tid]
            tmap = tmap_data.get(tid, {})
            out.write("\\t".join([
                tid,
                tmap.get("class_code", "NA"),
                tmap.get("ref_gene_id", "NA"),
                tmap.get("ref_id", "NA"),
                row["consensus"],
                str(row["votes"]),
                row["cpat"],
                row["feelnc"],
                row["plek"],
                str(row["cpat_score"]),
                str(row["feelnc_score"]),
                str(row["plek_score"]),
            ]) + "\\n")


def write_gtf(gtf_path, out_path, keep_ids):
    with open(gtf_path) as src, open(out_path, "w") as out:
        for line in src:
            if line.startswith("#"):
                continue
            if any(f'transcript_id "{tid}"' in line for tid in keep_ids):
                out.write(line)


def write_fasta(fasta_path, out_path, keep_ids):
    with open(out_path, "w") as out:
        for record in SeqIO.parse(fasta_path, "fasta"):
            if record.id in keep_ids:
                SeqIO.write(record, out, "fasta")


def main():
    cpat = read_cpat("$cpat_results")
    feelnc = read_simple_prediction("$feelnc_results", {"lncRNA", "Non-coding"})
    plek = read_simple_prediction("$plek_results", {"Non-coding"})
    tmap = read_tmap("$tmap")

    mode = "$task.ext.args" if "$task.ext.args" != "null" else "majority"
    if mode not in {"strict", "majority", "lenient"}:
        mode = "majority"

    prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "${meta4.id}"
    transcript_ids = set(cpat) | set(feelnc) | set(plek)
    predictions = {}

    for tid in transcript_ids:
        votes = []
        for source in (cpat, feelnc, plek):
            if tid in source:
                votes.append(source[tid][0])

        predictions[tid] = {
            "consensus": consensus(votes, mode),
            "votes": len(votes),
            "cpat": cpat.get(tid, ["NA", "NA"])[0],
            "feelnc": feelnc.get(tid, ["NA", "NA"])[0],
            "plek": plek.get(tid, ["NA", "NA"])[0],
            "cpat_score": cpat.get(tid, ["NA", "NA"])[1],
            "feelnc_score": feelnc.get(tid, ["NA", "NA"])[1],
            "plek_score": plek.get(tid, ["NA", "NA"])[1],
        }

    protein_coding_ids = {tid for tid, row in predictions.items() if row["consensus"] == "Protein-coding"}

    write_summary(f"{prefix}.prediction_summary.tsv", transcript_ids, predictions, tmap)
    write_summary(f"{prefix}.prediction_summary_protein_coding.tsv", protein_coding_ids, predictions, tmap)
    write_gtf("$gtf", f"{prefix}.final_protein_coding.gtf", protein_coding_ids)
    write_fasta("$fasta", f"{prefix}.final_protein_coding.fa", protein_coding_ids)

    noncoding_ids = transcript_ids - protein_coding_ids
    with open(f"{prefix}.coding_potential_report.txt", "w") as out:
        out.write(f"Consensus mode: {mode}\\n")
        out.write(f"CPAT cutoff: {CPAT_CUTOFF}\\n")
        out.write(f"Total transcripts analyzed: {len(transcript_ids)}\\n")
        out.write(f"Predicted protein-coding transcripts: {len(protein_coding_ids)}\\n")
        out.write(f"Predicted noncoding/no-consensus transcripts: {len(noncoding_ids)}\\n")

    with open("versions.yml", "w") as out:
        out.write("${task.process}:\\n")
        out.write(f"    python: {platform.python_version()}\\n")


if __name__ == "__main__":
    main()
