#!/usr/bin/env python3

import platform
import re
from collections import Counter


TRANSCRIPT_ID_RE = re.compile(r'transcript_id "([^"]+)"')


def parse_codes(codes):
    return {code.strip() for code in codes.split(",") if code.strip()}


def find_column(header, candidates, fallback):
    normalized = [field.lstrip("#") for field in header]
    for candidate in candidates:
        if candidate in normalized:
            return normalized.index(candidate)
    return fallback


def read_selected_transcripts(tmap_path, selected_codes, filtered_tmap_path):
    selected_ids = set()
    class_counts = Counter()
    selected_counts = Counter()

    with open(tmap_path) as tmap_in, open(filtered_tmap_path, "w") as tmap_out:
        header_line = tmap_in.readline()
        if not header_line:
            return selected_ids, class_counts, selected_counts

        tmap_out.write(header_line)
        header = header_line.rstrip("\\n").split("\\t")
        class_idx = find_column(header, ["class_code", "class"], 2)
        qry_idx = find_column(header, ["qry_id", "query_id", "qry_transcript_id"], 4)

        for line in tmap_in:
            if not line.strip():
                continue

            fields = line.rstrip("\\n").split("\\t")
            if len(fields) <= max(class_idx, qry_idx):
                continue

            class_code = fields[class_idx]
            query_id = fields[qry_idx]
            class_counts[class_code] += 1

            if class_code in selected_codes:
                selected_ids.add(query_id)
                selected_counts[class_code] += 1
                tmap_out.write(line)

    return selected_ids, class_counts, selected_counts


def filter_gtf(gtf_path, selected_ids, filtered_gtf_path):
    kept_lines = 0
    with open(gtf_path) as gtf_in, open(filtered_gtf_path, "w") as gtf_out:
        for line in gtf_in:
            if line.startswith("#"):
                #gtf_out.write(line)
                continue

            match = TRANSCRIPT_ID_RE.search(line)
            if match and match.group(1) in selected_ids:
                gtf_out.write(line)
                kept_lines += 1

    return kept_lines


def write_summary(summary_path, class_counts, selected_counts, selected_codes, selected_ids, kept_lines):
    with open(summary_path, "w") as out:
        out.write("metric\\tvalue\\n")
        out.write(f"selected_class_codes\\t{','.join(sorted(selected_codes))}\\n")
        out.write(f"selected_transcripts\t{len(selected_ids)}\\n")
        out.write(f"gtf_lines_kept\\t{kept_lines}\\n")
        out.write("\\nclass_code\\ttotal_transcripts\\tselected_transcripts\\n")
        for class_code in sorted(class_counts):
            out.write(
                f"{class_code}\\t{class_counts[class_code]}\\t{selected_counts.get(class_code, 0)}\\n"
            )


def write_versions():
    with open("versions.yml", "w") as out:
        out.write("${task.process}:\\n")
        out.write(f"    python: {platform.python_version()}\\n")


def run_filter(gtf, tmap, class_codes, prefix):
    selected_codes = parse_codes(class_codes)
    filtered_gtf = f"{prefix}.novel.gtf"
    filtered_tmap = f"{prefix}.novel.tmap"
    summary = f"{prefix}.class_summary.tsv"

    selected_ids, class_counts, selected_counts = read_selected_transcripts(
        tmap, selected_codes, filtered_tmap
    )
    kept_lines = filter_gtf(gtf, selected_ids, filtered_gtf)
    write_summary(summary, class_counts, selected_counts, selected_codes, selected_ids, kept_lines)
    write_versions()


if __name__ == "__main__":
    run_filter("$gtf", "$tmap", "$class_codes", "${prefix}")
