# Coding Potential Test Dataset

Synthetic dataset for validating the novel-transcript coding-potential branch.

Files:
- `novel_transcripts.fa`, `novel_transcripts.gtf`, `novel_transcripts.tmap`: 520 candidate novel transcripts.
- `cds.fa`: 520 coding sequences for CPAT positive training.
- `noncoding.fa`: 520 noncoding sequences for CPAT negative training.
- `mrna.fa`: 520 mRNA sequences for FEELnc coding reference.
- `genome.fa`, `reference.gtf`: matched synthetic reference files for testing the extraction modules.

Regenerate with:

```bash
python3 tests/data/coding_potential/generate_dataset.py
```
