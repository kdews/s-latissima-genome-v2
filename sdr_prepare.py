#!/usr/bin/env python3
"""Prepare SDR sample manifests and YGS read lists from genome metadata.

Run from the analysis directory containing assemblies/ and read_lists/.
Called by sdr_id.sbatch before any jobs are submitted.
"""

import argparse
import csv
import re
from pathlib import Path


def fail(message):
    raise SystemExit("Error: " + message)


def prepare(metadata, out, ygs_outdir, depth_outdir, plan_filename="sdr_plan.tsv"):
    """Validate metadata and read pairs, then write manifests and a job plan."""
    out.mkdir(parents=True, exist_ok=True)
    with metadata.open() as fh:
        rows = csv.reader(fh, delimiter='\t')
        header = next(rows)
        header[0] = header[0].lstrip('# ').lstrip('\ufeff')
        if len(set(header)) != len(header): fail('duplicate metadata column')
        required = ['Label', 'Sample_ID', 'Region', 'Sex', 'Ploidy',
                    'Analysis_Reads_Platform', 'Analysis_Reads_Method']
        if set(required) - set(header): fail('missing columns: ' + str(set(required)-set(header)))
        records = []
        for row in rows:
            if not row or row[0].startswith('#'): continue
            if len(row) != len(header): fail(f'expected {len(header)} columns: {row[0]}')
            records.append(dict(zip(header, row)))
    if not records: fail('metadata contains no assemblies')
    seen_labels = set()
    plan = []
    for row in records:
        label = row['Label']
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9._-]*', label): fail(f'invalid label: {label}')
        if label in seen_labels: fail(f'duplicate label: {label}')
        seen_labels.add(label)
        platform, method = row['Analysis_Reads_Platform'], row['Analysis_Reads_Method']
        if not platform and not method:
            print(f'{label}: no analysis reads configured; skipping.')
            continue
        if platform != 'illumina' or method not in ('script', 'sra'):
            fail(f'{label}: expected illumina and script/sra analysis-read method')
        sex = row['Sex'].upper()
        if sex not in ('U', 'V', 'UV'): fail(f'{label}: unsupported Sex={sex}')
        if not row['Sample_ID']: fail(f'{label}: missing Sample_ID')
        genome = Path('assemblies', label + '.fa').resolve()
        if not genome.is_file() or not genome.stat().st_size: fail(f'missing assembly: {genome}')
        samples, pairs, lists, seen_paths = {}, [], {}, set()
        for read_sex in ('female', 'male'):
            source = Path('read_lists', f'{label}_{read_sex}.list')
            if not source.is_file(): fail(f'missing read list: {source}')
            mates = {}
            paths = []
            for raw in source.read_text().splitlines():
                if not raw.strip() or raw.lstrip().startswith('#'): continue
                path = Path(raw).expanduser().resolve()  # relative to analysis working directory
                if not path.is_file() or not path.stat().st_size: fail(f'missing/empty read: {path}')
                if str(path) in seen_paths: fail(f'{label}: duplicate read: {path}')
                seen_paths.add(str(path))
                if any(c in str(path) for c in '\t\n\r'): fail(f'unsupported filename: {path}')
                name = path.name
                match = re.fullmatch(r'(.+)_R([12])((?:_.*)?\.(?:fastq|fq)(?:\.gz)?)', name)
                if not match: match = re.fullmatch(r'(.+)_([12])(\.(?:fastq|fq)(?:\.gz)?)', name)
                if not match: fail(f'cannot identify PE mate in {name}')
                # Include directory and suffix to distinguish lanes/libraries with the same prefix
                key = (str(path.parent), match[1], match[3])
                mate = int(match[2])
                if mate in mates.setdefault(key, {}): fail(f'duplicate mate: {path}')
                mates[key][mate] = str(path)
                paths.append(str(path))
            if not paths: fail(f'empty read list: {source}')
            lists[read_sex] = paths
            for key, mate_paths in mates.items():
                if set(mate_paths) != {1, 2}: fail(f'incomplete pair: {key}')
                sample = key[1].split('_')[0]
                if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9.-]*', sample): fail(f'unsafe sample ID: {sample}')
                if sample in samples and samples[sample] != read_sex: fail(f'sample in both sexes: {sample}')
                samples[sample] = read_sex
                # Local filename prefix: Sample_Library_RunID_Barcode_Plate_Lane
                fields = key[1].split('_')
                # For remote download, uses SRA ID
                library = fields[1] if len(fields) > 1 else sample
                pair_id = f'pair{len(pairs)+1:04d}'

                pairs.append([sample, read_sex, library, pair_id,
                              mate_paths[1], mate_paths[2]])
        depth = out / depth_outdir / label
        depth.mkdir(parents=True, exist_ok=True)
        for filename, header_fields, data in [
            ('samples.tsv', ['Sample', 'Sex'], list(samples.items())),
            ('read_pairs.tsv', 
             ['Sample', 'Sex', 'Library', 'Pair_ID', 'R1', 'R2'],
             pairs),
            ('assembly_metadata.tsv', header, [[row[col] for col in header]])]:
            with (depth / filename).open('w') as fh:
                writer = csv.writer(fh, delimiter='\t', lineterminator='\n')
                writer.writerow(header_fields); writer.writerows(data)
        for target in (['U', 'V'] if sex == 'UV' else [sex]):
            ygs = out / ygs_outdir / label / ('target_' + target)
            ygs.mkdir(parents=True, exist_ok=True)
            same = 'female' if target == 'U' else 'male'
            opposite = 'male' if target == 'U' else 'female'
            for filename, read_sex in [('ss.txt', same), ('os.txt', opposite)]:
                # Do not silently reuse cached k-mers for changed lists.
                dest = ygs / filename
                content = '\n'.join(lists[read_sex]) + '\n'
                if dest.exists() and dest.read_text() != content and any(ygs.glob('*.jf')):
                    fail(f'{ygs}: read list changed; move the previous direction directory before rerunning')
                dest.write_text(content)
            (ygs / 'assembly_metadata.tsv').write_text((depth / 'assembly_metadata.tsv').read_text())
        plan.append([label, str(genome), sex, str(len(samples))])
        print(f'{label}: Sample_ID={row["Sample_ID"]}; Region={row["Region"]}; Sex={sex}; '
              f'Ploidy={row["Ploidy"]}; {len(samples)} samples, {len(pairs)} read pairs')
    with (out / plan_filename).open('w') as fh:
        csv.writer(fh, delimiter='\t', lineterminator='\n').writerows(plan)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("metadata", type=Path, help="Genome metadata TSV")
    parser.add_argument("outdir", type=Path, help="Main SDR output directory")
    parser.add_argument("ygs_outdir", type=Path, help="YGS output subdirectory")
    parser.add_argument("depth_outdir", type=Path, help="Depth output subdirectory")
    parser.add_argument("plan_filename", nargs="?", default="sdr_plan.tsv",
                        help="Output plan filename (default: sdr_plan.tsv)")
    args = parser.parse_args()
    prepare(args.metadata, args.outdir, args.ygs_outdir,
            args.depth_outdir, args.plan_filename)


if __name__ == "__main__":
    main()
