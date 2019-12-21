import sys
import json

try:
    from aars_algorithms_fast import levenshtein_distance_c as levenshtein_distance
    from aars_algorithms_fast import align_c as align
except ImportError:
    print('using slow algorithms - this could take a while')
    from aars_algorithms_slow import levenshtein_distance_py as levenshtein_distance
    from aars_algorithms_slow import align_py as align

AA_LIST = ['gln', 'leu', 'glu', 'ile',
           'arg', 'met', 'val', 'cys',
           'trp', 'tyr']
KINGDOM_LIST = ['bact']


def main(filepath):
    with open(filepath + ".json") as f:
        db = json.load(f)
    fasta_dat = parse_fasta(filepath)
    for r1 in fasta_dat:
        for r2 in db:
            if r1['aa'] == r2['aa']:
                if r1['kingdom'] == r2['kingdom']:
                    if r1['pdb'] == r2['pdb']:
                        if r1['letter'] == r2['letter']:
                            if r1['genus'] == r2['genus']:
                                r1['ungapped'] = r2['aa_dat']
                                r1['aligned'] = align(
                                    r1['gapped'], r1['ungapped']
                                )

                                r1['nucleotide'] = r2['nuc_dat']
                                break
        else:
            print(r1)
            raise RuntimeError("no match found")
    with open(filepath + '.error_report.md', 'w') as err_fp:
        errs = 0
        for r in fasta_dat:
            n, n_rev = validate_translation(r['nucleotide'], r['ungapped'])
            if n >= 0.95:
                continue
            elif n < 0.95 and n_rev >= 0.95:
                errs += 1
                record_reverse_complement_translation(r, n, err_fp)
            else:
                errs += 1
                record_bad_translation(r, n, n_rev, err_fp)
        for r in fasta_dat:
            n = count_misalignments(r['ungapped'], r['aligned'])
            if n != 0:
                errs += 1
                record_bad_alignment(r, n, err_fp)
        if errs != 0:
            raise RuntimeError("{:d} Errors - see error report".format(errs))
        print(json.dumps(fasta_dat, indent=2))


def record_bad_translation(r, n, n_rev, fp):
    print('## {}\n'.format(pretty(r)), file=fp)
    print('This record seems to be incorrectly translated from', file=fp)
    print('nucleotides to amino acids. Variation up to 5.00% is', file=fp)
    print('tolerated by the script but the variation in this', file=fp)
    print('translation is {:.2f}%.\n'.format(100 - n*100), file=fp)
    print('Additionally, taking the reverse complement', file=fp)
    print('of the nucleotide sequence still results in a variation', file=fp)
    print('of {:.2f}%.\n'.format(100 - n_rev*100), file=fp)
    print_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_translation_errors(r, fp)


def record_reverse_complement_translation(r, n, fp):
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The translation of the nucleotides to', file=fp)
    print('amino acids seems to be the reverse complement.\n', file=fp)
    print_nucleotide_file_data(r, fp)
    print_reverse_complement_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_reverse_complement_translation_errors(r, fp)


def record_bad_alignment(r, n, fp):
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The amino acids provided in the FASTA alignment data', file=fp)
    print('do not align with amino acids in the source data file.\n', file=fp)
    print('There ', file=fp, end='')
    if n == 1:
        print('is 1 amino acid that does not align.\n', file=fp)
    else:
        print('are {:d} amino acids that do not align.\n'.format(n), file=fp)
    print('This could indicate a problem with the processing used to', file=fp)
    print('generate the data in the FASTA file.\n', file=fp)
    print_amino_acid_file_data(r, fp)
    print_alignment_error(r, fp)


def print_nucleotide_file_data(r, fp):
    print('### Original nucleotide data\n', file=fp)
    print('```text', file=fp)
    with open('data/{}.nuc'.format(pretty(r))) as f:
        for line in f:
            if line.strip() != '':
                print(line, file=fp, end='')
    print('\n```\n', file=fp)


def print_reverse_complement_nucleotide_file_data(r, fp, line_length=60):
    rev = reverse_seq(r['nucleotide'])
    print('### Reverse complement nucleotide data\n', file=fp)
    print('```text', file=fp)
    i = 0
    while rev[i:i+line_length] != "":
        print(rev[i:i+line_length], file=fp)
        i = i + line_length
    print('\n```\n', file=fp)


def print_amino_acid_file_data(r, fp):
    print('### Original amino acid data\n', file=fp)
    print('```text', file=fp)
    with open('data/{}.aa'.format(pretty(r))) as f:
        for line in f:
            if line.strip() != '':
                print(line, file=fp, end='')
    print('\n```\n', file=fp)


def print_translation_errors(r, fp, line_length=60):
    print('### Translation errors\n', file=fp)
    aa_trans = translate(r['nucleotide'])
    if len(r['ungapped']) != 0:
        bad = ['X' if x1 != x2 else ' ' for (
            x1, x2) in zip(r['ungapped'], aa_trans)]
        print('```text', file=fp)
        i = 0
        while bad[i:i+(line_length // 3)] != []:
            aas = ' ' + '  '.join(r['ungapped'][i:i+(line_length // 3)]) + ' '
            tns = ' ' + '  '.join(aa_trans[i:i+(line_length // 3)]) + ' '
            err = ' ' + '  '.join(bad[i:i+(line_length // 3)]) + ' '
            print("nuc: {}".format(r['nucleotide'][i:i+line_length]), file=fp)
            print("aas: {}".format(aas), file=fp)
            print("tns: {}".format(tns), file=fp)
            print("err: {}\n".format(err), file=fp)
            i = i + line_length
        print('```\n', file=fp)


def print_reverse_complement_translation_errors(r, fp, line_length=60):
    print('### Reverse complement translation errors\n', file=fp)
    rev = reverse_seq(r['nucleotide'])
    aa_trans = translate(rev)
    if len(r['ungapped']) != 0:
        bad = ['X' if x1 != x2 else ' ' for (
            x1, x2) in zip(r['ungapped'], aa_trans)]
        print('```text', file=fp)
        i = 0
        while bad[i:i+(line_length // 3)] != []:
            aas = ' ' + '  '.join(r['ungapped'][i:i+(line_length // 3)]) + ' '
            tns = ' ' + '  '.join(aa_trans[i:i+(line_length // 3)]) + ' '
            err = ' ' + '  '.join(bad[i:i+(line_length // 3)]) + ' '
            print("rev: {}".format(rev[i:i+line_length]), file=fp)
            print("aas: {}".format(aas), file=fp)
            print("tns: {}".format(tns), file=fp)
            print("err: {}\n".format(err), file=fp)
            i = i + line_length
        print('```\n', file=fp)


def print_alignment_error(r, fp, line_length=60):
    print('### Alignment errors\n', file=fp)
    bad = ['X' if c1 not in '-.*?' and c2 not in '*?' and c1 != c2 else ' '
           for c1, c2 in zip(r['aligned'], r['ungapped'])]
    print('```text', file=fp)
    i = 0
    while bad[i:i+line_length] != []:
        print("aas: {}".format(r['ungapped'][i:i+line_length]), file=fp)
        print("ali: {}".format(r['aligned'][i:i+line_length]), file=fp)
        print("err: {}\n".format(''.join(bad[i:i+line_length])), file=fp)
        i = i + line_length
    print('```\n', file=fp)


def parse_fasta(path):
    fasta_data = []
    gapped_sequence = None
    with open(path) as lines:
        for line in lines:
            if line[0] == '>':
                # store previous gapped sequence
                if gapped_sequence:
                    fasta_data[-1]['gapped'] = gapped_sequence
                gapped_sequence = ""
                xs = line[1:].strip().split('_')
                if not (xs and xs[0] and xs[0].lower() in AA_LIST):
                    raise RuntimeError(
                        "Amino Acid ({}) not recognized".format(
                            xs[0]
                        )
                    )
                aa = xs.pop(0).lower()
                if not (xs and xs[0] and xs[0].lower() in KINGDOM_LIST):
                    raise RuntimeError(
                        "Kingdom ({}) not recognized".format(
                            xs[0]
                        )
                    )
                kingdom = xs.pop(0).lower()
                pdb = None
                if xs and xs[0] and len(xs[0]) == 4:
                    pdb = xs.pop(0).lower()
                if not (xs and xs[0] and len(xs[0]) == 1):
                    raise RuntimeError(
                        "'{}' not recognized as first letter of genus".format(
                            xs[0]
                        )
                    )
                letter = xs.pop(0).upper()
                genus, num = '_'.join(xs).lower().split('/')
                fasta_data.append({
                    'aa': aa,
                    'kingdom': kingdom,
                    'pdb': pdb,
                    'letter': letter,
                    'genus': genus,
                    'num': num
                })
            else:
                gapped_sequence += line.strip()
    fasta_data[-1]['gapped'] = gapped_sequence
    return fasta_data


def count_misalignments(original, aligned):
    return sum([1 if c2 not in '-.*?' and c1 not in '*?' and c1 != c2 else 0
                for c1, c2 in zip(original, aligned)])


def pretty(r):
    return '_'.join(
        [x for x in [
            r['aa'],
            r['kingdom'],
            r['pdb'],
            r['letter'],
            r['genus'],
            # r['num']
        ] if x]
    )


def translate(nuc_seq):
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    aa_trans = ""
    if len(nuc_seq) % 3 == 0:
        for i in range(0, len(nuc_seq), 3):
            codon = nuc_seq[i:i + 3]
            try:
                aa_trans += table[codon]
            except KeyError:
                aa_trans += '?'
    return aa_trans


def validate_translation(nuc_seq, aa_seq):
    aa_trans = translate(nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    if len(aa_seq) == 0:
        return 0, 0
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n = (max(0, sum(x) - bad)) / len(aa_seq)

    # check reverse complement
    rev_nuc_seq = reverse_seq(nuc_seq)
    aa_trans = translate(rev_nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n_rev = (max(0, sum(x) - bad)) / len(aa_seq)
    return n, n_rev


def reverse_seq(nuc_seq):
    return "".join([alt_nuc(i) for i in reversed(nuc_seq)])


def alt_nuc(i):
    if i == 'A':
        return 'T'
    if i == 'T':
        return 'A'
    if i == 'C':
        return 'G'
    if i == 'G':
        return 'C'
    if i == 'a':
        return 't'
    if i == 't':
        return 'a'
    if i == 'c':
        return 'g'
    if i == 'g':
        return 'c'
    return i


if __name__ == "__main__":
    main(sys.argv[1])