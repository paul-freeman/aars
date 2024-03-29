import sys
import os
import json
import os.path
import argparse

import preprocessing

try:
    from aars_algorithms_fast import levenshtein_distance_c as levenshtein_distance
    from aars_algorithms_fast import align_c as align
except ImportError as e:
    print(e)
    print('using slow algorithms - this could take a while')
    from aars_algorithms_slow import levenshtein_distance_py as levenshtein_distance
    from aars_algorithms_slow import align_py as align

AA_LIST = ['ala', 'asn', 'asp', 'gln', 'leu', 'glu', 'gly', 'his', 'ile', 'lys',
           'arg', 'met', 'phe', 'pro', 'pyl', 'sep', 'ser', 'thr', 'val', 'cys', 'trp', 'tyr']
KINGDOM_LIST = ['bact', 'arch']


def main(filepath, regions, reverse_complement_errors=False):
    # parse regions
    reg = json.loads(regions)

    with open(filepath + ".json") as f:
        db = json.load(f)
    fasta_dat = parse_fasta(filepath, reg)
    accounted_fasta_dat = []
    missing_data_files = []
    for r1 in fasta_dat:
        for r2 in db:
            if r1['aa'] == r2['aa']:
                if r1['kingdom'] == r2['kingdom']:
                    if r1['regions'] == r2['regions']:
                        if r1['pdb'] == r2['pdb']:
                            if r1['letter'] == r2['letter']:
                                if r1['genus'] == r2['genus']:
                                    r1['ungapped'] = r2['aa_dat']
                                    r1['aligned'] = align(
                                        r1['gapped'], r1['ungapped']
                                    )
                                    r1['nucleotide'] = r2['nuc_dat']
                                    accounted_fasta_dat.append(r1)
                                    break
        else:
            missing_data_files.append(preprocessing.make_filename(r1))
            raise RuntimeError("no match found")
    fasta_dat = accounted_fasta_dat
    properly_translated_dat = []
    good_fasta_dat = []
    with open(filepath + '.error_report.md', 'w') as err_fp:
        with open(filepath + '.error_report_summary.txt', 'w') as err_fp_summary:
            for missing in missing_data_files:
                print('{} - missing data'.format(missing), file=err_fp_summary)
            errs = 0
            for r in fasta_dat:
                n, n_rev, n_comp, n_rev_comp, n_rev_comp_shift_left, n_rev_comp_shift_right = validate_translation(
                    r['nucleotide'], r['ungapped'])
                if n >= 0.95:
                    pass
                elif n_rev >= 0.95:
                    if reverse_complement_errors:
                        errs += 1
                        record_reverse_translation(
                            r, n, err_fp, err_fp_summary)
                elif n_comp >= 0.95:
                    r['nucleotide'] = complement(r['nucleotide'])
                    if reverse_complement_errors:
                        errs += 1
                        record_complement_translation(
                            r, n, err_fp, err_fp_summary)
                elif n_rev_comp >= 0.95:
                    r['nucleotide'] = reverse_complement(r['nucleotide'])
                    if reverse_complement_errors:
                        errs += 1
                        record_reverse_complement_translation(
                            r, n, err_fp, err_fp_summary)
                elif n_rev_comp_shift_left >= 0.95:
                    r['nucleotide'] = reverse_complement_shift_left(
                        r['nucleotide'])
                    if reverse_complement_errors:
                        errs += 1
                        record_reverse_complement_shift_left_translation(
                            r, n, err_fp, err_fp_summary)
                elif n_rev_comp_shift_right >= 0.95:
                    r['nucleotide'] = reverse_complement_shift_right(
                        r['nucleotide'])
                    if reverse_complement_errors:
                        errs += 1
                        record_reverse_complement_shift_right_translation(
                            r, n, err_fp, err_fp_summary)
                else:
                    errs += 1
                    record_bad_translation(
                        r, n, n_rev, n_comp, n_rev_comp, err_fp, err_fp_summary)
                    continue
                properly_translated_dat.append(r)
            for r in properly_translated_dat:
                n = count_misalignments(r['ungapped'], r['aligned'])
                if n > 0:
                    errs += 1
                    record_bad_alignment(r, n, err_fp, err_fp_summary)
                else:
                    good_fasta_dat.append(r)
            if errs != 0:
                print("{:d} Errors - see error report".format(errs))
    if errs == 0:
        os.remove(filepath + '.error_report.md')
        os.remove(filepath + '.error_report_summary.txt')
    # write_master_files(good_fasta_dat)
    # write_aars_regions(good_fasta_dat, filepath)
    # write_middle_base_regions(good_fasta_dat, filepath)
    write_middle_bases_with_numbers(good_fasta_dat, filepath)
    write_middle_bases_no_numbers(good_fasta_dat, filepath)
    write_whole_codon_with_numbers(good_fasta_dat, filepath)
    write_whole_codon_no_numbers(good_fasta_dat, filepath)
    write_amino_acids_with_numbers(good_fasta_dat, filepath)
    write_amino_acids_no_numbers(good_fasta_dat, filepath)
    # with open(filepath + '.csv', 'w') as csv_file:
    #     write_csv(fasta_dat, csv_file)
    os.remove(filepath + ".json")


def write_master_files(fasta_dat):
    dropbox_masters = os.path.join(os.path.expanduser(
        '~'), 'UniDropbox\\Dropbox\\With Paul Freeman\\master\\temp')
    try:
        os.makedirs(dropbox_masters)
    except FileExistsError:
        pass
    for r in fasta_dat:
        if r['regions']:
            continue
        filename = preprocessing.make_filename(r)
        for r2 in fasta_dat:
            if r['aa'] == r2['aa']:
                if r['kingdom'] == r2['kingdom']:
                    if r2['regions'] == True:
                        if r['pdb'] == r2['pdb']:
                            if r['letter'] == r2['letter']:
                                if r['genus'] == r2['genus']:
                                    regions = r2
                                    break
        else:
            regions = None
        with open(os.path.join(dropbox_masters, filename + '_master.txt'), 'w') as f:
            i = 0
            while i < len(r['nucleotide']):
                bad = ['X' if c1 not in '-.*?' and c2 not in '*?' and c1 != c2 else ' '
                       for c1, c2 in zip(r['aligned'], r['ungapped'])]
                if regions:
                    bad_regions = ['X' if c1 not in '-.*?' and c2 not in '*?' and c1 != c2 else ' '
                                   for c1, c2 in zip(regions['aligned'], regions['ungapped'])]
                print(''.join(['{:^3s}'.format(str((n+1) % 1000))
                               for n in range(i//3, (i+78)//3)]), file=f)
                print(r['nucleotide'][i:i+78], file=f)
                print(''.join(['{:^3s}'.format(n)
                               for n in r['ungapped'][i//3:(i+78)//3]]), file=f)
                print(''.join(['{:^3s}'.format(n)
                               for n in r['aligned'][i//3:(i+78)//3]]), file=f)
                print(''.join(['{:^3s}'.format(n)
                               for n in bad[i//3:(i+78)//3]]), file=f)
                if regions:
                    print(''.join(['{:^3s}'.format(n)
                                   for n in regions['aligned'][i//3:(i+78)//3]]), file=f)
                    print(''.join(['{:^3s}'.format(n)
                                   for n in bad_regions[i//3:(i+78)//3]]), file=f)
                print('', file=f)
                i += 78


def write_middle_base_regions(fasta_dat, filepath, max_width=10000):
    filepath = os.path.splitext(filepath)[0]
    with open(filepath + "_middle_base_regions.txt", 'w') as f:
        for r in fasta_dat:
            if r['regions']:
                continue
            filename = preprocessing.make_filename(r)
            for r2 in fasta_dat:
                if r['aa'] == r2['aa']:
                    if r['kingdom'] == r2['kingdom']:
                        if r2['regions'] == True:
                            if r['pdb'] == r2['pdb']:
                                if r['letter'] == r2['letter']:
                                    if r['genus'] == r2['genus']:
                                        regions = r2
                                        break
            else:
                err = 'could not find regions while writing middle base regions'.format(
                    filename)
                print(err)
                return
            middle_base = "".join([x if i % 3 == 1 else ' ' for i,
                                   x in enumerate(r['nucleotide'])])
            print("> {}".format(filename), file=f)
            i = 0
            while i < len(middle_base):
                aa_nums = ''.join(['{:^3s}'.format(str((n+1) % 1000))
                                   for n in range(i//3, (i+max_width)//3)])
                middles = middle_base[i:i+max_width]
                regs = ''.join(['{:^3s}'.format(n)
                                for n in regions['aligned'][i//3:(i+max_width)//3]])
                i += max_width
            aa_nums_final = ""
            middles_final = ""
            regs_final = ""
            i = 0
            gap = False
            while i < len(middle_base):
                if regs[i+1] == '-' and gap == False:
                    aa_nums_final += ' - '
                    middles_final += ' - '
                    regs_final += ' - '
                    gap = True
                if regs[i+1] != '-':
                    aa_nums_final += aa_nums[i:i+3]
                    middles_final += middles[i:i+3]
                    regs_final += regs[i:i+3]
                    gap = False
                i += 3
            if aa_nums_final[0:3] == ' - ':
                aa_nums_final = aa_nums_final[3:]
            if middles_final[0:3] == ' - ':
                middles_final = middles_final[3:]
            if aa_nums_final[-3:] == ' - ':
                aa_nums_final = aa_nums_final[:-3]
            if middles_final[-3:] == ' - ':
                middles_final = middles_final[:-3]

            print(aa_nums_final, file=f)
            print(middles_final, file=f)
            # print(regs_final, file=f)
            print('', file=f)


def write_middle_bases_no_numbers(fasta_dat, filepath):
    write_result_data(
        fasta_dat,
        filepath,
        "_middle_bases_no_numbers.txt",
        printMiddleBases=True
    )


def write_middle_bases_with_numbers(fasta_dat, filepath):
    write_result_data(
        fasta_dat,
        filepath,
        "_middle_bases_with_numbers.txt",
        printNumbers=True,
        printMiddleBases=True
    )


def write_whole_codon_no_numbers(fasta_dat, filepath):
    write_result_data(
        fasta_dat,
        filepath,
        "_whole_codon_no_numbers.txt",
        printCodons=True
    )


def write_whole_codon_with_numbers(fasta_dat, filepath):
    write_result_data(
        fasta_dat,
        filepath,
        "_whole_codon_with_numbers.txt",
        printNumbers=True,
        printCodons=True
    )


def write_amino_acids_no_numbers(fasta_dat, filepath):
    write_result_data(
        fasta_dat,
        filepath,
        "_amino_acids_no_numbers.txt",
        printAminoAcids=True
    )


def write_amino_acids_with_numbers(fasta_dat, filepath):
    write_result_data(
        fasta_dat,
        filepath,
        "_amino_acids_with_numbers.txt",
        printNumbers=True,
        printAminoAcids=True
    )


def write_result_data(fasta_dat, filepath, outputSuffix, printNumbers=False, printCodons=False, printMiddleBases=False, printAminoAcids=False):
    filepath = os.path.splitext(filepath)[0]
    with open(filepath + outputSuffix, 'w') as f:
        for r in fasta_dat:

            # print fasta name
            filename = preprocessing.make_filename(r)
            print("> {}".format(filename), file=f)

            # get middle bases of nucleotide
            codons = []
            middle_bases = []
            for i, x in enumerate(r['nucleotide']):
                if i % 3 == 1:
                    codons.append(r['nucleotide'][i-1:i+2])
                    middle_bases.append(x)

            # check lengths
            if len(middle_bases) != len(r['aligned']):
                raise RuntimeError("middle bases have length {}, but aligned amino acids have length {}".format(
                    len(middle_bases), len(r['aligned'])))

            # get amino acid
            xs = []
            for i, x in enumerate(middle_bases):
                xs.append({
                    "num": i+1,
                    "codon": codons[i],
                    "middle_base": x,
                    "aa": r['aligned'][i]
                })

            final_numbers = ''
            final_codons = ''
            final_middle_bases = ''
            final_amino_acids = ''
            i = 0
            for p in r['underscored']:
                if p == '_' or p == '-':
                    final_numbers += ' - '
                    final_codons += ' - '
                    final_middle_bases += ' - '
                    final_amino_acids += ' - '
                    continue
                while xs[i]['aa'] != p:
                    i += 1
                final_numbers += '{:^3s}'.format(str(xs[i]['num']))
                final_codons += '{:^3s}'.format(xs[i]['codon'])
                final_middle_bases += '{:^3s}'.format(xs[i]['middle_base'])
                final_amino_acids += '{:^3s}'.format(p)
            if printNumbers:
                print(final_numbers, file=f)
            if printCodons:
                print(final_codons, file=f)
            if printMiddleBases:
                print(final_middle_bases, file=f)
            if printAminoAcids:
                print(final_amino_acids, file=f)


def write_aars_regions(fasta_dat, filepath, max_width=10000):
    filepath = os.path.splitext(filepath)[0]
    with open(filepath + "_aars_regions.txt", 'w') as f:
        for r in fasta_dat:
            if r['regions']:
                continue
            filename = preprocessing.make_filename(r)
            for r2 in fasta_dat:
                if r['aa'] == r2['aa']:
                    if r['kingdom'] == r2['kingdom']:
                        if r2['regions'] == True:
                            if r['pdb'] == r2['pdb']:
                                if r['letter'] == r2['letter']:
                                    if r['genus'] == r2['genus']:
                                        regions = r2
                                        break
            else:
                err = 'could not find regions while writing aars regions'.format(
                    filename)
                print(err)
                return
            print("> {}".format(filename), file=f)
            i = 0
            while i < len(r['nucleotide']):
                aa_nums = ''.join(['{:^3s}'.format(str((n+1) % 1000))
                                   for n in range(i//3, (i+max_width)//3)])
                nucs = r['nucleotide'][i:i+max_width]
                regs = ''.join(['{:^3s}'.format(n)
                                for n in regions['aligned'][i//3:(i+max_width)//3]])
                i += max_width
            aa_nums_final = ""
            nuc_final = ""
            regs_final = ""
            i = 0
            gap = False
            while i < len(r['nucleotide']):
                if regs[i+1] == '-' and gap == False:
                    aa_nums_final += ' - '
                    nuc_final += ' - '
                    regs_final += ' - '
                    gap = True
                if regs[i+1] != '-':
                    aa_nums_final += aa_nums[i:i+3]
                    nuc_final += nucs[i:i+3]
                    regs_final += regs[i:i+3]
                    gap = False
                i += 3
            if aa_nums_final[0:3] == ' - ':
                aa_nums_final = aa_nums_final[3:]
            if nuc_final[0:3] == ' - ':
                nuc_final = nuc_final[3:]
            if regs_final[0:3] == ' - ':
                regs_final = regs_final[3:]
            if aa_nums_final[-3:] == ' - ':
                aa_nums_final = aa_nums_final[:-3]
            if nuc_final[-3:] == ' - ':
                nuc_final = nuc_final[:-3]
            if regs_final[-3:] == ' - ':
                regs_final = regs_final[:-3]

            print(aa_nums_final, file=f)
            print(nuc_final, file=f)
            print(regs_final, file=f)
            print('', file=f)


def write_csv(fasta_dat, csv_file):
    '''write CSV data'''
    header = 'Identifier,Amino Acid'
    counter = 1
    indices = [0 for _ in fasta_dat]
    strings = ['{},{}'.format(
        pretty(r), r['aa']) for r in fasta_dat]
    gap = True

    while any(i < len(r['aligned']) for i, r in zip(indices, fasta_dat)):
        header += ',{}'.format(counter)
        counter += 1

        # first, check if in a non-gap area with a removed section (i.e. a '.' character)
        dot = not gap and any(
            i < len(r['aligned']) and r['aligned'][i] == '.' for i, r in zip(indices, fasta_dat))

        for i, r in enumerate(fasta_dat):
            # check if at end of alignment
            if indices[i] >= len(r['aligned']):
                continue  # no need to write anything

            # get character
            c = r['aligned'][indices[i]]

            strings[i] += ','
            if gap and c == '-':  # only '-' alignments may advance
                strings[i] += '-'
                indices[i] += 1
                continue
            elif dot and c == '.':  # only '.' alignments may advance
                strings[i] += '.'
                indices[i] += 1
                continue
            elif not gap and c not in '.-':
                strings[i] += str(indices[i] + 1)
                indices[i] += 1
                continue

        # if all indices are not '-', we change to non-gap
        if all(i >= len(r['aligned']) or r['aligned'][i] != '-' for i, r in zip(indices, fasta_dat)):
            gap = False
        # if all indices are '-', we change to gap
        if all(i >= len(r['aligned']) or r['aligned'][i] == '-' for i, r in zip(indices, fasta_dat)):
            gap = True
        # otherwise, gap status stays the same

    csv_file.write(header + '\n')
    csv_file.write('\n'.join(strings) + '\n')


def record_bad_translation(r, n, n_rev, n_comp, n_rev_comp, fp, fp_summary):
    print('{} - bad translation ({:.2f}% incorrect)'.format(pretty(r),
                                                            100 - n*100), file=fp_summary)
    print('## {}\n'.format(pretty(r)), file=fp)
    print('This record seems to be incorrectly translated from', file=fp)
    print('nucleotides to amino acids. Variation up to 5.00% is', file=fp)
    print('tolerated by the script but the variation in this', file=fp)
    print('translation is {:.2f}%.\n'.format(100 - n*100), file=fp)
    print('Additionally, taking the reverse complement', file=fp)
    print('of the nucleotide sequence still results in a variation', file=fp)
    print('of {:.2f}%.\n'.format(100 - n_rev_comp*100), file=fp)
    print_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_translation_errors(r, fp)
    print_reverse_complement_translation_errors(r, fp)


def record_reverse_translation(r, n, fp, fp_summary):
    print('{} - reverse needed'.format(pretty(r)), file=fp_summary)
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The translation of the nucleotides to', file=fp)
    print('amino acids seems to be the reverse.\n', file=fp)
    print('This has been automatically been done, but the source', file=fp)
    print('data files should be updated.\n', file=fp)
    print_nucleotide_file_data(r, fp)
    print_reverse_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_reverse_translation_errors(r, fp)


def record_complement_translation(r, n, fp, fp_summary):
    print('{} - complement needed'.format(pretty(r)), file=fp_summary)
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The translation of the nucleotides to', file=fp)
    print('amino acids seems to be the complement.\n', file=fp)
    print('This has been automatically been done, but the source', file=fp)
    print('data files should be updated.\n', file=fp)
    print_nucleotide_file_data(r, fp)
    print_complement_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_complement_translation_errors(r, fp)


def record_reverse_complement_translation(r, n, fp, fp_summary):
    print('{} - reverse complement needed'.format(pretty(r)), file=fp_summary)
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The translation of the nucleotides to', file=fp)
    print('amino acids seems to be the reverse complement.\n', file=fp)
    print('This has been automatically been done, but the source', file=fp)
    print('data files should be updated.\n', file=fp)
    print_nucleotide_file_data(r, fp)
    print_reverse_complement_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_reverse_complement_translation_errors(r, fp)


def record_reverse_complement_shift_left_translation(r, n, fp, fp_summary):
    print('{} - reverse complement needed (left-shifted)'.format(pretty(r)), file=fp_summary)
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The translation of the nucleotides to', file=fp)
    print('amino acids seems to be the reverse complement.\n', file=fp)
    print('It also needs to be left-shifted by one place.\n', file=fp)
    print('This has been automatically been done, but the source', file=fp)
    print('data files should be updated.\n', file=fp)
    print_nucleotide_file_data(r, fp)
    print_reverse_complement_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_reverse_complement_translation_errors(r, fp)


def record_reverse_complement_shift_right_translation(r, n, fp, fp_summary):
    print('{} - reverse complement needed (right-shifted)'.format(pretty(r)), file=fp_summary)
    print('## {}\n'.format(pretty(r)), file=fp)
    print('The translation of the nucleotides to', file=fp)
    print('amino acids seems to be the reverse complement.\n', file=fp)
    print('It also needs to be right-shifted by one place.\n', file=fp)
    print('This has been automatically been done, but the source', file=fp)
    print('data files should be updated.\n', file=fp)
    print_nucleotide_file_data(r, fp)
    print_reverse_complement_nucleotide_file_data(r, fp)
    print_amino_acid_file_data(r, fp)
    print_reverse_complement_translation_errors(r, fp)


def record_bad_alignment(r, n, fp, fp_summary):
    print('{} - bad alignment ({} errors)'.format(pretty(r), n), file=fp_summary)
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


def print_reverse_nucleotide_file_data(r, fp, line_length=60):
    rev = reverse(r['nucleotide'])
    print('### Reverse nucleotide data\n', file=fp)
    print('```text', file=fp)
    i = 0
    while rev[i:i+line_length] != "":
        print(rev[i:i+line_length], file=fp)
        i = i + line_length
    print('\n```\n', file=fp)


def print_complement_nucleotide_file_data(r, fp, line_length=60):
    comp = complement(r['nucleotide'])
    print('### Complement nucleotide data\n', file=fp)
    print('```text', file=fp)
    i = 0
    while comp[i:i+line_length] != "":
        print(comp[i:i+line_length], file=fp)
        i = i + line_length
    print('\n```\n', file=fp)


def print_reverse_complement_nucleotide_file_data(r, fp, line_length=60):
    rev = reverse_complement(r['nucleotide'])
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
    aa_trans = r['nucleotide']
    if len(r['ungapped']) != 0:
        bad = ['X' if x1 != x2 else ' ' for (
            x1, x2) in zip(r['ungapped'], aa_trans)]
        print('```text', file=fp)
        i = 0
        while bad[i:i+(line_length // 3)] != []:
            aas = ' ' + '  '.join(r['ungapped'][i:i+(line_length // 3)]) + ' '
            tns = ' ' + '  '.join(aa_trans[i:i+(line_length // 3)]) + ' '
            err = ' ' + '  '.join(bad[i:i+(line_length // 3)]) + ' '
            print("nuc: {}".format(r['nucleotide']
                                   [i:i+line_length]), file=fp)
            print("aas: {}".format(aas), file=fp)
            print("tns: {}".format(tns), file=fp)
            print("err: {}\n".format(err), file=fp)
            i = i + line_length
        print('```\n', file=fp)


def print_reverse_translation_errors(r, fp, line_length=60):
    print('### Reverse translation errors\n', file=fp)
    rev = r['nucleotide']
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


def print_complement_translation_errors(r, fp, line_length=60):
    print('### Complement translation errors\n', file=fp)
    comp = r['nucleotide']
    aa_trans = translate(comp)
    if len(r['ungapped']) != 0:
        bad = ['X' if x1 != x2 else ' ' for (
            x1, x2) in zip(r['ungapped'], aa_trans)]
        print('```text', file=fp)
        i = 0
        while bad[i:i+(line_length // 3)] != []:
            aas = ' ' + '  '.join(r['ungapped'][i:i+(line_length // 3)]) + ' '
            tns = ' ' + '  '.join(aa_trans[i:i+(line_length // 3)]) + ' '
            err = ' ' + '  '.join(bad[i:i+(line_length // 3)]) + ' '
            print("com: {}".format(comp[i:i+line_length]), file=fp)
            print("aas: {}".format(aas), file=fp)
            print("tns: {}".format(tns), file=fp)
            print("err: {}\n".format(err), file=fp)
            i = i + line_length
        print('```\n', file=fp)


def print_reverse_complement_translation_errors(r, fp, line_length=60):
    print('### Reverse complement translation errors\n', file=fp)
    rev = r['nucleotide']
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


def parse_fasta(path, regionsList):
    regionsXs = []
    for rs in regionsList:
        if len(regionsXs) != 0:
            regionsXs.append('-')
        regionsXs.append("X" * rs)
    regionsXs = "".join(regionsXs)
    fasta_data = []
    gapped_sequence = None
    with open(path) as lines:
        for line in lines:
            if line[0] == '>':
                # store previous gapped sequence
                if gapped_sequence:
                    underscored = []
                    for i in range(len(regionsXs)):
                        if regionsXs[i] == '-':
                            underscored.append('-')
                        elif gapped_sequence[i] == '-':
                            underscored.append('_')
                        else:
                            underscored.append(gapped_sequence[i])
                    underscored = "".join(underscored)
                    fasta_data[-1]['underscored'] = underscored
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
                if (xs and xs[0] and xs[0].lower() == 'reg'):
                    xs.pop(0)
                    regions = True
                else:
                    regions = False
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
                try:
                    genus, num = '_'.join(xs).lower().split('/')
                    if genus[-4:] == "_aln":
                        genus = genus[:-4]
                except ValueError:
                    genus, num = xs[0].lower(), "0"
                fasta_data.append({
                    'aa': aa,
                    'kingdom': kingdom,
                    'regions': regions,
                    'pdb': pdb,
                    'letter': letter,
                    'genus': genus,
                    'num': num
                })
            else:
                gapped_sequence += line.strip()
    underscored = []
    for i in range(len(regionsXs)):
        if regionsXs[i] == '-':
            underscored.append('-')
        elif gapped_sequence[i] == '-':
            underscored.append('_')
        else:
            underscored.append(gapped_sequence[i])
    underscored = "".join(underscored)
    fasta_data[-1]['underscored'] = underscored
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

    # check reverse
    rev_nuc_seq = reverse(nuc_seq)
    aa_trans = translate(rev_nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n_rev = (max(0, sum(x) - bad)) / len(aa_seq)

    # check complement
    comp_nuc_seq = complement(nuc_seq)
    aa_trans = translate(comp_nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n_comp = (max(0, sum(x) - bad)) / len(aa_seq)

    # check reverse complement
    rev_comp_nuc_seq = reverse_complement(nuc_seq)
    aa_trans = translate(rev_comp_nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n_rev_comp = (max(0, sum(x) - bad)) / len(aa_seq)

    # check reverse complement (shift left)
    rev_comp_shift_left_nuc_seq = reverse_complement_shift_left(nuc_seq)
    aa_trans = translate(rev_comp_shift_left_nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n_rev_comp_shift_left = (max(0, sum(x) - bad)) / len(aa_seq)

    # check reverse complement (shift right)
    rev_comp_shift_right_nuc_seq = reverse_complement_shift_right(nuc_seq)
    aa_trans = translate(rev_comp_shift_right_nuc_seq)
    bad = abs(len(aa_seq) - len(aa_trans))
    x = [1 if x1 == x2 else 0 for (x1, x2) in zip(aa_seq, aa_trans)]
    n_rev_comp_shift_right = (max(0, sum(x) - bad)) / len(aa_seq)

    return n, n_rev, n_comp, n_rev_comp, n_rev_comp_shift_left, n_rev_comp_shift_right


def reverse(nuc_seq):
    return "".join([i for i in reversed(nuc_seq[3:] + "GTA")])


def complement(nuc_seq):
    return "".join([alt_nuc(i) for i in nuc_seq])


def reverse_complement(nuc_seq):
    return "".join([alt_nuc(i) for i in reversed(nuc_seq[3:] + "GTA")])


def reverse_complement_shift_left(nuc_seq):
    return "".join([alt_nuc(i) for i in reversed(nuc_seq[6:] + "GTAGTA")])


def reverse_complement_shift_right(nuc_seq):
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
    parser = argparse.ArgumentParser(description='Process a fasta file.')
    parser.add_argument('fasta_file', metavar='FILE',
                        help='the fasta file to process')
    parser.add_argument('--regions', required=True,
                        help='the regions to use with this file (ex. "[10,5,10,20]"')
    args = parser.parse_args()
    preprocessing.main(args.fasta_file)
    main(args.fasta_file, args.regions)
