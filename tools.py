import glob
import os.path

AA_LIST = ['gln', 'leu', 'glu', 'ile',
           'arg', 'met', 'val', 'cys', 'trp', 'tyr']
KINGDOM_LIST = ['bact']


def parse_fasta(path):
    fasta_data = []
    with open(path) as lines:
        for line in lines:
            if line[0] == '>':
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
    return fasta_data


def write_standardized_data(fasta_data):
    """Look through Alex's data and write file (if found) in standard format."""
    for ext in ['aa', 'nuc']:
        out_path = 'data/{}.{}'.format(make_filename(fasta_data), ext)
        if not os.path.exists(out_path):
            if fasta_data['pdb']:
                f = '{}*_' + ext
            else:
                f = '{}{}_{}_{}'.format(
                    fasta_data['letter'],
                    fasta_data['genus'],
                    fasta_data['aa'],
                    ext
                )

            g = glob.glob('data/**/*' + f, recursive=True)

            # SPECIAL CASE 1
            if fasta_data['genus'] == 'obscuriglobus':
                f = 'Gemmata_{}_{}'.format(
                    fasta_data['aa'], ext
                )
                g = glob.glob('data/**/*' + f, recursive=True)

            # SPECIAL CASE 2
            if fasta_data['aa'] == 'leu' and not g:
                f = '{}{}_{}ALPHA_{}'.format(
                    fasta_data['letter'],
                    fasta_data['genus'],
                    fasta_data['aa'],
                    ext
                )
                g = glob.glob('data/**/*' + f, recursive=True)

            # SPECIAL CASE 3
            if fasta_data['genus'] == 'asiaticus':
                f = 'CAmoebophilusAsiaticus_{}_{}'.format(
                    fasta_data['aa'],
                    ext
                )
                g = glob.glob('data/**/*' + f, recursive=True)

            if not g:
                # print('Missing data for: {}'.format(out_path))
                continue
            else:
                with open(g[0]) as f_in:
                    with open(out_path, 'w') as f_out:
                        for line in f_in:
                            f_out.write(line)


def write_binary_data(filename):
    for fasta_data in parse_fasta(filename):
        prefix = make_filename(fasta_data)
        aa_dat = read_fasta_file('data/{}.aa'.format(prefix))[2]
        nuc_dat = read_fasta_file('data/{}.nuc'.format(prefix))[2]
        if not aa_dat or not nuc_dat:
            print("Missing data for", prefix)
            continue
        if len(aa_dat) * 3 + 3 == len(nuc_dat):
            nuc_dat = nuc_dat[:-3]
        if len(aa_dat) * 3 != len(nuc_dat):
            print("Data incorrect length:", prefix, len(aa_dat), len(nuc_dat))


def make_filename(fasta_data):
    return '_'.join(
        [x for x in [
            fasta_data['aa'],
            fasta_data['kingdom'],
            fasta_data['pdb'],
            fasta_data['letter'],
            fasta_data['genus'],
            # fasta_data['num']
        ] if x]
    )


def read_fasta_file(path):
    """read the data from the fasta file"""
    if path is None:
        return None, None, None
    try:
        header, gi, dat = None, None, ''
        with open(path) as path_p:
            for next_dat in path_p.readlines():
                if next_dat[0] == '>':
                    header = next_dat.strip()
                    if next_dat[0:4] == '>gi|':
                        try:
                            gi = next_dat[4:].split('|')[0].split()[0]
                        except IndexError:
                            gi = None
                    else:
                        gi = None
                    continue
                dat += next_dat.strip()
        return header, gi, dat
    except:
        return None, None, None


def main(filename):
    for fasta_data in parse_fasta(filename):
        write_standardized_data(fasta_data)
    print("READING STANDARDIZED DATA")
    write_binary_data(filename)


if __name__ == "__main__":
    main('ALL2.fasta')
    # main('ALL3.fasta')
