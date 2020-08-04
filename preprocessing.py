import sys
import glob
import os.path
import json

AA_LIST = ['ala', 'asn', 'asp', 'gln', 'leu', 'glu', 'gly', 'his', 'ile', 'lys',
           'arg', 'met', 'phe', 'pro', 'pyl', 'sep', 'ser', 'thr', 'val', 'cys', 'trp', 'tyr']
KINGDOM_LIST = ['bact', 'arch']


def parse_fasta(path):
    fasta_data = []
    with open(path) as lines:
        for line in lines:
            if line[0] == '>':
                xs = line[1:].strip().split('_')
                if not (xs and xs[0] and xs[0].lower() in AA_LIST):
                    raise RuntimeError(
                        "Amino Acid ({}) not recognized in {}".format(
                            xs[0], line
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
                        "Kingdom ({}) not recognized in {}".format(
                            xs[0], line
                        )
                    )
                kingdom = xs.pop(0).lower()
                pdb = None
                if xs and xs[0] and len(xs[0]) == 4:
                    pdb = xs.pop(0).lower()
                if not (xs and xs[0] and len(xs[0]) == 1):
                    raise RuntimeError(
                        "'{}' not recognized as first letter of genus in {}".format(
                            xs[0], line
                        )
                    )
                letter = xs.pop(0).upper()
                try:
                    genus, num = '_'.join(xs).lower().split('/')
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
    return fasta_data


def search_data_folder(fasta_data, ext):
    """Glob the `data` folder looking for matches"""
    if fasta_data['pdb']:
        f1 = '{}*_'.format(fasta_data['pdb']) + ext
    else:
        f1 = '{}{}_{}_{}'.format(
            fasta_data['letter'],
            fasta_data['genus'],
            fasta_data['aa'],
            ext
        )
    return glob.glob('data/**/*' + f1, recursive=True)


def search_supplemental_folder(fasta_data, ext):
    """Glob the `supplemental` folder looking for matches"""
    if fasta_data['pdb']:
        f = '{}_*_{}_*_{}.fasta'.format(
            fasta_data['aa'],
            fasta_data['pdb'],
            ('aa' if ext == 'aa' else 'nuc')
        )
    else:
        f = '{}_{}_{}_{}_{}.fasta'.format(
            fasta_data['aa'],
            fasta_data['kingdom'],
            fasta_data['letter'],
            fasta_data['genus'],
            ('aa' if ext == 'aa' else 'nuc')
        )
    return glob.glob('data/supplemental/*' + f, recursive=False)


def search_downloads_folder(fasta_data, ext):
    """Search the Downloads folder looking for matches"""
    downloads = os.path.join(os.path.expanduser('~'), 'Downloads')
    if fasta_data['pdb']:
        f = '{}_*_{}_*_{}.fasta'.format(
            fasta_data['aa'],
            fasta_data['pdb'],
            ('aa' if ext == 'aa' else 'nuc')
        )
    else:
        f = '{}_{}_{}_{}_{}.fasta'.format(
            fasta_data['aa'],
            fasta_data['kingdom'],
            fasta_data['letter'],
            fasta_data['genus'],
            ('aa' if ext == 'aa' else 'nuc')
        )
    return glob.glob('{}/*'.format(downloads) + f, recursive=False)


def write_standardized_data(fasta_data):
    """Look through Alex's data and write file (if found) in standard format."""
    for ext in ['aa', 'nuc']:
        out_path = 'data/{}.{}'.format(make_filename(fasta_data), ext)
        if not os.path.exists(out_path):
            g1 = search_downloads_folder(fasta_data, ext)
            if not g1:
                g1 = search_supplemental_folder(fasta_data, ext)

            # possible location in Alex's data
            g2 = []
            if not g1:
                g2 = search_data_folder(fasta_data, ext)

                # SPECIAL CASE 1
                if fasta_data['genus'] == 'obscuriglobus':
                    f = 'Gemmata_{}_{}'.format(
                        fasta_data['aa'], ext
                    )
                    g2 = glob.glob('data/**/*' + f, recursive=True)

                # SPECIAL CASE 2
                if fasta_data['aa'] == 'leu' and not g2:
                    f = '{}{}_{}ALPHA_{}'.format(
                        fasta_data['letter'],
                        fasta_data['genus'],
                        fasta_data['aa'],
                        ext
                    )
                    g2 = glob.glob('data/**/*' + f, recursive=True)

                # SPECIAL CASE 3
                if fasta_data['genus'] == 'asiaticus':
                    f = 'CAmoebophilusAsiaticus_{}_{}'.format(
                        fasta_data['aa'],
                        ext
                    )
                    g2 = glob.glob('data/**/*' + f, recursive=True)

            g = g1 + g2
            if not g:
                # print('Missing data for: {}'.format(out_path))
                continue
            else:
                with open(g[0]) as f_in:
                    with open(out_path, 'w') as f_out:
                        for line in f_in:
                            f_out.write(line)


def write_binary_data(filename):
    dat = parse_fasta(filename)
    new_dat = []
    isMissingData = False
    for fasta_data in dat:
        prefix = make_filename(fasta_data)
        aa_file = 'data/{}.aa'.format(prefix)
        nuc_file = 'data/{}.nuc'.format(prefix)
        try:
            os.remove(aa_file + '.bad')
        except FileNotFoundError:
            pass
        try:
            os.remove(nuc_file + '.bad')
        except FileNotFoundError:
            pass
        aa_dat = read_fasta_file(aa_file)[2]
        nuc_dat = read_fasta_file(nuc_file)[2]
        if not aa_dat or not nuc_dat:
            if not aa_dat:
                print("Missing aa data for " + prefix)
            if not nuc_dat:
                print("Missing nuc data for " + prefix)
            isMissingData = True
            continue
        if len(aa_dat) * 3 + 3 == len(nuc_dat):
            nuc_dat = nuc_dat[:-3]
        elif len(aa_dat) * 3 != len(nuc_dat):
            err = "Data incorrect length: {} {} {} (expected {})".format(
                prefix,
                len(aa_dat),
                len(nuc_dat),
                len(aa_dat) * 3
            )
            try:
                os.rename(aa_file, aa_file + '.bad')
            except FileExistsError:
                pass
            try:
                os.rename(nuc_file, nuc_file + '.bad')
            except FileExistsError:
                pass
            print(err)
            isMissingData = True
            continue
            raise RuntimeError(err)
        fasta_data['aa_dat'] = aa_dat
        fasta_data['nuc_dat'] = nuc_dat
        new_dat.append(fasta_data)
    if isMissingData:
        pass
    with open(filename + '.json', 'w') as json_file:
        json.dump(new_dat, json_file, indent=2)


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
                if next_dat.strip() == '':
                    continue
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
                else:
                    dat += next_dat.strip()
        return header, gi, dat
    except FileNotFoundError:
        return None, None, None


def main(filename):
    for fasta_data in parse_fasta(filename):
        write_standardized_data(fasta_data)
    write_binary_data(filename)


if __name__ == "__main__":
    main(sys.argv[1])
