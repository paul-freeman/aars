import sys
import os
import glob

AA_LIST = ['gln', 'leu', 'glu', 'ile', 'lys',
           'arg', 'met', 'val', 'cys', 'trp', 'tyr']
KINGDOM_LIST = ['bact', 'arch']


def parse_fasta(path):
    for x in glob.glob(path + "/**/*nuc.fasta", recursive=True):
        name = None
        seq = ""
        filename = os.path.splitext(x)[0] + ".middle-base.fasta"
        with open(filename, "w") as f:
            with open(x) as lines:
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
                        if name is not None:
                            print(">{}".format(name), file=f)
                            print("".join([seq[i]
                                           for i in range(1, len(seq), 3)]), file=f)
                            seq = ""
                        name = make_filename({
                            'aa': aa,
                            'kingdom': kingdom,
                            'regions': regions,
                            'pdb': pdb,
                            'letter': letter,
                            'genus': genus,
                            'num': num
                        })
                    else:
                        seq += line.strip()
            print(">{}".format(name), file=f)
            print("".join([seq[i] for i in range(1, len(seq), 3)]), file=f)
    return


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


if __name__ == "__main__":
    parse_fasta(sys.argv[1])
