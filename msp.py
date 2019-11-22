import numpy as np
import pandas as pd
import re
import sys
import os

# Amino Acid letters
AAs = 'ACDEFGHIKLMNPQRSTVWY'

# Known immonium ions
imm = np.array(['AA', 'RA', 'RB', 'RC', 'RD', 'RE', 'RF', 'RG', 'RJ', 'NA', 'NB', 'DA', 'DB', 'CA', 'EA', 'QA', 'QB',
                'QC', 'QD', 'GA', 'HA', 'HB', 'HC', 'HD', 'HE', 'HF', 'IA', 'IB', 'IC', 'LA', 'LB', 'LC', 'KA', 'KB',
                'KC', 'KD', 'KE', 'KF', 'MA', 'MB', 'FA', 'FB', 'PA', 'SA', 'TA', 'WA', 'WB', 'WC', 'WD', 'WE', 'WF',
                'WG', 'YA', 'YB', 'YC', 'VA', 'VB', 'VC', 'VD', '0'])

#Make isotopes of an array of ions
def makeIsotopes(ions):
    iso = []
    # For each immonium ion, append the +i and +2i isotopes
    for im in ions:
        iso.append(im + '+i')
        iso.append(im + '+2i')

    ions = np.append(ions, np.array(iso))

    return set(ions)

# TODO: Use sklearn's OneHotEncoder
def OneHotEncode(aa):
    AA = aa.upper()
    hot = AAs.find(AA)
    hotV = np.zeros(20)
    if hot != -1:
        hotV[hot] = 1

    return hotV

def OneHotEncoder(pep):
    l = []
    for aa in pep:
        l.append(OneHotEncode(aa))

    return np.array(l)

# TODO: Charge is also important (separate data by charge)
# Get immonium ions from the file
# Sequence, charge, {Immonium Ions}

# Parse the MSP file and get the immonium ion dataset
def getImmDataset(fname=None, maxz=3, dict=None):

    seq = []
    X = []
    mz = []
    charge = []

    header = False
    dat = []
    with open(fname, 'r') as slib:
        for line in slib:
            # Check for empty line for spectrum end and header start
            if (line[0] == '\r' or line[0] == '#' or line[0] == '\n'):
                header = False
                if not dat:
                    dat = ['none']

                X.append(dat[:])
                dat.clear()
                continue

            # Check if spectrum header
            if (header == False):
                if (line[:7] != 'Comment'):
                    key, val = line.split(':')

                    if key == 'Name':
                        # Take out the sequence
                        nm,chg = val.split('/')
                        seq.append(nm[1:])

                        # Take out the charge
                        ch = chg.split('_')
                        charge.append(int(ch[0]))

                    if key == 'MW':
                        mz.append(float(val))

                    if key == 'Num peaks':
                        header = True
#                else:
#                    cc = line.find('Charge=')
#                    cz = line[cc+7:cc+8]
#                    charge.append(int(cz))

            # If not the header, then it is the peak list
            else:
                _, _, label = line.split('\t')

                label = label[1:-1]

                lbls = re.findall(r"(I[A-Z]+[+]*\d*i*)", label)
                lbls = list(map(lambda x: x[1:], lbls))

                # If anything to append then append
                if lbls:
                    for lab in lbls:
                        if lab not in dict:
                            dict.add(lab)
                            print ('Label Added:' + lab)
                    #print(lbls)
                    dat.extend(lbls)

    # For the last one
    if not dat:
        dat = ['none']

    X.append(dat[:])
    dat.clear()

    return pd.DataFrame(list(zip(seq, mz, charge, X)), columns=['seq', 'm/z', 'chg', 'ions'])

def main():
    # Check the Python version
    pyversion = float(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))
    if (pyversion < 3.5):
        print('\nERROR: This software requires Python 3.5 or later version')
        print('       Your version is Python' + str(pyversion) + '\n')
        exit(-1)

    filename = '/lclhome/mhase003/Data/sample.msp'
    # Convert to abs path
    filename = os.path.abspath(filename)

    # Check if MSP file exists
    if (os.path.exists(filename) == False):
        sys.exit(-1)

    dict = makeIsotopes(imm)
    data = getImmDataset(filename, 3, dict)

    print (data[:2])

if __name__ == "__main__":
    main()