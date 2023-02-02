'''
© 2022 Cambridge University
SPDX-FileCopyrightText: 2022 Cambridge University
SPDX-License-Identifier: GPL-3.0-or-later

Utility to generate CanRisk files for all combinations of gene tests
'''
import itertools
import random


def get_header():
    ''' Get a CanRisk file header with random risk factors and PRS. '''
    prs = ""
    c = random.choice([0, 1, 2, 3])
    if c == 1:
        prs = f"##PRS_BC=alpha={random.choice(['0.444', '0.442', '0.44'])}, zscore={round(random.uniform(-1.5, 1.5), 6)}\n"
    elif c == 2:
        prs = f"##PRS_BC=alpha={random.choice(['0.444', '0.442', '0.44'])}, zscore={round(random.uniform(-1.5, 1.5), 6)}\n"
        prs += f"##PRS_OC=alpha=0.223, zscore={round(random.uniform(-1.5, 1.5), 6)}\n"
    elif c == 3:
        prs = f"##PRS_OC=alpha=0.223, zscore={round(random.uniform(-1.5, 1.5), 6)}\n"

    menopause = ""
    if random.choice([0, 1]) == 1:
        menopause = f"##menopause={random.randint(38, 46)}\n"

    return (f"##CanRisk 2.0\n"
            f"##height={round(random.uniform(151.5, 171.5), 6)}\n"
            f"##bmi={random.randint(16,34)}\n"
            f"##birads={random.choice(['NA', 'a', 'b', 'c', 'd'])}\n"
            f"##menarche={random.randint(10,16)}\n"
            f"##alcohol={random.uniform(0.0, 30.0)}\n"
            f"##mht_use={random.choice(['NA', 'N', 'F','E','C'])}\n"
            f"##oc_use={random.choice(['NA', 'N', 'F:2', 'C:3'])}\n"
            f"##TL={random.choice(['NA', 'N', 'Y'])}\n"
            f"##endo={random.choice(['NA', 'N', 'Y'])}\n"
            f"{menopause}{prs}")

GENES = ["BRCA1", 'BRCA2', 'PALB2', 'ATM', 'CHEK2', 'BARD1', 'RAD51D', 'RAD51C', 'BRIP1']

x = ['0:0', 'S:N', 'S:P']                     # possible gene tests
c = [p for p in itertools.product(x, repeat=len(GENES))]    # generate all combinations for ER:PR:HER2:CK14:CK56

print(c)
print(len(c))

# generate CanRisk files for all combinations of BC1 pathology
for t in c:
    gts = '\t'.join(t)
    fn = '-'.join(t).replace(":", "") + ".canrisk2"

    with open("/home/tim/workspace/boadicea/boadicea/tests_selenium/pathtests/PPPNN.canrisk2") as file:
        lines = file.readlines()

    with open("/tmp/gentests/"+fn, 'a') as out:
        out.write(get_header())
        for line in lines:
            if '##FamID' in line or '##' not in line:
                if "1951" in line:      # assumes mother born in 1951
                    print(line.strip().split()[:-1-len(GENES)])
                    parts = line.strip().split()[:-1-len(GENES)]
                    out.write("\t".join(parts)+"\t"+gts+"\t0:0:0:0:0"+"\n")
                    print("\t".join(parts)+"\t"+gts+"\t0:0:0:0:0 :: "+fn)
                else:
                    out.write(line)
