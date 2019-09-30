import pandas as pd
from collections import Counter
import numpy as np

def getUniqueProbes(df):
    df['probe'] = df[0].str.split("_", n=1, expand=True)[1]
    probe_dict = Counter(df['probe'])
    uprobes = [x for x in probe_dict if probe_dict[x] == 1]
    return uprobes


def getprobes(filepath):
    df = pd.read_csv(filepath, sep=' ', header=None)
    df['probe'] = df[0].str.split("_", n=1, expand=True)[1]
    uniqprob = list(set(df['probe']))
    return uniqprob

def main():
    mer8_0file = './8mer0mismatch_hflu.txt'
    mer8_1file = './8mer1mismatch_hflu.txt'
    mer20_0file = './20mer0mismatch_hflu.txt'
    mer20_1file = './20mer1mismatch_hflu.txt'
    mer20_2file = '20mer2mismatch_hflu.txt'
    mer20_3file = '20mer3mismatch_hflu.txt'

    # mer8_0 = pd.read_csv(mer8_0file, sep=' ', header=None)
    # mer8probe_list = getUniqueProbes(mer8_0)

    mer20_0 = pd.read_csv(mer20_0file, sep=' ', header=None)
    mer20probe_list = getUniqueProbes(mer20_0)
    mer20_0.set_index('probe', inplace=True)

    mer20_1probs = getprobes(mer20_1file)
    mer20_2probs = getprobes(mer20_2file)
    mer20_3probs = getprobes(mer20_3file)
    mer8_0probs = getprobes(mer8_0file)
    mer8_1probs = getprobes(mer8_1file)
    allmis = mer20_1probs + mer20_2probs + mer20_3probs + mer8_0probs + mer8_1probs

    nomismat = np.setdiff1d(mer20probe_list, allmis)

    results = []
    for i in nomismat:
        line = mer20_0.loc[i]
        results.append(line)
    df = pd.DataFrame(results)
    df.to_csv('results.tsv', sep='\t')





if __name__ == '__main__':
    main()