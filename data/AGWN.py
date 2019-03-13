import numpy as np
import argparse

def main(args):
    print(args)
    dataset = np.random.normal(0, 1, size=(args.nconditions,1))
    for genes in range(1, args.ngenes):
        if (genes < args.nclusters):
            dataset = np.hstack((dataset, np.random.normal(0, 1, size=(args.nconditions,1))))
        else:
            dataset = np.hstack((dataset, dataset[:,[genes%args.nclusters]] + args.sigma * np.random.normal(0, 1, size=(args.nconditions, 1))))
    print(dataset)
    np.savetxt(args.outcsv, dataset, delimiter=",", fmt='%s')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gaussian blob.')
    parser.add_argument('--ngenes', type=int,
                        help='number of genes')
    parser.add_argument('--nclusters', type=int, 
                        help='number of clusters')
    parser.add_argument('--nconditions', type=int, 
                        help='number of conditions')
    parser.add_argument('--sigma', type=float, 
                        help='sigma')
    parser.add_argument('--outcsv', type=str, help='output csv directory')
    args = parser.parse_args()
    main(args)