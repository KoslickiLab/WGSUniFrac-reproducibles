import os
import wgsunifrac as wu
import argparse

#generates the 16S profiles and WGS profiles at the same time for all the settings of range and dissimilarity

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get raw input data for experiment 1.")
    parser.add_argument('-o', '--out_dir', type=str, help="Output directory.")
    parser.add_argument('-dm', '--distance_matrix', type=str, help="Distance matrix file.")
    parser.add_argument('-mf', '--mapping_file', type=str, help="Mapping file between otu and accession.")
    args = parser.parse_args()
    out_dir = args.out_dir
    dist_file = args.distance_matrix
    mapping_file = args.mapping_file

    distance_dict = wu.get_dist_dict(dist_file)
    otu_tax_dict = wu.get_dict_from_file(mapping_file, 0, 2)
    ranges = [200, 500, 1000, 5000, 10000, 15000, 20000]
    dissimilarity = [-1, 30000, 20000, 10000, 5000, 1000, 900, 800]
    range_out_dir = out_dir + "/testRange"
    diss_out_dir = out_dir + "/testDissimilarity"
    if not os.path.exists(range_out_dir):
        os.mkdir(range_out_dir)
    if not os.path.exists(diss_out_dir):
        os.mkdir(diss_out_dir)
    for r in ranges:
        for i in range(100):
            wu.run_one(distance_dict, otu_tax_dict, num_org=200, num_sample=25, Range=r, dissimilarity=-1, run=i,
                       out_dir=range_out_dir)
    for sim in dissimilarity:
        for i in range(100):
            wu.run_one(dist_dict=distance_dict, tax_dict=otu_tax_dict,
                    num_org=200, num_sample=25, Range=500, dissimilarity=sim, run=i, out_dir=diss_out_dir)
