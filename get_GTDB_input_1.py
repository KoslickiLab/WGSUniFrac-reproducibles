import wgsunifrac as wu
import argparse
import os
import shutil

def main():
    parser = argparse.ArgumentParser(description="Generate the input data for the first part of the GTDB experiments")
    parser.add_argument('-od', '--out_dir', type=str, help="Output directory.")

    cur_dir = os.getwd()
    args = parser.parse_args()
    out_dir = args.out_dir
    dist_matrix_file = 'data/bac120_valid_distance_matrix.txt'
    distance_dict = wu.get_dist_dict(dist_matrix_file)
    tax_dict = wu.parse_taxonomy_file('data/bac120_taxonomy_no_numeric.tsv')
    name_tax_dict = wu.get_name_taxid_dict('data/bac120_all_taxons_taxids_valid.txt')
    GTDBid_taxid_dict = wu.get_GTDBid_taxid_dict(tax_dict, name_tax_dict)

    ranges = [200, 400, 600, 800, 1000, 1500, 2000, 2500]
    dissimilarity = [-1, 4000, 3000, 2000, 1000, 800]
    for r in ranges: #fix dissimilarity to be 4000
        for i in range(100):
            exp_id = "r" + str(r) + "d4000" + "-" + str(i)
            exp_dir = out_dir + '/' + exp_id
            biom_file = exp_dir + "/distance_matrix.txt"
            if os.path.exists(exp_dir):
                print("directory exists")
                break
            os.mkdir(exp_dir)
            data_dict = wu.create_GTDB_data(distance_dict, 4000, r, 200, 25)
            updated_data = wu.create_GTDB_biom_table(exp_id, data_dict, biom_file, normalize=False)
            profile_dir = exp_dir + '/profiles'
            print(profile_dir)
            os.mkdir(profile_dir)
            for key, value in list(updated_data.items()):
                print(key)
                filename = str(key) + '.profile'
                wu.create_GTDB_profile(value, profile_dir, filename, tax_dict, name_tax_dict, GTDBid_taxid_dict)
    for d in dissimilarity:
        for i in range(100):
            exp_id = "r600" + "d" + str(d) + "-" + str(i)
            exp_dir = out_dir + '/' + exp_id
            biom_file = exp_dir + "/distance_matrix.txt"
            if os.path.exists(exp_dir):
                print("directory exists")
                break
            os.mkdir(exp_dir)
            data_dict = wu.create_GTDB_data(distance_dict, d, 600, 200, 25)
            updated_data = wu.create_GTDB_biom_table(exp_id, data_dict, biom_file, normalize=False)
            profile_dir = exp_dir + '/profiles'
            os.mkdir(profile_dir)
            for key, value in list(updated_data.items()):
                filename = str(key) + '.profile'
                wu.create_GTDB_profile(value, profile_dir, filename, tax_dict, name_tax_dict, GTDBid_taxid_dict)
    os.chdir(cur_dir)
    os.mkdir('results/GTDB/exp1/testRange')
    os.mkdir('results/GTDB/exp1/testDissimilarity')
    os.system('cp -r results/GTDB/exp1/r600d4000* results/GTDB/exp1/testDissimilarity/')
    os.system('mv results/GTDB/exp1/*d4000* results/GTDB/exp1/testRange/')
    os.system('mv results/GTDB/exp1/r600d* results/GTDB/exp1/testDissimilarity')


if __name__ == "__main__":
    main()
