import wgsunifrac as wu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Produces pairwise WGS UniFrac distance matrix with given user input.")
    parser.add_argument('-od', '--out_dir', type=str, help="A sub-directory 'abundance_file' will be created under this directory and abundance files deposited there.")
    parser.add_argument('-df', '--distance_file', type=str, help='Pairwise sorted distance matrix file.')
    parser.add_argument('-s', '--sample_num', type=int, help='Number of samples in each sample.')
    parser.add_argument('-o', '--org_num', type=int, help='Number of organisms in each environment')
    parser.add_argument('-e', '--env_num', type=int, help='Number of environments.')
    #parser.add_argument('-r', '--rnge', type=int, help='Range value')
    #parser.add_argument('-d', '--dissim', type=int, help='Dissimilarity value.')

    args = parser.parse_args()
    dist_file = args.distance_file
    sample_num = args.sample_num
    org_num = args.org_num
    out_dir = args.out_dir
    env_num = args.env_num
    #rnge = args.rnge
    #dissim = args.dissim
    distance_dict = wu.get_dist_dict(dist_file)
    ranges = [300, 500, 1000, 1500, 2000]
    dissimilarity = [3000, 2000, 1500, 1000, 500]
    range_out_dir = out_dir + "/testRange"
    diss_out_dir = out_dir + "/testDissimilarity"
    if not os.path.exists(range_out_dir):
        os.mkdir(range_out_dir)
    if not os.path.exists(diss_out_dir):
        os.mkdir(diss_out_dir)
    for r in ranges:
        for i in range(5):
            dir_name = range_out_dir +'/e'+ env_num + 'r' + r + 'd1500-' + str(i)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            wu.get_grinder_abundance_for_ogu(sample_num, org_num, dir_name, env_num, r, distance_dict, 1500)
    for dissim in dissimilarity:
        for i in range(5):
            dir_name = diss_out_dir + '/e' + env_num + 'r500d' + dissim + '-' + (i)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            wu.get_grinder_abundance_for_ogu(sample_num, org_num, dir_name, env_num, 500, distance_dict, dissim)



if __name__ == '__main__':
    main()
