import wgsunifrac as wu
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="Create abundance files for grinder amplicon according to user input")
    parser.add_argument('-s', '--sample_num', type=int, help='Number of samples for each environment.')
    parser.add_argument('-o', '--organism_num', type=int, help='Number of organisms per sample.')
    parser.add_argument('-od', '--out_dir', type=str, help="Output directory.")
    parser.add_argument('-e', '--env_num', type=int, help="Number of environments.", default=2)
    
    args = parser.parse_args()
    out_dir = args.out_dir
    env_num = args.env_num
    sample_num = args.sample_num
    org_num = args.organism_num
    
    ranges = [200, 500, 1000, 2000, 3000] #dissimilarity=4000
    dissimilarity = [1000, 2000, 3000, 4000, 5000, 6000] #range=500
    range_out_dir = out_dir + "/testRange-env" + str(env_num)
    diss_out_dir = out_dir + "/testDissimilarity-env" + str(env_num)
    if not os.path.exists(range_out_dir):
        os.mkdir(range_out_dir)
    if not os.path.exists(diss_out_dir):
        os.mkdir(diss_out_dir)
    for r in ranges:
        for i in range(5):
            dir_name = range_out_dir +'/e'+ str(env_num) + 'r' + str(r) + 'd4000-' + str(i)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            wu.get_grinder_abundances_for_both(sample_num, org_num, dir_name, env_num, r, 1500)
    for dissim in dissimilarity:
        for i in range(5):
            dir_name = diss_out_dir + '/e' + str(env_num) + 'r500d' + str(dissim) + '-' + str(i)
            if not os.path.exists(dir_name):
                os.mkdir(dir_name)
            wu.get_grinder_abundances_for_both(sample_num, org_num, dir_name, env_num, 500, dissim)






if __name__ == "__main__":
    main()
