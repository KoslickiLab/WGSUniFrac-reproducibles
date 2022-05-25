import wgsunifrac as wu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get pcoa plot of taxonomic profiles in a direcory.")
    parser.add_argument('-dir', '--dir', type=str, help="The directory where profiles are located.")
    parser.add_argument('-t', '--plt_title', type=str, help="Title of the pcoa plot.")
    parser.add_argument('-a', '--alpha', type=float, help="The factor for branch length function. -1, -0.5, 0, 1, 5")
    parser.add_argument('-by', '--by', type=str, help="For real data. choices: bodysites, study",default="bodysites")
    parser.add_argument('-type', '--data_type', type=str, help="real or simulated.")
    parser.add_argument('-s', '--save', type=str, help="File name to save file as.")    

    args = parser.parse_args()
    dir = args.dir
    plt_title = args.plt_title
    alpha = args.alpha
    dtype = args.data_type
    by = args.by

    if dtype == "simulated":
        sample_lst, dist_matrix, metadata = tu.pairwise_unifrac(dir, plt_title)
        wu.pairwise_unifrac(dir, plt_title, alpha)
    else:
        #real data
        metadata = wu.get_metadata_from_real_data_partial('data/body_sites12261.csv', dir)
        dist_matrix, sample_lst = wu.just_pairwise_unifrac(dir, alpha)
        wu.get_pcoa(dist_matrix,sample_lst,metadata,plt_title,args.save)


if __name__ == "__main__":
    main()
