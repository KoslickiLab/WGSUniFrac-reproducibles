import wgsunifrac as wu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get a combined dataframe from all experiments under a directory, saving in a file.")
    parser.add_argument('-d', '--dir', type=str, help="directory containing experiments")
    parser.add_argument('-s', '--save', type=str, help="File name of the saved file.")

    args = parser.parse_args()
    dir = args.dir
    save = args.save
    wu.get_ogu_vs_wgsunifrac_df(dir, save)


if __name__ == '__main__':
    main()
