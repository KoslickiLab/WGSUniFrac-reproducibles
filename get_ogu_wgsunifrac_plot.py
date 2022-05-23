import wgsunifrac as wu
import argparse

def main():
    parser = argparse.ArgumentParser(description="Get line plot for ogu vs. wgsunifrac experiment.")
    parser.add_argument('-f', '--file', type=str, help="Dataframe file.")
    parser.add_argument('-x', '--x', type=str, help="x axis.")
    parser.add_argument('-s', '--save', type=str, help="If wants to save df for future use, file name to save as")
    
    args = parser.parse_args()
    file = args.file
    x = args.x
    save = args.save
    wu.get_ogu_vs_wgsunifrac_plot(file, x, save)


if __name__ == '__main__':
    #df = tu.get_dataframe("data/testeverything")
    main()
