import taxunifrac as tu
import argparse
from ete3 import TreeNode
import copy

def main():
    parser = argparse.ArgumentParser(description="Get a pairwise distance matrix for a given file"
                                                 "containing otus and a tree.")
    parser.add_argument('-f', '--node_file', type=str, help='File containing otus in one of the columns.')
    parser.add_argument('-t','--tree', type=str, help='Tree file.')
    parser.add_argument('-c', '--col', type=int, help='Column in the file containing otus.')
    parser.add_argument('-o', '--output', type=str, help='Output file name.')
    parser.add_argument('-gtdb', '--gtdb', type=int, help='Get distance matrix from a GTDB tree? Yes:1, No:0')
    args = parser.parse_args()

    if args.gtdb == 0:
        file = args.node_file
        tree = TreeNode(args.tree, format=1, quoted_node_names=True)
        nodes = []
        col = args.col
        with open(file,'r') as f:
            for line in f.readlines():
                line = line.strip().split('\t')
                if len(line) > col and line[col][0] != '>':
                    nodes.append(line[col])
        real_nodes = list(map(lambda x: tree&x, nodes))
        print('nodes ready')
        with open(args.output, 'w+') as f:
            for node in nodes:
                print(node)
                (this_node, sorted_nodes) = tu._get_sorted_distance(tree&node, real_nodes)
                f.write("%s\t" % this_node)  # will be the first one anyway
                f.writelines("%s\t" % n for n in sorted_nodes)
                f.write("\n")
    else:
        #GTDB
        print("ready")
        tree = TreeNode(args.tree, format=1, quoted_node_names=True)
        node_lst = tree.get_leaves()
        node_lst = list(map(lambda x: x.name, node_lst))
        print("list of node names ready")
        valid_ids = tu.get_GTDB_valid_IDs('data/GTDB/bac120_all_taxons_taxids.txt', 'data/GTDB/bac120_taxonomy_no_numeric.tsv')
        node_lst = set(node_lst).intersection(set(valid_ids))
        print("intersecting completed")
        print("total nodes: ", len(node_lst))
        real_nodes = list(map(lambda x: tree&x, node_lst))
        print("list of real nodes ready")
        node_copy = list(map(lambda x: tree&x, node_lst))
        print("copied")
        with open(args.output, 'w+') as f:
            for node in real_nodes:
                print(node)
                (this_node, sorted_nodes) = tu._get_sorted_distance(node, node_copy)
                f.write("%s\t" % this_node)  # will be the first one anyway
                f.writelines("%s\t" % n for n in sorted_nodes)
                f.write("\n")


if __name__ == '__main__':
    main()
