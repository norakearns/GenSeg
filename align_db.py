# This builds a database of alignment values for each of the clusters
# If a new protein is found that hasn't been seen before, we add it
# If a match for an existing protein is found, but a new maximum alignment value is found, we update the entry
# If a cluster only has a single entry, there is no need to run alignment on it
# The db.json file that's created has entries like protein, segment number, max alignment pct, cluster number of max alignment, segment protein sequence
from Bio import pairwise2
from tinydb import TinyDB, Query, where
import os

def findpct_ident(X, Y):
    alignmentstring = pairwise2.align.globalxx(X, Y, one_alignment_only=1)
    alignment = alignmentstring[0][2]
    pct_ident = alignment / len(X)
    return (pct_ident)

db = TinyDB('cluster_db.json') # open or create the db
#db.truncate() # clear the database

relevant_path = "/Users/norakearns/Applications/ClusterFiles/"
included_extensions = ['fasta']
file_names = [fn for fn in os.listdir(relevant_path)
              if any(fn.endswith(ext) for ext in included_extensions)]

num_file_names = len(file_names)
print(num_file_names)

cluster_numbers = []
for f in file_names:
    clus_num = f.replace("cluster_", "")
    cluster_number = clus_num.replace(".fasta", "")
    cluster_numbers.append(cluster_number)

cluster_numbers.sort(key=int)

cluster_count = 1
for c in cluster_numbers:
    f = "cluster_" + c + ".fasta"
    print("Processing cluster number: " + str(cluster_count) + ", cluster_file: " + f)
    cluster_number = c

    fname = "/Users/norakearns/Applications/ClusterFiles/" + f
    my_file = open(fname, 'r')
    lines = my_file.read()
    my_file.close()
    # Import pairwise2 module
    lines1 = lines.split('\n')

    i_count = 1
    x_count = 0 # Offset of index into the list of protein subsequences
    y_count = 0 # Offset of index into the list of protein subsequences

    while i_count < len(lines1):
        if ((x_count % 10) == 0):
            print("Checking protein subseq " + str(x_count) + " against all others...")
        seg_i = lines1[i_count]
        description_line = lines1[i_count - 1].split()
        seg_info = description_line[0] # example: >sp|O43345|ZN208_HUMAN_S11
        seg_info_items = seg_info.split('|')
        prot = seg_info_items[1] # example: 043345
        segment_ident = seg_info_items[2] # example: ZN208_HUMAN_S11
        segment_number = segment_ident.split('_')[2] # example: S11

        j_count = 1
        y_count = 0
        highest_match_value = 0
        while j_count < len(lines1):
            if (i_count != j_count):
                seg_j = lines1[j_count]
                match_pct = findpct_ident(seg_i, seg_j)
                if ((match_pct > highest_match_value) and (match_pct < 1)):
                    highest_match_value = match_pct

            j_count += 2 # Offset into split overall fasta line array
            y_count += 1 # Offset into just the list of protein subsequences

        if (highest_match_value < 0.6):
            db.insert({ 'protein' : prot, 'seg_num' : segment_number, 'max' : highest_match_value, 'cluster' : cluster_number, 'segment' : seg_i })

        i_count += 2 # Offset into split overall fasta line array
        x_count += 1 # Offset into just the list of fasta subsequences

    cluster_count = cluster_count + 1
    #if (cluster_count > 10):
    #    break

#print("- database contents -")
#for item in db:
#     print(item)

# print("Building scatter plot of protein subseq indices vs alignment values for all other subseqs in cluster...")
# plt.ylabel('alignment identity score')
# plt.xlabel('position')
# plt.ylim([0,1.1])
# plt.scatter(X,Y)
# plt.show()



