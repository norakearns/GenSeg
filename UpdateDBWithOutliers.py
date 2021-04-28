from tinydb import TinyDB, Query, where
import os

fasta_file = open("/Users/davidkearns/Downloads/all_3_concat_fixed.fasta", 'r')
fasta_lines = fasta_file.read()
fasta_file.close()

clus_file = open("/Users/davidkearns/Downloads/all.clstr.sorted", 'r')
clus_lines = clus_file.read()
clus_file.close()

csv_file = open("AlignDB.csv", "a")

clusters = clus_lines.split('>Cluster')[1:-1] # Ignore the leading and trailing empty items, create array of cluster blocks
num_clusters = len(clusters)

i = 0
while i < num_clusters:
     cluster_block = clusters[i]
     cluster_number = cluster_block.split('\n')[0].strip()
     cluster_name = "cluster_" + str(cluster_number)
     cluster_fasta = ''
     entries = cluster_block.split('\n')[1:-1]
     #print("processing cluster: " + cluster_name )
     num_cluster_entries = len(entries)
     #if (i >= 12117):
     #    print("I think this is where the single entry clusters start")
     if (num_cluster_entries == 1):
           print("processing cluster: " + str(cluster_number))
           entry = entries[0] # there is only one line
           info = entry.split("|")
           prot_name = info[1]
           prot_detail = info[2]
           prot_sub = prot_detail.split(".")[0]
           segment_number = prot_sub.split('_')[2]  # example: S11
           prot_ID = ">sp|" + prot_name + "|" + prot_sub
           offset_into_fasta_char_array = fasta_lines.find(prot_ID)
           fasta_match_to_end_of_fasta_entries = fasta_lines[offset_into_fasta_char_array:]
           prot_seq = fasta_match_to_end_of_fasta_entries.split('\n')[1]

           item_str = prot_name + "," + str(segment_number) + "," + "-1" + "," + str(cluster_number) + "," + prot_seq + "\n"

           csv_file.write(item_str)

     i = i + 1

csv_file.close()