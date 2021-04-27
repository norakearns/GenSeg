fasta_file = open("/Users/norakearns/PLESA/all_3_concat_fixed.fasta", 'r')
fasta_lines = fasta_file.read()
fasta_file.close()

clus_file = open("/Users/norakearns/PLESA/all.clstr.sorted", 'r')
clus_lines = clus_file.read()
clus_file.close()

clusters = clus_lines.split('>Cluster')[1:-1] # Ignore the leading and trailing empty items, create array of cluster blocks
num_clusters = len(clusters)

i = 0
while i < num_clusters:
     cluster_block = clusters[i]
     cluster_name = "cluster_" + cluster_block.split('\n')[0].strip()
     cluster_fasta = ''
     entries = cluster_block.split('\n')[1:-1]
     print("processing cluster: " + cluster_name )
     num_cluster_entries = len(entries)
     if (i>= 12117):
     if (num_cluster_entries == 1):
          for entry in entries:
               info = entry.split("|")
               prot_name = info[1]
               prot_detail = info[2]
               prot_sub = prot_detail.split(".")[0]
               prot_ID = ">sp|" + prot_name + "|" + prot_sub
               fasta_lines_index = fasta_lines.find(prot_ID)
               new_string = fasta_lines[fasta_lines_index:]
               new_string1 = new_string.split('\n>sp')
               new_string1 = new_string1[0] + '\n'
               cluster_fasta = cluster_fasta + new_string1

          print(thing_I_want)
          cluster_fname = "/Users/norakearns/PLESA/ClusterFiles/" + cluster_name + ".fasta"
          cluster_file = open(cluster_fname, 'w')
          cluster_file.write(cluster_fasta)
          cluster_file.close()
          print("Wrote cluster file " + cluster_fname )
     else:
          print("Cluster " + cluster_name + " only has a single entry, no point in generating a cluster file")

     i = i + 1