import pandas as pd
from timeit import default_timer as timer
import os
import urllib
import re
from urllib import request
from urllib import error

# files_location will differ for users and kinds of computers - Mac vs Windows, etc.
files_location = "/Users/norakearns/Applications/HugeExcel/"
df = pd.read_csv(files_location + "Big_spreadsheet.csv", low_memory=False, dtype=str)
entry_array = df['Protein'].array
# row = df.loc[df['Protein'] == 'O60779']
# num_subseqs = row['# Subseqs']
# num_subseqs_cnt = int(num_subseqs.tolist()[0])
# print(num_subseqs_cnt)

fat_fasta = open("/Users/norakearns/Applications/HugeExcel/fatboy.fasta", "w")  # write mode

# skipped_entries = []
failed_query_entries = []
protein_processed_count = 1
max_count = 1 # stop at this many proteins in big spreadsheet
item_count = 0
# max_segments_count = 0
# all_series = []
# # DEBUG - uncomment if you just want to look at a single entry

for one_protein in entry_array:
#     print("Processing protein #", protein_processed_count, " :", one_protein)
#     protein_processed_count = protein_processed_count + 1
     URL = 'https://www.uniprot.org/uniprot/' + one_protein + '.fasta'
     protein_info_file_loc = files_location + one_protein + '.fasta'
     try:
         urllib.request.urlretrieve(URL, protein_info_file_loc)
     except urllib.error.HTTPError as e:
         print('    The server could not fulfill the request for ', URL)
         print('    Error code: ', e.code)
         failed_query_entries.append(URL)
         if (len(failed_query_entries) > 10):
             print("FAIL: ongoing server issues, terminating run")
             exit(-1)
         continue
     except urllib.error.URLError as e:
         print('    We failed to reach a server when requesting ', URL)
         print('    Reason: ', e.reason)
         failed_query_entries.append(URL)
         if (len(failed_query_entries) > 10):
             print("FAIL: ongoing server issues, terminating run")
             exit(-1)
         continue
     except ConnectionResetError as e:
         print('    Connection reset by peer when requesting ', URL)
         print('    Reason: ', e.reason)
         failed_query_entries.append(URL)
         if (len(failed_query_entries) > 10):
             print("FAIL: ongoing server issues, terminating run")
             exit(-1)
         continue
     else:
         print('    Successfully downloaded ', URL)

     my_file = open(protein_info_file_loc, "r")
     # read the whole file into a single string
     wholefile = my_file.read()
     my_file.close()
     # once we've read the file in, delete it so we don't fill up the hard drive with temp files
     os.remove(protein_info_file_loc)
#
     #just get the first line from the text file
     proteinseq = wholefile.splitlines()
     labeled_seq = (proteinseq[0]) # this is the first line of the fasta file
     #print(labeled_seq)

     # print the protein subseqs for one item
     row = df.loc[df["Protein"] == one_protein]
     num_subseqs = row['# Subseqs']
     num_subcells = int(num_subseqs.tolist()[0])
     loop_cnt = 0
     while (loop_cnt < num_subcells):
        col_prefix = "S" + str(loop_cnt)
        one_cell = row[col_prefix].tolist()[0]
        one_cell = one_cell.replace("<","")
        one_cell = one_cell.replace(">","")
        fat_fasta.write(labeled_seq + "\n")
        fat_fasta.write(col_prefix + "\n")
        fat_fasta.write(one_cell + "\n")
        loop_cnt = loop_cnt + 1

     print("    Protein " + one_protein + " information written to /Users/norakearns/Applications/HugeExcel/fatboy.fasta")

     item_count = item_count + 1
     if (item_count == max_count):
         print("*** reached " + str(max_count) + " entries, stopping ***")
         break

fat_fasta.close()