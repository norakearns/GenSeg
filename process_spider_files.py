import pandas as pd
import statistics
from timeit import default_timer as timer
import os
import subprocess
import urllib
import re
from urllib import request
from urllib import error
import xlsxwriter

# files_location will differ for users and kinds of computers - Mac vs Windows, etc.
files_location = "/Users/norakearns/downloads/no_markov_fastas/"
chars_and_attr_list = []  # this is a two-dimensional data structure

# at each index position there is a tuple consisting of a single letter and an array
# for example, the entry at index position 6 might be ['a',['DOMAIN', 'STRAND, 'HELIX']]

# Build up the chars_and_attr_list for a fasta string. Each entry in the list contains a single char and an empty attr array
def build_default_attr_structure_for_fasta(fasta_string):
    str_len = len(fasta_string)
    char_index = 0
    while (char_index < str_len):
        initial_entry = []
        attribute_list = []
        initial_entry.append(fasta_string[char_index])  # the letter
        initial_entry.append(attribute_list)  # the type - by default there is nothing set
        chars_and_attr_list.append(initial_entry)
        char_index = char_index + 1
    # print(chars_and_attr_list)

# For a range of indices in a fasta string, append the attributes appropriate for the protein structure
def add_attr_entries_for_zone(start_index, end_index, protein_structure):
    if (start_index < 0) or (end_index > len(chars_and_attr_list) - 1):
        print("ERROR: user provided invalid start or end indices:", start_index, ",", end_index)
        exit(-1)

    ind = start_index
    while (ind <= end_index):
        # chars_and_attr_list[i] is the letter, attribute array entry, for example ['A', ["DOMAIN", "HELIX"]] at index i
        # the letter itself is at chars_and_attr_list[i][0]
        attributes_for_letter = chars_and_attr_list[ind][1] # for example ["DOMAIN", "HELIX"]
        if (protein_structure not in attributes_for_letter):
            # If the attribute you are about to add is not already present in the attribute list, add it
            chars_and_attr_list[ind][1].append(protein_structure)
        ind = ind + 1

def calc_average_segment_length(array_seggos):
    num_segs = len(array_seggos)
    seg_counter = 0
    sample = []
    while (seg_counter < num_segs):
        one_seg = array_seggos[seg_counter]
        seg_length = len(one_seg) - one_seg.count("<")- one_seg.count(">")
        sample.append(seg_length)
        seg_counter += 1
    mean_length = statistics.mean(sample)

    if num_segs > 1:
        return(mean_length, statistics.stdev(sample, mean_length))
    else:
        return(mean_length, 0)

start = timer()
df = pd.read_excel(files_location + 'spider_entries.xlsx', col=1)
entry_array = df["Entry"].tolist()
print(entry_array)


failed_query_entries = []
protein_processed_count = 1
item_count = 0
max_segments_count = 0
all_series = []
# DEBUG - uncomment if you just want to look at a single entry
#entry_array=['Q6AI39']
for one_protein in entry_array:
    one_protein = one_protein.rstrip()
    #one_protein = one_protein.strip(".i1") #remove extension
    print("Processing protein #", protein_processed_count, " :", one_protein)
    protein_processed_count = protein_processed_count + 1
    protein_info_file_loc = files_location + one_protein + '.txt'
    i1_file_loc = files_location + one_protein + '.i1'

    if (os.path.isfile(protein_info_file_loc) == False):
        print('   The protein does not have a txt companion file')
        continue

    my_file = open(protein_info_file_loc, 'r')
    list = my_file.read()
    my_file.close()
    listarray = list.split('\n')
    seq = listarray[1:]
    actual_sequence = ''.join(seq)

    build_default_attr_structure_for_fasta(actual_sequence)

    #get just the family and domains info
    i1file = open(i1_file_loc, 'r')
    i1list = i1file.read()
    i1file.close()
    # print(i1list)
    i1listarray = i1list.split('\n')
    withoutheader = i1listarray[1:]

    i1seq = []
    for i in withoutheader[:-1]:
        i1seq.append(i[0])
    seq_len = len(i1seq)
    char_counter = 0
    while char_counter < seq_len:
        domain_name = i1seq[char_counter]
        range_start = char_counter
        range_end = char_counter

        if (domain_name != "C"): #coil is default, we don't mark those
            add_attr_entries_for_zone(range_start,range_end,domain_name)

        char_counter = char_counter + 1

    print(chars_and_attr_list)


    min_cut_length = 100
    segments = ''
    array_of_segments = []
    row_count = 2 # Start at row 2 after legend in row 1
    len_list = len(chars_and_attr_list)
    in_a_domain = False
    i = 0
    while (i < len(chars_and_attr_list)):
        existing_entry = chars_and_attr_list[i]
        existing_entry_char = existing_entry[0]
        existing_entry_attrs = existing_entry[1]
        if (len(existing_entry_attrs) > 0): # Letter belongs to at least one domain
            if (in_a_domain == False):
                # This is the beginning of at least one domain. Mark with a '<'
                segments += "<"
                in_a_domain = True
            segments += existing_entry_char  # Add letter to segment
        else:
            if (in_a_domain == True): # Letter does not belong to any domain
                # This is the end of all domains for this part of the protein. Mark with a '>'
                segments += ">"
                in_a_domain = False
            if (len(segments)> min_cut_length): # transition from marked to unmarked
                # print("cut at letter: ", existing_entry_char, " at index ", i)
                array_of_segments.append(segments) #finish up the previous segment
                segments = '' #make segment empty
                segments += existing_entry_char #add new character as the first character in a new segment
            else:
                segments += existing_entry_char #if it hasn't reached 100 yet, keep appending
        i = i + 1

    # Write whatever is left over at the end to the tail of the previous item unless over length of 300
    if (len(array_of_segments)>0):
        last_segment = array_of_segments[-1]
        last_segment_length = len(last_segment)
        if (last_segment_length < 300):
            # append the fragment to the previous segment if overall len < 300
            array_of_segments[-1] += segments
        else:
            # go ahead and create another line
            array_of_segments.append(segments)
    else:
        array_of_segments.append(segments)

    # Keep track of the max number of segments in any protein entry. At the end, when we create the pandas
    # dataframe we'll need to pad shorter entries so that each row has the same number of segments, even
    # if some are blank
    len_of_array_of_segments = len(array_of_segments)
    if (len_of_array_of_segments > max_segments_count):
        max_segments_count = len_of_array_of_segments
    avg, stdev = calc_average_segment_length(array_of_segments)
    series_w_segments = [one_protein, len(actual_sequence), actual_sequence, len(array_of_segments), avg, stdev]
    num_segments = len(array_of_segments)
    seg_cnt = 0
    while (seg_cnt < num_segments):
        series_w_segments.append(array_of_segments[seg_cnt])
        seg_cnt = seg_cnt + 1

    new_series = pd.Series(series_w_segments)

    all_series.append(new_series)

    # empty out our chars and attr list for the next block
    chars_and_attr_list = []

print("Building data structures to contain all results")

# at this point we have processed all protein entries. build the final data structures.
series_labels = ["Protein", "Length", "Sequence", "# Subseqs", "Avg Segment Length", "Stdev Segment Length"]
segment_labels = []
seg_cnt = 0
while (seg_cnt < max_segments_count):
    segment_labels.append("S" + str(seg_cnt))
    seg_cnt = seg_cnt + 1

# Figure out the total number of columns needed to accommodate any of the rows
max_cols_needed = len(series_labels) + max_segments_count

# Build up the column labels for the whole table
seg_labels_len = len(segment_labels)
seg_label_cnt = 0
while (seg_label_cnt < seg_labels_len):
    series_labels.append(segment_labels[seg_label_cnt])
    seg_label_cnt = seg_label_cnt + 1

# Go through all of the row series items, adding as many blank column entries as needed
# so that all row entries have the same number of columns, even if some (farthest to right) cols are blank
buffered_series_array = []
for one_series in all_series:
    one_series_len = len(one_series)
    if (one_series_len < max_cols_needed):
        cols_to_add = max_cols_needed - one_series_len
        col_cnt = 0
        while (col_cnt < cols_to_add):
            empty_series = pd.Series([""])
            one_series = one_series.append(empty_series)
            col_cnt = col_cnt + 1
    buffered_series_array.append(one_series)

df_list_array = []
series_cnt = len(buffered_series_array)
counter = 0
while (counter < series_cnt):
    a_series = buffered_series_array[counter]
    a_list = a_series.values.tolist()
    df_list_array.append(a_list)
    counter = counter + 1

df = pd.DataFrame(df_list_array, columns=series_labels)

excel_file = files_location + "Protein_segments_nomarkov.xlsx"
print("Writing final data structure to Excel file at ", excel_file)
df.to_excel(excel_file, index=True)
import pandas as pd


#if (len(skipped_entries) > 0):
    #print("These proteins were skipped as they did not have requested domains")
    #print(skipped_entries)

#skipped_wrkbk = files_location + 'skipped_entries.xlsx'
#df_new = pd.DataFrame(skipped_entries)
#df_new.to_excel(skipped_wrkbk, index=True)

end = timer()
print("Total time in seconds:", end-start )

