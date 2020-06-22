import xlsxwriter
import pandas as pd
import urllib
import re
import xlsxwriter
from urllib import request

chars_and_attr_list = []  # this is a two-dimensional data structure

# at each index position there is a tuple consisting of a single letter and an array
# for example, the entry at index position 6 might be ['a',['DOMAIN', 'STRAND, 'HELIX']]

# Build up the chars_and_attr_list for a fasta string. Each entry in the list contains a single char and an empty attr array
def build_default_attr_structure_for_fasta(fasta_string):
    for i in range(len(fasta_string)):
        initial_entry = []
        attribute_list = []
        initial_entry.append(fasta_string[i])  # the letter
        initial_entry.append(attribute_list)  # the type - by default there is nothing set
        chars_and_attr_list.append(initial_entry)


# For a range of indices in a fasta string, append the attributes appropriate for the protein structure
def add_attr_entries_for_zone(start_index, end_index, protein_structure):
    if (start_index < 0) or (end_index > len(chars_and_attr_list) - 1):
        print("ERROR: user provided invalid start or end indices:", start_index, ",", end_index)
        return

    i = start_index
    while (i <= end_index):
        # chars_and_attr_list[i] is the letter, attribute array entry, for example ['A', ["DOMAIN", "HELIX"]] at index i
        # the letter itself is at chars_and_attr_list[i][0]
        attributes_for_letter = chars_and_attr_list[i][1] # for example ["DOMAIN", "HELIX"]
        if (protein_structure not in attributes_for_letter):
            # If the attribute you are about to add is not already present in the attribute list, add it
            chars_and_attr_list[i][1].append(protein_structure)
        i = i + 1


ID = ['P01833']

for i in ID:
    URL = 'https://www.uniprot.org/uniprot/' + i + '.txt'
    urllib.request.urlretrieve(URL, '\\Users\\Nora Kearns\\Downloads\\' + i + '.txt')

#download the text file
my_file = open("C:\\Users\\Nora Kearns\\Downloads\\" + i + '.txt', "r")
wholefile = my_file.read()
my_file.close()

#just get the protein sequence from the text file
proteinseq = wholefile.split('SQ   ')
labeled_seq = (proteinseq[1])

#get just the sequence
protseq1 = wholefile.split("CRC64;")
justthesequence = protseq1[1]

#print(justthesequence)
justtheseq = justthesequence.split('\n')
#print(justtheseq)
separator = ""
want = separator.join(justtheseq)
actual_sequence = "".join(want.split())
# This is the FASTA string
print(actual_sequence)

build_default_attr_structure_for_fasta(actual_sequence)

#get just the family and domains info
allthestuffsplitatFT = wholefile.split("Transmembrane helix")
justthestuffwithFT = allthestuffsplitatFT[1]
justFT = justthestuffwithFT.split("SQ   ")
FTs = justFT[0]
print(FTs)
arrayofFTs = FTs.split('/n')
#print(arrayofFTs)

#find the strands, domains, helices, and turns with the labeled numbers
matches = re.findall('FT\s+(STRAND|DOMAIN|HELIX|TURN)\s+(\d+..\d+)',FTs)
Big_array = []
for onematch in matches:

    # example: ('DOMAIN', '19..120') - onematch[0] is string for match type (DOMAIN)
    print(str(onematch))

    # example: ['19','120'] - range[0] and range[1] are strings for start and end indices
    range = onematch[1].split('..')
    #print(range)

    # print the sequence from the range given
    # example: DOMAIN AHTLMN
    each_item = (onematch[0] + ' ' + (actual_sequence[int(range[0]):int(range[1])]))

    add_attr_entries_for_zone(int(range[0]),int(range[1]),onematch[0])

print(chars_and_attr_list)

workbook = xlsxwriter.Workbook('C:\\Users\\Nora Kearns\\Downloads\\split_strings.xlsx')
worksheet = workbook.add_worksheet()

# Set up some formats to use.
bold = workbook.add_format({'bold': True})
italic = workbook.add_format({'italic': True})
black = workbook.add_format({'color': 'black'}) # LINKER
red = workbook.add_format({'color': 'red'}) # HELIX
blue = workbook.add_format({'color': 'blue'}) # STRAND
green = workbook.add_format({'color': 'green'}) # DOMAIN
navy = workbook.add_format({'color': 'navy'})  # TURN
magenta = workbook.add_format({'color': 'magenta'}) # MIX
purple = workbook.add_format({'color': 'purple'})

worksheet.set_column('A:A', 60) #make the cell 60 characters wide
col_heading_segments = [italic, 'COLOR LEGEND: ', black, 'LINKER  ', red, 'HELIX  ', green, 'DOMAIN  ', blue, 'STRAND  ', purple, 'TURN  ', magenta, 'MULTI  ']
worksheet.write_rich_string('A1', *col_heading_segments) #write the heading "color legend" in the first cell

segments = []
in_marked_area = False
row_count = 2 # Start at row 2 after legend in row 1
len_list = len(chars_and_attr_list)
i = 0
while (i < len(chars_and_attr_list)):
    existing_entry = chars_and_attr_list[i]
    existing_entry_char = existing_entry[0]
    existing_entry_attrs = existing_entry[1]
    if (len(existing_entry_attrs) > 1): # letter belongs to more than one class - use magenta for color
        segments.append(magenta)
        segments.append(existing_entry_char)  # APPEND THE LETTER TO WHAT WILL BE WRITTEN OUT
        in_marked_area = True
    elif ("HELIX" in existing_entry_attrs):
        segments.append(red)
        segments.append(existing_entry_char)  # APPEND THE LETTER TO WHAT WILL BE WRITTEN OUT
        in_marked_area = True
    elif ("DOMAIN" in existing_entry_attrs):
        segments.append(green)
        segments.append(existing_entry_char)  # APPEND THE LETTER TO WHAT WILL BE WRITTEN OUT
        in_marked_area = True
    elif ("STRAND" in existing_entry_attrs):
        segments.append(blue)
        segments.append(existing_entry_char)  # APPEND THE LETTER TO WHAT WILL BE WRITTEN OUT
        in_marked_area = True
    elif ("TURN" in existing_entry_attrs):
        segments.append(purple)
        segments.append(existing_entry_char)  # APPEND THE LETTER TO WHAT WILL BE WRITTEN OUT
        in_marked_area = True
    elif (len(existing_entry_attrs) == 0):
        if (in_marked_area == 1): # transition from marked to unmarked
            print("cut at letter: ", existing_entry_char, " at index ", i)
            row_ident = "A" + str(row_count) #creating a new row for a new segment, terminating at the cut line
            worksheet.write_rich_string(row_ident, *segments)
            row_count = row_count + 1
            in_marked_area = 0
            segments = []
            segments.append(black)
            segments.append(existing_entry_char)
        else:
            segments.append(black)
            segments.append(existing_entry_char)
    i = i + 1

# Write whatever is left over at the end, if anything, adding a blank to the end in case of a single char by itself
segments.append(" ")
row_ident = "A" + str(row_count)
worksheet.write_rich_string(row_ident, *segments)

workbook.close()

#     making_items_array = each_item.split(' ')
#     Big_array.append(making_items_array)
#     print(making_items_array)
# print(Big_array)
#
# #create a workbook and add a worksheet
# workbook = xlsxwriter.Workbook('Protein_domains.xlsx')
# worksheet = workbook.add_worksheet()
#
# #start from first cell, rows and columns are zero indexed
# row = 0
# col = 0
#
# #Iterate of the data and write it out row by row
# for i in Big_array:
#     worksheet.write(row, col, i[0])
#     worksheet.write(row, col+1, i[1])
#     row += 1
#
# workbook.close()
#