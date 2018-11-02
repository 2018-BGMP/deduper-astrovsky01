import argparse
import numpy as np
import re
pcr_store = dict() #Stores non-pcr duplicate strands by UMI
headers = []
umifile = ""
def info():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help = "SAM file to be read in",
        required = True, type=str)
    parser.add_argument('-p', '--paired', help = "Optional: defines the import file as paired if given, s for single end, p for paired end. Default is s",
        required = False, action="store_true")
    parser.add_argument('-u', '--umifile', help = "Optional file with list of UMIs, unset if randomers",
        required = False, type=str, default = "ASDF")
    return parser.parse_args()
args = info()
file = args.file
paired = args.paired
umifile = args.umifile
out_file = str(file) + "_dedup"
provided = 0
if paired == True:
    exit("Error: Unable to use paired end data right now.")
if umifile != "ASDF":
    print("File provided")
    provided = 1
    with open(umifile, 'r') as umi:
        for x in umi:
            x = x.strip()
            pcr_store[x] = []
else:
    print("No umi file given")

def linesep(sam_line):
    """Creates array from line in samfile"""
    sam_line = sam_line.strip()
    x = re.split('\t|:', sam_line)
    #Order of array: [UMI, bitflag, chromosome, start position, CIGAR, line]
    curr_line = [x[7][-8:], x[8], x[9], x[10], x[12], sam_line]
    return(curr_line)

def pos_adj(position, cigar, bitflag):
    """Adjusts position for indels, splicing, etc."""
    position = int(position)
    new = re.split("([0-9]{1,10}[A-Z])", cigar) #Breaks cigar into each part
    new = list(filter(lambda a: a != "", new)) #Removes blanks due to re.splt
    change = 0
    if (int(bitflag) & 16) == 16:
        reverse = "R"
    else:
        reverse = "F"
    if reverse == "R":
        if "S" in new[-1]:
            change += int(new[-1][:-1])
        for x in range(0,len(new)):
            if "D" in new[x]:
                change += int(new[x][:-1]) #Adds in deletions to adjust start position of reverse strand
            if "M" in new[x]:
                change += int(new[x][:-1]) # Adds in mapped to adjust start position of reverse strand
            if "N" in new[x]:
                change += int(new[x][:-1]) # Adds in alt. splicing event to adjust start position of reverse strand
        adjusted = position + change
    else:
        if "S" in new[0]:
            change += int(new[0][:-1])
        adjusted = position - change
    return(adjusted, reverse)

def new_pcr_test(umi, chromosome, position, direction):
    """Compares current sam string to all lines already stored for that UMI"""
    for value in pcr_store[umi]:
        if chromosome == value[0] and position == value[1] and direction == value[2]:
            print("PCR dup")
            return(False)
    return(True)

with open(file, 'r') as f:
    for x in f:
        if x[0] == "@":
            headers.append(x)
            continue
        sam_line = x #Read in SAM string
        cur_line = linesep(sam_line) #Break into constituent parts
        umi = cur_line[0]
        cur_line = cur_line[1:]
        adjusted, reverse = pos_adj(cur_line[2], cur_line[3], cur_line[0])
        cur_line = [cur_line[1], adjusted, reverse, sam_line.strip()] #creates list of chromosome, adjusted position, direction, total sam line
        if provided == 1: #Given umi file
            if umi not in pcr_store.keys(): #Novel umi not in file
                print("Non-valid umi:", umi)
                continue
            if umi in pcr_store: #string from a given umi
                if new_pcr_test(umi, cur_line[0], cur_line[1], cur_line[2]) == True: #If novel read for an umi
                    pcr_store[umi].append(cur_line)
        if provided == 0: #No umi file given
            if umi not in pcr_store.keys(): #Novel umi
                pcr_store[umi] = []
                pcr_store[umi].append(cur_line) #add new dictionary call for specific umi
            elif umi in pcr_store.keys(): #Previously seen UMI
                if new_pcr_test(umi, cur_line[0], cur_line[1], cur_line[2]) == True: #Novel read of already seen UMI
                    pcr_store[umi].append(cur_line)
with open(out_file, "w") as out:
    for x in headers:
        out.write(x)
    for x in pcr_store:
        for n in pcr_store[x]:
            out.write(n[-1])
            out.write('\n')
