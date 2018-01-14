#####################################################################################################################
############################################ Set up Argparse ########################################################
#####################################################################################################################

import argparse

parser = argparse.ArgumentParser(description="DeDuper: Removes all PCR duplicates from a sorted SAM file. Identifies the chromosome, start position adjusted by soft clipping in the cigar strand, strand, and UMI for each read and uses those parameters to classify a duplicate.")
parser.add_argument("-f", "--file", help="Required parameter. Please pass a SORTED SAM file", required=True, type=str)
parser.add_argument("-p", "--paired", help="Optional parameter designates the file contains paired end reads. If paired end, pass the parameter 'True'", required=False, type=bool, default=False)
parser.add_argument("-u", "--umis", help="Optional parameter designates the file passed contains a list of valid UMIs. The file specified must be structured one UMI per line. Not specifying this parameter resorts to the program handling 'randomers'. In the case of randomers, randomers that include 'N' will be excluded from the output file.", required=False, type=str, default='')
parser.add_argument("-read_kept", "--read_kept", 
                    help="Optional parameter to specify the read to keep in the event of duplicates. Specifying 'Q' returns the read with the highest mapping quality. Default, by not specifying the parameter, returns the first read encountered in a pair/group of duplicates.", 
                    required=False, type=str, default='')

args = parser.parse_args() #sets arguments as call-able 
#print("Running on file:", args.filepath)

f = args.file
u = args.umis
paired = args.paired
duplicate_flag = args.read_kept

# f = "./Dataset1_sortSample.sam"
# u = "./KnownUMIs.txt"
# u = ""
# paired = False
# duplicate_flag = 'Q'

#####################################################################################################################
############################################ Define functions #######################################################
#####################################################################################################################

def valid_umisList(u):
    '''Creates and returns a list of valid UMIs from UMI file passed in parameter argument'''
    
    umis_list = [] #initiatilize dictionary of known UMIs
    with open(u, 'r') as valid_umiFile:
        umis_list = open(u).read().splitlines() #save each UMI as a new element in the list of known UMIs while removing new line characters/tabs
        return umis_list    

        
def bit_checker(bit, p):
    '''Takes the bitwise flag and checks it for strandedness. Assumes read is 
    mapped otherwise returns "None" and data are single-end. Returns "+" or "-" depending on strand.'''
    
    if (bit & 4) == 4: #if bit is equal to 4 the read is unmapped and therefore we dont want to continue with that read
        return None
    strand = "+"
    if (bit & 16) == 16: #changes strand if bit 16 is set
        strand = "-"            
    return strand 

    if p == True:
        if (bit & 40) == 40 and (bit & 80) != 80:
            strandedness = 1
        elif (bit & 40) != 40 and (bit & 80) == 80:
            strandedness = 2
        return strandedness
    
    
def valid_UMI_checker(currentUMI, umisList):
    '''When passed a single UMI extracted from a single SAM file read (self.umi) and a list of known UMIs, 
    function returns True if the read has a valid UMI, confirmed in the list. Otherwise, returns False. '''
    
    if currentUMI in umisList:
        return True
    else:
        return False

    
def adjust_softClipping(cigar):
    '''Takes in a single line's cigar parameter. The function evaluates the cigar for soft clipping present in
    the beginning of the read, specifically the first two indexes. Returns the number of positions soft clipped as
    an integer'''
    amt_adjs = 0
    if cigar[1]=="S":
        amt_adjs=cigar[0]        
    elif cigar[2]=="S":
        amt_adjs=cigar[:2]
        
    return(int(amt_adjs))    


#####################################################################################################################
####################################### Beginning of Main Arguments #################################################
#####################################################################################################################

#If an UMI list is specified via the user
if u != "": 
    #Create a list of valid UMIs 
    umisList = valid_umisList(u)
else:
    #No UMI list is specified, create an empty UMI list
    umisList = []

# If the -p parameter is passed by the user
if paired == True:
    raise TypeError("Program not equipped to handle paired-end reads at this time.")

#Customize output file being written to be "<file name passed by user>_deduped"
newFile = f.split("/")[-1].split(".")[1]
string_deDuped = "".join([newFile, "_deduped"])   
outputFile = f.replace(newFile, string_deDuped)

#Initialize counters, variables, and dictionaries
invalidUMI_ctr = 0
validUMI_ctr = 0
duplicate_ctr = 0
current_chromosome = ""
referenceDict = {} #Key: (umi, adj_pos, strand) as tuple
                   #Value: raw line
umiN_purgatory = []


with open(f, "r") as fh, open(outputFile, "w+") as output:
    i = 0
    # Initialize header list meant to contain all SAM file header lines
    header = []
    
    for line in fh:  
        #Bypass all headers that aren't SAM entries, but save header lines for new deduped SAM file later
        if line.startswith('@'): 
            header.append(line)
            i+=1
        
        # We've found a SAM file read of interest
        if line.startswith('NS5'): 
            i+=1
            # Extract all fields of interest from that line
            umi = line.split("\t")[0].split(':')[7]
            bit = int(line.split("\t")[1])
            strand = bit_checker(bit, paired)
            chromosome = line.split("\t")[2]
            cigar = str(line.split("\t")[5])
            start_pos = int(line.split("\t")[3])
            adj_pos = start_pos - adjust_softClipping(cigar)
            mapQ = line.split("\t")[4]
            
            #If the user passed a file of known UMIs
            if umisList != []: 
                #Check if the UMI for the current line is valid compared to the UMI file passed (Function returns True if valid)
                umiStatus = valid_UMI_checker(umi, umisList) 
                if umiStatus == True:  
                    validUMI_ctr+=1
                    #If the line being read has a new chromosome than previous
                    if chromosome != current_chromosome: 
                        # Write all valid reads tracked in dictionary to output file
                        for value in referenceDict.values():
                            output.write(value)
                        
                        #Set the new chromosome
                        current_chromosome = chromosome 
                        #Reset dictionary to save memory and run program more effciently: part 1
                        referenceDict = referenceDict.clear() 
                        #Reset dictionary: part 2 (need this step because clear() makes dictionary Nonetype)
                        referenceDict = {} 
                    
                    #Create "key" of type 'tuple' for dictionary entry
                    tupleKey = (umi, adj_pos, strand) 
                    
                    #If the dictionary is empty, meaning a new chromosome was reached
                    if any(referenceDict) == False: 
                        #Add the first line read pertaining to the new chromosome in to the dictionary
                        referenceDict[tupleKey] = line 
                    
                    #If the dictionary is not empty, meaning the current line still has the same chromosome
                    if referenceDict != {}:
                        #Check if the -read_kept flag indicated specifies keeping the highest mapping quality duplicate
                        if duplicate_flag == 'Q':
                            #If a line with the same (umi, position, strand) exists as a key in the dictionary already
                            if tupleKey in referenceDict:
                                duplicate_ctr+=1 
                                #Meaning two lines were  exactly equivalent, the line in the dictionary suffices
                                if line in referenceDict[tupleKey]: 
                                    continue
                                else: 
                                    #Extract the mapping quality field out of the line in the dictionary entry
                                    dict_mapQ = referenceDict[tupleKey].split("\t")[4]
                                    #If mapping quality in dictionary entry is equal or higher than current line's, keep dictionary entry
                                    if dict_mapQ >= mapQ:
                                        pass
                                    else:
                                        #Current line has higher mapping quality, replace current dictionary entry
                                        referenceDict[tupleKey] = line
                            else:
                                #Not a duplicate line, add to dictionary
                                referenceDict[tupleKey] = line 
                                                
                        
                        #Duplicate flag not specified, -read_kept defaults to first read kept
                        else: 
                            if tupleKey in referenceDict: 
                                duplicate_ctr+=1 
                                #Meaning UMI, position, and strand of a line to an existing line already read were equivalent 
                                if line in referenceDict[tupleKey]: 
                                    continue
                            else:
                                #Not a duplicate line, add to dictionary
                                referenceDict[tupleKey] = line 
                
                else:
                    # List of valid UMIs were specified and the umi of the current line was invalid
                    invalidUMI_ctr+=1
 


            #No list of valid UMIs were provided
            else: 
                #If current line's chromosome is different from 'current chromosome' being processed
                if chromosome != current_chromosome:
                    # Write all valid reads tracked in dictionary to output file
                    for value in referenceDict.values():
                        output.write(value)
                    #Set the new chromosome to 'current chromosome'
                    current_chromosome = chromosome
                    #Reset dictionary for new chromosome and to save memory/program runtime
                    referenceDict = referenceDict.clear()
                    referenceDict = {}
                
                #Create "key" of type tuple for dictionary entry    
                tupleKey = (umi, adj_pos, strand) 
                
                #If the dictionary is empty, meaning a new chromosome was reached
                if any(referenceDict) == False: 
                    #Add the first line with a new chromosome to the dictionary
                    referenceDict[tupleKey] = line 
                    
                # Line being read still has same chromosome as previous    
                if referenceDict != {}:
                    #Check if the -read_kept flag indicated specifies keeping the highest mapping quality duplicate
                    if duplicate_flag == 'Q':
                        if tupleKey in referenceDict: 
                            duplicate_ctr+=1 
                            #Meaning UMI, position, and strand of a line to an existing line were equivalent, current entry in dictionary suffices
                            if line in referenceDict[tupleKey]: 
                                continue
                            else: 
                                #Extract the mapping quality field out of the line in the dictionary entry
                                dict_mapQ = referenceDict[tupleKey].split("\t")[4]
                                #If mapping quality in dictionary entry is equal or higher than current line's, keep dictionary entry
                                if dict_mapQ >= mapQ:
                                    pass
                                else:
                                    # Add current line to dictionary to replace one with lower mapping quality
                                    referenceDict[tupleKey] = line
                        else:
                            #Not a duplicate line, add to dictionary
                            referenceDict[tupleKey] = line 
                            
                    #Duplicate kept flag not specified, defaults to first read being kept
                    else: 
                        if tupleKey in referenceDict: 
                            duplicate_ctr+=1   
                            #Meaning UMI, position, and strand of a line to an existing line were equivalent 
                            if line in referenceDict[tupleKey]: 
                                continue
                        else:
                            #If current UMI contains an 'N', toss the read, the randomer is invalid
                            if "N" in umi:
                                umiN_purgatory.append(line)
                                invalidUMI_ctr+=1
                            else:
                                #Not a duplicate line, add to dictionary
                                referenceDict[tupleKey] = line 


  

    #Write header lines to output SAM file
    if header != []: 
        output.write(''.join(header))
    #In the special case the whole SAM file has only reads on one chromosome
    if referenceDict != {}: 
        for value in referenceDict.values():
            output.write(value)               


#####################################################################################################################
############################### Print counters for debugging and reassurance ########################################
#####################################################################################################################

print("Invalid UMIs:\t\t", invalidUMI_ctr)           
print("Valid UMIs:\t\t", validUMI_ctr) 
print("Header lines:\t\t", len(header))
print("Total lines:\t\t", i)
print("Duplicates identified:\t", duplicate_ctr) 
#print("UMI purgatory:\n", umiN_purgatory)
    
