#!/usr/bin/env python

import os
from Bio import Entrez
import xml.etree.ElementTree as ET
from optparse import OptionParser
import pandas as pd

usage= """
This script will take a list of numeric NCBI taxon IDs and return the taxonomic lineage for those taxa. 
\n The lineage file returned can be used as input into METATRYP for LCA analysis. 
\n Input files should be .csv files with the following headers:
\n Idividual genomes: Taxon_fasta_file_name, Taxon_name,  NCBI_taxon_id
\n Metagenomes: Taxon_fasta_file_name, ORF_name, Taxon_name, NCBI_taxon_id
\n An e-mail address is also required as input to access NCBI API tools. 
usage: %prog [-f FILE] [-e STR] [-o STR] """
parser = OptionParser(usage=usage, version="%prog 0.1")

parser.add_option("-f", dest="taxaID_file",
                  help="Specify the input file of numeric taxon IDs.")
parser.add_option("-e", dest="email",
                  help="Specify a valid email to be used to access NCBI API tools.")
parser.add_option("-o", dest="out_name",
                  help="(OPTIONAL) Specify a preferred output file name for the .csv output file.")

(options, args) = parser.parse_args()
#Makes sure all mandatory options appear
mandatories = ["taxaID_file", "email"]
for m in mandatories:
    if not options.__dict__[m]:
        print("A mandatory option is missing!\n")
        parser.print_help()
        exit(-1)

#Opens the blast output info file and the new output file
taxaParse = os.path.abspath(options.taxaID_file)
df = pd.read_csv(taxaParse, low_memory=False)

if options.out_name is None:
    outputFileName = "NCBI_taxonomic_lineage.csv"
else:
    outputFileName = str(options.out_name) + ".csv"

email_api = str(options.email)

taxaList = df["NCBI_taxon_id"].astype('str').str.replace(' ','').to_list()
stringTaxaID = ','.join(list(set(taxaList)))

#Need to specify an e-mail to get data from NCBI... they limit how much you can hit their servers
Entrez.email = email_api

#Pull taxonomy information and output as xml file
handle = Entrez.efetch(db="taxonomy", id=stringTaxaID, format="xml")

#Create a new file to store all parsed taxonomic lineage information
taxaKeysInfo = open("NCBI_taxon_temp.xml", "w")

#Read the XML output... need to create a new xml file, as the xml parsing package ElementTree requires in in a file
###This might be an unnecessary step
readXMLFile = handle.readlines()

for line in readXMLFile:
    taxaKeysInfo.write(line)

taxaKeysInfo.close()

#Use ElementTree to read XML and parse desired information
### The Entrez package has an xml parser, but not sure how flexible it is... 
tree = ET.parse("NCBI_taxon_temp.xml") #temp file will be deleted at end of program
root = tree.getroot()

#Create a file which will store all the parsed taxonomy info
newFile = open(outputFileName, "w")

totalTaxa = len(root) 
taxonEntry = 0
taxaLineageDict = {}
lineageDoesNotExist = []

#while taxonEntry < totalTaxa:
while taxonEntry < totalTaxa:
    taxonKey = root[taxonEntry][0].text
    orgName = str(root[taxonEntry][1].text)
    lineageTree = root[taxonEntry].find('LineageEx')
    kingdomLoc = None
    phylumLoc = None
    classLoc = None
    orderLoc = None
    familyLoc = None
    genusLoc = None
    speciesLoc = None
    rankLevel = root[taxonEntry].find('Rank').text
    

    if lineageTree != None:    
        index = 0
        for node in lineageTree.findall('.//Taxon/Rank'): #This finds the location index of the lineage levels
           # print node.tag, node.text, index
            if node.text == "superkingdom":
                kingdomLoc = index
            elif node.text == "phylum":
                phylumLoc = index
            elif node.text == "class":
                classLoc = index
            elif node.text == "order":
                orderLoc = index
            elif node.text == "family":
                familyLoc = index
            elif node.text == "genus":
                genusLoc = index
            elif node.text == "species":   
                speciesLoc = index  #Using species in lineage for species

            index += 1

        taxonEntry += 1

        Group = lineageTree[0][1].text

        taxaValuesList = [Group]

        lineageLevels = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom', 'Group']
        taxaList = [kingdomLoc, phylumLoc, classLoc, orderLoc, familyLoc, genusLoc, speciesLoc]
        subList = []
        firstOccurrence = True
        

        #Check to see if taxaList is filled with 'None'
        allSame = taxaList.count(taxaList[0]) == len(taxaList)
        if allSame == True:
            taxaList = [None]

        counter = 8
        for item in taxaList:
            counter = counter - 1

            if taxaList[0] is None:

                if rankLevel == 'superkingdom':
                    UnclassX = ("Unclassified " + str(orgName))
                    newList = [orgName] + [UnclassX] * 6
                    taxaValuesList = [taxaValuesList[0]] + newList
                    subList.append(orgName)

                else:
                    taxaValuesList = ["Unknown"] * 8

            elif item is None:

                if firstOccurrence == True:
                    firstOccurrence = False #reset
                    if str(rankLevel) == str(lineageLevels[counter - 1]):
                        
                        taxaValuesList.append(orgName)
                        subList.append(orgName)

                    else:
                        UnclassX = ("Unclassified " + str(taxaValuesList[-1]))
                        taxaValuesList.append(UnclassX)
                        subList.append(UnclassX)
                
                else:

                    if str(rankLevel) == str(lineageLevels[counter - 1]):
                        taxaValuesList.append(orgName)
                        subList.append(orgName)
                    elif str(subList[-1]).find("Unclassified") >= 0:
                        UnclassX = (str(subList[-1]))
                        taxaValuesList.append(UnclassX)
                        subList.append(UnclassX)
                    else:
                        UnclassX = ("Unclassified " + str(taxaValuesList[-1]))
                        taxaValuesList.append(UnclassX)
                        subList.append(UnclassX)

            #Pulls this lineage level information for the specified index location
            else:
                subList.append(lineageTree[item][1].text)
                taxonomyLevel = lineageTree[item][1].text
                taxaValuesList.append(taxonomyLevel)

        firstOccurrence = True
        counter = 8
        taxaValuesList.append(orgName)
        taxaLineageDict[str(taxonKey)] = taxaValuesList

    else:
        lineageDoesNotExist.append(taxonKey)
        taxaLineageDict[str(taxonKey)] = ["Unknown"] * 9
        taxonEntry += 1

taxaLineageDict[str('-1')] = ["Unknown"] * 9
#Create dictionary where the NCBI taxonomy ID = key; lineage = values
allKeys = taxaLineageDict.keys()

if len(lineageDoesNotExist) > 0:
    print("The following taxon IDs were not found in the NCBI Taxonomy database API:")
    print(list(set(lineageDoesNotExist)))


#Convert dictionary of lineages to dataframe and merge with input file info
df_lineage = pd.DataFrame.from_dict(taxaLineageDict, orient='index', columns= \
                                    ["Group","Kingdom","Phylum", \
                                     "Class","Order","Family","Genus","Species","NCBI_taxon_name"])
df_lineage['NCBI_taxon_id'] = df_lineage.index
df_lineage['NCBI_taxon_id'] = df_lineage['NCBI_taxon_id'].astype(str)
df['NCBI_taxon_id'] = df['NCBI_taxon_id'].astype(str)
df_out = pd.merge(df, df_lineage, how='left', on='NCBI_taxon_id')


df_out.to_csv(outputFileName, index=False)

#Remove temp xml lineage file from NCBI
os.remove("NCBI_taxon_temp.xml")
