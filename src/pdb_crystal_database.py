import sys, json, pickle, traceback, os
from time import sleep
from collections import OrderedDict
from pathlib import Path
import xml.etree.ElementTree as etree

# Make sure directories exist
if not os.path.exists("Structures"):
    os.makedirs("Structures")

# Base project directory
BASE_DIR = Path("/home/ruth/Protein_Crystallization_Data_Extraction")

# Subdirectories
INPUT_DIR = BASE_DIR / "Input"
STRUCTURE_DIR = BASE_DIR / "Structures"

# Input Files
COMPOUND_DICTIONARY_FILE = INPUT_DIR / "compound_dictionary.json"
SMILES_DICTIONARY_FILE = INPUT_DIR / "smiles_dictionary.json"

STOP_WORDS_FILE = INPUT_DIR / "stop_words.json"
UNKNOWN_LIST_FILE = INPUT_DIR / "unknown_list.json"
LOWERCASE_REPLACEMENT_FILE = INPUT_DIR / "replacementLowercase.json"
SENSITIVE_REPLACEMENT_FILE = INPUT_DIR / "replacementSensitive.json"
MIXTURES_FILE = INPUT_DIR / "mixture_compounds.json"
CLASSIFICATION_FILE = INPUT_DIR / "classification_dictionary.json"

STRUCTURES_FILE = STRUCTURE_DIR / "structures.pkl" # The database file. This is the main output of this code, and is required for the rest of the pipeline to run.
SENSIBLE_STRUCTURES_FILE = STRUCTURE_DIR / "sensible_structures.pkl" # 

def loadJson(jsonFile, object_pairs_hook=None):
    """Loads a json file, returns whatever python object is in the file
    Optional arguments are the same as in json.load"""
    with open(jsonFile, "r") as inFile:
        return json.load(inFile, object_pairs_hook=object_pairs_hook)
    
# Load lists and dictionaries from files
print("Loading input files...")
try:
    unknownList = loadJson(UNKNOWN_LIST_FILE)
    compoundDictionary = loadJson(COMPOUND_DICTIONARY_FILE)
    smilesDictionary = loadJson(SMILES_DICTIONARY_FILE)
    lowercaseReplacement = loadJson(LOWERCASE_REPLACEMENT_FILE, object_pairs_hook=OrderedDict)
    sensitiveReplacement = loadJson(SENSITIVE_REPLACEMENT_FILE, object_pairs_hook=OrderedDict)
    mixturesDictionary = loadJson(MIXTURES_FILE)
    classificationDictionary = loadJson(CLASSIFICATION_FILE)
    STOP_WORDS = set(loadJson(STOP_WORDS_FILE)) - {'m'} - {'am'}
except FileNotFoundError as notFound:
    print(notFound)
    print("ERROR: The file {} cannot be found. Verify that it is in the proper directory.".format(notFound.filename))

# Global variables set by exportOutputFiles()
sensibleStructureList = []
nonsensibleStructureList = []

# Compounds that contain numbers (eg jeffamine 600)
NUMBERED_COMPOUNDS = ["jeffamine", "propoxylate", "polypropylene", "ndsb"]

def printList(list, delimiter=", "): # string
    """Takes a list and returns a string deliminated by a specific substring"""
    string = ""
    count = 1
    for element in list:
        string += element
        if (count != len(list)):
            string += delimiter
        count += 1
    return string

def listToFile(list, filename): # void
    """Writes a list to a text file, with each element being a line in the file"""
    with open(filename, "w") as outFile:
        for element in list:
            outFile.write(str(element)+"\n")

def getKey(compound): # string
    """Takes a string name of a compound and removes certain characters to make a simplified dictionary key"""
    charactersToRemove = [" ", "(", ")", "-", ",", "/", "*", "_"]
    for char in charactersToRemove:
        compound = compound.replace(char, "")
    return compound.lower()

class Structure:
    """
    Stores crystallization + compound information for an X-ray structure.
    """

    # Legacy filtering keywords (optional future use)
    RESEVOIR_INDICATOR_WORDS = {
        "reservoir", "resevoir", "crystallization",
        "precipitant", "well", "solution", "cocktail"
    }

    def __init__(self, pdbid, details, compounds ):

        self.pdbid = pdbid
        self.details = details
        # ⚠️ FIX: do NOT overwrite compounds
        self.compounds = compounds if compounds is not None else []

        # Ensure only X-ray diffraction entries are kept
        if not self._is_xray():
            raise ValueError(f"{pdbid} is not X-ray diffraction")

    # ─────────────────────────────────────────────
    # X-ray filter
    # ─────────────────────────────────────────────
    def _is_xray(self):
        return (
            self.method is not None and
            "x-ray diffraction" in self.method.lower()
        )

    # ─────────────────────────────────────────────
    # String representation
    # ─────────────────────────────────────────────
    def __str__(self):
        return f"""
        PDB: {self.pdbid}
        Compounds: {self.compounds}
        Details: {self.details}
        """

    def parseDetails(self, detailsString=None, debug=False):
        """Parses the details string and returns a list of compounds, followed by concentration, or None if conc. is not found
        Format of compounds = ['compound name', '100', 'another compound name', '45%']
        By default (if detailsString=None), the input will be self.details, but
        detailsString can be set to another string instead
        If debug=True, then this function will print the words after every step
        Also sets the compounds list of the structure
        Returns None if failed or if details are unavailible
        """
        global nltk

        # Make sure NLTK is imported
        importNLTK()

        if (detailsString==None):
            details = self.details
        else:
            details = detailsString

        if details == None:
            return None

        if debug:
            print("Debugging PDB {}...\n".format(self.pdbid))
            print("Raw details:\n"+details+"\n")

        # Add spaces after commas before numbers (eg "sodium acetate,4% PEG4K" --> "sodium acetate, 4% PEG4K")
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(1, len(details)-1):
                if details[j] == "," and details[j-1].isalpha() and isNumber(details[j+1]):
                    # Insert a space
                    details = details[:j+1] + " " + details[j+1:]
                    runAgain = True
                    break

        if debug:
            print("Dealing with lack of space after commas:\n"+details+"\n")

        details = wordReplacement(details, sensitiveReplacement, debug=debug) # Replace with case sensitivity
        details = wordReplacement(details.lower(), lowercaseReplacement, debug=debug) # Replace without case sensitivity

        if debug:
            print("Word replacement:\n"+details+"\n")

        # Remove all text after 'cryo'
        cryoLocation = details.find("cryo")
        if cryoLocation != -1:
            details = details[:cryoLocation]

        # Removed all text after 'soaked'
        soakedLocation = details.find("soaked")
        if soakedLocation != -1:
            details = details[:soakedLocation]

        # Remove all text before "reservoir solution" to include only resevoir soltuion
        reservoirLocation = details.find("reservoir solution")
        if reservoirLocation != -1:
            details = details[reservoirLocation:]
            # Remove protein solution if it comes after reservoir solution
            proteinLocation = details.find("protein solution")
            if proteinLocation != -1:
                details = details[:proteinLocation]

        if debug:
            print("Isolate resevoir solution:\n"+details+"\n")


        words = nltk.word_tokenize(details)

        if debug:
            print("Tokenized words:\n"+str(words)+"\n")

        if debug:
            print("Remove protein solution:\n"+str(words)+"\n")

        # Fix commas between numbers (turn 6000,10 in to 6000 , 10)
        runAgain = True
        while (runAgain):
            runAgain = False
            for j in range(0, len(words)):
                if "," in words[j] and isNumber(words[j].replace("%","").replace("m","").replace("k","")):
                    numbers = words[j].split(",")
                    if len(numbers) == 2 and (len(numbers[1]) != 3 or len(numbers[0]) > 3 or "." in numbers[1]) and numbers[0] != "" and numbers[1] != "":
                        runAgain = True
                        words[j] = numbers[0]
                        words.insert(j+1, ",")
                        words.insert(j+2, numbers[1])
                        if len(numbers[1]) > 2 and numbers[1][-2:] == "mm":
                            words[j+2] = numbers[1][:-2]
                            words.insert(j+3, "mm")
                        elif numbers[1][-1:] == "m" or numbers[1][-1:] == "%":
                            words[j+2] = numbers[1][:-1]
                            words.insert(j+3, numbers[1][-1])
                        break

        if debug:
            print("Fix commas:\n"+str(words)+"\n")

        # Process micromolar (uM) concentrations, add to number like: "35 um" --> "35um"
        runAgain = True
        while (runAgain):
            runAgain = False
            for j in range(1, len(words)):
                if words[j] == "um" and isNumber(words[j-1]):
                    words[j-1] = words[j-1]+"um"
                    del words[j]
                    runAgain = True
                    break

        if debug:
            print("Process uM concentrations:\n"+str(words)+"\n")

        # Remove pH and temperature
        runAgain = True
        while (runAgain): # pH
            runAgain = False
            for j in range(0, len(words)):
                if j+1 < len(words) and words[j] == "ph" and isNumber(words[j+1]):
                    del words[j:j+2]
                    runAgain = True
                    break
                elif len(words[j]) > 2 and words[j][:2] == "ph" and isNumber(words[j][2:].replace("=","").replace(" ","")):
                    del words[j]
                    runAgain = True
                    break
                elif j+2 < len(words) and words[j] == "ph" and words[j+1] == "=" and isNumber(words[j+2]):
                    del words[j:j+3]
                    runAgain = True
                    break
        runAgain = True
        while (runAgain): # Temperature
            runAgain = False
            for j in range(1, len(words)):
                if (words[j][-1] == 'k' and len(words[j]) > 3 and isNumber(words[j][:-1])):
                    del words[j]
                    runAgain = True
                    break
                elif (words[j][-1] == 'k' and len(words[j]) > 1 and isNumber(words[j][:-1]) and words[j-1] == "temperature"):
                    del words[j]
                    runAgain = True
                    break
                elif words[j] == 'k' and isNumber(words[j-1]):
                    temp = float(words[j-1].replace(",",""))
                    if temp > 200 and temp < 400:
                        del words[j-1:j+1]
                        runAgain = True
                        break

        if debug:
            print("Remove pH and temperature:\n"+str(words)+"\n")

        # Process compounds with numbers in them
        for j in range(0, len(words)-1):
            if words[j] in NUMBERED_COMPOUNDS and isNumber(words[j+1]):
                words[j] = words[j] + " " + words[j+1]
                del words[j+1]
                break

        # Process PEG and MPEG compounds
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)):
                if (words[j][:3] == "peg" or words[j][:4] == "mpeg"):
                    pegNumber = ""
                    if j+1<len(words) and words[j+1] in STOP_WORDS:
                        runAgain = True
                        del words[j+1]
                        break
                    elif len(words[j]) > 3 and words[j][:3] == "peg": # If the number is already attached, eg PEG5000
                        if j+1<len(words) and words[j+1] == "mme": # Ex. PEG3000 MME
                            pegNumber = "MME " + words[j][3:].replace(",","").replace("-", " ")
                            del words[j+1]
                        else:
                            pegNumber = words[j][3:].replace("-", " ").replace(",","")
                    elif len(words[j]) > 4 and words[j][:4] == "mpeg": #MPEG4000
                        pegNumber = words[j][4:].replace("-", " ").replace(",","")
                    elif j+2 < len(words) and (words[j+1] == "mme" or words[j+1] == "monomethyl") and words[j+2] not in STOP_WORDS: # ex. PEG MME 5000
                        if isNumber(words[j+2].replace("k","000")):
                            pegNumber = "MME " + words[j+2].replace(",","")
                        else:
                            pegNumber = "MME"
                        del words[j+2]
                        del words[j+1]
                    elif j+2 < len(words) and words[j+2] == "mme" and words[j+1] not in STOP_WORDS: # ex. PEG 5000 MME
                        if isNumber(words[j+1].replace("k","000")):
                            pegNumber = "MME " + words[j+1].replace(",","")
                        else:
                            pegNumber = "MME"
                        del words[j+2]
                        del words[j+1]
                    elif j+1<len(words) and words[j+1] not in STOP_WORDS: # Ex PEG 5000
                        if isNumber(words[j+1].replace("k","000")):
                            pegNumber = words[j+1].replace(",","")
                            del words[j+1]
                    if len(pegNumber) > 3 and pegNumber[-3:] == "mme":
                        pegNumber = "MME " + pegNumber[:-3]
                    pegNumber = pegNumber.replace("000k","000").replace("k","000").replace("-","")
                    if words[j][:4] == "mpeg":
                        words[j] = "MPEG "+pegNumber
                    else:
                        words[j] = "PEG "+pegNumber
                    runAgain = True
                    break

        if debug:
            print("Process PEG and MPEG compounds:\n"+str(words)+"\n")

        # Remove terminal periods from numbers
        for j in range(0, len(words)):
            if isNumber(words[j]) and words[j][-1] == ".":
                words[j] = words[j][:-1]

        if debug:
            print("Remove terminal periods from numbers:\n"+str(words)+"\n")

        # Put together numbers separated by "-"
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(1, len(words)-1):
                if words[j] == "-" and isNumber(words[j-1]):
                    if isNumber(words[j+1]) or (words[j+1][-1] == "m" and isNumber(words[j+1][:-1])) or (len(words[j+1]) > 2 and words[j+1][-2] == "mm" and isNumber(words[j+1][:-2])):
                        words[j-1] = printList(words[j-1:j+2], "")
                        del words[j:j+2]
                        runAgain = True
                        break

        if debug:
            print("Combine numbers separated by '-':\n"+str(words)+"\n")

        # Average ranges of numbers
        for j in range(0, len(words)):
            if '-' in words[j]:
                if averageNumberString(words[j]) != None: # If there is no letter in front of the numbers
                    words[j] = str(averageNumberString(words[j]))
                elif averageNumberString(words[j][:-1]) != None: # If there is one letter in front of the numbers (e.g. 'm' or '%')
                    words.insert(j+1, words[j][-1])
                    words[j] = str(averageNumberString(words[j][:-1]))
                elif averageNumberString(words[j][:-2]) != None: # If there are two letters in front of the numbers (i.e. Mm)
                    words.insert(j+1, words[j][-2:])
                    words[j] = str(averageNumberString(words[j][:-2]))

        # Remove temperature again
        runAgain = True
        while (runAgain): # Temperature
            runAgain = False
            for j in range(1, len(words)):
                if (words[j][-1] == 'k' and len(words[j]) > 3 and isNumber(words[j][:-1])):
                    del words[j]
                    runAgain = True
                    break
                elif (words[j][-1] == 'k' and len(words[j]) > 1 and isNumber(words[j][:-1]) and words[j-1] == "temperature"):
                    del words[j]
                    runAgain = True
                    break
                elif words[j] == 'k' and isNumber(words[j-1]):
                    temp = float(words[j-1].replace(",",""))
                    if temp > 200 and temp < 400:
                        del words[j-1:j+1]
                        runAgain = True
                        break

        if debug:
            print("Average ranges of numbers and remove temperature again:\n"+str(words)+"\n")

        if len(words) == 0:
            return None

        # Process concentrations
        runAgain = True
        while(runAgain):
            for j in range(0, len(words)):
                if words[j][-1] == "m" and isNumber(words[j][:-1]): # Units are M
                    words[j] = str(float(words[j][:-1].replace("%","").replace(",",""))*1000)
                    runAgain = True
                    break
                elif words[j][-2:] == "mm" and isNumber(words[j][:-2]): # Units are mM
                    words[j] = words[j][:-2]
                    runAgain = True
                    break
                elif words[j][-2:] == "um" and isNumber(words[j][:-2]): # Units are uM
                    words[j] = str(float(words[j][:-2].replace("%","").replace(",",""))/1000)
                    runAgain = True
                    break
                elif isNumber(words[j]) and j+1 < len(words): # Conc. is separated from number by space, e.g. "0.1 M"
                    if words[j+1] == "m":
                        del words[j+1]
                        words[j] = str(float(words[j].replace(",",""))*1000)
                        runAgain = True
                        break
                    elif words[j+1] == "mm":
                        del words[j+1]
                        words[j] = str(float(words[j].replace(",","")))
                        runAgain = True
                        break
                    elif words[j+1] == "%":
                        if j+2 < len(words) and (words[j+2] == "w/v" or words[j+2] == "v/v"):
                            units = "% " + words[j+2]
                            del words[j+2]
                        else:
                            units = "%"
                        del words[j+1]
                        words[j] = str(words[j]+units)
                        runAgain = True
                        break
                else:
                    runAgain = False

        if debug:
            print("Process concentrations:\n"+str(words)+"\n")

        # Remove excess punctuation
        for j in range(0, len(words)):
            if words[j] in {",", "."}:
                if all(word in {",", "."} for word in words[j:]):
                    del words[j:]
                    break

        if debug:
            print("Remove excess punctuation:\n"+str(words)+"\n")

        # Combline parts of chemical compounds
        runAgain = True
        while runAgain:
            runAgain = False
            for j in range(0, len(words)-1):
                if isCompound(words[j]) and isCompound(words[j+1]) and words[j][:4] != "PEG " and words[j+1][:4] != "PEG " and words[j][:4] != "MPEG " and words[j+1][:4] != "MPEG ":
                    k = j
                    chemicalCompound = []
                    while(k < len(words) and isCompound(words[k])):
                        chemicalCompound.append(words[k])
                        k += 1
                    del words[j+1:k]
                    words[j] = printList(chemicalCompound, " ")
                    runAgain = True
                    break

        if debug:
            print("Combine parts of chemical compounds:\n"+str(words)+"\n")

        # Remove periods from compounds
        for j in range(0, len(words)):
            if isCompound(words[j]):
                words[j] = words[j].replace(".","")

        if debug:
            print("Remove periods from compounds:\n"+str(words)+"\n")

        # Combine w/v or v/v with non-percent nubers, and assume it's supposed to be percent
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)-1):
                if isNumber(words[j]) and (words[j+1] == "w/v" or words[j+1] == "v/v"):
                    words[j] = words[j] + "% " + words[j+1]
                    del words[j+1]
                    runAgain = True
                    break

        # If percent concentration is before and w/v is after
        runAgain = True
        while(runAgain):
            runAgain = False
            for j in range(0, len(words)-2):
                if isPercent(words[j]) and isCompound(words[j+1]) and (words[j+2] == "w/v" or words[j+2] == "v/v"):
                    words[j] = words[j] + " " + words[j+2]
                    del words[j+2]
                    runAgain = True
                    break

        if debug:
            print("Parse w/v and v/v:\n"+str(words)+"\n")

        # Determine if concentration or compounds come first
        concentrationBeforeCompound = True
        for j in range(0, len(words)-1):
            if isConcentraton(words[j]) and isCompound(words[j+1]) and getKey(words[j+1]) in compoundDictionary:
                concentrationBeforeCompound = True
                break
            elif isConcentraton(words[j+1]) and isCompound(words[j]) and getKey(words[j]) in compoundDictionary:
                concentrationBeforeCompound = False
                break

        # Remove stop words
        words = [w for w in words if w not in STOP_WORDS]

        if debug:
            print("Remove stop words:\n"+str(words)+"\n")

        # Add compounds to list
        compounds = []
        if len(words) > 0 and words[0] != "ERROR":
            runAgain = True
            while(runAgain):
                runAgain = False
                for j in range(0, len(words)):
                    if isCompound(words[j]):
                        if concentrationBeforeCompound and j > 0 and isConcentraton(words[j-1]): # If a concentration is found (before compound)
                            compound = words[j]
                            concentration = words[j-1]
                            # If w/v or v/v comes after compound
                            if j+1<len(words) and "%" in concentration and (words[j+1] == "w/v" or words[j+1] == "v/v"):
                                concentration = concentration + " " + words[j+1]
                                del words[j-1:j+2]
                            else:
                                del words[j-1:j+1]
                        elif not concentrationBeforeCompound and j+1 < len(words) and isConcentraton(words[j+1]): # if a concentration is found (after compound)
                            compound = words[j]
                            concentration = words[j+1]
                            del words[j:j+2]
                        else: # No concentration found
                            compound = words[j]
                            concentration = None
                            del words[j]
                        compounds.append(compound)
                        compounds.append(concentration)
                        runAgain = True
                        break
        self.compounds = compounds

        if debug:
            print("Extract compounds:\n"+str(compounds)+"\n")

        return compounds

    def standardizeNames(self):
        """Standardizes the name of an individual structure"""
        runAgain = True
        while(runAgain):
            runAgain = False
            for i in range(0, len(self.compounds), 2):
                if getKey(self.compounds[i]) in compoundDictionary:
                    self.compounds[i] = compoundDictionary[getKey(self.compounds[i])]

                    # Parse multiple compounds delimited with " / "
                    if " / " in self.compounds[i]:
                        multipleCompounds = self.compounds[i].split(" / ")
                        concentration = self.compounds[i+1]
                        for c in multipleCompounds:
                            if c in set(compoundDictionary.values()):
                                self.compounds.append(c)
                                self.compounds.append(concentration)
                            else:
                                print("ERROR: Compound {} was found in multiple compound {} but is not in dictionary.".format(c, self.compounds[i]))
                                input = input("Add the key {} : {} to the dictionary? (y/n): ".format(getKey(c), c))
                                if input == "y":
                                    compoundDictionary[getKey(c)] = c
                                    print("Added key to dictionary")
                                    self.compounds.append(c)
                                    self.compounds.append(concentration)
                        del self.compounds[i:i+2]
                        runAgain = True
                        break

                    # Parse mixture compounds, that is,
                    # names of mixtures which refer to multiple chemical compounds
                    if self.compounds[i] in mixturesDictionary:
                        mixtureDict = mixturesDictionary[self.compounds[i]] # Dictionary of the compounds in the mixture
                        concentration = self.compounds[i+1]
                        if concentration != None:
                            findPercent = concentration.find("%")
                            percentUnits = ""
                            if findPercent != -1:
                                percentUnits = concentration[findPercent:] # Everything after and including '%'
                                concentration = concentration[:findPercent] # Everything berfore "%"
                        for c in mixtureDict:
                            if c in set(compoundDictionary.values()):
                                self.compounds.append(c)
                                if concentration != None:
                                    if percentUnits == "":
                                        newConcentration = str(float(concentration)*mixtureDict[c])
                                    else:
                                        # If concentration is a percent, then calculate millimolar concentration
                                        newConcentration = str(float(concentration)*mixtureDict[c]*10)
                                else:
                                    newConcentration = None
                                self.compounds.append(newConcentration)
                            else:
                                print("ERROR: Compound {} was found in mixture compound {} but is not in dictionary.".format(c, self.compounds[i]))
                                input = input("Add the key {} : {} to the dictionary? (y/n): ".format(getKey(c), c))
                                if input == "y":
                                    compoundDictionary[getKey(c)] = c
                                    print("Added key to dictionary")
                                    self.compounds.append(c)
                                    if concentration != None:
                                        if percentUnits == "":
                                            newConcentration = str(float(concentration)*mixtureDict[c])
                                        else:
                                            # If concentration is a percent, then calculate millimolar concentration
                                            newConcentration = str(float(concentration)*mixtureDict[c]*10)
                                    else:
                                        newConcentration = None
                                    self.compounds.append(newConcentration)
                        del self.compounds[i:i+2]
                        runAgain = True
                        break

                    # Remove duplicate compounds
                    if self.compounds[i] in self.compounds[:i:2]: # If compound has appeared before
                        locationOfFirstCompound = self.compounds.index(self.compounds[i], 0, i)
                        if self.compounds[i+1] == None and self.compounds[locationOfFirstCompound+1] != None: # If only the first compound has a concentrations
                            del self.compounds[i:i+2] # Remove the second compound
                            runAgain = True
                            break
                        else:
                            del self.compounds[locationOfFirstCompound:locationOfFirstCompound+2] # Remove first compound
                            runAgain = True
                            break

    def hasUnknown(self): # boolean
        """returns True if the structure has a compound in the unknownList"""
        for compound in self.compounds[::2]:
            if compound in unknownList:
                return True

    def printError(self, errorMessage, errorObject):
        """Prints an error message regarding a certain structure.
        errorMessage is a string which describes when the error occrued
        errorObject is the exception object which is caught"""
        print("--------------------\nERROR: {} for PDB ID {}.\n".format(errorMessage, self.pdbid))
        traceback.print_tb(errorObject.__traceback__)
        print(str(errorObject)+"\n")
        print("Crystallization Details: {}\n\nCompounds: {}\n--------------------\n".format(self.details, self.compounds))

def parseAllDetails(structureList, structureFile=None, searchString=None):
    """Reparses all of the details for a list of structures
    Should be called when the parseDetails function has been modified
    If search string is not None, then it will only parse Structures
    which have the search string in their details, making it faster
    If structureFile is specified, the structure list is saved to that file
    """
    global nltk
    # Make sure NLTK is imported
    if 'nltk' not in sys.modules:
        importNLTK()

    print("Parsing details with search string: \"{}\"...".format(searchString))
    count = 1

    if searchString != None:
        structureList = [structure for structure in structureList if searchString.lower() in structure.details.lower()]
        if structureList == []:
            print("No structures found with search string '{}'".format(searchString))

    print("Parsing details of {} structures...".format(len(structureList)))
    for structure in structureList:
        if count % 10000 == 0:
            print("Parsing structure {} of {}...".format(count, len(structureList)))
        try:
            structure.parseDetails()
        except Exception as e:
            structure.printError("Unable to parse details", e)
        count += 1

    if structureFile != None:
        print("Writing to structure file {}...".format(structureFile))
        writeStructures(structureList, structureFile)

def standardizeAllNames(structureList, structureFile=None): # void
    """Standardizes names of a list of compounds based on the compound Dictionary
    Also parses dictionary values which represent multiple compounds (e.g. acetic acid / sodium acetate)
    Also parses mixture compounds, which are mixtures of multiple compounds (e.g Molecular Dimensions Buffers)
    If structureFile is specified, the structure list is saved to that file
    """
    for structure in structureList:
        try:
            structure.standardizeNames()
        except Exception as e:
            structure.printError("Unable to standardize compound names", e)
    if structureFile != None:
        print("Writing to structure file {}...".format(structureFile))
        writeStructures(structureList, structureFile)

def importNLTK():
    global nltk
    try:
        import nltk
    except ModuleNotFoundError:
        print("NLTK not installed")
        sys.exit()

def downloadPunkt():
    global nltk
    try:
        import nltk
        try:
            nltk.data.find('tokenizers/punkt')
        except LookupError:
            nltk.download('punkt')

        # NEW FIX (important for newer NLTK versions)
        try:
            nltk.data.find('tokenizers/punkt_tab')
        except LookupError:
            nltk.download('punkt_tab')

    except ModuleNotFoundError:
        print("NLTK not installed")
        sys.exit()

def getSensibleStructures(structureList):
    """Returns a list of structures that make sense
    That is, all of their compounds are found in the dictionary"""
    print("Searching for sensible structures...")
    sensibleStructures = []
    for structure in structureList:
        isSensible = True
        if structure.compounds == []:
            isSensible = False
        for compound in structure.compounds[::2]:
            if compound not in set(compoundDictionary.values()):
                isSensible = False
                break
        if isSensible:
            sensibleStructures.append(structure)
    return sensibleStructures

def getUnknowns(structureList):
    """Takes a structure list, and returns a list of all structures which have an unknown compound"""
    print("Finding unknown compounds...")
    return [structure for structure in structureList if structure.hasUnknown()]


def getDatabaseSubset(structureList, pdbidList, sensibleOnly=True, structureFile=None):
    """Takes a list of PDB IDs and returns a list of structure files associated with them
    If sensibleOnly is True, then only 'sensible' structures will be included
    structureFile is an optional filename argument, if it is specified, then it will
    export the subset of structure objects to the structure file specified"""
    print("\nFetching database subset to export to {}...".format(structureFile))
    structureSubsetList = []
    pdbidList = list(set(pdbidList))
    for pdbid in pdbidList:
        structure = getStructure(structureList, pdbid)
        if structure != None:
            structureSubsetList.append(structure)

    print("Fetched {} structures for databse subset".format(len(structureSubsetList)))

    if sensibleOnly:
        structureSubsetList = getSensibleStructures(structureSubsetList)
        print("Fetched {} sensible structures from database subset".format(len(structureSubsetList)))

    if structureFile != None:
        print("Exporting structure subset to {}".format(structureFile))
        writeStructures(structureSubsetList, structureFile)

    return structureSubsetList

def getStructure(structureList, pdbid): # Structure
    """Returns a specific structure in a list based on its pdbid"""
    for structure in structureList:
        if pdbid == structure.pdbid:
            return structure
    return None

def isNumber(s): # Boolean
    try:
        float(s.replace(",",""))
        return True
    except ValueError:
        return False # Boolean

def isConcentraton(s):
    return isNumber(s) or isPercent(s) # Boolean

def isPercent(s): # Boolean
    if (s[-3:] == "w/v" or s[-3:] == "v/v") and "%" in s:
        return True
    elif len(s) > 1 and s[-1] == "%":
        try:
            float(s[:-1])
            return True
        except ValueError:
            return False
    else:
        return False # Boolean

def isCompound(word): # boolean
    """Takes a string and returns true if the string is part of a chemical compound
    Returns false if the string is a comma or a number"""
    isPunctuation = word in [",", ".", ":", ";", "-"]
    return (not isConcentraton(word) and not isPunctuation and not word == "ERROR" and word not in STOP_WORDS)

def averageNumberString(numberRange): # float
    """Takes a string such as 20-25 and averages the two (or more) numbers (ex. "20-25" --> "22.5")
    Returns None if the split doesn't produce all numbers (ex. 20%-25%)
    """
    stringSplit = numberRange.split('-')
    if all(isNumber(number) for number in stringSplit):
        numbers = [float(i.replace(",","")) for i in stringSplit]
        return sum(numbers) / len(numbers)
    else:
        return None

def wordReplacement(s, replacementDictionary, debug=False):
    """Replaces certain strings to make parsing easier
    replacementDictionary: an OrderedDict dictionary which controls the replacements to make
    If debug is True, then every step is printed"""
    details = s
    ignoredSequence = "$*"
    for key, value in replacementDictionary.items():
        details = details.replace(key.replace(ignoredSequence, ""), value)
        if debug:
            print("Replacing '{}' with '{}':".format(key.replace(ignoredSequence, ""),value))
            print(details+"\n")
    return details

def loadStructures(structureFile=STRUCTURES_FILE): # list
    """Returns a list of structures from the pickled structure file"""
    print("Loading structures from file {}...".format(structureFile))
    with open(structureFile, "rb") as f:
        return pickle.load(open(structureFile, "rb"))

def writeStructures(structureList, structureFile, count=0):
    """Writes a list of structures to a pickle file
    count keeps track of how many times the function has had to wait to write"""
    if count > 5:
        print("ERROR: Permission denied {} times when trying to write structures to {}".format(count-1, structureFile))
        return None
    try:
        with open(structureFile, "wb") as f:
            pickle.dump(structureList, f)
            if count > 0:
                print("Successfully wrote structures")
            return True
    except KeyboardInterrupt:
        print("Keyboard interrupt detected. Writing to file...")
        with open(structureFile, "wb") as f:
            pickle.dump(structureList, f)
        sys.exit()
    except PermissionError:
        if count == 0:
            print("Permission denied when attempting to write structures. Waiting and trying again...")
        sleep(3)
        nextTry = writeStructures(structureList, structureFile, count=count+1)
        return nextTry

if __name__ == "__main__":
    structureList = loadStructures(STRUCTURES_FILE)
    parseAllDetails(structureList, searchString=None, structureFile=STRUCTURES_FILE)
    standardizeAllNames(structureList, structureFile=STRUCTURES_FILE)
   