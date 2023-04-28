import pandas as pd

def toCSV(data, nameOfFile = "../data/data.csv"):
    df = pd.read_csv(nameOfFile, 
                     sep = ';',
                     keep_default_na = False,
                     na_values = ['N/A'],
                     header = 0,
                     dtype = {'Sub-goal section': str, 'SNR': "Int64", 'R': "Int64", 'N': 'Int64'},
                     encoding='latin-1')
    Tcontinuum = []
    Tcontinuum_comment = []
    Tcore = []
    Tcore_comment = [] 
    for line in data:
        Tcontinuum.append(line[0])
        Tcontinuum_comment.append(line[1])
        Tcore.append(line[2])
        Tcore_comment.append(line[3])
    df["T continuum"] = Tcontinuum
    df["T continuum comment"] = Tcontinuum_comment
    df["T core"] = Tcore
    df["T core comment"] = Tcore_comment
    nameOfFile = nameOfFile[:-4] + "with_transmissions.csv"
    df.to_csv(nameOfFile, sep = ';', index = False)
    print("CSV file saved in path: " + nameOfFile)

def toTransmissions(data):
    df = pd.DataFrame(data, columns = ["T continuum", "T continuum comment", "T core", "T core comment"])
    df.to_csv("./transmissions.csv", sep = ';', index = False)
    print("CSV file saved in path: " + "./transmissions.csv")
    