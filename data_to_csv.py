import pandas as pd

def toCSV(data, nameOfFile = "../data/20230516_OP.csv", exportRoute = "../data/data_with_transmissions.csv"):
    df = pd.read_csv(nameOfFile, 
                     sep = ';',
                     keep_default_na = False,
                     na_values = ['N/A'],
                     header = 0,
                     dtype = {'Wavelength (nm)': "Float64", 'Sub-goal section': str, 'SNR': "Int64", 'R': "Int64", 'N': 'Int64'},
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
    df.to_csv(exportRoute, sep = ';', index = False)
    print("CSV file saved in path: " + exportRoute)

def toTransmissions(data):
    df = pd.DataFrame(data, columns = ["T continuum", "T continuum comment", "T core", "T core comment"])
    df.to_csv("./transmissions.csv", sep = ';', index = False)
    print("CSV file saved in path: " + "./transmissions.csv")
    