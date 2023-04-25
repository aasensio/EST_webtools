import pandas as pd

def toCSV(data, nameOfFile = "./measures_with_transmissions.csv"):
    df = pd.DataFrame(data, columns =
                      ["Wavelength (nm)",
                       "required d (\")",
                       "SIS Instrument",
                       "R",
                       "BP (nm)",
                       "SNR",
                       "tint (s)",
                       "Polarimetry",
                       "T continuum",
                       "T continuum comment"])
    df.to_csv(nameOfFile, sep = ';', index = False)
    print("CSV file saved in path: " + nameOfFile)
    