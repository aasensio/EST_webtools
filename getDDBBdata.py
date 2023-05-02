import pandas as pd
import numpy as np

def getDDBBdata(nameOfFile = "./measures.csv"):
    try:
        open(nameOfFile, 'r')
    except IOError:
        print("File does not exist")
        return
    
    df = pd.read_csv(nameOfFile, 
                     sep = ';',
                     keep_default_na = True,
                     na_values = 'NA',
                     header = 0,
                     dtype = {'SNR': "Int64", 'R': "Int64", "BP (nm)": "Float64", "tint (s)": "Float64"},
                     encoding='latin-1')

    df_transmissions = df[["Wavelength (nm)", "required d (\")", "SIS Instrument", "R", "BP (nm)", "SNR", "tint (s)", "Polarimetry"]]
    df_transmissions = df_transmissions.replace({np.nan: None})

    wavelengths = []
    spatialResolutions = []
    instruments = []
    spectralResolutions = []
    bandpass = []
    snr = []
    integrationTimes = []
    polarimetrys = []

    for (index, row) in df_transmissions.iterrows():
        wavelengths.append(row["Wavelength (nm)"])
        spatialResolutions.append(row["required d (\")"])
        instruments.append(row["SIS Instrument"])
        spectralResolutions.append(row["R"])
        bandpass.append(row["BP (nm)"])
        snr.append(row["SNR"])
        integrationTimes.append(row["tint (s)"])
        polarimetrys.append(row["Polarimetry"])

    joinedData = [[]]
    for (index, row) in df_transmissions.iterrows():
        joinedData.append([row["Wavelength (nm)"], 
                           row["required d (\")"],
                           row["SIS Instrument"],
                           row["R"],
                           row["BP (nm)"],
                           row["SNR"],
                           row["tint (s)"],
                           row["Polarimetry"]])
    joinedData = joinedData[1:]
    return joinedData