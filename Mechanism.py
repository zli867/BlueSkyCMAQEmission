# convert unit
convert_to_g_s = 907184.74 / 3600


def convert_to_mole_s(molecule_weight):
    # convert to g/s -> moles/s
    return convert_to_g_s / molecule_weight


# CB6
CB6_species_mapping = {
    'CO': [(1 * convert_to_mole_s(28), 'CO')],
    'SO2': [(1 * convert_to_mole_s(64), 'SO2')],
    'NH3': [(1 * convert_to_mole_s(17), 'NH3')],
    # SERA
    'PMC': [(1.2480483 * convert_to_g_s, 'PM2.5'), (-1 * convert_to_g_s, 'PM2.5')],
    # SERA
    'NO2': [(0.35 * convert_to_mole_s(46), 'NOx')],
    # SERA
    'NO': [(0.65 * convert_to_mole_s(46), 'NOx')],
    'ALD2_PRIMARY': [(0.026318 * convert_to_mole_s(44.0526), 'VOC')],
    'FORM_PRIMARY': [(0.115995 * convert_to_mole_s(30.026), 'VOC')],
    'SOAALK': [(0.021027 * convert_to_mole_s(80.5347), 'VOC')],
    'ACET': [(0.029243 * convert_to_mole_s(58.0791), 'VOC')],
    'ALD2': [(0.026318 * convert_to_mole_s(44.0526), 'VOC')],
    'ALDX': [(0.081322 * convert_to_mole_s(40.7486), 'VOC')],
    'BENZ': [(0.025761 * convert_to_mole_s(78.1118), 'VOC')],
    'CH4': [(0.247726 * convert_to_mole_s(16.0425), 'VOC')],
    'ETH': [(0.059042 * convert_to_mole_s(28.0532), 'VOC')],
    'ETHA': [(0.03161 * convert_to_mole_s(30.069), 'VOC')],
    'ETHY': [(0.014204 * convert_to_mole_s(26.0373), 'VOC')],
    'ETOH': [(0.000949 * convert_to_mole_s(46.0684), 'VOC')],
    'FORM': [(0.115995 * convert_to_mole_s(30.026), 'VOC')],
    'IOLE': [(0.012975 * convert_to_mole_s(55.6133), 'VOC')],
    'ISOP': [(0.00527 * convert_to_mole_s(68.117), 'VOC')],
    'KET': [(0.007112 * convert_to_mole_s(17.8646), 'VOC')],
    'MEOH': [(0.105412 * convert_to_mole_s(32.0419), 'VOC')],
    'OLE': [(0.092184 * convert_to_mole_s(31.8777), 'VOC')],
    'PAR': [(0.200242 * convert_to_mole_s(21.0484), 'VOC')],
    'PRPA': [(0.013175 * convert_to_mole_s(44.0956), 'VOC')],
    'TERP': [(0.011595 * convert_to_mole_s(136.234), 'VOC')],
    'TOL': [(0.03537 * convert_to_mole_s(93.0188), 'VOC')],
    'UNR': [(0.269171 * convert_to_mole_s(36.5914), 'VOC')],
    'XYLMN': [(0.007905 * convert_to_mole_s(106.165), 'VOC')],
    'PAL': [(4.6 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PCA': [(7.2 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PCL': [(0.00239 * convert_to_g_s, 'PM2.5')],
    'PEC': [(0.1093 * convert_to_g_s, 'PM2.5')],
    'PFE': [(4.45 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PK': [(0.00135 * convert_to_g_s, 'PM2.5')],
    'PMN': [(1.1 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PMOTHR': [(0.0125 * convert_to_g_s, 'PM2.5')],
    'PNA': [(0.00135 * convert_to_g_s, 'PM2.5')],
    'PNCOM': [(0.3513 * convert_to_g_s, 'PM2.5')],
    'PNH4': [(0.00341 * convert_to_g_s, 'PM2.5')],
    'PNO3': [(0.0107 * convert_to_g_s, 'PM2.5')],
    'POC': [(0.5019 * convert_to_g_s, 'PM2.5')],
    'PSI': [(1.0 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PSO4': [(0.0033 * convert_to_g_s, 'PM2.5')],
    'PTI': [(6.7 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')]
}


# SAPRC07
SAPRC07_species_mapping = {
    'CO': [(1 * convert_to_mole_s(28), 'CO')],
    # SERA
    'NO2': [(0.35 * convert_to_mole_s(46), 'NOx')],
    'NO': [(0.65 * convert_to_mole_s(46), 'NOx')],
    'NH3': [(1 * convert_to_mole_s(17), 'NH3')],
    'SO2': [(1 * convert_to_mole_s(64), 'SO2')],
    'ACYE': [(0.084 * convert_to_mole_s(26.0373), 'VOC')],
    'ALK1': [(0.1047 * convert_to_mole_s(30.069), 'VOC')],
    'ALK2': [(0.0035 * convert_to_mole_s(44.0956), 'VOC')],
    'ALK3': [(0.004108 * convert_to_mole_s(60.8495), 'VOC')],
    'ALK4': [(0.0255 * convert_to_mole_s(101.5224), 'VOC')],
    'ALK5': [(0.2443 * convert_to_mole_s(166.6683), 'VOC')],
    'APIN': [(0.0143 * convert_to_mole_s(208.3353), 'VOC')],
    'ARO1': [(0.0263 * convert_to_mole_s(166.6883), 'VOC')],
    'ARO2': [(0.422 * convert_to_mole_s(187.5018), 'VOC')],
    'B124': [(0.002927 * convert_to_mole_s(187.5018), 'VOC')],
    'BDE13': [(0.0052 * convert_to_mole_s(54.0904), 'VOC')],
    'CH4': [(0.0982 * convert_to_mole_s(16.0425), 'VOC')],
    'CRES': [(0.002397 * convert_to_mole_s(145.8347), 'VOC')],
    'ETHE': [(0.1911 * convert_to_mole_s(28.0532), 'VOC')],
    'MVK': [(0.0003175 * convert_to_mole_s(83.3341), 'VOC')],
    'NROG': [(0.0206 * convert_to_mole_s(41.6671), 'VOC')],
    'OLE1': [(0.0451 * convert_to_mole_s(74.1464), 'VOC')],
    'OLE2': [(0.0163 * convert_to_mole_s(83.8603), 'VOC')],
    'PRD2': [(0.0174 * convert_to_mole_s(125.0012), 'VOC')],
    'PRPE': [(0.0393 * convert_to_mole_s(42.0797), 'VOC')],
    'TERP': [(0.0123 * convert_to_mole_s(208.3353), 'VOC')],
    'PAL': [(4.6 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PCA': [(7.2 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PCL': [(0.00239 * convert_to_g_s, 'PM2.5')],
    'PEC': [(0.109 * convert_to_g_s, 'PM2.5')],
    'PFE': [(4.45 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PK': [(0.00135 * convert_to_g_s, 'PM2.5')],
    'PMFINE': [(0.375 * convert_to_g_s, 'PM2.5')],
    'PMN': [(1.1 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PMOTHR': [(0.012995 * convert_to_g_s, 'PM2.5')],
    'PNA': [(0.00135 * convert_to_g_s, 'PM2.5')],
    'PNCOM': [(0.351 * convert_to_g_s, 'PM2.5')],
    'PNH4': [(0.00341 * convert_to_g_s, 'PM2.5')],
    'PNO3': [(0.0107 * convert_to_g_s, 'PM2.5')],
    'POC': [(0.502 * convert_to_g_s, 'PM2.5')],
    'PSI': [(1.0 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    'PSO4': [(0.0033 * convert_to_g_s, 'PM2.5')],
    'PTI': [(6.7 * (10 ** (-4)) * convert_to_g_s, 'PM2.5')],
    # SERA
    'PMC': [(1.2480483 * convert_to_g_s, 'PM2.5'), (-1 * convert_to_g_s, 'PM2.5')]
}