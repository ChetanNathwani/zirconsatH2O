import numpy as np
import pandas as pd
import periodictable as pt

def BA_rutile_solubility(mole_oxide, t):
    # Calculate TiO2 activity using the rutile solubility model of Borisov & Aranovich (2020; Chemical Geology)
    t = t + 273.15
    
    # Recalculate to mole fraction oxide
    mole_oxide = mole_oxide.divide(mole_oxide.sum(axis = 1), axis = 0)

    logXTiO2 = ( -1417*mole_oxide['SiO2']**2 
                 -18377*mole_oxide['Al2O3']**2 
                 -7532*mole_oxide['FeO']**2 
                 -3686*mole_oxide['MgO']**2 
                 -4338*mole_oxide['SiO2'] *
                 mole_oxide['Al2O3'] 
                 -3593*mole_oxide['SiO2'] * 
                 mole_oxide['CaO'] - 3433 *
                 mole_oxide['SiO2']*mole_oxide['Na2O'] 
                - 9886*mole_oxide['Al2O3']*
                mole_oxide['MgO']+26304*mole_oxide['Al2O3']*mole_oxide['CaO']+
                23456*mole_oxide['MgO']*mole_oxide['Na2O']-2156
               )/t+1
    
    XTiO2 = 10**logXTiO2
    act_coeff = 1/XTiO2
    aTiO2 = act_coeff * mole_oxide['TiO2']
    return(aTiO2)

def calc_opt_bas(df):
    
    # Method is from Mills (1993) appendix
    
    elemental_opt = {'K2O':1.4, 
                 'Na2O':1.15,
                 'CaO':1.0,
                 'MgO':0.78,
                 'MnO':1.0,
                 'FeO':1.0,
                 'Al2O3':0.60,
                 'SiO2':0.48,
                 'P2O5':0.33}
    
    # df is a pandas dataframe
    
    if 'Fe2O3' in df.columns:
        df['FeO'] = df['FeO'] + df['Fe2O3']/1.1111
    
    SiO2_mol = df['SiO2']/pt.formula("SiO2").mass
    K2O_mol = df['K2O']/pt.formula("K2O").mass
    CaO_mol = df['CaO']/pt.formula("CaO").mass
    MgO_mol = df['MgO']/pt.formula("MgO").mass
    MnO_mol = df['MnO']/pt.formula("MnO").mass
    FeO_mol = df['FeO']/pt.formula("FeO").mass
    Al2O3_mol = df['Al2O3']/pt.formula("Al2O3").mass
    
    total_mol = (SiO2_mol+K2O_mol+CaO_mol+MgO_mol+
                 MnO_mol+FeO_mol+Al2O3_mol)
      
    term1 = (
        ((SiO2_mol/total_mol)*2*elemental_opt['SiO2']) +
        ((K2O_mol/total_mol)*1*elemental_opt['K2O']) +
        ((CaO_mol/total_mol)*1*elemental_opt['CaO']) +
        ((MgO_mol/total_mol)*1*elemental_opt['MgO']) +
        ((MnO_mol/total_mol)*1*elemental_opt['MnO'])+
        ((FeO_mol/total_mol)*1*elemental_opt['FeO']) +
        ((Al2O3_mol/total_mol)*3*elemental_opt['Al2O3'])
         )
    
    term2 = (
        ((SiO2_mol/total_mol)*2) +
        ((K2O_mol/total_mol)*1) +
        ((CaO_mol/total_mol)*1) +
        ((MgO_mol/total_mol)*1) +
        ((MnO_mol/total_mol)*1)+
        ((FeO_mol/total_mol)*1) +
        ((Al2O3_mol/total_mol)*3)
         )
        
    opt_bas = term1/term2

    return(opt_bas)

def fo2buffer(T, P, delta, buff):
    #REFERENCE: Reviews in Mineralogy Volume 25
    #T is in C
    #P is in MPa
    #delta is the delta value from NNO or QFM
    #buff - text indicating NNO/FMQ/QFM

    T += 273.15
    P *= 10.0

    if buff in ['FMQ','QFM']:
        FO2 = 10**((-25096.3/T)+8.735+(0.110*(P-1)/T)+delta)
    else:
        FO2 = 10**((-24930/T)+9.36+(0.046*(P-1)/T)+delta)

    return FO2

def calc_Fe2O3_FeO(comp_dict, t, p, buffer, delta):
    # Calculates new composition with the fo2 (and thus Fe valence) taken into account for a new temp or pressure
    # T in K, P in Mpa, composition as a dictionary, buffer as NNO/FMQ, delta as delta from the buffer
    t = t + 273.15
    fo2 = fo2buffer(t, p, delta, buffer)
    mm = {'SiO2': 60.08, 'TiO2': 79.866, 'Al2O3': 101.96, 'Fe2O3': 159.69, 'FeO': 71.844, 'MnO': 70.9374, 'MgO': 40.3044, 'CaO': 56.0774, 'Na2O':
          61.9789, 'K2O': 94.2, 'P2O5': 283.89, 'H2O': 18.01528, 'CO2': 44.01, 'S': 32.065, 'Cl': 35.453, 'F': 18.9984}
    P = p*1e6
    param = [0.196, 11492.0, -6.675, -2.243, -1.828, 3.201, 5.854, 6.215, -3.36, 1673, -7.01e-7, -1.54e-10, 3.85e-17]
    mol = {k: comp_dict[k]/mm[k] for k in mm.keys() & comp_dict} # Divide composition by molar mass to get moles
    mol_frac = {k:v/sum(mol.values(), 0.0) for k,v in mol.items()} # Divide by sum to get mole fraction

    # Use Kress & Carmichael Eq 7 to calculate Fe2O3/FeO ratio
    Fe2O3FeO = np.exp((param[0]*np.log(fo2))+(param[1]/t)+param[2]+
                      ((param[3]*mol_frac['Al2O3'])+(param[4]*mol_frac['FeO'])+
                       (param[5]*mol_frac['CaO'])+(param[6]*mol_frac['Na2O'])+
                       (param[7]*mol_frac['K2O']))+
                      param[8]*(1-param[9]/t-np.log(t/param[9])+param[10]*(p/t)+param[11]*((t-param[9])*(p/t)+param[12]*(p**2)/t))
                     ) 

    FeOtotal = ((comp_dict['FeO'] + comp_dict['Fe2O3']/1.1111)/mm['FeO'])/sum(mol.values(), 0.0) # Calculate mol fraction FeO total
    XFe2O3 = Fe2O3FeO * FeOtotal/(2*Fe2O3FeO+1) # Mole Fe2O3 in the melt
    XFeO = FeOtotal/(1+2*Fe2O3FeO) # Mole FeO in the melt
    molFe3 = XFe2O3*2
    molFe2 = XFeO
    molFe3_Fetot_ratio = molFe3/(molFe2+molFe3)
    mol_frac['FeO'] = XFeO
    mol_frac['Fe2O3'] = XFe2O3
    xmw = {k: mol_frac[k]*mm[k] for k in mm.keys() & mol_frac} # Divide composition by molar mass to get moles
    wt = {k:v*(100/sum(xmw.values(), 0.0)) for k,v in xmw.items()} # Multiply to get the normalised wt%
    new_comp = comp_dict
    new_comp['FeO'] = wt['FeO']
    new_comp['Fe2O3'] = wt['Fe2O3']
    return(new_comp)
