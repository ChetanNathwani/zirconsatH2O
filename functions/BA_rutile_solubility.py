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
