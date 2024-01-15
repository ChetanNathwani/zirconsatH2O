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
