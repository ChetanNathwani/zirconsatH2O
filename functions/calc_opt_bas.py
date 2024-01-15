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
