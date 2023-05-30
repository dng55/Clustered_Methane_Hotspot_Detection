def yhat_predict(observation,model_coeff):
    """
    Applying MLR Model:
    yhat = ß0 + ß1*month.f6 + ß2*month.f7 + ß3*month.f8 + ß4*month.f9 +ß5*NDVI + ß6*NDWI + ß7*NDMI + ß8*temp 
           + [ß9,ß10,ß11,ß12]*month.f6,f7,f8,f9*NDVI + [ß13,ß14,ß15,ß16]*month.f6,f7,f8,f9*NDWI + [ß17,ß18,ß19,ß20]*month.f6,f7,f8,f9*NDMI + [ß21,ß22,ß23,ß24]*month.[f6,f7,f8,f9]*temp
    Inputs:
        model = dataframe of landsat 8 map to be used as x-data
        model_coeff = MLR model coefficients
    """
    import numpy as np
    x_var = model_coeff['x_var']
    fixed = []
    interaction = []
    for idx,x_name in enumerate(x_var):
        # x variable coefficient
        b = model_coeff['x'][idx] 
        
        # Fixed effects
        if ".f" not in x_name and x_name != '(Intercept)': 
            X = observation[x_name] # All observation data of this x variable
            fixed.append(X*b)

    # Month factor
    if any("month" in string for string in model_coeff['x_var']): # If month.f exists as a predictor in model
        f6 = model_coeff[model_coeff['x_var'] == 'month.f6'].index[0] # month.f6 index
        f7 = model_coeff[model_coeff['x_var'] == 'month.f7'].index[0] # month.f7 index
        f8 = model_coeff[model_coeff['x_var'] == 'month.f8'].index[0] # month.f8 index
        f9 = model_coeff[model_coeff['x_var'] == 'month.f9'].index[0] # month.f9 index
        factor_map = {5:0,6:model_coeff['x'][f6],7:model_coeff['x'][f7],8:model_coeff['x'][f8],9:model_coeff['x'][f9]}
    else:
        factor_map = {5:0,6:0,7:0,8:0,9:0}
    factor1 = [factor_map[item] for item in observation['month.f']]
    factor = factor1
    # Year factor
    factor_map = {2021:0}
    if any("year" in string for string in model_coeff['x_var']): # If year.f exists as a predictor in model
        for idx,term in enumerate(model_coeff['x_var']):
            if "year" in term and ":" not in term:
                year = int(term.split('.f')[1])
                factor_map[year] = model_coeff['x'][idx]
        factor2 = [factor_map[item] for item in observation['year.f']]
        factor = np.add(factor1,factor2)

    # All interaction terms (interact1 = month.f, interact2 = year.f)
    interact1_terms = [term for string in model_coeff['x_var'] for term in string.split() if ":" in term and "month" in term]
    interact2_terms = [term for string in model_coeff['x_var'] for term in string.split() if ":" in term and "year" in term]
    # All interaction terms' index (1 = month.f, 2 = year.f)
    interact1_terms_idx = [(index) for index, string in enumerate(model_coeff['x_var']) 
                                 if ':' in string and 'month' in string]
    interact2_terms_idx = [(index) for index, string in enumerate(model_coeff['x_var']) 
                                 if ':' in string and 'year' in string]
    
    # Populating dict that can call interaction coefficient by the term name.
    interact_dict = {'month.f':{},'year.f':{}}
    # populating ß for month.f interaction
    for i,term in enumerate(interact1_terms):
        interact_dict['month.f'][term]=model_coeff['x'][interact1_terms_idx[i]]
    # populating ß for year.f interaction
    for i,term in enumerate(interact2_terms):
        interact_dict['year.f'][term]=model_coeff['x'][interact2_terms_idx[i]]
    # All fixed terms that get interacted with.
    interact1_fixed = np.unique([string.split(':')[1] for string in model_coeff['x_var'] if ":" in string and "month" in string])
    interact2_fixed = np.unique([string.split(':')[1] for string in model_coeff['x_var'] if ":" in string and "year" in string])
    
    # Iterating through each i in xi
    interaction = []
    # 1st factor: month.f
    for i,month in enumerate(observation['month.f']):
        total_term = 0
        # Month factor
        if month != 5: # if not intercept month
            for fixed_term in interact1_fixed:
                coeff_name = 'month.f'+str(month)+':'+fixed_term
                current_term = interact_dict['month.f'][coeff_name] * observation[fixed_term][i]
                total_term += current_term 
        # 2nd factor: year.f if it exists as a predictor
        if 'year.f' in observation.keys():
            year = observation['year.f'][i]
            if year != min(observation['year.f']): #if not intercept year
                
                for fixed_term in interact2_fixed:
                    coeff_name = 'year.f'+str(year)+':'+fixed_term
                    current_term = interact_dict['year.f'][coeff_name] * observation[fixed_term][i]
                    total_term += current_term
                
        interaction.append(total_term)
    
    # Summing
    summed_fixed = [sum(element) for element in zip(*fixed)]

    yhat_pred = [sum(element)+model_coeff['x'][0] for element in zip(summed_fixed,factor,interaction)]

    return yhat_pred

def yhat_predict_1factor(observation,model_coeff):
    """
    Applying MLR Model:
    yhat = ß0 + ß1*month.f6 + ß2*month.f7 + ß3*month.f8 + ß4*month.f9 +ß5*NDVI + ß6*NDWI + ß7*NDMI + ß8*temp 
           + [ß9,ß10,ß11,ß12]*month.f6,f7,f8,f9*NDVI + [ß13,ß14,ß15,ß16]*month.f6,f7,f8,f9*NDWI + [ß17,ß18,ß19,ß20]*month.f6,f7,f8,f9*NDMI + [ß21,ß22,ß23,ß24]*month.[f6,f7,f8,f9]*temp
    Inputs:
        model = dataframe of landsat 8 map to be used as x-data
        model_coeff = MLR model coefficients
    """
    x_var = model_coeff['x_var']
    fixed = []
    interaction = []
    for idx,x_name in enumerate(x_var):
        # x variable coefficient
        b = model_coeff['x'][idx] 
        
        # Fixed effects
        if "month" not in x_name and x_name != '(Intercept)': 
            X = observation[x_name] # All observation data of this x variable
            fixed.append(X*b)

    # Factor
    if any("month" in string for string in model_coeff['x_var']): # If month.f exists as a predictor
        f6 = model_coeff[model_coeff['x_var'] == 'month.f6'].index[0] # month.f6 index
        f7 = model_coeff[model_coeff['x_var'] == 'month.f7'].index[0] # month.f7 index
        f8 = model_coeff[model_coeff['x_var'] == 'month.f8'].index[0] # month.f8 index
        f9 = model_coeff[model_coeff['x_var'] == 'month.f9'].index[0] # month.f9 index
        factor_map = {5:0,6:model_coeff['x'][f6],7:model_coeff['x'][f7],8:model_coeff['x'][f8],9:model_coeff['x'][f9]}
    else:
        factor_map = {5:0,6:0,7:0,8:0,9:0}
    factor = [factor_map[item] for item in observation['month.f']]

    # Interaction - Can't get away with not doing a for loop
    # All interaction terms
    interact_terms = [term for string in model_coeff['x_var'] for term in string.split() if ":" in term]
    # All interaction terms' index
    interact_terms_idx = [(index) for index, string in enumerate(model_coeff['x_var']) for b, char in enumerate(string) if char == ":"]
    # Populating dict that can call interaction coefficient by the term name.
    interact_dict = {}
    for i,term in enumerate(interact_terms):
        interact_dict[term]=model_coeff['x'][interact_terms_idx[i]]
    # All fixed terms that get interacted with.
    interact_fixed = np.unique([term.split(':')[1] for string in model_coeff['x_var'] for term in string.split() if ":" in term])
    
    # Iterating through each i in xi
    interaction = []
    for i,month in enumerate(observation['month.f']):
        total_term = 0
        if month != 5:
            for fixed_term in interact_fixed:
                coeff_name = 'month.f'+str(month)+':'+fixed_term
                current_term = interact_dict[coeff_name] * observation[fixed_term][i]
                total_term += current_term
        interaction.append(total_term)
    
    # Summing
    summed_fixed = [sum(element) for element in zip(*fixed)]

    yhat_pred = [sum(element)+model_coeff['x'][0] for element in zip(summed_fixed,factor,interaction)]

    return yhat_pred

