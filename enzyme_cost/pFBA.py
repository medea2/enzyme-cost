""" run pFBA  """

import gurobi as grb
import cobra
import numpy as np
import rev_reaction
import enzyme
import pickle
import sys
import math

# put functions in separate files
def create_lp(cobra_model, parameters={}, attributes={}):
    """Creates a Gurobi LP model from a COBRA model."""

    lp = grb.Model(cobra_model.id)

    # Add variables and get constraints
    constr = {m: [] for m in cobra_model.metabolites if m._reaction}
    for r in cobra_model.reactions:
        x = lp.addVar(r.lower_bound, r.upper_bound, r.objective_coefficient,
                      'continuous', r.id)
        for m, c in r.metabolites.items():
            constr[m].append((c, x))

    lp.update()

    # Add constraints
    for m in constr:
        c = grb.LinExpr(constr[m])
        lp.addConstr(c, '=', m._bound, m.id)

    # Set parameters and attributes
    lp.setParam('OutputFlag', False)
    for k, v in parameters.items():
        lp.setParam(k, v)
    for k, v in attributes.items():
        lp.setAttr(k, v)

    lp.update()

    return lp

def make_irreversible(model, subset=set()):
    """Makes reactions in subset positive and irreversible."""

    model.update()

    if not subset:
        subset = set(model.getVars())
    else:
        subset = set(subset)

    to_reverse = []

    for r in list(subset):
        if r.ub > 0:
            # Make existing reaction irreversible
            r.lb = np.max([0, r.lb])
        if r.lb < 0:
            if r.ub > 0:
                # Create new reverse reaction
                r_rev = model.addVar(np.max([0, -r.ub]), -r.lb, r.obj, r.VType,
                                     r.varName + '_rev')
                to_reverse.append((r, r_rev))
            else:
                # Reverse existing reaction
                r.varName += '_rev'
                r.lb, r.ub = np.max([0, -r.ub]), -r.lb
                for c in model.getConstrs():
                    model.chgCoeff(c, r, -model.getCoeff(c, r))

    model.update()

    for (r, r_rev) in to_reverse:
        for c in model.getConstrs():
            model.chgCoeff(c, r_rev, -model.getCoeff(c, r))
            subset.add(r_rev)

    model.update()

    return subset

def pfba(model, rct_cost_per_flux):
    irr_model = model.copy()

    model = create_lp(irr_model)
    make_irreversible(model)
    model.setParam('OutputFlag', True)
    model.setParam('OptimalityTol', 1e-9)
    model.setParam('FeasibilityTol', 1e-9)

    # Minimize internal microbial fluxes
    if rct_cost_per_flux is None:
        inr = [v for v in model.getVars() if model.getCol(v).size() > 1]
        model.setObjective(grb.LinExpr(np.ones(len(inr)), inr))
    else:
        vars = list()
        vals = list()
        for key, val in rct_cost_per_flux.items():
            
            if not np.isinf(val):
                test = model.getVarByName(key)
                
                if test is not None:
                    # had to add this because some where nan & gurobi couldn't handle it
                    if math.isnan(val) is not True: 
                        vars.append(test)
                        vals.append(val)
                
                test_rev = model.getVarByName(key + '_rev')
                if test_rev is not None:
                    if math.isnan(val) is not True:# had to add this because some where nan
                        vars.append(test_rev)
                        vals.append(val)
        
        
        model.setObjective(grb.LinExpr(vals, vars))
       
    model.ModelSense = grb.GRB.MINIMIZE
    model.setParam('OutputFlag', True)
    model.optimize()

    return model

def extract_solution(model):
    sol_dict = {v.VarName: v.x for v in model.getVars()}
    # unify reversed fluxes
    removals = []
    for key in sol_dict:
        if key.endswith("_rev"):
            orig_key = key[:-4]
            if orig_key in sol_dict:
                sol_dict[orig_key] -= sol_dict[key]
                removals.append(key)

    return sol_dict

def get_constrained_solutions(model, rct_cost_per_flux=None):
    # run with unconstrained model
    unc_pfba = pfba(model, rct_cost_per_flux)
    unc_solution = extract_solution(unc_pfba)
    r_list = [r for r in model.reactions if 'c' in r.compartments and r.reversibility and not r.boundary]
    variability = cobra.flux_analysis.variability.flux_variability_analysis(model, reaction_list=r_list,
                                                                            fraction_of_optimum=0)
    bool_series = (variability.minimum < 0) & (variability.maximum > 0)

    # get out all reversible reactions and save their names as a list
    true_reversible = variability[bool_series]
    reactionNames = true_reversible.index.values
    listOfReactionNames = list(reactionNames)

    # for reaction in rev_reactions:
    attributes = ["lower_bound", "upper_bound"]
    solutions = dict()
    # for reaction_name in rev_reactions:
    for reaction_name in listOfReactionNames:
        reaction = model.reactions.get_by_id(reaction_name)
        solutions[reaction_name] = dict()
        for attribute in attributes:
            bound = getattr(reaction, attribute)
            setattr(reaction, attribute, 0)
            model.repair()
            unc_pfba = pfba(model, rct_cost_per_flux)
            setattr(reaction, attribute, bound)
            solutions[reaction_name][attribute] = extract_solution(unc_pfba)
            solutions[reaction_name][attribute + "_cost"] = unc_pfba.ObjVal

    return solutions

def main():
     
    # define directory
    xml_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/Bacterial-RBA-models/Escherichia-coli-K12-WT/'
    #
    # read in the model    
    model = cobra.io.read_sbml_model(xml_dir + "/data/iJO1366.xml")
    solution = model.optimize()
    # create a dictionary containing the reaction names and the corresponding cost per unit flux values
    # get reversible reactions
    with open('save_rev_rcts.p', 'rb') as pickle_file:
         rev_reactions = pickle.load(pickle_file)
    # get list of enzymes
    with open('save_enzymes.p', 'rb') as pickle_file:
        enzymes = pickle.load(pickle_file)
    
    for enz in enzymes:
        enzymes[enz].reaction_name = enzymes[enz].reaction_name[2:] # to get rid of 'R_'
    
    # calculate cost per unit flux
    # create empty dictionary and get all reaction ids
    rct_cost_per_flux = {}
    reaction_names = [r.id for r in list(model.reactions)]
    for reaction in reaction_names:
        # initialize variable for turnover number
        rct_tn = 0
        
        for enz in enzymes:
            # get the right enzyme for the reaction
            if enzymes[enz].reaction_name == reaction:
                rct_enzcomp = enzymes[enz].cost
                rct_tn = enzymes[enz].turnover
                try:
                    rct_cost_per_flux[reaction] = rct_enzcomp / rct_tn   #new = CORRECT formula
                # in case original flux = 0
                except ZeroDivisionError:
                    rct_cost_per_flux[reaction] = np.inf
    
    print(rct_cost_per_flux)
    
 #   # create random numbers for cost per flux (to check concept before having the real numbers)
 #   rct_cost_per_flux = {}
 #   for r in list(model.reactions):
 #       rct_cost_per_flux[r.id] = np.random.rand()
 #       if rct_cost_per_flux[r.id] > 0.9:
 #           rct_cost_per_flux[r.id] = np.inf
            
    
    # biomass
    biomass_rxn = model.reactions.get_by_id("BIOMASS_Ec_iJO1366_WT_53p95M")
    # biomass_rxn = model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM") #biomass rct. for E.Coli core model
    biomass_rxn.lower_bound = 0.001
    model.repair()
   
    # calculate pfba and weighted pfba solutions
   # pfba_solutions = get_constrained_solutions(model)
    weighted_solutions = get_constrained_solutions(model, rct_cost_per_flux=rct_cost_per_flux)
   
    # store the solutions as pickle files
   #  pickle.dump(pfba_solutions, open( "pfba_solutions.p", "wb" ))
    pickle.dump(weighted_solutions, open( "weighted_solutions.p", "wb" ))

    print("hurray!")

main()

    
    
    
    