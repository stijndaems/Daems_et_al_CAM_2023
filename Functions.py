#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def generateMultiphaseLeafModel(model,number_of_phases=12,
                                met2accumulate=["STARCH_p","MAL_v","aMAL_v","CIT_v","aCIT_v", "SUCROSE_v",
                                                "PHOSPHO_ENOL_PYRUVATE_c","PYRUVATE_c","ADP_D_GLUCOSE_p",
                                                "FRUCTOSE_6P_p","GLC_1_P_c","GLC_6_P_c","GLYCEROL_3P_c","CIS_ACONITATE_c"],
                                verbose=False):
    # Generate a multiphase model with 12 copies of core model
    from cobra.core import Model, Reaction, Metabolite
    ModelF = Model("Final")
    ModelF.solver="glpk"
    for i in range(1,number_of_phases+1):
        model1 = model.copy()
        model1.reactions
        for rxn in model1.reactions:
            rxn.id = rxn.id + str(i)
        model1.metabolites
        for met in model1.metabolites:
            met.id = met.id + str(i)
        ModelF = ModelF.merge(model1)

        setMet1 = set()
        for met in model1.metabolites:
            setMet1.add(met.id)
        setMet2 = set()
        for met in ModelF.metabolites:
            setMet2.add(met.id)
        if verbose:
            print("What's in model1 but not in ModelF")
            print(str(setMet1 - setMet2))

        for met in setMet1 - setMet2:
            met = model1.metabolites.get_by_id(met)
            metCopy = met.copy()
            ModelF.add_metabolites([metCopy])
    
    # Add linker reactions to connect the different model phases together
    from cobra.core import Reaction
    for met in met2accumulate:
        for i in range(1,number_of_phases+1):
            rxn = Reaction(met+str(i)+"_accumulation")
            if i != number_of_phases:
                #reference metabolite from phase 1
                met1 = ModelF.metabolites.get_by_id(met+str(i))
                #reference metabolite from phase 2
                met2 = ModelF.metabolites.get_by_id(met+str(i+1))
                #add met1 and met2 to reaction such that one is consumed and other is produced
            else:
                #reference metabolite from phase 1
                met1 = ModelF.metabolites.get_by_id(met+str(i))
                #reference metabolite from phase 2
                met2 = ModelF.metabolites.get_by_id(met+str(1))
                #add met1 and met2 to reaction such that one is consumed and other is produced
            rxn.add_metabolites({met1:-1, met2:1})
            rxn.lower_bound = 0
            #check reaction is correct
            if verbose:
                print(rxn.reaction)
            #add reaction to model
            ModelF.add_reactions([rxn])
    
    #turn off all sugar uptake
    for phase in range(1,13):
        ModelF.reactions.get_by_id("Sucrose_tx"+str(phase)).lower_bound = 0
        ModelF.reactions.get_by_id("Sucrose_tx"+str(phase)).upper_bound = 0
        ModelF.reactions.get_by_id("GLC_tx"+str(phase)).lower_bound = 0
        ModelF.reactions.get_by_id("GLC_tx"+str(phase)).upper_bound = 0
    
    #making sure no other reaction is set as objective
    for rxn in ModelF.reactions:
        rxn.objective_coefficient = 0
    
    #add reaction to drain day and night phloem from the metabolic system in 3:1 ratio and set that as the objective of the system
    PhlRxn = Reaction("Diel_phloem_export")
    PhlRxn.lower_bound = 0
    PhlRxn.upper_bound = 1000
    ModelF.add_reactions([PhlRxn])
    PhlRxn.objective_coefficient = 1

    #assuming 12h daylenght
    sunset = 6

    #introduce metabolite to represent night-time phloem
    for i in range(1,13):
        met1 = Metabolite("Phloem_e"+str(i))
        rxn = ModelF.reactions.get_by_id("Phloem_output_tx"+str(i))
        rxn.add_metabolites({met1:1})
        if i>sunset:
            PhlRxn.add_metabolites({met1:-0.25})
        else:
            PhlRxn.add_metabolites({met1:-0.75})
    
    return ModelF

def generateMetabolitePlots(ModelF,sol,data,keyDict):
    df= data
    import numpy as np
    for met_id in keyDict.keys():
        scale = "umols"
        if keyDict[met_id]==["",]:
            continue
        xlist = range(1,13)
        pred_list = list()
        tot_list = list()
        data_list= list()
        error_max = list()
        error_min = list()
        lb_list=list()
        ub_list=list()
        met = met_id
        for i in range(1,13):
            tot = 0
            lb_tot = 0
            ub_tot = 0
            for rxn in keyDict[met]:
                rxn = rxn + str(i)+"_accumulation"
                tot = tot+float(sol.fluxes[rxn])
                lb_tot = lb_tot+ModelF.reactions.get_by_id(rxn).lower_bound
                ub_tot = ub_tot+ModelF.reactions.get_by_id(rxn).upper_bound
            lb_list.append(lb_tot)
            ub_list.append(ub_tot)
            val = float(list(df[df["Hour"]==((i*2)+8)%24][met+" average"])[0])
            data_list.append(val)
            pred_list.append(tot)
            stdev = float(list(df[df["Hour"]==((i*2)+8)%24][met+" STDEV"])[0])
            error_max.append(val+stdev)
            error_min.append(val-stdev)
        if max(data_list)<0.01:
            scale = "nmols"
            pred_list = [m * 1000 for m in pred_list]
            data_list = [m * 1000 for m in data_list]
            error_min = [m * 1000 for m in error_min]
            error_max = [m * 1000 for m in error_max]
            lb_list = [m * 1000 for m in lb_list]
            ub_list = [m * 1000 for m in ub_list]
        import matplotlib.pyplot as plt
        plt.plot(xlist,pred_list,label="predicted data", color='red',zorder=1)
        plt.scatter(xlist,data_list,label="experimental data",color="black")
        for i in range(0,len(xlist)):
            plt.plot([xlist[i],xlist[i]],[error_max[i],error_min[i]],color="black")
        plt.plot(xlist,lb_list, color='black', linestyle='--', linewidth=0.75)
        plt.plot(xlist,ub_list, color='black', linestyle='--', linewidth=0.75)
        plt.legend()
        plt.xlabel("Phase")
        if scale=="nmols":
            plt.ylabel("flux ("+r'$ nmol.m^{-2}.s^{-1}$'+")")
        else:
            plt.ylabel("flux ("+r'$ µmol.m^{-2}.s^{-1}$'+")")
        start=6.5
        end=12
        step=0.5
        xlist = [i for i in np.arange(start, end+step, step)]
        plt.fill_between(xlist,65,facecolor='grey',alpha=0.3)
        plt.title(met)
        plt.ylim(0,max(max(pred_list),max(error_max))*1.1)
        plt.savefig(met, dpi=1200)
        plt.show()
    return

def checkProtonFluxes(ModelF3,tag=""):
    import pandas as pd
    #summarizing top 5 proton producing and consuming reactions during daytime and nighttime relative to the cytosol
    if tag!="":
        tag="_"+str(tag)
    
    #Daytime
    outfile = "DaytimeProducing"+tag+".csv"
    fout = open(outfile,"w")
    for i in range(1,7):
        met = ModelF3.metabolites.get_by_id("PROTON_c"+str(i))
        for rxn in met.reactions:
            if round(rxn.flux,3)!=0 and rxn.flux*rxn.metabolites[met]>0:
                fout.write(rxn.id+","+str(((rxn.flux*rxn.metabolites[met])*60*60*2)/1000)+"\n") #convert per sec to per 2h (60*60*2) --> mmol/m-2/2h
                            
    fout.close()
    
    df = pd.read_csv(outfile, header=None)
    df.rename(columns={0: 'RxnID', 1: 'DaytimeFlux'}, inplace=True)
    df['RxnID'] = df['RxnID'].map(lambda x: str(x)[:-1])
    df.to_csv(outfile, index=False)
    dff = df.groupby(['RxnID']).DaytimeFlux.sum().reset_index()
    dff.sort_values("DaytimeFlux", axis=0, ascending=False, inplace=True, na_position='first') 
    dff.to_csv(outfile, index=True)
    
    outfile = "DaytimeConsuming"+tag+".csv"
    fout = open(outfile,"w")
    for i in range(1,7):
        met = ModelF3.metabolites.get_by_id("PROTON_c"+str(i))
        for rxn in met.reactions:
            if round(rxn.flux,3)!=0 and rxn.flux*rxn.metabolites[met]<0:
                fout.write(rxn.id+","+str(((rxn.flux*rxn.metabolites[met])*60*60*2)/1000)+"\n") 
                            
    fout.close()
    
    df = pd.read_csv(outfile, header=None)
    df.rename(columns={0: 'RxnID', 1: 'DaytimeFlux'}, inplace=True)
    df['RxnID'] = df['RxnID'].map(lambda x: str(x)[:-1])
    df.to_csv(outfile, index=False)
    dff = df.groupby(['RxnID']).DaytimeFlux.sum().reset_index()
    dff.sort_values("DaytimeFlux", axis=0, ascending=True, inplace=True, na_position='first') 
    dff.to_csv(outfile, index=True)
    
    #Nighttime
    outfile = "NighttimeProducing"+tag+".csv"
    fout = open(outfile,"w")
    for i in range(7,13):
        met = ModelF3.metabolites.get_by_id("PROTON_c"+str(i))
        for rxn in met.reactions:
            if round(rxn.flux,3)!=0 and rxn.flux*rxn.metabolites[met]>0:
                fout.write(rxn.id+","+str(((rxn.flux*rxn.metabolites[met])*60*60*2)/1000)+"\n")
                            
    fout.close()
    
    df = pd.read_csv(outfile, header=None)
    df.rename(columns={0: 'RxnID', 1: 'NighttimeFlux'}, inplace=True)
    df['RxnID'] = df['RxnID'].map(lambda x: str(x)[:-1])
    df['RxnID'] = df['RxnID'].map(lambda x: x.rstrip('1'))
    df.to_csv(outfile, index=False)
    dff = df.groupby(['RxnID']).NighttimeFlux.sum().reset_index()
    dff.sort_values("NighttimeFlux", axis=0, ascending=False, inplace=True, na_position='first') 
    dff.to_csv(outfile, index=True)

    outfile = "NighttimeConsuming"+tag+".csv"
    fout = open(outfile,"w")
    for i in range(7,13):
        met = ModelF3.metabolites.get_by_id("PROTON_c"+str(i))
        for rxn in met.reactions:
            if round(rxn.flux,3)!=0 and rxn.flux*rxn.metabolites[met]<0:
                fout.write(rxn.id+","+str(((rxn.flux*rxn.metabolites[met])*60*60*2)/1000)+"\n")
                            
    fout.close()
    
    df = pd.read_csv(outfile, header=None)
    df.rename(columns={0: 'RxnID', 1: 'NighttimeFlux'}, inplace=True)
    df['RxnID'] = df['RxnID'].map(lambda x: str(x)[:-1])
    df['RxnID'] = df['RxnID'].map(lambda x: x.rstrip('1'))
    df.to_csv(outfile, index=False)
    dff = df.groupby(['RxnID']).NighttimeFlux.sum().reset_index()
    dff.sort_values("NighttimeFlux", axis=0, ascending=True, inplace=True, na_position='first') 
    dff.to_csv(outfile, index=True)
    return

def customFVA(ModelF3,rxnlist = ["Photon_tx",]):
    opt = ModelF3.slim_optimize()
    if ModelF3.solver.status!='optimal':
        print("FBA solution is not optimal")
        return    
    
    from cobra.flux_analysis import pfba
    SoF = pfba(ModelF3)
    backup = ModelF3.copy()
    
    objectives = list()
    for rxn in ModelF3.reactions:
        if rxn.objective_coefficient!=0:
            objectives.append(rxn.id)
    
    
    fva_dict = dict()
    if len(objectives)!=1:
        print("Cannot run this simple function for multiobjective models")
        return
    else:
        rxn = ModelF3.reactions.get_by_id(objectives[0])
        rxn.objective_coefficient = 0
        rxn.upper_bound = opt
        rxn.lower_bound = opt
        backup2 = ModelF3.copy()
        tempModel = backup2.copy()

        from cobra.flux_analysis.parsimonious import add_pfba
        prob = tempModel.problem
        with tempModel:
            add_pfba(tempModel, fraction_of_optimum=0)
            ub = tempModel.slim_optimize(error_value=None)
            flux_sum = prob.Variable("flux_sum", ub=1 * ub)
            flux_sum_constraint = prob.Constraint(
                tempModel.solver.objective.expression - flux_sum,
                lb=0,
                ub=0,
                name="flux_sum_constraint",
            )
        tempModel.add_cons_vars([flux_sum, flux_sum_constraint])
        
        for rxnID in rxnlist:
            tempModel.reactions.get_by_id(rxnID).objective_coefficient = 1
            maxim = round(tempModel.slim_optimize(),6)
            tempModel.reactions.get_by_id(rxnID).objective_coefficient = -1
            minim = round(-1*tempModel.slim_optimize(),6)
            tempModel.reactions.get_by_id(rxnID).objective_coefficient = 0
            fva_dict[rxnID]=(maxim, minim)
    return fva_dict

