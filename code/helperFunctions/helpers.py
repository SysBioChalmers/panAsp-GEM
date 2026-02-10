import medusa
import cobra
from cobra.flux_analysis.loopless import loopless_solution
import pandas as pd
import numpy as np
from pathlib import Path
from cobra.io import read_sbml_model
from medusa.flux_analysis import flux_balance
from copy import deepcopy
from pickle import load

from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# my_functions.py
def add(a, b):
    return a + b

def reactionToComp(model, reaction_id, new_compartment):
    reaction = model.reactions.get_by_id(reaction_id)
    for metabolite in list(reaction.metabolites):  # Use list() to safely iterate
            # Check if the metabolite is in the current compartment
            if metabolite.compartment != new_compartment:
                # Create the new metabolite ID by replacing the current compartment with the new one
                new_metabolite_id = metabolite.id[:-3] + '[' + new_compartment + ']'  # Remove the last 3 characters (compartment) and add the new one
                
                try:
                    # Check if the new metabolite exists in the model
                    new_metabolite = model.metabolites.get_by_id(new_metabolite_id)
                except KeyError:
                    # If the metabolite does not exist, create it
                    new_metabolite = cobra.Metabolite(id=new_metabolite_id, compartment=new_compartment, name=metabolite.name)
                    model.add_metabolites([new_metabolite])
                    print(f"Created new metabolite: {new_metabolite_id}")
                
                # Get the coefficient (stoichiometry) of the metabolite in the reaction
                coefficient = reaction.metabolites[metabolite]
                
                # Add the new metabolite to the reaction with the same coefficient
                reaction.add_metabolites({new_metabolite: coefficient})
                
                # Remove the old metabolite (old compartment) from the reaction
                # Achieved by adding the same metabolite to the other side of the equation
                reaction.add_metabolites({metabolite: coefficient * -1})
    return model

def geneToComp(model, gene_id, new_compartment):
     
     # Get the reactions catalyzed by the gene
    reactions_by_gene = []
    for reaction in model.reactions:
        for gene in reaction.genes:
            if gene.id == gene_id:  # Check by the gene ID
                reactions_by_gene.append(reaction)
                break  # No need to check further once we find the match

    for reaction in reactions_by_gene:
        reaction_id = reaction.id
        model = reactionToComp(model, reaction_id, new_compartment)
    
    return model