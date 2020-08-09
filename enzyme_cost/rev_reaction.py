""" definition of the class enzyme """

class rev_reaction:
    """ class to store the modes """
    def __init__(self, reaction_name, original_flux, cost_rev_flux, cost_delta, gr_reversed, cost_per_gr, duplicates, orig_direction):

        self.reaction_name = reaction_name
        self.original_flux = original_flux
        self.cost_rev_flux = cost_rev_flux
        self.cost_delta = cost_delta
        self.gr_reversed = gr_reversed
        self.cost_per_gr = cost_per_gr
        # duplicates = the number of reactions that where taken together for this one reaction, INCLUDING the original one
        self.duplicates = duplicates
        # store the direction of the original flux
        self.orig_direction = orig_direction
        
    