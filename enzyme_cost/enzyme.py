""" definition of the class enzyme """

class enzyme:
    """ class to store the modes """
    def __init__(self, reaction_name, macrocomponents, cost, EC_number, turnover):
       
        # macrocomponents = dictionary including the species name/ id and the corresponding stoichiometry in this specific enzyme
        self.reaction_name = reaction_name
        self.macrocomponents = macrocomponents
        
        
        # these are in general not known by the time an enzyme is created initially
        # cost = amount of amino acids of which the enzyme is composed of (1 AA = cost of 1)
        self.cost = 0
        self.EC_number = 0
        self.turnover = 0 
        
        
        
         
        