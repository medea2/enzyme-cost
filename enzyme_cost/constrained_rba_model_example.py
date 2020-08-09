"""Solve standard RBA problem."""

# python 2/3 compatibility
from __future__ import division, print_function
from rba.xml.targets import TargetGroup, TargetReaction

# global imports
import sys
import os.path

# package imports
import rba


def main():
    if len(sys.argv) < 2:
        xml_dir = '/Users/medeafux/Desktop/ETH/Masterarbeit/Bacterial-RBA-models/Escherichia-coli-K12-WT'
    else:
        xml_dir = sys.argv[1]

    if len(sys.argv) >= 3:
        output_dir = sys.argv[2]
    else:
        output_dir = xml_dir

    # load model, build matrices and solve
    print('Model building from XML files ...')
    model = rba.RbaModel.from_xml(xml_dir)

    print('Starting iterative RBA resolution (default model)...')
    results = model.solve()
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(output_dir)
    
    
### co
    print('Starting iterative RBA resolution (constrained model)...')
    # create a  new TargetGroup object = group of target fluxes
    target_group = TargetGroup('mode_targets')
    # create a new TargetReaction = Association of a target value with a chemical reaction (use id of wanted reaction)
    target_reaction = TargetReaction('R_GHMT2r')
    # set upper bound of target reaction to zero (in this case)
    target_reaction.upper_bound = 'zero'
    # add the target flux value of the target reaction to the list of reaction fluxes in the target group object
    target_group.reaction_fluxes.append(target_reaction)
    # add the target reaction flux to the model (why is this step needed?)
    model.targets.target_groups.append(target_group)
    results = model.solve()
    print('Optimal growth rate is {}.'.format(results.mu_opt))
    results.write(output_dir)
    

if __name__ == '__main__':
    main()
