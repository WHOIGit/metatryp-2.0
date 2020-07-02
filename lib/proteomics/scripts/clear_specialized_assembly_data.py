"""
name: clear_specialized_assembly_data.py

usage: clear_specialized_assembly_data.py [--taxon-ids-file=taxon_ids_file] [--taxon-ids= id1 id2 ..]


description: This script clears data in the redundancy DB for a given set of taxon ids.
Assumptions:
    - The redundancy db has already been created and is writeable.
"""

"""
Imports and setup.
"""

from proteomics.services.clear_specialized_assembly_data import ClearSpecializedAssemblyDataTask
import argparse
import logging


"""
Process arguments.
"""
argparser = argparse.ArgumentParser(description=(
    'Clear data in the DB for the given Specialized Assemblies.'))
argparser.add_argument('--specialized-assembly-names-file', help=(
    'File containing a list of specialized assembly ids to clear, one id per line'))
argparser.add_argument('--specialized-assembly-names', nargs='*', help=(
    'List of specialized assembly ids to clear.'))
"""
Main method.
"""
def main():
    args = argparser.parse_args()

    logger = logging.getLogger('clear_specialized_assemblies_assemblies')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    if args.specialized_assembly_names_file:
        with open(args.specialized_assembly_names_file) as f:
            specialized_assembly_ids = f.readlines()
    else:
        specialized_assembly_ids = args.specialized_assembly_names

    if not specialized_assembly_ids:
        argparser.error("Provide specialized_assembly_ids via the --specialized-assembly-ids option, or the" 
                        " --specialized_assembly-ids-file option")

    # Confirm deletion w/ the user.
    print("You are about to delete the following taxons:\n")
    print("\n".join(specialized_assembly_ids), "\n")
    confirmation = input("Type 'yes' and hit enter if this is really "
                             "what you want to do: ")
    if confirmation != 'yes':
        logger.info("You did not enter 'yes', quitting. Nothing has been done.")
        exit()

    # Run the clearing task.
    task = ClearSpecializedAssemblyDataTask(
        logger=logger,
        specialized_assembly_ids=specialized_assembly_ids,
    )
    task.run()
    logger.info("Done.")

if __name__ == '__main__':
    main()
