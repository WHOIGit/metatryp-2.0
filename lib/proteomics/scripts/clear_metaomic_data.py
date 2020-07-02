"""
name: clear_metaomic_data.py

usage: clear_metaomic_data.py [--taxon-ids-file=taxon_ids_file] [--taxon-ids= id1 id2 ..]


description: This script clears data in the redundancy DB for a given set of taxon ids.
Assumptions:
    - The redundancy db has already been created and is writeable.
"""

"""
Imports and setup.
"""

from proteomics.services.clear_metaomic_data import ClearMetaomicDataTask
import argparse
import logging


"""
Process arguments.
"""
argparser = argparse.ArgumentParser(description=(
    'Clear data in the DB for the given Meta-omic Assemblies.'))
argparser.add_argument('--metaomic-ids-file', help=(
    'File containing a list of Meta-omic Assembly ids to clear, one id per line'))
argparser.add_argument('--metaomic-ids', nargs='*', help=(
    'List of meta-omic assembly ids to clear.'))
"""
Main method.
"""
def main():
    args = argparser.parse_args()

    logger = logging.getLogger('clear_metaomic_assemblies')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)

    if args.metaomic_ids_file:
        with open(args.metaomic_ids_file) as f:
            metaomic_ids = f.readlines()
    else:
        metaomic_ids = args.metaomic_ids

    if not metaomic_ids:
        argparser.error("Provide metaomic_ids via the --metaomic-ids option, or the" 
                        " --metaomic-ids-file option")

    # Confirm deletion w/ the user.
    print("You are about to delete the following taxons:\n")
    print("\n".join(metaomic_ids), "\n")
    confirmation = input("Type 'yes' and hit enter if this is really "
                             "what you want to do: ")
    if confirmation != 'yes':
        logger.info("You did not enter 'yes', quitting. Nothing has been done.")
        exit()

    # Run the clearing task.
    task = ClearMetaomicDataTask(
        logger=logger,
        metaomic_ids=metaomic_ids,
    )
    task.run()
    logger.info("Done.")

if __name__ == '__main__':
    main()
