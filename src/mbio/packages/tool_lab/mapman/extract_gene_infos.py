from mapman_packages import *
import os

species_dir = "/mnt/lustre/users/sanger-dev/app/database/Tool_lab/mapman/species"
species_dir = "/mnt/lustre/users/sanger-dev/app/database/Tool_lab/mapman_e44/mapman"
all_species = os.listdir(species_dir)
for species in all_species:
    single_species_dir  = os.path.join(species_dir,species)
    all_versions = os.listdir(single_species_dir)
    for version in all_versions:
        single_version_species_dir = os.path.join(single_species_dir,version)
        single_version_species_file = os.path.join(single_version_species_dir,os.listdir(single_version_species_dir)[0])
        mapping_file_name,g2b,b2g = extract_bind_id_to_genes(single_version_species_file)
        single_species_infos = dict()
        print("species :{}".format(species))
        print("version :{}".format(version))
        print("g2b :{}".format(str(g2b)))
        print("b2g :{}".format(str(b2g)))
        single_species_infos["species"] = species
        single_species_infos["version"] = version
        single_species_infos["g2b"] = g2b
        single_species_infos["b2g"] = b2g
        with open(os.path.join(single_version_species_dir,"species_infos"),"w") as f:
            json.dump(single_species_infos, f, sort_keys=True, indent=4,default=set_default)

