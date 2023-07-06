import csv
from collections import defaultdict

# Map each strain to its subspecies and abbreviation
strain_subspecies_mapping = {}
with open('strain_subspecies_ab.csv', 'r') as file:
    reader = csv.reader(file, delimiter=';')
    next(reader)  # Skip the header row
    for row in reader:
        if len(row) < 3:
            continue  # Skip rows with fewer than 3 elements
        strain, subspecies, abbreviation = [item.strip() for item in row]  # Strip leading/trailing spaces
        strain_subspecies_mapping[strain] = (subspecies, abbreviation)

# Map each identifier to a dictionary of {abb_subspec}_{strain} and counts
identifier_mapping = defaultdict(lambda: defaultdict(int))
with open('all_strains_lookup.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        identifier = row[0]
        strains = row[1].split(",")  # Split strains by comma
        for strain in strains:
            if strain in strain_subspecies_mapping:
                subspecies, abbreviation = strain_subspecies_mapping[strain]
                subspecies_strain = f"{abbreviation}_{strain}"
                identifier_mapping[identifier][subspecies_strain] += 1  # Increase the count for this subspecies

# Generate a list of all unique {abb_subspec}_{strain}
subspecies_strain_list = set()
for strain_dict in identifier_mapping.values():
    subspecies_strain_list.update(strain_dict.keys())

subspecies_strain_list = sorted(list(subspecies_strain_list))

# Write the result to a CSV file
with open('identifier_to_strain.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    # Write a header row with all {abb_subspec}_{strain}
    writer.writerow(['Identifier'] + subspecies_strain_list)
    for identifier, strain_dict in identifier_mapping.items():
        row = [identifier]
        for subspecies_strain in subspecies_strain_list:
            row.append(strain_dict.get(subspecies_strain, 0))  # If there is no count for a subspecies, put a 0
        writer.writerow(row)
