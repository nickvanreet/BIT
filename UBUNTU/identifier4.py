import csv
from collections import defaultdict

# Map each strain to its subspecies
strain_subspecies_mapping = {}
with open('strain_subspecies_ab.csv', 'r') as file:
    reader = csv.reader(file, delimiter=';')
    next(reader)  # Skip the header row
    for row in reader:
        if len(row) < 2:
            continue  # Skip rows with fewer than 2 elements
        strain, subspecies = [item.strip() for item in row]  # Strip leading/trailing spaces
        strain_subspecies_mapping[strain] = subspecies

# Generate a set of all unique subspecies
subspecies_set = set(strain_subspecies_mapping.values())

# Map each identifier to a dictionary of subspecies and counts
identifier_subspecies_mapping = defaultdict(lambda: defaultdict(int))
with open('all_strains_lookup.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        identifier = row[0]
        strains = row[1].split(",")  # Split strains by comma
        for strain in strains:
            if strain in strain_subspecies_mapping:
                subspecies = strain_subspecies_mapping[strain]
                identifier_subspecies_mapping[identifier][subspecies] += 1  # Increase the count for this subspecies

# Write the result to a CSV file
with open('identifier_to_subspecies.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    # Write a header row with all subspecies
    writer.writerow(['Identifier'] + sorted(list(subspecies_set)))
    for identifier, subspecies_dict in identifier_subspecies_mapping.items():
        row = [identifier]
        for subspecies in sorted(list(subspecies_set)):  # Write a row for each subspecies count
            row.append(subspecies_dict.get(subspecies, 0))  # If there is no count for a subspecies, put a 0
        writer.writerow(row)
