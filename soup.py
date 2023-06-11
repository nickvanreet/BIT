# importing the libraries
from bs4 import BeautifulSoup
import requests

url="https://tandem.bu.edu/trf/output/tmpve8meldr.2.7.7.80.10.50.500.1.txt.html#1--788,177,4.4,177,4"

# Make a GET request to fetch the raw HTML content
html_content = requests.get(url).text

# Create a BeautifulSoup object
soup = BeautifulSoup(html_content, 'html.parser')

# Find the <pre> tag containing the sequences
pre_tag = soup.find('pre')

# Extract the text within <pre> tag
pre_text = pre_tag.get_text()

# Extract the full sequence
sequence_start = pre_text.index('Indices:') + len('Indices:')
sequence_end = pre_text.index('\n', sequence_start)
sequence_lines = pre_text[sequence_start:sequence_end].strip().split('\n')
sequence = ''.join(line.split()[-1] for line in sequence_lines)

# Write the full sequence to a file
with open('full.fasta', 'w') as full_file:
    full_file.write('>Full Sequence\n')
    full_file.write(sequence)

print('Full sequence saved successfully.')