import json
import csv

# open json file small_dataset/ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl
# convert this to a csv table. 
# Define the input and output file paths
input_file = 'small_dataset/ncbi_dataset/ncbi_dataset/data/assembly_data_report.jsonl'
output_file = 'small_dataset/ncbi_dataset/ncbi_dataset/data/assembly_data_report.csv'

def get_minimum_values(data):
    keys_to_fetch = [
        'accession',
        'assemblyInfo > biosample > accession',
        'assemblyInfo > biosample > strain',
        'annotationInfo > pipeline',
        'annotationInfo > softwareVersion',
        'annotationInfo > releasedate',
        'annotationInfo > stats > geneCounts > nonCoding',
        'annotationInfo > stats > geneCounts > proteinCoding',
        'annotationInfo > stats > geneCounts > psuedogene',
        'annotationInfo > stats > geneCounts > total',
        'assemblyInfo > bioprojectAccession',
        'assemblyStats > contigN50',
        'assemblyStats > gcPercent',
        'assemblyStats > numberOfContigs',
        'assemblyStats > totalSequenceLength',
        'averageNucleotideIdentity > bestAniMatch > ani',
        'averageNucleotideIdentity > bestAniMatch > organismName',
        'averageNucleotideIdentity > bestAniMatch > assembly',
        'averageNucleotideIdentity > bestAniMatch > assemblyCoverage',
        'checkmInfo > checkmMarkerSet',
        'checkmInfo > checkmMarkerSetRank',
        'checkmInfo > species',
        'checkmInfo > completeness',
        'checkmInfo > checkmVersion',
        'organism > organismName',
        'organism > taxid',
    ]

    def get_nested_value(data, key):
        if key == 'checkmInfo > species':
            if not data.get('checkmInfo', {}).get('species'):
                return data.get('checkmInfo', {}).get('checkmSpeciesTaxId', {})
            return data.get('checkmInfo', {}).get('species', {})
        keys = key.split(' > ')
        for k in keys:
            data = data.get(k, {})
        return data if data else None
    def cleanup_key(key):
        if len(key.split(' > ')) > 2:
            # only use first and last key
            return key.split(' > ')[0] + '_' + key.split(' > ')[-1]
        return key.replace(' > ', '_').replace(' ', '_').lower()
    

    return {cleanup_key(key): get_nested_value(data, key) for key in keys_to_fetch}



# Open the input JSONL file and the output CSV file
with open(input_file, 'r') as jsonl_file, open(output_file, 'w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    
    # Write the header to the CSV file
    header_written = False
    
    # Process each line in the JSONL file
    for line in jsonl_file:
        data = json.loads(line)
        values = get_minimum_values(data)
        # Write the header based on the keys of the first JSON object
        if not header_written:
            header = values.keys()
            csv_writer.writerow(header)
            header_written = True
        
        # Write the data to the CSV file
        csv_writer.writerow(values.values())