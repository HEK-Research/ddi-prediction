
import os # Portable way of using operating system dependent functionality
# Documentation: https://docs.python.org/3/library/os.html
import csv # Implements classes to read and write tabular data in CSV format
# Documentation: https://docs.python.org/3/library/csv.html
import gzip # Reads and writes gzip-format files,
# Documentation: https://docs.python.org/3/library/gzip.html
import zipfile # Reads and writes zip-format files
# Documentation: https://docs.python.org/3/library/zipfile.html
import collections # Implements specialized container datatypes providing alternatives to Python's general purpose built-in containers
# Documentation: https://docs.python.org/3/library/collections.html
import json #  Jsons is a library that allows you to serialize your plain old Python objects to readable json (dicts or strings) and deserialize them back
# Documentation: https://docs.python.org/3/library/json.html
import xml.etree.ElementTree as ET # Built in python library for parsing xml files
# Documentation:
# XML: Extensible Markup Language
import pandas as pd # Python library offering high-performance, easy-to-use data structures and data analysis tools
# Documentation: https://pandas.pydata.org/docs/user_guide/index.html

from dataclasses import dataclass, asdict # Better way to represent structured data?

"""
The parsing logic of this code:

1) Open a gzip XML file
2) Parse the root element (<drugbank>)
3) Loop through each <drug> child
4) For EACH drug, loop through FIELD_DEFINITIONS
5) Use th expath + extractor to find and extract the value
6) Build Drug object with all extractd fields.

"""



FIELD_DEFINITIONS = { # Simply add a field defintion so there is no need to adjust parsing function
    'name': {
        'xpath': '{ns}name',
        'extractor': 'text' 
    },
    'smiles': {
        'xpath': "{ns}calculated-properties/{ns}property",
        'extractor': 'property_value',  # custom logic
        'filter_key': 'SMILES' #only extracts this value from the properties 
    },
    'groups': {
        'xpath': '{ns}groups/{ns}group',
        'extractor': 'list'
    }, 
    'drug_type': {
        'xpath': '{ns}',  # This is a placeholder, will use attribute on root element
        'extractor': 'attribute',
        'attribute_name': 'type'    
    },
    'atc_code': {
        'xpath': '{ns}atc-codes/{ns}atc-code',
        'extractor': 'atc_code'
    }

}

# Dataclass = the blueprint for a house
# Dictionary of Drug Objects = the neighborhood of houses built from the blueprint
# FIELD_DEFINIITIONS = configurable intractions on how to extract each field

@dataclass # A lightweight class designed to hold data. Automatically generates boilerplate code, IDE can automatically check if has extracted the wrong type 
class Drug:
    drugbank_id: str # A drug will always need a drugbank ID
    name: str # A drug will always need a name
    smiles: str = None
    inchi: str = None
    groups: list = None
    drug_type: str = None
    atc_code: dict = None

xml_path = r"C:\Users\ashto\OneDrive - Eastern Connecticut State University\Project 5. Data\raw\drugbank_all_full_database_V5.1.14.xml.zip" # MUST CHANGE PATH DEPENDING ON WHERE ACCESSING XML ZIP FILE
# Maybe make an user interface to be able to input own xml file, and how many drugs it wants to extract

    
class DrugBankParser: # Main class which will encapsulate the workflow
    
    def __init__(self, xml_path):
        self.xml_path = xml_path
        self.ns = '{http://www.drugbank.ca}'
    
    def parse_drugs_from_xml(self, limit=None, required_fields=None):
        """Parse all drugs from the DrugBank XML file and return dictionary of Drug objects
        
        Args:
            limit: Maximum number of drugs to parse (None = all)
            required_fields: List of field names that must be present (None = all fields required)
                            Example: ['name', 'smiles', 'groups']
        """
        # Any global initializations
        drugs = {}
        count = 0
        drugs_skipped_count = 0

        # Default: require all fields if not specified
        if required_fields is None:
            required_fields = ['drugbank_id', 'name', 'smiles', 'groups', 'drug_type', 'atc_code']
       
        # Open the xml file - handle both .zip and .gz formats
        try:
            if self.xml_path.endswith('.zip'):
                # For .zip files, find the XML file inside
                with zipfile.ZipFile(self.xml_path, 'r') as zip_file:
                    # Find the first .xml file in the zip
                    xml_files = [f for f in zip_file.namelist() if f.endswith('.xml')]
                    if not xml_files:
                        print("ERROR: No XML file found in zip archive")
                        return drugs
                    xml_filename = xml_files[0]
                    print(f"Found XML file: {xml_filename}")
                    with zip_file.open(xml_filename) as xml_file:
                        tree = ET.parse(xml_file)
                        root = tree.getroot()
            else:
                # For .gz files
                with gzip.open(self.xml_path, 'rb') as xml_file:
                    tree = ET.parse(xml_file)
                    root = tree.getroot()
        except Exception as e:
            print(f"ERROR opening file: {e}")
            return drugs
       
        # Loop through each drug element
        for drug_element in root:
            # Stop if we've reached the limit
            if limit and count >= limit:
                print(f"Reached limit of {limit} drugs")
                break
            
            # Extract the primary drugbank ID
            drug_id = drug_element.findtext(f"{self.ns}drugbank-id[@primary='true']")
            if not drug_id:
                continue
            
            # Initialize the drug data dictionary with the ID
            drug_data = {'drugbank_id': drug_id}
            count += 1
            
            # Loop through FIELD_DEFINITIONS and extract each field
            for field_name, field_config in FIELD_DEFINITIONS.items():
                extractor_type = field_config.get('extractor')
                
                # Call the appropriate extractor function based on the type
                if extractor_type == 'text':
                    value = self.extract_text_value(drug_element, field_config)
                elif extractor_type == 'property_value':
                    value = self.extract_filtered_property(drug_element, field_config)
                elif extractor_type == 'list':
                    value = self.extract_element_list(drug_element, field_config)
                elif extractor_type == 'attribute':
                    value = self.extract_attribute(drug_element, field_config)
                elif extractor_type == 'atc_code':
                    value = self.extract_atc_code(drug_element, field_config)
                else:
                    value = None # Error handling so it keeps the value at None
                
                # Add the extracted value to the drug data dictionary
                drug_data[field_name] = value
            
            # Check if all required fields are present
            if any(drug_data.get(field) is None or drug_data.get(field) == "" for field in required_fields):
                print(f"Skipping {drug_id} - missing required fields")
                drugs_skipped_count += 1
                continue
            
            # Create a Drug object with all the extracted data
            try:
                drug = Drug(**drug_data)
                drugs[drug_id] = drug
            except TypeError as e:
                print(f"Error creating Drug object for {drug_id}: {e}")
                continue
            
        return drugs
    
    def extract_text_value(self, element, config):
        """Extract simple text value from an element specified by xpath"""
        xpath = config.get('xpath', '').format(ns=self.ns)
        found_element = element.find(xpath)
        if found_element is not None:
            return found_element.text
        return None
    
    def extract_filtered_property(self, element, config):
        """Search through properties, filter by key, and return the value"""
        xpath = config.get('xpath', '').format(ns=self.ns)
        filter_key = config.get('filter_key')
        
        # Find all property elements
        for prop_element in element.findall(xpath):
            # Look for the 'kind' child element that matches our filter_key
            kind_elem = prop_element.find(f"{self.ns}kind")
            if kind_elem is not None and kind_elem.text == filter_key:
                # Once found, get the 'value' child element
                value_elem = prop_element.find(f"{self.ns}value")
                if value_elem is not None:
                    return value_elem.text
        
        return None
    
    def extract_element_list(self, element, config):
        """Gather multiple elements into a list of text values"""
        xpath = config.get('xpath', '').format(ns=self.ns)
        elements = element.findall(xpath)
        # Return list of text from each element, filtering out None values and stripping whitespace
        return [elem.text.strip() for elem in elements if elem.text and elem.text.strip()] if elements else []
    
    def extract_attribute(self, element, config):
        """Get attribute value from element"""
        attribute_name = config.get('attribute_name')
        if attribute_name:
            return element.get(attribute_name)
        return None

    def extract_atc_code(self, element, config):
        """Extract ATC code with code attribute and all level descriptions
        
        Returns a dict like:
        {
            'code': 'B01AE02',
            'levels': ['Direct thrombin inhibitors', 'ANTITHROMBOTIC AGENTS', ...]
        }
        """
        xpath = config.get('xpath', '').format(ns=self.ns)
        atc_element = element.find(xpath)
        
        if atc_element is None:
            return None
        
        # Get the code attribute
        code = atc_element.get('code')
        
        # Get all level text values
        levels = []
        for level_elem in atc_element.findall(f"{self.ns}level"):
            if level_elem.text:
                levels.append(level_elem.text.strip())
        
        # Return dict with code and levels
        return {
            'code': code,
            'levels': levels
        } if code else None
        
        
# Example usage:
if __name__ == "__main__":
    # Create parser instance
    parser = DrugBankParser(xml_path)
    
    # Parse all drugs from XML - specify which fields are required
    # Example: only require name and smiles
    drugs_dict = parser.parse_drugs_from_xml(
        limit=15, 
        required_fields=['name', 'atc_code']
    )
    
    # Or use default (all fields required):
    # drugs_dict = parser.parse_drugs_from_xml(limit=15)
    
    print(f"\nSuccessfully parsed {len(drugs_dict)} drugs")
    
    # Access individual drugs
    for drug_id, drug in list(drugs_dict.items()):
        print(f"\nID: {drug.drugbank_id}")
        print(f"Name: {drug.name}")
        print(f"Drug Type: {drug.drug_type}")
        print(f"Groups: {drug.groups}")
        print(f'ATC: {drug.atc_code}')
    
