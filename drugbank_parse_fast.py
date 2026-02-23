"""
Fast DrugBank XML Parser - Optimized for speed
Improvements:
1. Uses iterparse() for streaming (not loading entire XML into memory)
2. Single-pass extraction (no loop through FIELD_DEFINITIONS for each drug)
3. Returns dictionaries directly (faster than dataclass creation)
4. Clears elements immediately to reduce memory footprint
5. Optional lxml fallback for even faster parsing
"""

# Import the standard Python XML ElementTree parser
import xml.etree.ElementTree as ET
# Import zipfile module to handle .zip compressed XML files
import zipfile
# Import gzip module to handle .gz compressed XML files
import gzip
# Import type hints for better code documentation (Dict, List, Optional are type annotations)
from typing import Dict, List, Optional
import time  # For timing the parsing process (optional, for performance measurement)
import pandas as pd  # For saving parsed data to CSV (optional, for output handling)

# Try to import faster C-based XML parsing libraries (optional but recommended)
try:
    # Attempt to import cElementTree (faster C implementation of ElementTree)
    import xml.etree.cElementTree as ET
    USING_CXML = True  # Flag to track if C version was loaded
except ImportError:
    # If cElementTree not available, try lxml (much faster, requires pip install)
    try:
        from lxml import etree as ET_lxml
        USE_LXML = True  # Flag to track if lxml was loaded
    except ImportError:
        # If neither available, falls back to standard ElementTree
        USE_LXML = False
        USING_CXML = False



# Main parser class optimized for speed using streaming XML parsing
class FastDrugBankParser:
    """Fast parser using iterparse for streaming XML processing"""
    
    # Constructor - initializes the parser with the path to the XML file
    def __init__(self, xml_path):
        # Store the path to the XML file (can be .zip, .gz, or plain .xml)
        self.xml_path = xml_path
        # Store the DrugBank XML namespace - all DrugBank elements are prefixed with this
        self.ns = '{http://www.drugbank.ca}'
    
    # Main parsing method - reads the XML file and extracts drug information
    def parse_drugs_from_xml(self, limit: Optional[int] = None, 
                            required_fields: Optional[List[str]] = None) -> Dict:
        """Parse drugs using streaming iterparse - much faster for large files
        
        Args:
            limit: Maximum number of drugs to parse (None = all)
            required_fields: List of field names that must be present
        
        Returns:
            Dictionary of drug dictionaries {drugbank_id: drug_data}
        """
        # If no required fields specified, set defaults (all important fields)
        if required_fields is None:
            required_fields = ['drugbank_id', 'name', 'smiles', 'groups', 'drug_type', 'atc_codes']
        
        # Dictionary to store all parsed drugs (key=drugbank_id, value=drug_data)
        drugs = {}
        # Counter for successful drug parses
        count = 0
        # Counter for drugs that were skipped (missing required fields)
        drugs_skipped = 0
        
        try:
            # Check if the file is a zip archive
            if self.xml_path.endswith('.zip'):
                # Open the zip file in read mode
                with zipfile.ZipFile(self.xml_path, 'r') as zip_file:
                    # Get list of all XML files inside the zip
                    xml_files = [f for f in zip_file.namelist() if f.endswith('.xml')]
                    # Check if any XML files were found
                    if not xml_files:
                        print("ERROR: No XML file found in zip archive")
                        return drugs
                    # Get the first XML file name (usually only one)
                    xml_filename = xml_files[0]
                    print(f"Found XML file: {xml_filename}")
                    # Open the XML file from inside the zip
                    with zip_file.open(xml_filename) as xml_file:
                        # Loop through each drug yielded by the generator function
                        for drug_data, skip_reason in self._iterparse_drugs(xml_file, required_fields, limit):
                            # If drug_data is None, the drug was skipped (missing fields)
                            if drug_data is None:
                                drugs_skipped += 1
                                continue
                            # Add the drug to our dictionary using its drugbank_id as key
                            drugs[drug_data['drugbank_id']] = drug_data
                            # Increment successful parse counter
                            count += 1
                            # Stop if we've reached the limit
                            if limit and count >= limit:
                                break
            elif self.xml_path.endswith('.gz'):
                # File is a .gz compressed file
                # Open the gzip file in binary read mode
                with gzip.open(self.xml_path, 'rb') as xml_file:
                    # Loop through each drug yielded by the generator function
                    for drug_data, skip_reason in self._iterparse_drugs(xml_file, required_fields, limit):
                        # If drug_data is None, the drug was skipped (missing fields)
                        if drug_data is None:
                            drugs_skipped += 1
                            continue
                        # Add the drug to our dictionary using its drugbank_id as key
                        drugs[drug_data['drugbank_id']] = drug_data
                        # Increment successful parse counter
                        count += 1
                        # Stop if we've reached the limit
                        if limit and count >= limit:
                            break
            else:
                # File is a plain .xml file (uncompressed)
                # Open the XML file in binary read mode
                with open(self.xml_path, 'rb') as xml_file:
                    # Loop through each drug yielded by the generator function
                    for drug_data, skip_reason in self._iterparse_drugs(xml_file, required_fields, limit):
                        # If drug_data is None, the drug was skipped (missing fields)
                        if drug_data is None:
                            drugs_skipped += 1
                            continue
                        # Add the drug to our dictionary using its drugbank_id as key
                        drugs[drug_data['drugbank_id']] = drug_data
                        # Increment successful parse counter
                        count += 1
                        # Stop if we've reached the limit
                        if limit and count >= limit:
                            break

                
        
        # Catch any errors that occur during file opening or parsing
        except Exception as e:
            print(f"ERROR opening/parsing file: {e}")
            return drugs
        
        # Print summary of parsing results
        print(f"Parsed {count} drugs, skipped {drugs_skipped}")
        # Return the dictionary of all parsed drugs
        return drugs
    
    # Generator function that yields drugs one at a time (memory efficient!)
    def _iterparse_drugs(self, xml_file, required_fields, limit):
        """Generator that yields drug data using iterparse for memory efficiency"""
        # Use iterparse() to stream through the XML file one element at a time
        # events=['end'] means we process elements when they close (complete)
        for event, elem in ET.iterparse(xml_file, events=['end']):
            # Check if this element is a 'drug' element (skip everything else)
            if elem.tag != f'{self.ns}drug':
                continue
            
            # Extract all fields for this drug (returns a dictionary)
            drug_data = self._extract_drug_data(elem, required_fields)
            
            # Check if all required fields have values (not None or empty string)
            if drug_data and not any(drug_data.get(field) is None or drug_data.get(field) == "" 
                                     for field in required_fields):
                # If all required fields present, yield the drug data and None (no skip reason)
                yield drug_data, None
            else:
                # If missing fields, yield None and the skip reason
                yield None, "missing_fields"
            
            # **CRITICAL FOR MEMORY EFFICIENCY**: Delete the element from memory after processing
            # This prevents the entire XML tree from accumulating in RAM
            elem.clear()
    
    # Extract all drug fields from a single drug element
    def _extract_drug_data(self, drug_element, required_fields=None) -> Optional[Dict]:
        """Extract all drug fields in a single efficient pass
        
        Instead of looping through FIELD_DEFINITIONS multiple times,
        this extracts all fields at once.
        """
        # Find the primary drugbank ID (attribute primary='true' filters to the main ID)
        drug_id = drug_element.findtext(f"{self.ns}drugbank-id[@primary='true']")
        # If no primary ID found, skip this drug (return None)
        if not drug_id:
            return None
        
        # Create a new dictionary to hold all drug information
        drug_data = {'drugbank_id': drug_id}
        
        # Extract the drug name (simple text field)
        drug_data['name'] = self._get_text(drug_element, f'{self.ns}name')
        # Extract SMILES (filtered from calculated-properties - only gets SMILES type)
        drug_data['smiles'] = self._get_filtered_property(drug_element, 'SMILES')
        # Extract InChI (filtered from calculated-properties - only gets InChI type)
        drug_data['inchi'] = self._get_filtered_property(drug_element, 'InChI')
        # Extract list of groups (multiple elements converted to a list)
        drug_data['groups'] = self._get_list(drug_element, f'{self.ns}groups/{self.ns}group')
        # Extract drug type (get from XML attribute 'type')
        drug_data['drug_type'] = drug_element.get('type')
        # Extract ATC code information (special parsing for nested structure)
        drug_data['atc_codes'] = self._get_atc_code(drug_element)
        
        # Return the complete drug data dictionary
        return drug_data
    
    # Helper method: extract a simple text value from an XML element
    def _get_text(self, element, xpath: str) -> Optional[str]:
        """Get text from element by xpath"""
        # Use findtext() to get the text content of an element by its path
        found = element.findtext(xpath)
        # Return the text with whitespace stripped, or None if not found/empty
        return found if found and found.strip() else None
    
    # Helper method: find a property by type and return its value
    def _get_filtered_property(self, element, property_name: str) -> Optional[str]:
        """Extract property value by filtering on 'kind' field"""
        # Loop through all property elements in the calculated-properties section
        for prop in element.findall(f"{self.ns}calculated-properties/{self.ns}property"):
            # Get the 'kind' field which tells us what type of property this is
            kind = prop.findtext(f"{self.ns}kind")
            # Check if this property's kind matches what we're looking for (e.g., 'SMILES')
            if kind == property_name:
                # Get the 'value' field which contains the actual data
                value = prop.findtext(f"{self.ns}value")
                # Return the value with whitespace stripped, or None if empty
                return value.strip() if value else None
        # If we didn't find the property, return None
        return None
    
    # Helper method: extract multiple elements into a list
    def _get_list(self, element, xpath: str) -> Optional[List[str]]:
        """Get list of text values from elements"""
        # Use findall() to get all matching elements, extract text from each
        # Only include elements that have non-empty text
        items = [elem.text.strip() for elem in element.findall(xpath) 
                 if elem.text and elem.text.strip()]
        # Return the list if it has items, otherwise return None
        return items if items else None
    
    # Helper method: extract all ATC codes (multiple can exist per drug)
    def _get_atc_code(self, element) -> Optional[List[str]]:
        """Extract all ATC codes"""
        # Find all atc-code elements in the atc-codes section
        atc_elems = element.findall(f"{self.ns}atc-codes/{self.ns}atc-code")
        # If no ATC codes found, return None
        if not atc_elems:
            return None
        
        atc_codes = []
        # Loop through each ATC code element
        for atc_elem in atc_elems:
            # Get the 'code' attribute (e.g., 'A01AB01')
            code = atc_elem.get('code')
            # Skip if no code found
            if not code:
                continue
            # Add this ATC code to the list
            atc_codes.append(code)
        
        # Return the list of ATC codes, or None if empty
        return atc_codes if atc_codes else None



# Example usage / main execution block
if __name__ == "__main__":
    # Set the path to the DrugBank XML file (change this to your actual file path)
    xml_path = r"C:\Users\ashto\OneDrive\DDI codes\DrugBank Parsing Code\drugbank_all_full_database_V5.1.14.xml\drugbank_full_database_V5.1.14.xml"
    
    start = time.perf_counter()  # Start timing the parsing process


    # Create an instance of the fast parser with the xml file path
    parser = FastDrugBankParser(xml_path)
    # Parse drugs from the XML file (limit to 100 drugs for this example)
    drugs = parser.parse_drugs_from_xml(limit = 100)

    df = pd.DataFrame.from_dict(drugs, orient='index')  # Convert the drugs dictionary to a DataFrame

    try: 
        df.to_csv(r"C:\Users\ashto\OneDrive\DDI codes\DrugBank Parsing Code\drugbank_parsed.csv", index=False)
    except Exception as e:
        print(f"Error saving to CSV: {e}")

    end = time.perf_counter()  # End timing
    print(f"\nParsing completed in {end - start:.2f} seconds")  # Print how long parsing took


    # ALTERNATIVE: Use the VERY fast parser (with lxml if available)
    # Uncomment the lines below to use lxml instead:
    # parser = VeryFastDrugBankParser(xml_path)
    # drugs = parser.parse_drugs_from_xml(limit=50)
    
    # Print how many drugs were successfully parsed
    print(f"\nSuccessfully parsed {len(drugs)} drugs")
    
    # Loop through the first 5 drugs and print their information
    """for drug_id, drug in list(drugs.items())[:5]:
        # Print the drug's ID
        print(f"\nID: {drug_id}")
        # Print the drug's name (use .get() to safely handle missing values)
        print(f"Name: {drug.get('name')}")
        # Print the drug's type (e.g., 'small molecule', 'protein')
        print(f"Drug Type: {drug.get('drug_type')}")
        # Print the drug's groups (therapeutic categories)
        print(f"Groups: {drug.get('groups')}")
        # Print the drug's ATC code information (includes code and hierarchy levels)
        print(f'ATC: {drug.get("atc_code")}')"""

