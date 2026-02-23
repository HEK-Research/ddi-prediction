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
import re  # For regex pattern matching in interaction descriptions

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
import re # For regular expression matching in interaction descriptions

import random




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
                            required_fields: Optional[List[str]] = None, require_interaction: Optional[bool] = True) -> Dict[str, Dict]:
        """Parse drugs using streaming iterparse - much faster for large files
        
        Args:
            limit: Maximum number of drugs to parse (None = all)
            required_fields: List of field names that must be present
            require_interaction: If True, only include drugs with interactions
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
                        for drug_data, skip_reason in self._iterparse_drugs(xml_file, required_fields, limit, require_interaction=require_interaction):
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
                    for drug_data, skip_reason in self._iterparse_drugs(xml_file, required_fields, limit, require_interaction=require_interaction):
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
                    for drug_data, skip_reason in self._iterparse_drugs(xml_file, required_fields, limit, require_interaction=require_interaction):
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
    def _iterparse_drugs(self, xml_file, required_fields, limit, require_interaction=True):
        """Generator that yields drug data using iterparse for memory efficiency"""
        # Use iterparse() to stream through the XML file one element at a time
        # events=['end'] means we process elements when they close (complete)
        for event, elem in ET.iterparse(xml_file, events=['end']):
            # Check if this element is a 'drug' element (skip everything else)
            if elem.tag != f'{self.ns}drug':
                continue
            # If require_interaction is True, check if this drug has any interactions before processing
            if require_interaction:
                interactions_elem = elem.find(f"{self.ns}drug-interactions")
                if interactions_elem is None or len(interactions_elem.findall(f"{self.ns}drug-interaction")) == 0:
                    # If no interactions found, skip this drug
                    yield None, "no_interactions"
                    elem.clear()  # Clear element to save memory before continuing
                    continue
            # Extract all fields for this drug (returns a dictionary)
            drug_data = self._extract_drug_data(elem, required_fields)
            
            # Check if all required fields have values (not None)
            if drug_data and not any(drug_data.get(field) is None 
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
    
    # Random sampling method - separate from iterparse
    def random_sample(self, drugs: Dict[str, Dict], sample_size: Optional[int] = None, 
                     sample_rate: Optional[float] = None, num_samples: Optional[int] = None, send_to_csv: Optional[bool] = True) -> Dict[str, Dict]:
        """
        Randomly sample drugs from a parsed drugs dictionary
        
        Args:
            drugs: Dictionary of drugs {drugbank_id: drug_data}
            sample_size: Number of drugs to randomly select (e.g., 100)
            sample_rate: Fraction of drugs to select (e.g., 0.1 for 10%)
                        Only used if sample_size is None
            num_samples: Number of independent random samples to generate
            send_to_csv: If True, save samples to CSV; if False, return first sample
        
        Returns:
            Dictionary of randomly sampled drugs (or first sample if multiple generated)
            
        Examples:
            # Get 100 random drugs
            sample = parser.random_sample(drugs, sample_size=100)
            
            # Get 10% of drugs
            sample = parser.random_sample(drugs, sample_rate=0.1)
            
            # Get 10 different random samples of 100 drugs each
            parser.random_sample(drugs, sample_size=100, num_samples=10, send_to_csv=True)
        """

        #file_directory = input('Specify file directory: ') ... for now just hardcoding to raw data directory
        file_directory = r'C:\Users\ashto\ddi-prediction\data\raw'

        if not drugs: # ensures we have the correct dictionary structure and it's not empty
            print("Warning: No drugs available to sample from. Either incorrect structure or empty dictionary.")
            return {}
        if num_samples is None:
            num_samples = 1
        
        # Validate that at least one sampling method is provided
        if sample_size is None and sample_rate is None:
            print("Warning: Neither sample_size nor sample_rate provided. Returning all drugs.")
            if not send_to_csv:
                return drugs
            else:
                self.random_sample_csv(drugs, 0)
                return {}
        
        # Generate multiple samples
        for sample_num in range(num_samples):
            # Determine actual sample size
            if sample_size is not None:
                actual_sample_size = min(sample_size, len(drugs))
            else:
                # sample_rate specified
                sample_rate_clamped = max(0.0, min(1.0, sample_rate)) # Clamp sample_rate to [0, 1]
                actual_sample_size = max(1, int(len(drugs) * sample_rate_clamped)) # Ensure at least 1 drug is sampled if sample_rate > 0
            
            # Get random sample of keys from the drugs dict
            sampled_keys = random.sample(list(drugs.keys()), actual_sample_size)
            # Create new dict with only sampled drugs
            sampled_drugs = {key: drugs[key] for key in sampled_keys}
            
            # Save to CSV or return first sample
            if send_to_csv:
                self.random_sample_csv(sampled_drugs, sample_num, actual_sample_size, file_directory)
            elif sample_num == 0:
                return sampled_drugs
        
        return {} if send_to_csv else sampled_drugs
        
    def random_sample_csv(self, drugs, sample_number, actual_sample_size, file_directory): 
        """Save a drug sample to CSV file"""
        df = pd.DataFrame.from_dict(drugs, orient='index')
        df.to_csv(f'{file_directory}/drug_random_sample_#{sample_number}_sample_size_{actual_sample_size}.csv', index=True)
        print(f"Saved sample {sample_number} with {len(df)} drugs")
        


### UNDER CONSTRUCTION ### 

class ExtractNegativeInteractions: 
    
    def __init__(self, xml_path: str):
        """Initialize the extractor with path to DrugBank XML file"""
        self.xml_path = xml_path
        self.ns = '{http://www.drugbank.ca}'
    
    def extract_interactions(self, confirmed_only: bool = True, negative_only: bool = True, limit: Optional[int] = None) -> pd.DataFrame:
        """
        Extract drug-drug interactions from DrugBank
        
        Args:
            confirmed_only: If True, only extract interactions marked as "confirmed"
            negative_only: If True, only extract interactions with negative/harmful effects
            limit: If set, limits the number of drugs to process
        
        Returns:
            DataFrame with columns: [drug_1_id, drug_1_name, drug_2_id, drug_2_name, description]
        """
        interactions = [] # List to hold all extracted interactions
        count = 0 # Counter for how many drugs have been processed (for progress tracking)
        drugs_with_interactions = 0
        
        try:
            # Open the XML file (handling zip, gz, or plain xml) with helper method
            print(f"Opening file: {self.xml_path}")
            xml_file = self.open_xml_file()
            print("File opened successfully")
            
            # Parse interactions using iterparse for memory efficiency
            for event, elem in ET.iterparse(xml_file, events=['end']):
                # Process only 'drug' elements
                if elem.tag != f'{self.ns}drug':
                    elem.clear()
                    continue
                
                # Extract drug ID
                drug_id = elem.findtext(f"{self.ns}drugbank-id[@primary='true']")
                drug_name = elem.findtext(f"{self.ns}name")
                
                if not drug_id: # Saving memory by skipping drugs without ID immediately
                    elem.clear()
                    continue
                
                # Extract interactions for this drug using helper method
                drug_interactions = self.extract_drug_interactions(elem, drug_id, drug_name, confirmed_only=confirmed_only, negative_only=negative_only)
                
                if drug_interactions:
                    drugs_with_interactions += 1
                    interactions.extend(drug_interactions)
                
                count += 1
               
                if limit and count >= limit: # If limit is set to a value, exit the parsing loop after processing that many drugs
                    break
                if count % 1000 == 0:
                    filter_type = "negative" if negative_only else "all"
                    print(f"Processed {count} drugs, found {len(interactions)} {filter_type} interactions so far...")
                
                # Clear element to save memory
                elem.clear()
            
            xml_file.close() if hasattr(xml_file, 'close') else None
            print(f"File parsing completed")
        
        except Exception as e:
            print(f"ERROR processing file: {e}")
            import traceback
            traceback.print_exc()
            return pd.DataFrame()
        
        # Create DataFrame from interactions
        if interactions:
            df = pd.DataFrame(interactions)
            print(f"\nTotal processed: {count} drugs")
            print(f"Drugs with interactions: {drugs_with_interactions}")
            print(f"Total interactions extracted: {len(df)}")
            return df
        else:
            print(f"No interactions found after processing {count} drugs")
            print(f"Drugs with interactions: {drugs_with_interactions}")
            return pd.DataFrame()
    
    def open_xml_file(self):
        """Open XML file handling zip, gz, or plain xml formats"""
        if self.xml_path.endswith('.zip'):
            zip_file = zipfile.ZipFile(self.xml_path, 'r')
            xml_files = [f for f in zip_file.namelist() if f.endswith('.xml')]
            if not xml_files:
                raise ValueError("No XML file found in zip archive")
            return zip_file.open(xml_files[0])
        elif self.xml_path.endswith('.gz'):
            return gzip.open(self.xml_path, 'rb')
        else:
            return open(self.xml_path, 'rb')
    
    def extract_drug_interactions(self, 
                                  drug_elem, 
                                  drug_id: str, 
                                  drug_name: str,
                                  confirmed_only: bool = True, 
                                  negative_only: bool = True) -> List[Dict]:
        """
        Extract all interactions for a single drug
        
        Args:
            drug_elem: XML element for the drug
            drug_id: DrugBank ID of the current drug
            drug_name: Name of the current drug
            confirmed_only: Only extract confirmed interactions (must have description)
            negative_only: Only extract interactions with negative/harmful descriptions
        Returns:
            List of interaction dictionaries
        """
        interactions = [] # List to hold interactions for this drug
        
        # Keywords indicating negative/bad interactions
        bad_keywords = [
            r"increase.*risk", "adverse", "toxicity", "bleeding", "seizure",
            "contraindicated", "avoid", "serious", "severe", "dangerous",
            "reduced effect", "decreased effect", "inhibit", "toxin",
            "warn", "caution", "monitor", "harmful", "side effect",
            "increase", "risk", "reduce", "decrease"
        ]
        
        # Find all drug-interactions elements
        drug_interactions_elem = drug_elem.find(f"{self.ns}drug-interactions")
        if drug_interactions_elem is None:
            return interactions
        
        # Extract each drug-interaction
        for interaction_elem in drug_interactions_elem.findall(f"{self.ns}drug-interaction"):
            # Get the interacting drug's ID
            drug_2_id = interaction_elem.findtext(f"{self.ns}drugbank-id")
            drug_2_name = interaction_elem.findtext(f"{self.ns}name")
            
            if not drug_2_id:
                continue
            
            # Get interaction description
            description = interaction_elem.findtext(f"{self.ns}description")
            
            # Filter by confirmed status if requested
            if confirmed_only:
                # Must have a description to be confirmed
                if not description:
                    continue
            
            # Filter by negative/harmful interactions if requested
            if negative_only:
                # Only include if description contains bad keywords
                if not self.is_bad_interaction(description, keywords=bad_keywords):
                    continue
            
            # Avoid duplicate pairs (only add if drug_1_id < drug_2_id lexicographically)
            # This prevents listing same pair twice in opposite order
            if drug_id > drug_2_id:
                continue
            
            interactions.append({
                'drug_1_id': drug_id,
                'drug_1_name': drug_name,
                'drug_2_id': drug_2_id,
                'drug_2_name': drug_2_name,
                'description': description,
            })
        
        return interactions

    @staticmethod
    def is_bad_interaction(description: str, keywords: List[str]) -> bool:
        """
        Check if an interaction description contains keywords indicating negative interaction
        
        DrugBank descriptions of bad interactions typically contain words like:
        - "increase risk", "adverse", "toxicity", "bleeding", "seizure"
        - "contraindicated", "avoid", "serious", "severe"
        - "reduced effect", "inhibit", "warn", "caution"
        
        Args:
            description: The interaction description text
            keywords: List of keywords (can use regex patterns) to search for
            
        Returns:
            True if description contains bad keywords, False otherwise
        """
        if not description:
            return False
        
        description_lower = description.lower()
        
        # Check if any keywords match in the description
        for keyword in keywords:
            # Try regex match first (for patterns like "increase.*risk")
            try:
                if re.search(keyword, description_lower):
                    return True
            except:
                # Fall back to simple substring match if regex fails
                if keyword.lower() in description_lower:
                    return True
        
        return False
# Example usage / main execution block
if __name__ == "__main__":
    # Set the path to the DrugBank XML file (change this to your actual file path)
    # Using forward slashes to avoid unicode escape issues with Windows paths
    xml_path = r"C:\Users\ashto\OneDrive\DDI codes\DrugBank Parsing Code\drugbank_all_full_database_V5.1.14.xml\drugbank_full_database_V5.1.14.xml"
    
    drug_parser = FastDrugBankParser(xml_path)
    print("Starting drug parsing...")
    start_time = time.time()
    drugs = drug_parser.parse_drugs_from_xml(limit=None, require_interaction=True)
    end_time = time.time()
    print(f"Finished parsing {len(drugs)} drugs in {end_time - start_time:.2f} seconds")

    print('starting random sampling...')
    sample = drug_parser.random_sample(drugs, sample_size=100, num_samples=10, send_to_csv=True)
    print(f"Random sample of {len(sample)} drugs created and saved to CSV")
