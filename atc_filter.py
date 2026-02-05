#!/usr/bin/env python
# coding: utf-8

"""
ATC Classification Filter Module
=================================

This module provides utilities to filter DrugBank compounds based on 
ATC (Anatomical Therapeutic Chemical) classification codes.

ATC codes follow a hierarchical structure:
- Level 1 (1 letter): Anatomical Main Group (e.g., A = Alimentary tract and metabolism)
- Level 2 (2 digits): Therapeutic Subgroup
- Level 3 (1 letter): Pharmacological Subgroup
- Level 4 (1 letter): Chemical Subgroup
- Level 5 (2 digits): Chemical substance

Example ATC codes:
- C09AA05: Lisinopril (Cardiovascular > ACE inhibitors)
- A10BA02: Metformin (Alimentary tract > Diabetes agents)
- N06AB06: Sertraline (Nervous system > Serotonin reuptake inhibitors)
"""

import json
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union


class ATCFilter:
    def __init__(self, df, atc_column): # Can change this later if column name is different

        self.df = df.copy() # Ensures original df is not altered
        self.atc_column = atc_column  # Column name containing ATC codes specified by the user

        # Validation: Check if the column exists
        if atc_column not in self.df.columns: # Check if the specified ATC column exists
            available = list(self.df.columns) # List of available columns
            raise ValueError( # Raise an error if the column is not found
                f"Column '{atc_column}' not found. Available columns: {available}" 
            )
        
        print(f"ATCFilter initialized with {len(df)} drugs")
        print(f"Using ATC column: '{atc_column}'")

     
    def filter_main_group(self):
        """ Filter drugs by their Anatomical Main Group (first letter of ATC code). """

        groups = {}  # Dictionary to store a list for each letter
        
        for idx, row in self.df.iterrows(): # Iterate over each row in the dataframe
            atc_code = row[self.atc_column] # Get the ATC code for the current row
            first_letter = atc_code[0]  # Get first letter of ATC code Anatomical Main Group
            
            # If this letter hasn't been seen before, create a new list
            if first_letter not in groups: 
                groups[first_letter] = [] # Create a new list for this main group
            
            # Add the drug to the appropriate list
            groups[first_letter].append(row.to_dict()) # Add the drug to the appropriate list
        
        return groups 


    def filter_therapeutic_subgroup(self):

        """ Filter drugs by their Therapeutic Subgroup (second and third characters of ATC code). """
       
       
        groups = {} 

        for idx, row in self.df.iterrows(): # Iterate over each row in the dataframe
            atc_code = row[self.atc_column] # Get the ATC code for the current row
            thera_code = atc_code[1:3]  # Get the two digit code for Theraputic Subgroup
            
         # If this letter hasn't been seen before, create a new list
            if thera_code not in groups:
                groups[thera_code] = []
            
        # Add the drug to the appropriate list
            groups[thera_code].append(row.to_dict())
        
        return groups 

        
    def filter_pharmacological_subgroup(self):

        """ Filter drugs by their Pharmacological Subgroup (fourth character of ATC code). """

        groups = {} 

        for idx, row in self.df.iterrows(): # Iterate over each row in the dataframe
            atc_code = row[self.atc_column] # Get the ATC code for the current row
            pharma_code = atc_code[3]  # Get 4th character of ATC code
            
            # If this letter hasn't been seen before, create a new list
            if pharma_code not in groups: # If this letter hasn't been seen before, create a new list
                groups[pharma_code] = [] # Create a new list for this pharmacological subgroup
            
            # Add the drug to the appropriate list
            groups[pharma_code].append(row.to_dict())
        
        return groups 



    def filter_chemical_subgroup(self):

        """ Filter drugs by their Chemical Subgroup (fifth character of ATC code). """
        
        groups = {} 

        for idx, row in self.df.iterrows(): # Iterate over each row in the dataframe
            atc_code = row[self.atc_column] # Get the ATC code for the current row
            chem_code = atc_code[4]  # Get 5th character of ATC code
            
            # If this letter hasn't been seen before, create a new list
            if chem_code not in groups: # If this letter hasn't been seen before, create a new list
                groups[chem_code] = [] # Create a new list for this chemical subgroup
            
            # Add the drug to the appropriate list
            groups[chem_code].append(row.to_dict()) # Add the drug to the appropriate list
        
        return groups 

    def filter_by_level(self, level, code): # Filter drugs by a specific ATC level and code
        
        """Filter drugs by a specific ATC level and code.
        
            level: 1, 2, 3, or 4 (ATC hierarchy level)
            code: The prefix value to filter by (e.g., 'A' for level 1, 'A10' for level 2)
        
        Returns:
            Filtered dataframe with matching drugs
        """
        # Error handling for user misinput
        if level not in [1, 2, 3, 4]:
            raise ValueError("Level must be 1, 2, 3, or 4")
        if not code:
            raise ValueError("Code cannot be empty.")
        return self.df[self.df[self.atc_column].astype(str).str.startswith(code)].copy()



# User Functions
# Make a quick user interface where a user can give a dataframe, what they want to filter the ATC codes by, and then it returns a dataframe of the filtered ATC codes with their drug bank names

def load_dataframe_from_csv():
    """Load a dataframe from a CSV file given by the user."""

    path = input("Enter CSV file path: ").strip().strip('"').strip("'")
    if not path:
        raise ValueError("CSV path cannot be empty.")
    return pd.read_csv(path), path

def get_atc_column(df):

    """Get the ATC column name from the user."""

    col = input("Enter the ATC column name: ").strip()
    if col not in df.columns:
        raise ValueError(f"Column '{col}' not found. Available columns: {list(df.columns)}")
    return col

def get_level():

    """Get the ATC hierarchy level from the user."""

    print("\nChoose ATC hierarchy level:")
    print("1 = Level 1 (1 letter)  e.g., A")
    print("2 = Level 2 (2 digits)  e.g., A10")
    print("3 = Level 3 (1 letter)  e.g., A10B")
    print("4 = Level 4 (1 letter)  e.g., A10BA")
    level = input("Enter level (1-4): ").strip()
    if level not in {"1", "2", "3", "4"}:
        raise ValueError("Level must be 1, 2, 3, or 4.")
    return int(level)

def build_prefix(level):

    """Build the ATC code prefix based on the selected level."""

    if level == 1:
        return input("Enter Level 1 code (1 letter, e.g., A): ").strip().upper()
    if level == 2:
        return input("Enter Level 2 code (3 chars, e.g., A10): ").strip().upper()
    if level == 3:
        return input("Enter Level 3 code (4 chars, e.g., A10B): ").strip().upper()
    if level == 4:
        return input("Enter Level 4 code (5 chars, e.g., A10BA): ").strip().upper()



def main():

    """Main function to run the ATC code filter."""

    print('Welcome to the ATC_Code_Filter :)')
    df, src_path = load_dataframe_from_csv()
    atc_column = get_atc_column(df)
    level = get_level()
    prefix = build_prefix(level)

    # Use ATCFilter class to filter by level and prefix
    filter_obj = ATCFilter(df, atc_column)
    result_df = filter_obj.filter_by_level(level, prefix)

    out_path = input("Enter output CSV path (or press Enter for default): ").strip()
    if not out_path:
        out_path = src_path.replace(".csv", f"_filtered_{prefix}.csv")

    result_df.to_csv(out_path, index=False)
    print(f"Saved {len(result_df)} rows to: {out_path}")

if __name__ == "__main__":

    """Run the ATC code filter."""

    main()


