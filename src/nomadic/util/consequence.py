import re
import warnings
from dataclasses import dataclass


@dataclass
class Consequence:
    """
    Encapsulate consequence information from  `bcftools csq`
    
    """
    
    csq: str
    target_gene: str
    transcript: str
    biotype: str
        
    strand: str=""
    aa_change: str=""
    nt_change: str=""
        
    @classmethod
    def from_string(cls, csq_string: str):
        """
        Parse from the output string
        
        What to do if "double"?
        Or @
        Or *
        
        """
        if csq_string == ".": # intergenic
            return cls(".", ".", ".", ".")
        
        if csq_string.startswith("@"): # compound variety, recorded elsewhere
            return cls(".", ".", ".", f"compound{csq_string}")
        
        consequences = csq_string.split(",")
        if len(consequences) > 1:
            warnings.warn(f"Found multiple consequences of variant: {csq_string}! Keeping only first.")
        
        fields = consequences[0].split("|")
        assert len(fields) >= 4, f"Failed for {csq_string}"
        assert len(fields) <= 7, f"Failed for {csq_string}" 
        
        return cls(*fields)
    
    
    def get_concise_aa_change(self) -> str:
        """
        Get a more concise encoding of the amino acid change
        
        """
        
        if not self.aa_change:
            return self.aa_change
        
        AAs = "ARNDCEQGHILKMFPSTWYV"
        stop = "\*"
        match = re.match(f"([0-9]+)([{AAs}|{stop}]+)>[0-9]+([{AAs}|{stop}]+)", self.aa_change)
        
        if match is None:
            if self.csq != "synonymous":
                warnings.warn(f"Unable to parse AA change for: {self}. Returning as-is.")
            return self.aa_change
        
        pos, from_aa, to_aa = match.groups()
        
        return f"{from_aa}{pos}{to_aa}"
    

